##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline eqtl
===========================

Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import CGATCore.Experiment as E
from CGATCore import Pipeline as P
import re

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS["projectsrc"] = os.path.dirname(__file__)
#for key, value in PARAMS.iteritems():
#    print "%s:\t%s" % (key,value)


# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.

PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    'genesets',
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


# ---------------------------------------------------
# Specific pipeline tasks

@follows(mkdir("readgroups.dir"))
@transform("input_files.dir/*.bam",formatter(),r"readgroups.dir/{basename[0]}.readgroups.bam")
def add_read_groups(infile, outfile):
    platform = PARAMS["platform"]
    groupsample = P.snip(os.path.basename(infile), ".bam")
    statement = '''java -Xmx8G -jar /shared/sudlab1/General/apps/bio/picard-tools-1.135/picard.jar
                   AddOrReplaceReadGroups
                   I=%(infile)s
                   O=%(outfile)s
                   RGLB=lib1
                   RGPL=%(platform)s
                   RGPU=unit1
                   RGSM=%(groupsample)s'''

    job_memory = "16G"
    P.run(statement)



@follows(mkdir("deduped.dir"))
@transform(add_read_groups,
           regex(r"readgroups.dir/(.+).readgroups.bam"),
           r"deduped.dir/\1.bam")
def dedup_bams(infile, outfile):
    '''Use MarkDuplicates to mark dupliceate reads'''
    job_memory = "16G"

    tempfile=P.snip(outfile, ".bam") + ".temp.bam"   
    metrics=P.snip(outfile, ".bam") + ".metrics.tsv"
    temporary = PARAMS["tmpdir"]
    statement = '''MarkDuplicates I=%(infile)s
                                  O=%(tempfile)s
                                  M=%(metrics)s
                                  TMP_DIR=%(temporary)s > %(outfile)s.log;

                                samtools view 
                                -F 1024
                                -b
                                %(tempfile)s
                                > %(outfile)s;
                  
                                rm -r %(tempfile)s;

                                samtools index %(outfile)s'''
    P.run(statement)


@follows(mkdir("split.dir"))
@transform(dedup_bams,regex(r"deduped.dir/(.+).bam"),r"split.dir/\1.split.bam")
def splitbams(infile,outfile):
    '''use GATK splitNcigar to split reads into exon segements'''
    fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"
    fastamap = PARAMS["mapfasta"]
    drctry= PARAMS["tmpdir"]
    statement = '''java -Xmx10G -Djava.io.tmpdir=%(drctry)s -jar /shared/sudlab1/General/git_repositories/GATK_file/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar 
                   -T SplitNCigarReads 
                   -R %(fastamap)s
                   -I %(infile)s 
                   -o %(outfile)s 
                   -rf ReassignOneMappingQuality 
                   -RMQF 255 
                   -RMQT 60 
                   -U ALLOW_N_CIGAR_READS ''' 
                                      

   
    job_memory = "32G"
    P.run(statement)

#@follows(mkdir("BaseRecalibration.dir"))
#@transform(splitbams,regex(r"split.dir/(.+).split.bam"),r"BaseRecalibration.dir/\1.recal.bam")
#def baserecal(infile,outfile):
#    fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"
#    fastamap = PARAMS["mapfasta"]
#    drctry= PARAMS["tmpdir"]
#    statement = '''java -Xmx10G -Djava.io.tmpdir=%(drctry)s -jar ~/Downloads/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar 
#                   -T BaseRecalibrator
#                   -R %(fastamap)s
#                   -I %(infile)s 
#                   -o %(outfile)s
#                   '''
#  
#    job_memory = "12G"
#    P.run(statement)
             

#fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"

@follows(mkdir("Variantcalls.dir"))
@transform(splitbams,regex(r"split.dir/(.+).split.bam"),r"Variantcalls.dir/\1.vcf.gz")
def variantcalling(infile,outfile):
    fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"
    fastamap = PARAMS["mapfasta"]
    drctry= PARAMS["tmpdir"]
    tempfile=P.snip(outfile,".gz")
    statement = '''java -Xmx10G -Djava.io.tmpdir=%(drctry)s -jar /shared/sudlab1/General/git_repositories/GATK_file/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar 
                   -T HaplotypeCaller
                   -R %(fastamap)s 
                   -I %(infile)s
                   -dontUseSoftClippedBases 
                   -stand_call_conf 10.0
                   --dbsnp /shared/sudlab1/General/projects/Sumeet/dbSNP/All_20180418_chr.vcf.gz
                   -o %(tempfile)s &&
                   bgzip -c %(tempfile)s > %(outfile)s &&
                   tabix -p vcf %(outfile)s
                   '''  

    job_memory = "16G"
    P.run(statement)         

###################################

@follows(variantcalling,mkdir("phased.dir"))
@transform("Variantcalls.dir/*.vcf.gz", regex(r"Variantcalls.dir/(.+).vcf.gz"),
            add_inputs(r"split.dir/\1.split.bam"),
            r"phased.dir/\1.vcf.gz")
def phasevariants(infiles, outfile):
    pass
    vcf, bam = infiles
    fastamap = PARAMS["mapfasta"]
    drctry= PARAMS["tmpdir"]
    tempfile=P.snip(outfile,".gz")
    statement = ''' java -Xmx10g -Djava.io.tmpdir=%(drctry)s -jar /shared/sudlab1/General/git_repositories/GATK_file/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
                          -T ReadBackedPhasing
                          -R %(fastamap)s
                          -I %(bam)s
                          --variant %(vcf)s
                          -L %(vcf)s
                          -o %(outfile)s
                          --phaseQualityThresh 20.0
                           '''

    job_memory = "16G"
    P.run(statement)

#######################################################
@follows(phasevariants,mkdir("Genotype_vcf.dir"))
@transform("phased.dir/*.vcf.gz",regex(r"phased.dir/(.+).vcf.gz"),
           add_inputs(r"split.dir/\1.split.bam"),r"Genotype_vcf.dir/\1.vcf.gz")


def gvcf(infiles,outfile):
    vcf,bam = infiles
    geno = P.snip(outfile,".gz")
    fastamap = PARAMS["mapfasta"]
    statement = '''freebayes
                   -f %(fastamap)s
                   -@ %(vcf)s %(bam)s > %(geno)s &&
                    bgzip -c %(geno)s > %(outfile)s &&
                    tabix -p vcf %(outfile)s '''
    job_memory = "16G"
    P.run(statement)


#-------------------------------------------------
@follows (mkdir("dbsnpid.dir"))
@transform(gvcf,regex(r"Genotype_vcf.dir/(.+).vcf.gz"),r"dbsnpid.dir/\1.vcf.gz")

def filters(infile,outfile):    
    vcf = P.snip(outfile,".gz")
    statement = ''' bcftools 
                    annotate
                    -a /shared/sudlab1/General/projects/Sumeet/dbSNP/All_20180418_chr.vcf.gz
                    -c ID %(infile)s > %(vcf)s &&
                    bgzip -c %(vcf)s > %(outfile)s &&
                    tabix -p vcf %(outfile)s '''
    P.run(statement)
#---------------------------------------------------

#------------------------------------------------

@follows (filters,mkdir("BiallicSNPs.dir"))
@transform(filters,regex(r"dbsnpid.dir/(.+).vcf.gz"),r"BiallicSNPs.dir/\1.vcf.gz")

def readquality (infile,outfile):
    vcf  = P.snip (outfile,".gz")
    statement = '''bcftools view 
                   -m2 -M2 -v snps %(infile)s > %(vcf)s && 
                    bgzip -c %(vcf)s > %(outfile)s &&
                    tabix -p vcf %(outfile)s '''
    P.run(statement)


#######################################

@follows (readquality,mkdir ("merge.dir"))
@merge ("BiallicSNPs.dir/*.gz","merge.dir/merged.vcf")

def merge (infiles,outfile):
    infile = infiles
    inputs = "BiallicSNPs.dir/*.gz"
    statement = '''bcftools merge %(inputs)s > %(outfile)s'''
    P.run (statement)
 
#-----------------------------------------------------
@follows (merge,mkdir("FilterGenotypecalls.dir"))
@transform("merge.dir/*.vcf",regex(r"merge.dir/merged.vcf"),r"FilterGenotypecalls.dir/FilterGenotypecalls.recode.vcf")

def msinggenotype (infile,outfile):
    vcf = P.snip(outfile,".recode.vcf")
    statement =  ''' vcftools --vcf %(infile)s --max-missing 0.5 --recode --recode-INFO-all --out %(vcf)s'''
    P.run (statement)



#-----------------------------------------
#@follows(merge,mkdir("correcting_genotype.dir"))
#@transform("merge.dir/*", regex(r"merge.dir/merged.vcf"),r"correcting_genotype.dir/corrected_genotype.vcf")

#def correctedgenotype(infile,outfile):
#    job_memory="6G"
#    statement= '''bcftools +missing2ref %(infile)s  > %(outfile)s'''

#    P.run(statement)
   
#------------------------------------------

@follows(msinggenotype,mkdir("FilterMAF.dir"))
@transform("FilterGenotypecalls.dir/*", regex(r"FilterGenotypecalls.dir/FilterGenotypecalls.recode.vcf"),r"FilterMAF.dir/MAF_filtered.recode.vcf")

def MAF(infile,outfile):
    vcf = P.snip (outfile,".recode.vcf")
    statement= ''' vcftools --vcf %(infile)s --maf 0.01 --recode --recode-INFO-all --out %(vcf)s'''
    P.run(statement)

#-----------------------------------------

@follows(MAF,mkdir("FilterHWE.dir"))
@transform("FilterMAF.dir/*", regex(r"FilterMAF.dir/MAF_filtered.recode.vcf"),r"FilterHWE.dir/HWE_filter.recode.vcf")

def HWE (infile,outfile):
    output = P.snip(outfile, ".recode.vcf")
    statement= '''vcftools --vcf %(infile)s --hwe 0.0001 --recode --out %(output)s'''
    P.run(statement)

#---------------------------------------
@follows (HWE, mkdir("Genotype.dir"))
@transform("FilterHWE.dir/*",regex(r"FilterHWE.dir/HWE_filter.recode.vcf"),["Genotype.dir/sample_names.txt", "Genotype.dir/genotype_info.txt"])

def genotype (infile,outfiles):
    names,genotype = outfiles
    statement = ''' bcftools query -l %(infile)s > %(names)s &&
                    bcftools query -f '%%CHROM \\t%%POS \\t%%ID [\\t%%GT]\\n' %(infile)s > %(genotype)s '''

    P.run(statement)
 

# ---------------------------------------------------
@follows (genotype,mkdir("Eqtl_input.dir"))
@transform ("Genotype.dir/*", regex("Genotype.dir/sample_names.txt"),
            add_inputs(r"Genotype.dir/genotype_info.txt"),
            ["Eqtl_input.dir/SNP_info.txt", "Eqtl_input.dir/genotype_info.txt"])

def eqtl_input (infiles,outfiles):
    job_memory = "12G"
    job_threads = 2
    names,geno = infiles
    SNP,genotype = outfiles 
    path_to_script = os.path.dirname(__file__)
    statement = '''Rscript %(path_to_script)s/Input_eqtl.R %(names)s %(geno)s %(SNP)s %(genotype)s 
                     ''' 
    P.run (statement)

#--------------------------------------------------
@follows (eqtl_input, mkdir ("eQTL_results.dir"))
@transform("Eqtl_input.dir/*", regex(r"Eqtl_input.dir/genotype_info.txt"),
            add_inputs(r"Genome_cordinates.dir/*.tsv","Eqtl_input.dir/SNP_info.txt"),
            ["eQTL_results.dir/trans_eqtl.txt","eQTL_results.dir/cis_eqtl.txt"])

def eqtl_results(infiles,outfiles):
    job_memory = "25G"
    job_threads = 3 
    genoptye,gene_location,SNP = infiles
    trans,cis = outfiles
    genotpe = "Eqtl_input.dir/genotype_info.txt"
    expression = PARAMS["Expression_dataset"]
    cov = PARAMS["Covrts"]
    expression_data_path = PARAMS["Expression_dataset"]
    path_to_script = os.path.dirname(__file__)
    statement = '''Rscript %(path_to_script)s/eQTL.R %(genoptye)s %(SNP)s %(expression)s %(expression_data_path)s %(gene_location)s %(cov)s %(trans)s %(cis)s '''
    P.run (statement)

#----------------------------------------------


# Generic pipeline tasks
@follows(add_read_groups,dedup_bams,splitbams,variantcalling,phasevariants,gvcf,filters,readquality,msinggenotype,MAF,HWE,genotype,eqtl_input,eqtl_results)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
