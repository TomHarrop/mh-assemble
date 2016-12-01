#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

################################
# ASW genome assembly pipeline #
################################

################
# Requirements #
################

import tompltools
import tompytools
import ruffus
import os


#############
# Functions #
#############

# sanity check list of pairs
def pair_sanity_check(pairs_list):
    tompytools.generate_message('Sanity checking pair:')
    print("R1 file:\t%s\nR2 file:\t%s\n\n" % (pairs_list[0], pairs_list[1]))
    if not os.path.isfile(pairs_list[0]):
        raise FileNotFoundError('R1 file\n\t' + pairs_list[0] +
                                '\nnot found')
    if not os.path.isfile(pairs_list[1]):
        raise FileNotFoundError('R2 file\n\t' + pairs_list[1] +
                                '\nnot found')


# create pairs from sanity-checked dictionary
def create_list_of_pairs(pairs_dict):
    return(list([x, pairs_dict[x]] for x in pairs_dict))


############
# Pipeline #
############

def main():

    #########
    # SETUP #
    #########

    # test function for checking input/output passed to job_script and parsing
    # by src/sh/io_parser
    test_job_function = tompltools.generate_job_function(
        job_script='src/sh/io_parser',
        job_name='test')

    # parse email etc. here?
    parser = ruffus.cmdline.get_argparse(
        description='ASW genome assembly pipeline.')
    # parser.add_argument('--email', '-e',
    #                     help='Logon email address for JGI',
    #                     type=str,
    #                     dest='jgi_logon')
    # parser.add_argument('--password', '-p',
    #                     help='JGI password',
    #                     type=str,
    #                     dest='jgi_password')
    options = parser.parse_args()
    # jgi_logon = options.jgi_logon
    # jgi_password = options.jgi_password

    # initialise pipeline
    main_pipeline = ruffus.Pipeline.pipelines['main']

    # find fastq.gz files
    dir_listing = [x[0] for x in os.walk(top='data', followlinks=True)]
    fastq_file_list = []
    for directory in dir_listing:
        file_list = os.scandir(directory)
        fastq_file_list.append([x.path for x in file_list
                               if (x.name.endswith('fastq.gz')
                                   or x.name.endswith('.fastq'))
                               and x.is_file()])

    fastq_files = list(tompytools.flatten_list(fastq_file_list))

    # extract only MH gDNA fastq data, i.e.
    # 2125-06-11-1 = MH PE
    # 2125-06-06-1 = ASW MP
    active_fq_files = [x for x in fastq_files
                       if ('2125-06-11-1' in x
                           or '2125-06-06-1' in x)]

    # load files into ruffus 
    raw_fq_files = main_pipeline.originate(
        name='raw_fq_files',
        task_func=os.path.isfile,
        output=active_fq_files)

    # merge libraries
    merged_fq_files = main_pipeline.collate(
        name='merged_fq_files',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/merge_fq',
            job_name='merge_fq'),
        input=raw_fq_files,
        filter=ruffus.formatter(
            r'data/NZGL02125/.*/'
            '[^-]+-(?P<LIB>[^_]+).+_R(?P<RN>\d)_.*.fastq.gz'),
        output=[r'output/fq_merged/{LIB[0]}_R{RN[0]}_merged.fastq.gz'])

    # make pairs and send to cutadapt
    pe_trimmed = main_pipeline.collate(
        name='pe_trimmed',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/cutadapt_pe',
            job_name='pe_trimmed'),
        input=merged_fq_files,
        filter=ruffus.regex(
            r'output/fq_merged/2125-06-11-1_R(\d)_merged.fastq.gz'),
        output=['output/cutadapt/pe/2125-06-11-1_R1_trimmed.fastq.gz',
                'output/cutadapt/pe/2125-06-11-1_R2_trimmed.fastq.gz'])

    mp_trimmed = main_pipeline.collate(
        name='mp_trimmed',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/cutadapt_mp',
            job_name='mp_trimmed'),
        input=merged_fq_files,
        filter=ruffus.regex(
            r'output/fq_merged/2125-06-06-1_R(\d)_merged.fastq.gz'),
        output=['output/cutadapt/mp/2125-06-06-1_R1_trimmed.fastq.gz',
                'output/cutadapt/mp/2125-06-06-1_R2_trimmed.fastq.gz'])

    # decontaminate PhiX (other?) sequences
    decon = main_pipeline.transform(
        name='decon',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/decon',
            job_name='decon',
            ntasks=1,
            cpus_per_task=4),
        input=[pe_trimmed, mp_trimmed],
        filter=ruffus.formatter(
            r'output/cutadapt/(?P<LT>\w{2})/'
             '(?P<LIB>[^_]+)_R(?P<RN>\d)_trimmed.fastq.gz',
            r'output/cutadapt/(?P<LT>\w{2})/'
             '(?P<LIB>[^_]+)_R(?P<RN>\d)_trimmed.fastq.gz'),
        output=['output/bbduk/{LT[0]}/'
                '{LIB[0]}_R{RN[0]}_decon.fastq.gz',
                'output/bbduk/{LT[1]}/'
                '{LIB[1]}_R{RN[1]}_decon.fastq.gz'])

    # run fastqc on trimmed libraries
    fastqc = main_pipeline.subdivide(
        name='fastqc',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/fastqc',
            job_name='fastqc',
            ntasks=1,
            cpus_per_task=2),
        input=decon,
        filter=ruffus.formatter(
            r'output/bbduk/\w{2}/'
             '(?P<LIB>[^_]+)_R(?P<RN>\d)_decon.fastq.gz',
            r'output/bbduk/\w{2}/'
             '(?P<LIB>[^_]+)_R(?P<RN>\d)_decon.fastq.gz'),
        output=['output/fastqc/{LIB[0]}_R{RN[0]}_decon_fastqc.html',
                'output/fastqc/{LIB[1]}_R{RN[1]}_decon_fastqc.html'])

    # reverse complement mp reads
    mp_revcomp = main_pipeline.transform(
        name='mp_revcomp',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/mp_revcomp',
            job_name='mp_revcomp',
            ntasks=2,
            cpus_per_task=1),
        input=decon,
        filter=ruffus.formatter(
            r'output/bbduk/mp/'
             '(?P<LIB>[^_]+)_R(?P<RN>\d)_decon.fastq.gz',
            r'output/bbduk/mp/'
             '(?P<LIB>[^_]+)_R(?P<RN>\d)_decon.fastq.gz'),
        output=['output/revcomp/{LIB[0]}_R{RN[0]}_rc.fastq.gz',
                'output/revcomp/{LIB[1]}_R{RN[1]}_rc.fastq.gz'])

    # prepare files with velveth
    hash_files = main_pipeline.collate(
        name='hash_files',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/hash_files',
            job_name='hash_files',
            ntasks=1,
            cpus_per_task=8),
        input=decon,
        add_inputs=ruffus.add_inputs(ruffus.output_from(mp_revcomp)),
        filter=ruffus.regex(
            r'output/bbduk/pe/2125-06-11-1_R(\d)_decon.fastq.gz'),
        output=['output/velveth/Sequences'])



    ###################
    # RUFFUS COMMANDS #
    ###################

    # print the flowchart
    ruffus.pipeline_printout_graph(
        "ruffus/flowchart.pdf", "pdf",
        pipeline_name="ASW genome assembly pipeline")

    # run the pipeline
    ruffus.cmdline.run(options, multithread=8)

if __name__ == "__main__":
    main()
