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
    # 2125-06-06-1 = MH MP
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


    # make pairs and send to cutadapt for trimming external adaptors
    trim_cutadapt = main_pipeline.collate(
        name='trim_cutadapt',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/cutadapt_pe',
            job_name='cutadapt'),
        input=merged_fq_files,
        filter=ruffus.formatter(
            r'.+/(?P<LIB>[^_]+)_R(?P<RN>\d)_merged.fastq.gz'),
        output=[['output/cutadapt/{LIB[0]}_R1_trimmed.fastq.gz',
                'output/cutadapt/{LIB[0]}_R2_trimmed.fastq.gz']])

    # send trimmed reads to splitnextera
    mp_splitnextera = main_pipeline.subdivide(
        name='splitnextera',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/splitnextera',
            job_name='splitnextera'),
        input=trim_cutadapt,
        filter=ruffus.regex(
            r'.+?/2125-06-06-1_R(?P<RN>\d)_trimmed.fastq.gz'),
        output=['output/splitnextera/2125-06-06-1.pe.fastq.gz',
                'output/splitnextera/2125-06-06-1.se.fastq.gz',
                'output/splitnextera/2125-06-06-1.mp.fastq.gz',
                'output/splitnextera/2125-06-06-1.unknown.fastq.gz'])

    # decontaminate PhiX (other?) sequences
    decon_mp = main_pipeline.transform(
        name='decon_mp',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/decon',
            job_name='decon_mp'),
        input=mp_splitnextera,
        filter=ruffus.formatter(
            r'.+/2125-06-06-1\.(?P<VL>[^.]+)\.fastq.gz'),
        output=['output/decon/2125-06-06-1_{VL[0]}.fastq.gz'])

    decon_pe = main_pipeline.transform(
        name='decon_pe',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/decon',
            job_name='decon_pe'),
        input=trim_cutadapt,
        filter=ruffus.regex(
            r'.+?/2125-06-11-1_R(?P<RN>\d)_trimmed.fastq.gz'),
        output=[r'output/decon/2125-06-11-1.fastq.gz'])

    decon = [decon_mp, decon_pe]

    # digital normalisation and error correction w/ bbnorm
    bbnorm = main_pipeline.subdivide(
        name='bbnorm',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/bbnorm',
            job_name='bbnorm',
            mem_per_cpu=7000,
            cpus_per_task=8),
        input=decon,
        filter=ruffus.formatter(r'.+/(?P<LN>[^(_|.)]+)(?P<VL>_?\w*).fastq.gz'),
        output=[r'output/bbnorm/{LN[0]}{VL[0]}.fastq.gz'])

    # run fastqc on decontaminated and normalised libraries
    main_pipeline.transform(
        name='fastqc',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/fastqc',
            job_name='fastqc',
            cpus_per_task=1),
        input=bbnorm,
        filter=ruffus.formatter(r'.+/(?P<LN>[^(_|.)]+)(?P<VL>_?\w*).fastq.gz'),
        output=[r'output/fastqc/{LN[0]}{VL[0]}_fastqc.html'])

    # select files for hashing
    # velveth_input_files = [x.path for x in os.scandir('output/khmer')
    #                        if (x.name.endswith('fastq.gz')
    #                        or x.name.endswith('.fastq'))
    #                        and x.is_file()
    #                        and (('proper' in x.name and 'se' not in x.name)
    #                             or ('se_orphans' in x.name))]

    # # # prepare files with velveth
    # # # set threads for velvet to 1 !!!
    # min_kmer = 71
    # max_kmer = 87
    # step = 8
    # kmer_lengths = [x for x in range(min_kmer, max_kmer + 1, step)]
    # velveth_output = list(
    #     tompytools.flatten_list(
    #         [('output/velveth_' + str(x) + '/Sequences')
    #          for x in kmer_lengths]))
    # velveth = main_pipeline.merge(
    #     name='velveth',
    #     task_func=test_job_function,
    #     input=velveth_input_files,
    #     output=velveth_output)\
    #     .follows(diginorm)

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
