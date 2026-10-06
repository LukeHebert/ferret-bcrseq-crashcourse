#!/usr/bin/env python3
"""
Takes one forward and one reverse read file from Illumina platform paired-end
sequencing of B cell receptors and:
1) generates fastQC .html files to manually evaluate quality metrics of sequencing data
2) trims library prep adapter sequences & trims poor quality 3' ends of reads
3) merges the paired forward & reverse reads

The primary output file "basename.assembled.fastq" will be written to the same 
directory as the input R1 & R2 files, and other outputs & a log file will be 
saved to a new 'trim_merge' (quality control) subdirectory.

Use `python trim_merge.py --help` for instructions
"""

import datetime
import argparse
import subprocess
import os
import shutil

def parse_arguments():
    """ Parse command-line arguments """
    parser = argparse.ArgumentParser(
        description='Generates .html files for quality assessment, trims adapters and low quality 3 prime read ends, and merges paired forward & reverse reads.')
    parser.add_argument('R1', metavar='forward_reads', type=str, 
        help='R1 fastq.gz file')
    parser.add_argument('R2', metavar='reverse_reads', type=str, 
        help='R2 fastq.gz file')
    parser.add_argument('--notrim', action='store_true', 
                        help='If set, skips the adapter & quality score trimming step before merging paired reads.')
    parser.add_argument('--threads', type=int, default=min(os.cpu_count() or 1, 2),
                        help='Worker threads for cutadapt and VSEARCH (default: 2 or available CPUs).')
    parser.add_argument('--fastqc', default='fastqc', help='FastQC executable or PATH name.')
    parser.add_argument('--cutadapt', default='cutadapt', help='cutadapt executable or PATH name.')
    parser.add_argument('--vsearch', default='vsearch', help='VSEARCH executable or PATH name.')
    return parser.parse_args()

def resolve_executable(value):
    """Return an explicit executable path, with a useful portable error."""
    found = value if os.path.exists(value) else shutil.which(value)
    if not found:
        raise FileNotFoundError(f"Required executable was not found: {value}")
    return found

def run_logged(command, log_pathway):
    """Run a command and append its combined output to the pipeline log."""
    with open(log_pathway, 'a') as log:
        log.write("$ " + subprocess.list2cmdline(command) + "\\n")
        log.flush()
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)

def run_fastqc(files, output_dir, log_pathway, fastqc):
    """ Runs FastQC on a list of files and outputs the reports to the specified directory """
    command = [fastqc, '-o', output_dir, *files]
    run_logged(command, log_pathway)
    return subprocess.list2cmdline(command)

def call_cutadapt(r1, r2, output_dir, log_pathway, cutadapt, threads):
    """ Trims adapters and low-quality ends from reads using cutadapt """
    # These are the TruSeq standard adapter sequences
    adapter1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
    adapter2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    quality_cutoff = '20'  # Quality score cutoff for 3' ends
    trimmed_r1 = os.path.join(output_dir, os.path.basename(r1).replace('.fastq.gz', '_trimmed.fastq.gz'))
    trimmed_r2 = os.path.join(output_dir, os.path.basename(r2).replace('.fastq.gz', '_trimmed.fastq.gz'))
    # A minimum length of 50 prevents the merger from receiving empty reads. Use
    # pair-filter=any so a pair is removed when either mate falls below that
    # limit; merging requires two non-empty synchronized mate files.
    command = [cutadapt, '-q', quality_cutoff, '-a', adapter1, '-A', adapter2,
               '-o', trimmed_r1, '-p', trimmed_r2, r1, r2, '-m', '50',
               '--pair-filter=any', '-j', str(threads)]
    run_logged(command, log_pathway)
    return trimmed_r1, trimmed_r2, subprocess.list2cmdline(command)

def call_vsearch(r1, r2, output_dir, log_pathway, vsearch, threads):
    """Merge overlapping forward/reverse pairs with VSEARCH."""
    base = os.path.basename(r1).replace('_trimmed.fastq.gz', '').replace('_R1','')
    assembled = os.path.join(output_dir, f"{base}.assembled.fastq")
    command = [vsearch, '--fastq_mergepairs', r1, '--reverse', r2,
               '--fastqout', assembled, '--fastq_minovlen', '10',
               '--fastq_maxmergelen', '700', '--fastq_minmergelen', '50',
               '--threads', str(threads)]
    run_logged(command, log_pathway)
    return subprocess.list2cmdline(command)

def time_passed(start_time):
    """ Makes a human-readable string of elapsed time from start_time to 
    when this function is called"""
    elapsed_time = datetime.datetime.now() - start_time
    elapsed_hr = int(elapsed_time.total_seconds() // 3600)
    elapsed_min = int((elapsed_time.total_seconds() % 3600) // 60)
    elapsed_sec = int(elapsed_time.total_seconds() % 60)
    return f"{elapsed_hr:02}:{elapsed_min:02}:{elapsed_sec:02}"

def main(args):
    """ Main entry point of the module """
    
    if args.threads < 1:
        raise ValueError('--threads must be at least 1')
    fastqc = resolve_executable(args.fastqc)
    vsearch = resolve_executable(args.vsearch)
    cutadapt = None if args.notrim else resolve_executable(args.cutadapt)

    # Create log file with current date and time in its name
    now = datetime.datetime.now()
    input_dir = os.path.dirname(args.R1)
    out_dir = os.path.join(input_dir, "trim_merge")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    log_name = f'log_QC-trimming-merging_{now.strftime("%Y-%m-%d_%H-%M-%S")}.txt'
    log_path = os.path.join(out_dir, log_name)

    # Run FastQC on initial input files
    qc_args = run_fastqc([args.R1, args.R2], out_dir, log_path, fastqc)
    
    if not args.notrim:
        # Trim adapters & poor quality read ends
        trimmed_r1, trimmed_r2, trim_args = call_cutadapt(
            args.R1, args.R2, out_dir, log_path, cutadapt, args.threads)
    else:
        # Or use original files if user specifies to skip trimming
        # In this case, copy the original files to the qc directory to keep generated files together
        trimmed_r1 = os.path.join(out_dir, os.path.basename(args.R1))
        trimmed_r2 = os.path.join(out_dir, os.path.basename(args.R2))
        shutil.copy2(args.R1, trimmed_r1)
        shutil.copy2(args.R2, trimmed_r2)
        trim_args = "Trimming skipped due to --notrim flag"
    
    # Merge overlapping mates into one read for IgBLAST annotation.
    merge_args = call_vsearch(trimmed_r1, trimmed_r2, out_dir, log_path, vsearch, args.threads)

    # Determine base name for the assembled file
    base = os.path.basename(trimmed_r1).replace('_trimmed.fastq.gz', '').replace('_R1','')
    # The assembled file is originally created in the 'trim_merge' directory
    old_assembled_file = os.path.join(out_dir, f"{base}.assembled.fastq")
    # Move the assembled file to the input directory
    assembled_file = os.path.join(input_dir, f"{base}.assembled.fastq")
    shutil.move(old_assembled_file, assembled_file)
    
    # Run FastQC on the output assembled file
    run_fastqc([assembled_file], out_dir, log_path, fastqc)
    
    # Append the tool-calling command and the time it took to log file
    elapsed_str = time_passed(now)
    with open(log_path, 'a') as log:
        log.write(f'\n\nQUALITY CHECK COMMAND:\n{qc_args}'
                  f'\n\nADAPTER & PHRED TRIMMING COMMAND:\n{trim_args}'
                  f'\n\nMERGE COMMAND:\n{merge_args}'
                  f'\n\nTOTAL SCRIPT RUN TIME\n(HR:MIN:SEC):\n{elapsed_str}')

if __name__ == "__main__":
    """ This is executed when run from the command line """
    args = parse_arguments()
    main(args)
