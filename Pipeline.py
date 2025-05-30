#!/usr/bin/env python3
#The goal of this script is to to build a pipeline for Dr. Harrington to analyze the data and identify the mutations that resulted in the different color molds across 50 patients

import os
import sys
import argparse
import shutil
import subprocess
import glob
import pysam

script_dir = os.path.dirname(os.path.realpath(__file__))
necessary_scripts_path = os.path.join(script_dir, 'necessary_scripts')
sys.path.append(necessary_scripts_path)

from parseFastq import ParseFastQ

#Demultiplex the pooled fastq
def trim_read(sequence, quality):
    for i in range(len(quality) - 2, -1, -1):
        if quality[i] in 'DF' and quality[i + 1] in 'DF':
            sequence = sequence[:i]
            quality = quality[:i]
            break
    return sequence, quality

def demultiplex_and_trim(fastq_file, clinical_data_file, output_dir):
    barcodes = {}
    with open(clinical_data_file, 'r') as file:
        for line in file:
            name, color, barcode = line.strip().split('\t')
            barcodes[barcode] = name
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    fastq_parser = ParseFastQ(fastq_file)
    
    # Process each record in the FASTQ file and extract the barcode and trim beginning and ends of reads
    for seqHeader, seqStr, qualHeader, qualStr in fastq_parser:
        barcode = seqStr[:5]
        if barcode in barcodes:
            trimmed_sequence, trimmed_quality = trim_read(seqStr[5:], qualStr[5:])
            output_path = os.path.join(output_dir, f"{barcodes[barcode]}_trimmed.fastq")
            with open(output_path, 'a') as output_file:
                output_file.write(f"{seqHeader}\n{trimmed_sequence}\n{qualHeader}\n{trimmed_quality}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Demultiplex, trim, and write FASTQ records to per-sample files.")
    parser.add_argument("-f", "--fastq", required=True, help="Path to the input pooled FASTQ file")
    parser.add_argument("-c", "--clinical", required=True, help="Path to the clinical data file")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for demultiplexed FASTQ files")
    args = parser.parse_args()
    
    demultiplex_and_trim(args.fastq, args.clinical, args.output_dir)

#Perform alignment on each FASTQ to reference sequence
def align_fastq_to_reference(reference_path, fastq_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    sample_name = os.path.basename(fastq_path).replace('_trimmed.fastq', '')
    sam_file_path = os.path.join(output_dir, f"{sample_name}.sam")
    bwa_command = f"bwa mem {reference_path} {fastq_path} > {sam_file_path}"
    try:
        subprocess.run(bwa_command, shell=True, check=True)
        print(f"Alignment completed for {sample_name}. SAM file created at {sam_file_path}.") #Debug
    except subprocess.CalledProcessError as e:
        print(f"Error during alignment for {sample_name}: {e}") #Debug

reference_genome = 'dgorgon_reference.fa'
index_command = f"bwa index {reference_genome}"
subprocess.run(index_command, shell=True)
fastq_directory = 'fastqs'
fastq_files = glob.glob(f"{fastq_directory}/*_trimmed.fastq")
output_directory = 'sam_files'
for fastq_file in fastq_files:
    align_fastq_to_reference(reference_genome, fastq_file, output_directory)

#Convert the samfiles to bamfiles
def sam_to_sorted_bam(sam_file_path, sam_files_dir, bam_files_dir):
    base_name = os.path.basename(sam_file_path).replace('.sam', '')
    bam_file_path = os.path.join(bam_files_dir, f"{base_name}.bam")
    sorted_bam_file_path = os.path.join(bam_files_dir, f"{base_name}.sorted.bam")
    subprocess.run(f"samtools view -bS {sam_file_path} > {bam_file_path}", shell=True, check=True)
    subprocess.run(f"samtools sort -m 100M -o {sorted_bam_file_path} {bam_file_path}", shell=True, check=True)
    subprocess.run(f"samtools index {sorted_bam_file_path}", shell=True, check=True)
    
#Remove the prior .sam files and the .bam files
    os.remove(sam_file_path)
    os.remove(bam_file_path)
    
sam_files_directory = 'sam_files'
bam_files_directory = 'bam_files'
os.makedirs(bam_files_directory, exist_ok=True)
sam_files = glob.glob(f"{sam_files_directory}/*.sam")
for sam_file in sam_files:
    sam_to_sorted_bam(sam_file, sam_files_directory, bam_files_directory)


#Discover variants in each sorted bam file.
def pileup(bam_file_path):
    samfile = pysam.AlignmentFile(bam_file_path, "rb")
    processed_reads = set()
    for pileupcolumn in samfile.pileup():
        ntdict = {}
        for pileupread in pileupcolumn.pileups:
            query_name = pileupread.alignment.query_name 
            if query_name not in processed_reads:
                processed_reads.add(query_name)
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    ntdict[base] = ntdict.get(base, 0) + 1
        total_reads = sum(ntdict.values())
        for base, count in ntdict.items():
            frequency = (count / total_reads) * 100
            if frequency != 100 and frequency > 0:
                print(f"Base {base} at position {pileupcolumn.pos} has frequency {frequency}%")

    samfile.close()

if __name__ == "__main__":
    bam_files_dir = "bam_files"
    for bam_file in os.listdir(bam_files_dir):
        if bam_file.endswith(".sorted.bam"):
            bam_file_path = os.path.join(bam_files_dir, bam_file)
            print(f"Processing {bam_file}...")
            pileup(bam_file_path)

#Create a report that outputs nucleotide position & mutation reponsible for each color of the mold and print out the number of sequence that were used for each sample
def parse_clinical_data(clinical_data_path):
    sample_colors = {}
    with open(clinical_data_path, 'r') as file:
        for line in file:
            name, color, _ = line.strip().split('\t')
            sample_colors[name] = color
    return sample_colors

def analyze_pileup(bam_file_path):
   
    samfile = pysam.AlignmentFile(bam_file_path, "rb")
    mutations = []
    for pileupcolumn in samfile.pileup():
        ntdict = {}
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                ntdict[base] = ntdict.get(base, 0) + 1
        total_reads = sum(ntdict.values())
        for base, count in ntdict.items():
            percentage = (count / total_reads) * 100
            if percentage > 0 and percentage < 100:
                mutations.append((pileupcolumn.pos, base, percentage))
                break
    samfile.close()
    return mutations

def count_sequences(bam_file_path):
    if bam_file_path.endswith('.sorted.bam'):
        samfile = pysam.AlignmentFile(bam_file_path, "rb")
        sequence_count = samfile.count()
        samfile.close()
        return sequence_count
    else:
        return 0  #if file is not a .sorted.bam file do not consider

def create_report(bam_files_dir, clinical_data_path, output_report_path):
    sample_colors = parse_clinical_data(clinical_data_path)
    with open(output_report_path, 'w') as report_file:
        for sample_name, color in sample_colors.items():
            bam_file_name = f"{sample_name}.sorted.bam"
            bam_file_path = os.path.join(bam_files_dir, bam_file_name)
            if os.path.exists(bam_file_path):
                mutations = analyze_pileup(bam_file_path)
                sequence_count = count_sequences(bam_file_path)
                for pos, base, percentage in mutations:
                    report_line = f"Sample {sample_name} had a {color} mold, {sequence_count} reads, and had {percentage:.2f}% of reads at position {pos} had the mutation {base}.\n"
                    report_file.write(report_line)

clinical_data_path = 'harrington_clinical_data.txt'
bam_files_dir = 'bam_files'
output_report_path = 'report.txt'

create_report(bam_files_dir, clinical_data_path, output_report_path)
