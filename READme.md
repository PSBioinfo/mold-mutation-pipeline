
# Mold Mutation Detection Pipeline

**Author:** PSBioinfo   
**Objective:** Build a Python pipeline to identify mutations responsible for different mold colors in ear samples of 50 patients.

---

## Overview

This pipeline performs the following steps:

1. **Demultiplexes** pooled FASTQ files by barcodes from clinical metadata.
2. **Trims** low-quality reads.
3. **Aligns** reads to a reference genome using BWA.
4. **Converts & sorts** alignment files using SAMtools.
5. **Detects mutations** using pileup and basic variant analysis.
6. **Generates a report** mapping mutations to mold color and read counts per sample.

---

##  Project Structure

mold-mutation-pipeline/
├── data/ # Input files (mock FASTQ and metadata)
├── scripts/ # Python scripts (Pipeline and parser)
├── output/ # Output report file
├── fastqs/ # Auto-generated trimmed FASTQ files
├── sam_files/ # Auto-generated SAM files
├── bam_files/ # Auto-generated sorted BAM files
├── requirements.txt # Python dependencies (pysam)
├── README.md # This project description
└── LICENSE # MIT License

## Sample Data


Mock data files are provided for demonstration purposes:

- `data/sample.fastq`: Two short FASTQ reads with barcodes
- `data/sample_clinical_data.txt`: Clinical data file with sample names, mold colors, and barcodes

These files simulate the input format expected by the pipeline and allow the code to be tested without real patient data.

## How to Run

From the root directory, execute:
python3 scripts/Pipeline.py \
  -f data/sample.fastq \
  -c data/sample_clinical_data.txt \
  -o fastqs

This will create intermediate files and generate a file report.txt in the output/ directory

## Example Output

Sample Patient01 had a green mold, 1452 reads, and had 48.62% of reads at position 125 with the mutation T.
Sample Patient02 had a blue mold, 1380 reads, and had 36.44% of reads at position 147 with the mutation A.

## Notes

- The reference genome (dgorgon_reference.fa) should be placed in the project directory and indexed with bwa index before alignment.
- Intermediate files (fastqs/, sam_files/, bam_files/) are auto-created during runtime.
- The script expects barcode matching to be exact.

## Dependencies

Install Python dependencies:

```bash
pip install -r requirements.txt


