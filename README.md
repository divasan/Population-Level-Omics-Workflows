## ATAC-seq

ATAC-seq data is processed using the [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) with Caper and SLURM.

### Files

- `ATAC-seq/samplesheet.csv` – Sample names, biological replicates, and paired-end FASTQ paths.
- `ATAC-seq/generate_json.py` – Generates the ENCODE input JSON files from the samplesheet.
- `ATAC-seq/run_atac.sh` – Submits the samples to the ENCODE pipeline using Caper and SLURM.

### Requirements

- ENCODE ATAC-seq pipeline
- Caper, Java, Python 3, and Apptainer
- SLURM scheduler
- Genome reference files and the corresponding `hg38.tsv` for hg38 analysis

### Running the pipeline

**1.** Clone the ENCODE ATAC-seq pipeline:

```bash
git clone https://github.com/ENCODE-DCC/atac-seq-pipeline.git
```

**2.** Set up Caper with the SLURM backend and prepare the genome reference files following the [ENCODE documentation](https://github.com/ENCODE-DCC/atac-seq-pipeline).

**3.** Edit `ATAC-seq/samplesheet.csv` with the sample names, replicate numbers, and full paths to the FASTQ files. Samples with multiple biological replicates should have separate rows with the same sample name.

**4.** Update the paths at the beginning of `run_atac.sh`:

- `PIPELINE_DIR` – Directory containing the ENCODE `atac.wdl` file.
- `GENOME_TSV` – Full path to the genome reference TSV file (`hg38.tsv`).
- `BASE_DIR` – Directory where the results will be saved.

**5.** Submit the pipeline from the `ATAC-seq` directory on the HPC:

```bash
sbatch run_atac.sh
```

The script generates one input JSON per sample, groups biological replicates within the same experiment, and submits the workflows through Caper. Results and logs are saved in the specified output directory.

### Downstream analysis

Differential accessibility is analyzed using [PyDESeq2](https://github.com/owkin/PyDESeq2), followed by motif enrichment with [HOMER](http://homer.ucsd.edu/homer/) and GO enrichment with [rGREAT](https://bioconductor.org/packages/rGREAT/).

### Downstream files

| File | Description |
|---|---|
| [`differential_accessibility_pyDESeq2.ipynb`](ATAC-seq/downstream/differential_accessibility_pyDESeq2.ipynb) | Merges IDR peaks, counts fragments from BAM files, performs PCA and differential accessibility analysis using PyDESeq2, and identifies DARs using adjusted P-value < 0.05 and \|log2FC\| ≥ 0.5. |
| [`run_HOMER.sh`](ATAC-seq/downstream/run_HOMER.sh) | Runs known transcription factor motif enrichment and genomic peak annotation for regions with increased and decreased accessibility. |
| [`run_GREAT_GO_BP.R`](ATAC-seq/downstream/GREAT/run_GREAT_GO_BP.R) | Runs GO Biological Process enrichment using rGREAT and identifies region-to-gene associations for DARs in both directions. |
| [`submit_GREAT.sh`](ATAC-seq/downstream/GREAT/submit_GREAT.sh) | Example SLURM script for submitting the GREAT analysis on an HPC cluster. |

### Running downstream analyses

**1.** Update the sample information, input paths, and conditions in the PyDESeq2 notebook.

**2.** Prepare the DAR BED files, update the paths in `run_HOMER.sh`, and submit:

```bash
sbatch run_HOMER.sh
```

**3.** Update the paths in `run_GREAT_GO_BP.R` and the environment settings in `submit_GREAT.sh`, then submit:

```bash
sbatch submit_GREAT.sh
```

## CUT&Tag / CUT&RUN

CUT&Tag and CUT&RUN data are processed using [nf-core/cutandrun v3.2.2](https://nf-co.re/cutandrun/3.2.2/) with Nextflow and Apptainer on a SLURM cluster.

### Files

- `CUTandRUN-CUTandTAG/samplesheet.csv` – Contains sample names, biological replicates, paired-end FASTQ paths, and IgG controls.
- `CUTandRUN-CUTandTAG/run_cutandrun.slurm` – SLURM script for running nf-core/cutandrun.

### Requirements

- nf-core/cutandrun v3.2.2
- Nextflow 24.10.0
- Java 17 and Apptainer
- SLURM scheduler
- hg38 reference FASTA (`hg38.fa`)
- Bowtie2 index files for hg38
- GENCODE gene annotation (GTF)
- ENCODE hg38 blacklist (`hg38-blacklist.v2.bed.gz`)

### Running the pipeline

**1.** Install Nextflow and nf-core/cutandrun v3.2.2 following the [nf-core documentation](https://nf-co.re/cutandrun/3.2.2/).

**2.** Edit `samplesheet.csv` with the sample names, replicate numbers, paired-end FASTQ paths, and corresponding IgG controls.

**3.** Update the paths at the beginning of `run_cutandrun.slurm`:

- `PIPELINE` – Directory containing `main.nf` from nf-core/cutandrun v3.2.2.
- `NEXTFLOW` – Full path to the Nextflow executable.
- `BASE_DIR` – Directory for results and intermediate files.
- `FASTA` – Path to the hg38 reference FASTA.
- `BOWTIE2_DIR` – Directory containing the hg38 Bowtie2 index files.
- `GTF` – Path to the GENCODE gene annotation.
- `BLACKLIST` – Path to the hg38 blacklist file.

**4.** Adjust the SLURM resources and partition for your HPC, then submit from the `CUTandRUN-CUTandTAG` directory:

```bash
sbatch run_cutandrun.slurm
```

The pipeline uses Bowtie2 for alignment, SEACR and MACS2 for peak calling, and CPM normalization. IgG controls are specified in the samplesheet.

The analysis runs within a single SLURM allocation using Nextflow's local executor. Results are saved in the specified output directory.


## RNA-seq

RNA-seq data are processed using [nf-core/rnaseq v3.18.0](https://nf-co.re/rnaseq/3.18.0/) with Nextflow and Apptainer on a SLURM cluster.

### Files

- `RNA-seq/samplesheet.csv` – Contains sample names, paired-end FASTQ paths, and strandedness information.
- `RNA-seq/run_rnaseq.slurm` – SLURM script for running nf-core/rnaseq.

### Requirements

- nf-core/rnaseq v3.18.0
- Nextflow 24.10.0
- Java 17 and Apptainer
- SLURM scheduler
- hg38 reference FASTA (`hg38.fa`)
- GENCODE v50 gene annotation (`gencode.v50.annotation.gtf`)

### Running the pipeline

**1.** Install Nextflow and nf-core/rnaseq v3.18.0 following the [nf-core documentation](https://nf-co.re/rnaseq/3.18.0/).

**2.** Edit `samplesheet.csv` with the sample names, paired-end FASTQ paths, and strandedness. Each biological replicate should be listed as a separate sample. Use `auto` for automatic strandedness inference.

**3.** Update the paths at the beginning of `run_rnaseq.slurm`:

- `PIPELINE` – Directory containing `main.nf` from nf-core/rnaseq v3.18.0.
- `NEXTFLOW` – Full path to the Nextflow executable.
- `BASE_DIR` – Directory for results and intermediate files.
- `FASTA` – Path to the hg38 reference FASTA.
- `GTF` – Path to the GENCODE v50 annotation file.

**4.** Adjust the SLURM resources and partition for your HPC, then submit from the `RNA-seq` directory:

```bash
sbatch run_rnaseq.slurm
```

The pipeline uses STAR and RSEM for alignment and quantification. Alignment intermediates and generated reference indices are saved.

The analysis runs within a single SLURM allocation using Nextflow's local executor. Results are saved in the specified output directory.
