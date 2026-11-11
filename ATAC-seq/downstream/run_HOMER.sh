#!/usr/bin/env bash
#SBATCH --job-name=HOMER_DAR
#SBATCH --output=HOMER_%j.out
#SBATCH --error=HOMER_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=standard

set -euo pipefail

# Load HOMER on the cluster.
module load homer/5.1

# Set the paths to the two DAR BED files.
# These should contain regions with lower or higher accessibility in condition B.
LOWER_BED="/path/to/lower_in_condition_B.bed"
HIGHER_BED="/path/to/higher_in_condition_B.bed"
OUT_DIR="/path/to/homer_results"
GENOME="mm10"

mkdir -p "${OUT_DIR}"

# Check that the input files are available.
for bed in "${LOWER_BED}" "${HIGHER_BED}"; do
    if [[ ! -s "${bed}" ]]; then
        echo "Missing or empty BED file: ${bed}" >&2
        exit 1
    fi
done

# Use the full DAR intervals for known-motif enrichment.
# -nomotif skips de novo motif discovery; HOMER chooses the background.
for group in lower_in_condition_B higher_in_condition_B; do
    if [[ "${group}" == "lower_in_condition_B" ]]; then
        bed="${LOWER_BED}"
    else
        bed="${HIGHER_BED}"
    fi

    echo "Running HOMER for ${group}"
    motif_dir="${OUT_DIR}/${group}/known_motifs"
    mkdir -p "${motif_dir}"

    findMotifsGenome.pl "${bed}" "${GENOME}" "${motif_dir}" \
        -size given \
        -mask \
        -nomotif \
        -p "${SLURM_CPUS_PER_TASK:-4}"

    # Annotate the same DARs with nearby genes and genomic features.
    annotatePeaks.pl "${bed}" "${GENOME}" \
        > "${OUT_DIR}/${group}_annotation.tsv"

done

echo "HOMER analysis complete."
