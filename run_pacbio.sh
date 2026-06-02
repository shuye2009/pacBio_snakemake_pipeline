#!/bin/bash
#SBATCH --job-name=pacbio_pipeline
#SBATCH --account=hakemgroup
#SBATCH --partition=all
#SBATCH --time=2-1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Pipeline directory
PIPELINE_DIR="/cluster/home/t128737uhn/snakemake_pipelines/pacBio_snakemake_pipeline"

# Target region type for methbat (e.g., enhancer, cgi, centromere, repeat)
# Override via: sbatch run_pacbio.sh --config target=cgi
TARGET="${TARGET:-enhancer}"

# Change to pipeline directory
cd "$PIPELINE_DIR"

# Create slurm output directory if it doesn't exist
mkdir -p slurm_out

# Load conda
source $(conda info --base)/etc/profile.d/conda.sh

# Activate snakemake environment (adjust if needed)
module load snakemake/7.3.8

# Clear any existing locks
snakemake \
    --unlock \
    --snakefile "$PIPELINE_DIR/Snakefile" 

# Run snakemake with SLURM cluster execution
snakemake \
    --configfile "$PIPELINE_DIR/config/config.yaml" \
    --config target="$TARGET" \
    --cluster-config "$PIPELINE_DIR/config/cluster.json" \
    --snakefile "$PIPELINE_DIR/Snakefile" \
    --cluster "sbatch \
        --account={cluster.account} \
        --partition={cluster.partition} \
        --job-name={cluster.job-name} \
        --time={cluster.time} \
        --nodes={cluster.nodes} \
        --mem={cluster.memory} \
        --ntasks-per-node={cluster.ntasks-per-node} \
        \$(if [ -n \"{cluster.gres}\" ]; then echo \"--gres={cluster.gres}\"; fi) \
        --chdir={cluster.chdir} \
        --output={cluster.output} \
        --error={cluster.error}" \
    --jobs 80 \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    "$@"

module unload snakemake/7.3.8
