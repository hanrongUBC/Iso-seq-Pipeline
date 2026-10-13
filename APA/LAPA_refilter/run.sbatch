#!/bin/bash
#SBATCH --account=def-kladams
#SBATCH --job-name=LAPA_refilter
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=LAPA_refilter_%j.out
#SBATCH --error=LAPA_refilter_%j.err
# Re-filter the finished LAPA run (no LAPA re-run) and test differential
# poly(A) site usage with all replicates. Submit from this directory:
#   sbatch run.sbatch
set -eo pipefail
HERE="$SLURM_SUBMIT_DIR"
APA_ROOT=/home/hanrong/scratch/hanrong_LR_stress/APA
LAPA_RUN="$APA_ROOT/LAPA_Westar/runs/lapa_end_flnc_job63782941_20261009_023104_n3qqA0"
# Must be the FASTA LAPA itself used (see $LAPA_RUN/command.sh), with its .fai.
FASTA="$APA_ROOT/Westar_Gennome/Brassica_napus.Westar.v0.genome.fa"
SAMPLES="$APA_ROOT/APA_software_server_20261009/samples.tsv"
OUT="$HERE/runs/job${SLURM_JOB_ID}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /home/hanrong/scratch/conda_envs/movapa_westar   # has pyBigWig, DEXSeq, data.table
mkdir -p "$OUT"
cp "$HERE/lapa_refilter.py" "$HERE/run_dexseq_apa.R" "$OUT/"

python "$HERE/lapa_refilter.py" --lapa-results "$LAPA_RUN/results" --fasta "$FASTA" \
    --samples "$SAMPLES" --out "$OUT/refilter" \
    --min-count 3 --min-reps 2 --a10 7 --a20 12 --a-run 6
Rscript --vanilla "$HERE/run_dexseq_apa.R" "$OUT/refilter/count_matrix.tsv" "$SAMPLES" \
    "$OUT/DEXSeq" "${SLURM_CPUS_PER_TASK:-1}" Control
date > "$OUT/COMPLETED.txt"
