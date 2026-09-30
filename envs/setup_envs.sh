#!/usr/bin/env bash
# Create the four conda environments for pipeline_prostate.py.
#
#   bash envs/setup_envs.sh                  # environments + ASCAT (small downloads)
#   bash envs/setup_envs.sh --vep-cache DIR  # also the VEP GRCh38 cache (~20 GB)
#   bash envs/setup_envs.sh --sig-genome     # also SigProfiler's GRCh38 reference (~3 GB)
#
# Re-running is safe: existing environments are updated, not recreated.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

# conda >= 23.10 uses the fast libmamba solver by default; mamba also works
CONDA=$(command -v mamba || command -v conda)

conda config --set channel_priority strict

for env in wgs manta ascat vep; do
  if conda env list | awk '{print $1}' | grep -qx "$env"; then
    echo "== updating $env"
    "$CONDA" env update -n "$env" -f "$HERE/$env.yml" --prune
  else
    echo "== creating $env"
    "$CONDA" env create -f "$HERE/$env.yml"
  fi
done

echo "== ASCAT from GitHub"
conda run --no-capture-output -n ascat Rscript -e \
  'if (!requireNamespace("ASCAT", quietly = TRUE))
     remotes::install_github("VanLoo-lab/ascat/ASCAT", upgrade = "never");
   library(ASCAT); cat("ASCAT", as.character(packageVersion("ASCAT")), "\n")'

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vep-cache)
      CACHE=${2:?--vep-cache needs a directory}; shift
      mkdir -p "$CACHE"
      conda run --no-capture-output -n vep vep_install -a cf -s homo_sapiens -y GRCh38 \
        -c "$CACHE" --CONVERT
      echo "set VEP_CACHE = Path(\"$CACHE\") in pipeline_prostate.py"
      ;;
    --sig-genome)
      conda run --no-capture-output -n wgs python -c \
        'from SigProfilerMatrixGenerator import install as g; g.install("GRCh38")'
      ;;
  esac
  shift
done

echo "== versions"
conda run -n wgs   gatk --version 2>/dev/null | head -1
conda run -n wgs   samtools --version | head -1
conda run -n wgs   bcftools --version | head -1
conda run -n wgs   mosdepth --version
conda run -n manta configManta.py --version
conda run -n vep   vep --help 2>/dev/null | grep -m1 -i 'version' || true
echo "done: conda activate wgs"
