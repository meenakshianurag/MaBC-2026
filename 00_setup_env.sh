#!/usr/bin/env bash
# ============================================================================
# pySCENIC environment + cisTarget resources
# Run once. Everything lands in ./resources/
# ============================================================================
set -euo pipefail

# ---- 1. conda environment -------------------------------------------------
# pySCENIC 0.12.1 needs python 3.10 and pinned numpy/pandas; a clean env avoids
# the dependency conflicts that are common with newer numpy.
conda create -y -n pyscenic python=3.10
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pyscenic
pip install pyscenic==0.12.1

python -c "import pyscenic; print('pySCENIC', pyscenic.__version__)"

# ---- 2. cisTarget resources (human hg38) ----------------------------------
mkdir -p resources && cd resources

# TF list
curl -O https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt

# motif -> TF annotation table (v10nr_clust)
curl -O https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

# ranking database: 10 kb up/down of TSS, gene-based (~1.1 GB)
curl -O https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

ls -lh
cd ..
echo "Setup done."
