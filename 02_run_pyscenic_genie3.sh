#!/usr/bin/env bash
# ============================================================================
# pySCENIC on the Callari bulk cohort, using GENIE3 for network inference.
#
# GENIE3 is the original SCENIC network algorithm and was developed for
# microarray / bulk expression data, which is what the Callari cohort is
# (custom cDNA arrays, 37 MBC + 53 FBC).
#
#   Step 1  pyscenic grn   --method genie3   -> adjacencies (TF -> target)
#   Step 2  pyscenic ctx                     -> motif-pruned regulons
#   Step 3  pyscenic aucell                  -> per-sample regulon activity
#
# Runtime: step 1 is the slow one (~2-5 h for 7.8k genes x 90 samples on 8
# cores); steps 2-3 take minutes.
# ============================================================================
set -euo pipefail

NW=${NW:-8}                       # workers; override with  NW=16 ./02_...sh
SEED=777

EXPR=Input/callari_expr.csv
RES=resources
TFS=$RES/allTFs_hg38.txt
DB=$RES/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
ANNO=$RES/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

OUT=Output_genie3
mkdir -p $OUT
for f in "$EXPR" "$TFS" "$DB" "$ANNO"; do
  [ -f "$f" ] || { echo "missing: $f"; exit 1; }
done

# ---- Step 1: GRN inference with GENIE3 ------------------------------------
echo "[1/3] GRN inference (GENIE3) ..."
pyscenic grn "$EXPR" "$TFS" \
    --method genie3 \
    -o $OUT/adjacencies_genie3.tsv \
    --num_workers $NW \
    --seed $SEED

# ---- Step 2: cisTarget motif pruning -> regulons --------------------------
# NOTE: --mask_dropouts is a single-cell option (it ignores zero counts when
# computing correlations). Bulk microarray data has no dropouts, so it is
# deliberately NOT used here.
echo "[2/3] Motif pruning (cisTarget) ..."
pyscenic ctx $OUT/adjacencies_genie3.tsv "$DB" \
    --annotations_fname "$ANNO" \
    --expression_mtx_fname "$EXPR" \
    --mode custom_multiprocessing \
    -o $OUT/regulons_genie3.csv \
    --num_workers $NW

# ---- Step 3: AUCell per-sample regulon activity ---------------------------
echo "[3/3] AUCell scoring ..."
pyscenic aucell "$EXPR" $OUT/regulons_genie3.csv \
    -o $OUT/auc_genie3.csv \
    --num_workers $NW \
    --seed $SEED

echo
echo "Done. Outputs in $OUT/:"
ls -lh $OUT
echo "Next: knit 03_regulon_analysis.Rmd"
