#!/bin/bash
#SBATCH --job-name=HP4_align_tree
#SBATCH --output=HP4_align_tree_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00


# --- SETUP ENVIRONMENT (Do this once) ---
CONDA_PATH="/home/mbata001/envs/miniconda3/etc/profile.d/conda.sh"
source $CONDA_PATH
conda activate hybphaser_v1

# --- DEFINE PATHS ---
BASE_DIR="/home/mbata001/seqs/target/hybphaser/all/HP4-non-hybrids"
INPUT_DIR="$BASE_DIR/merged_loci_consensus"
ALIGNED_DIR="$BASE_DIR/loci_aligned"
DISCARD_DIR="$BASE_DIR/loci_discarded"
GENETREE_DIR="$BASE_DIR/gene_trees"
MIN_TAXA=4
ASTRAL_JAR="/home/mbata001/envs/miniconda3/envs/hybphaser_v1/share/astral-tree-5.7.8-1/astral.5.7.8.jar"


mkdir -p "$ALIGNED_DIR"
mkdir -p "$DISCARD_DIR"
mkdir -p "$GENETREE_DIR"

################################### ALIGNING ##################################

echo "--- Starting MAFFT Alignment (Parallel) ---"

# OPTIMIZATION: Run 16 MAFFT jobs at once, 1 thread each. Much faster for gene lists.
find "$INPUT_DIR" -name "*.fasta" -print0 | xargs -0 -P 16 -I {} bash -c '
    INPUT_FILE="{}"
    ALIGNED_DIR=$1
    BASE_NAME=$(basename "$INPUT_FILE")
    OUTPUT_FILE="$ALIGNED_DIR/$BASE_NAME"

    if [ -s "$INPUT_FILE" ]; then
        if [ ! -f "$OUTPUT_FILE" ]; then
            # --auto automatically selects strategy (FFT-NS-2 for large, L-INS-i for small)
            mafft --thread 1 --auto --quiet "$INPUT_FILE" > "$OUTPUT_FILE"
            echo "Aligned: $BASE_NAME"
        else
            echo "Alignment for $BASE_NAME previously done; skipping"
        fi
    fi
' _ "$ALIGNED_DIR"

echo "Alignment complete."

################################### FILTERING #################################

echo "--- Starting Filtering (Min Taxa: $MIN_TAXA) ---"

removed_count=0
kept_count=0

for file in "$ALIGNED_DIR"/*.fasta; do
    [ -e "$file" ] || continue
    taxa_count=$(grep -c ">" "$file")

    if [ "$taxa_count" -lt "$MIN_TAXA" ]; then
        mv "$file" "$DISCARD_DIR/"
        ((removed_count++))
    else
        ((kept_count++))
    fi
done

echo "Filtered. Kept: $kept_count | Discarded: $removed_count"

################################### GENE TREES #################################

echo "Starting IQ-TREE 2 (Parallel)..."

# Run IQ-TREE 2

find "$ALIGNED_DIR" -name "*.fasta" -print0 | xargs -0 -P 16 -I {} bash -c '
    INPUT_FILE="{}"
    OUT_DIR=$1
    BASE_NAME=$(basename "$INPUT_FILE" .fasta)
    OUT_PREFIX="$OUT_DIR/$BASE_NAME"

    # Only run if output does not exist (resume capability)
    if [ ! -f "${OUT_PREFIX}.treefile" ]; then
         iqtree2 -s "$INPUT_FILE" -m MFP -B 1000 -T 1 --quiet --prefix "$OUT_PREFIX"
    fi
' _ "$GENETREE_DIR"


echo "Concatenating trees..."
# Combine all resulting trees
cat "$GENETREE_DIR"/*.treefile > "$GENETREE_DIR/all_gene_trees.newick"


##################################### ASTRAL ##################################

cd "$GENETREE_DIR"

echo "Collapsing low support branches (BS < 10)..."
nw_ed all_gene_trees.newick 'i & b < 10' o > all_gene_trees_collapsed.tre

echo "Running ASTRAL..."

# INCREASED RAM TO 16G
java -Xmx16g -jar "$ASTRAL_JAR" -i all_gene_trees_collapsed.tre -o astral_all_genes_phased.tre 2> astral.log

echo "Done."