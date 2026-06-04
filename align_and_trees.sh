#!/bin/bash
#SBATCH --job-name=align_tree
#SBATCH --output=align_tree_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00


#align_and_trees.sh
#Usage: align_and_trees.sh -j path/to/astral.jar -d path/to/base_directory/


# --- DEFINE PATHS AND PARAMETERS ---
MIN_TAXA=4
MIN_SUPPORT=50
ASTRAL_JAR=""
BASE_DIR=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -j|--jar) ASTRAL_JAR="$2"; shift ;;
        -d|--dir) BASE_DIR="$2"; shift ;;
        -samples) MIN_TAXA="$2"; shift ;;
        -support) MIN_SUPPORT="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# --- PARSE ARGUMENTS ---
if [ -z "$ASTRAL_JAR" ]; then
    echo "No path to astral jar file provided. Rerun using 'align_and_trees.sh -j path/to/astral.jar -d path/to/base_directory/'"
    exit 1
fi

#Check if BASE_DIR was provided
if [ -z "$BASE_DIR" ]; then
    echo "No directory provided. Using current directory: $(pwd)"
    BASE_DIR=$(pwd)
fi
INPUT_DIR="$BASE_DIR/merged_loci"
ALIGNED_DIR="$BASE_DIR/loci_aligned"
DISCARD_DIR="$BASE_DIR/loci_discarded"
GENETREE_DIR="$BASE_DIR/gene_trees"


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

echo "Starting IQ-TREE 2"

find "$ALIGNED_DIR" -name "*.fasta" -print0 | xargs -0 -P 16 -I {} bash -c '
    INPUT_FILE="{}"
    OUT_DIR=$1
    BASE_NAME=$(basename "$INPUT_FILE" .fasta)
    OUT_PREFIX="$OUT_DIR/$BASE_NAME"

    # Only run if output does not exist (resume capability)
    if [ ! -f "${OUT_PREFIX}.treefile" ]; then
         iqtree2 -s "$INPUT_FILE" -m MFP -B 1000 -T 1 --quiet --prefix "$OUT_PREFIX"
         echo "Made $BASE_NAME tree"
    else
        echo "$BASE_NAME previously done; skipping"
    fi
' _ "$GENETREE_DIR"


echo "Concatenating trees..."
# Combine all resulting trees
cat "$GENETREE_DIR"/*.treefile > "$GENETREE_DIR/all_gene_trees.newick"


##################################### ASTRAL ##################################

cd "$GENETREE_DIR"

echo "Collapsing low support branches"
nw_ed all_gene_trees.newick 'i & b < '"$MIN_SUPPORT" o > all_gene_trees_collapsed.tre
echo "Running ASTRAL..."

# -Xmx16g to increase RAM to 16G
java -Xmx16g -jar "$ASTRAL_JAR" -i all_gene_trees_collapsed.tre -o astral_all_genes.tre 2> astral.log

echo "Done."