#!/bin/bash
# Generate test datasets for CRAM performance benchmarking
#
# This script converts BAM files to CRAM and creates test datasets
# of various sizes for benchmarking purposes.

set -e

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Generate test datasets for CRAM benchmarking.

Required:
  -b, --bam FILE          Input BAM file
  -r, --ref FILE          Reference FASTA file
  -o, --output DIR        Output directory

Optional:
  -s, --sizes SIZE_LIST   Comma-separated sizes: small,medium,large (default: all)
  -c, --contig NAME       Specific contig to extract (default: chr20 or smallest available)
  --small-size MB         Small dataset size in MB (default: 1)
  --medium-size MB        Medium dataset size in MB (default: 100)
  --large-size MB         Large dataset size in MB (default: 1000)
  -t, --threads N         Number of threads (default: 4)
  -h, --help              Show this help message

Examples:
  # Generate all test sizes
  $0 -b input.bam -r ref.fa -o test_data/

  # Generate only small and medium
  $0 -b input.bam -r ref.fa -o test_data/ -s small,medium

  # Custom sizes and contig
  $0 -b input.bam -r ref.fa -o test_data/ -c chr22 --small-size 5
EOF
}

# Default parameters
SIZES="small,medium,large"
SMALL_SIZE_MB=1
MEDIUM_SIZE_MB=100
LARGE_SIZE_MB=1000
THREADS=4
CONTIG=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--bam)
            BAM_FILE="$2"
            shift 2
            ;;
        -r|--ref)
            REF_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -s|--sizes)
            SIZES="$2"
            shift 2
            ;;
        -c|--contig)
            CONTIG="$2"
            shift 2
            ;;
        --small-size)
            SMALL_SIZE_MB="$2"
            shift 2
            ;;
        --medium-size)
            MEDIUM_SIZE_MB="$2"
            shift 2
            ;;
        --large-size)
            LARGE_SIZE_MB="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$BAM_FILE" ]] || [[ -z "$REF_FILE" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments"
    usage
    exit 1
fi

# Validate input files
if [[ ! -f "$BAM_FILE" ]]; then
    echo "Error: BAM file not found: $BAM_FILE"
    exit 1
fi

if [[ ! -f "$REF_FILE" ]]; then
    echo "Error: Reference file not found: $REF_FILE"
    exit 1
fi

# Check for required tools
for tool in samtools bcftools; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool not found in PATH"
        exit 1
    fi
done

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=================================================="
echo "CRAM Test Data Generation"
echo "=================================================="
echo "Input BAM:    $BAM_FILE"
echo "Reference:    $REF_FILE"
echo "Output:       $OUTPUT_DIR"
echo "Sizes:        $SIZES"
echo "Threads:      $THREADS"
echo "=================================================="
echo ""

# Function to estimate file size
get_file_size_mb() {
    local file=$1
    local size_bytes=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
    echo $(( size_bytes / 1024 / 1024 ))
}

# Function to get available contigs
get_contigs() {
    samtools view -H "$BAM_FILE" | grep '^@SQ' | cut -f2 | cut -d: -f2
}

# Function to get contig length
get_contig_length() {
    local ctg=$1
    samtools view -H "$BAM_FILE" | grep "^@SQ.*SN:${ctg}" | sed 's/.*LN://' | cut -f1
}

# Select contig if not specified
if [[ -z "$CONTIG" ]]; then
    echo "Auto-selecting contig..."

    # Try common small contigs first
    for try_ctg in chr20 chr22 chrM 20 22 MT; do
        if get_contigs | grep -q "^${try_ctg}$"; then
            CONTIG=$try_ctg
            echo "Selected contig: $CONTIG"
            break
        fi
    done

    # If no common contig found, use the smallest
    if [[ -z "$CONTIG" ]]; then
        echo "Analyzing contig sizes..."
        CONTIG=$(get_contigs | head -1)
        min_length=$(get_contig_length "$CONTIG")

        for ctg in $(get_contigs); do
            length=$(get_contig_length "$ctg")
            if [[ $length -lt $min_length ]]; then
                CONTIG=$ctg
                min_length=$length
            fi
        done
        echo "Selected smallest contig: $CONTIG (length: $min_length bp)"
    fi
fi

# Verify contig exists
if ! get_contigs | grep -q "^${CONTIG}$"; then
    echo "Error: Contig '$CONTIG' not found in BAM file"
    echo "Available contigs:"
    get_contigs | head -10
    exit 1
fi

CONTIG_LENGTH=$(get_contig_length "$CONTIG")
echo "Using contig: $CONTIG (length: $CONTIG_LENGTH bp)"
echo ""

# Function to create test dataset
create_test_dataset() {
    local size_name=$1
    local target_size_mb=$2
    local output_prefix="${OUTPUT_DIR}/${size_name}"

    echo "=================================================="
    echo "Creating ${size_name} dataset (target: ${target_size_mb} MB)"
    echo "=================================================="

    # Calculate region to extract based on target size
    # Rough estimate: 100x coverage ~ 100 MB/Mb for BAM
    local coverage_factor=100
    local target_region_mb=$(echo "scale=2; $target_size_mb / $coverage_factor" | bc)
    local target_region_bp=$(echo "$target_region_mb * 1000000" | bc | cut -d. -f1)

    # Cap at contig length
    if [[ $target_region_bp -gt $CONTIG_LENGTH ]]; then
        target_region_bp=$CONTIG_LENGTH
        echo "Warning: Target size exceeds contig length, using full contig"
    fi

    local region="${CONTIG}:1-${target_region_bp}"
    echo "Extracting region: $region"

    # Extract BAM region
    echo "[1/4] Extracting BAM region..."
    samtools view -@ $THREADS -b -h "$BAM_FILE" "$region" > "${output_prefix}.bam"
    samtools index -@ $THREADS "${output_prefix}.bam"

    local bam_size_mb=$(get_file_size_mb "${output_prefix}.bam")
    echo "  BAM size: ${bam_size_mb} MB"

    # Convert to CRAM
    echo "[2/4] Converting to CRAM..."
    samtools view -@ $THREADS -C -T "$REF_FILE" -o "${output_prefix}.cram" "${output_prefix}.bam"
    samtools index -@ $THREADS "${output_prefix}.cram"

    local cram_size_mb=$(get_file_size_mb "${output_prefix}.cram")
    echo "  CRAM size: ${cram_size_mb} MB"
    echo "  Compression ratio: $(echo "scale=2; $bam_size_mb / $cram_size_mb" | bc)x"

    # Verify CRAM integrity
    echo "[3/4] Verifying CRAM integrity..."
    local bam_count=$(samtools view -c "${output_prefix}.bam")
    local cram_count=$(samtools view -c -T "$REF_FILE" "${output_prefix}.cram")

    if [[ "$bam_count" != "$cram_count" ]]; then
        echo "Error: Read count mismatch (BAM: $bam_count, CRAM: $cram_count)"
        exit 1
    fi
    echo "  Verified: $cram_count reads"

    # Create summary file
    echo "[4/4] Creating summary..."
    cat > "${output_prefix}_info.txt" << SUMMARY
Dataset: ${size_name}
Region: ${region}
Reads: ${bam_count}
BAM Size: ${bam_size_mb} MB
CRAM Size: ${cram_size_mb} MB
Compression Ratio: $(echo "scale=2; $bam_size_mb / $cram_size_mb" | bc)x
Created: $(date)
SUMMARY

    echo "Completed: ${size_name} dataset"
    echo ""
}

# Generate requested datasets
IFS=',' read -ra SIZE_ARRAY <<< "$SIZES"
for size in "${SIZE_ARRAY[@]}"; do
    case $size in
        small)
            create_test_dataset "small" $SMALL_SIZE_MB
            ;;
        medium)
            create_test_dataset "medium" $MEDIUM_SIZE_MB
            ;;
        large)
            create_test_dataset "large" $LARGE_SIZE_MB
            ;;
        *)
            echo "Warning: Unknown size '$size', skipping"
            ;;
    esac
done

# Create reference subset for the contig
echo "=================================================="
echo "Extracting reference subset"
echo "=================================================="
samtools faidx "$REF_FILE" "$CONTIG" > "${OUTPUT_DIR}/ref_${CONTIG}.fa"
samtools faidx "${OUTPUT_DIR}/ref_${CONTIG}.fa"
echo "Reference subset: ${OUTPUT_DIR}/ref_${CONTIG}.fa"
echo ""

# Summary
echo "=================================================="
echo "SUMMARY"
echo "=================================================="
echo "Generated datasets in: $OUTPUT_DIR"
echo ""
ls -lh "$OUTPUT_DIR"/*.{bam,cram,fa} 2>/dev/null | awk '{print $9, $5}'
echo ""
echo "To run benchmarks:"
echo "  python benchmarks/benchmark_bam_vs_cram.py \\"
echo "    --bam ${OUTPUT_DIR}/small.bam \\"
echo "    --cram ${OUTPUT_DIR}/small.cram \\"
echo "    --ref ${OUTPUT_DIR}/ref_${CONTIG}.fa \\"
echo "    --truth_vcf <your_truth.vcf> \\"
echo "    --output benchmark_results/"
echo ""
echo "Done!"
