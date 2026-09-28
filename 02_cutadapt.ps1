# ===============================
# RNA-seq Preprocessing Pipeline
# Windows PowerShell Version (Kallisto)
# ===============================

$THREADS = 8

$PROJECT_DIR = Get-Location

$RAW_DIR      = "$PROJECT_DIR\data\raw_fastq"
$TRIM_DIR     = "$PROJECT_DIR\data\trimmed_fastq"
$QC_DIR       = "$PROJECT_DIR\qc"
$REF_DIR      = "$PROJECT_DIR\reference"
$META_FILE    = "$PROJECT_DIR\metadata\sample_metadata.csv"
$KALLISTO_DIR = "$PROJECT_DIR\kallisto"

$ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
$ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Ensure directories exist
New-Item -ItemType Directory -Force -Path $TRIM_DIR | Out-Null
New-Item -ItemType Directory -Force -Path $QC_DIR | Out-Null
New-Item -ItemType Directory -Force -Path $KALLISTO_DIR | Out-Null

Write-Host "Loading metadata..."

$metadata = Import-Csv $META_FILE

foreach ($sample in $metadata) {

    $sampleID = $sample.SampleID
    $condition = $sample.Condition

    Write-Host "Processing sample: $sampleID ($condition)"

    # Input FASTQs
    $R1_IN = "$RAW_DIR\${sampleID}_R1.fastq.gz"
    $R2_IN = "$RAW_DIR\${sampleID}_R2.fastq.gz"

    # Output FASTQs
    $R1_OUT = "$TRIM_DIR\${sampleID}_R1.trimmed.fastq.gz"
    $R2_OUT = "$TRIM_DIR\${sampleID}_R2.trimmed.fastq.gz"

    # ---------------------------
    # 1. CUTADAPT
    # ---------------------------

    Write-Host "  Cutadapt trimming..."

    cutadapt `
        -a $ADAPTER_R1 `
        -A $ADAPTER_R2 `
        -o $R1_OUT `
        -p $R2_OUT `
        $R1_IN $R2_IN `
        --minimum-length 20 `
        --cores $THREADS
}
Write-Host "Preprocessing complete!"

