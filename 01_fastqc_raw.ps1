# ===============================
# RNA-seq Preprocessing Pipeline
# Windows PowerShell Version (Kallisto)
# ===============================

$THREADS = 8

$PROJECT_DIR = Get-Location

$RAW_DIR      = "$PROJECT_DIR\data\raw_fastq"
$TRIM_DIR     = "$PROJECT_DIR\data\trimmed_fastq"
$QC_DIR_RAW   = "$PROJECT_DIR\qc\raw"
$QC_DIR_TRIM  = "$PROJECT_DIR\qc\trim"
$REF_DIR      = "$PROJECT_DIR\reference"
$META_FILE    = "$PROJECT_DIR\metadata\sample_metadata.csv"
$KALLISTO_DIR = "$PROJECT_DIR\kallisto"

$ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
$ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Ensure directories exist
New-Item -ItemType Directory -Force -Path $TRIM_DIR | Out-Null
New-Item -ItemType Directory -Force -Path $QC_DIR_RAW | Out-Null
New-Item -ItemType Directory -Force -Path $QC_DIR_TRIM | Out-Null
New-Item -ItemType Directory -Force -Path $KALLISTO_DIR | Out-Null


# ---------------------------
# 2. FASTQC
# ---------------------------

Write-Host "  FastQC..."

perl  .\fastqc\fastqc `
    $RAW_DIR\*.fastq.gz `
    -o $QC_DIR_RAW `
    --threads $THREADS

Write-Host "Preprocessing complete!"