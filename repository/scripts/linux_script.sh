#!/bin/bash
#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vt-pid@vt.edu # Change to whichever email you would like to receive job updates
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --output=full_pipeline_%j.out

#cd to directory (NEEDS TO CHANGE)
cd /PATH/TO/REPOSITORY/DIRECTORY

#Set Conda Environment
module load Miniconda3
module load Biopython

source ~/.bashrc
conda activate indiv_prod_mg

set -euo pipefail

###############################################
# Input files/directorys/databases (Expected script outputs for both Linux and RStudio are located in expected_outputs folder)
###############################################

INPUT_DIR="/projects/intro2gds/I2GDS2025/G4_Viruses/mitchell/individual_project/inputs/sequences" #DO NOT CHANGE
OUTPUT_DIR="outputs/kraken_output"
REPORT_DIR="outputs/reports"
EXTRACT_DIR="outputs/extracted"
SPADES_DIR="outputs/spades_out"

KRAKEN_DB="/projects/intro2gds/I2GDS2025/G4_Viruses/mitchell/individual_project/database/k2_cust_viral_db" #DO NOT CHANGE
KRAKEN_TOOLS="scripts/extract_kraken_reads.py"

THREADS=16

# Make log file for debugging
LOG_DIR="logs"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/pipeline_steps_$(date +%Y%m%d_%H%M%S).log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR" "$EXTRACT_DIR" "$SPADES_DIR"

shopt -s nullglob

log "Scanning for R1 files in $INPUT_DIR"
ls -lh "$INPUT_DIR"/*_R1_*.fastq* 2>/dev/null || {
    log "[ERROR] No R1 FASTQ files found — exiting"
    exit 1
}

# MAIN LOOP: PER SAMPLE
# Split samples by common naming across R1 and R2 files
for R1 in "$INPUT_DIR"/*_R1_*.fastq*; do

    R2="${R1/_R1_/_R2_}"
    if [[ ! -e "$R2" ]]; then
        log "[WARN] Missing R2 for: $R1"
        continue
    fi

    sample="$(basename "$R1" | sed -E 's/_R1_.*$//')"
    log "============================================"
    log "Processing sample: $sample"
    log "============================================"

    ###############################################
    # STEP 1 — RUN KRAKEN2
    ###############################################
    # Run k2 classify for each sample R1 and R2 
    kraken_file="$OUTPUT_DIR/${sample}.kraken"
    k2 classify \
        --db "$KRAKEN_DB" \
        --paired "$R1" "$R2" \
        --threads "$THREADS" \
        --output "$kraken_file"

    log "[INFO] Kraken2 classification finished for $sample"

    ###############################################
    # STEP 2 — GET TOP 3 TAXIDS WITH READ COUNTS
    ###############################################
    # Pull out top three read counts from kraken2 classification to extract and run assembly on
    top_taxids=$(awk '$1=="C"{count[$3]++} END{for(t in count){print count[t], t}}' "$kraken_file" \
                | sort -k1,1nr | head -n 3)

    if [[ -z "$top_taxids" ]]; then
        log "[WARN] No classified reads found for $sample — skipping"
        continue
    fi

    taxid_array=()
    read_count_array=()
    while read -r reads taxid; do
        taxid_array+=("$taxid")
        read_count_array+=("$reads")
        log "Most abundant classified taxid = $taxid (reads = $reads)"
    done <<< "$top_taxids"

    ###############################################
    # STEP 3 — EXTRACT READS FOR ALL TOP TAXIDS
    ###############################################
    # Extract the top three taxid reads from each sample using the extract_kraken_reads.py script from KrakenTools
    log "[INFO] Extracting reads for top taxids..."
    for i in "${!taxid_array[@]}"; do
        taxid="${taxid_array[$i]}"
        read_count="${read_count_array[$i]}"

        out_R1="$EXTRACT_DIR/${sample}_${taxid}_R1.fastq"
        out_R2="$EXTRACT_DIR/${sample}_${taxid}_R2.fastq"

        python3 "$KRAKEN_TOOLS" \
            -k "$kraken_file" \
            -s1 "$R1" \
            -s2 "$R2" \
            --taxid "$taxid" \
            -o "$out_R1" \
            -o2 "$out_R2" \
            --fastq-output

        log "[INFO] Extracted reads for taxid=$taxid (reads=$read_count)"
    done
    
    ###############################################
    # STEP 3.5 — GZIP EXTRACTED READS
    ###############################################
    log "[INFO] Gzipping extracted FASTQ files..."

    for taxid in "${taxid_array[@]}"; do
        for mate in R1 R2; do
            fq="$EXTRACT_DIR/${sample}_${taxid}_${mate}.fastq"
            fq_gz="${fq}.gz"

            if [[ -s "$fq" && ! -s "$fq_gz" ]]; then
                log "[INFO] Gzipping $(basename "$fq")"
                gzip -f "$fq"
            fi
        done
    done
    
    ###############################################
    # STEP 4 — RUN SPADES FOR ALL EXTRACTED READ SETS
    ###############################################
    log "[INFO] Running SPAdes for all top taxids..."

    for taxid in "${taxid_array[@]}"; do
        out_R1="$EXTRACT_DIR/${sample}_${taxid}_R1.fastq.gz"
        out_R2="$EXTRACT_DIR/${sample}_${taxid}_R2.fastq.gz"
        sample_spades="$SPADES_DIR/${sample}_${taxid}"
        mkdir -p "$sample_spades"

        if [[ ! -s "$out_R1" || ! -s "$out_R2" ]]; then
            log "[WARN] Missing gzipped reads for sample=$sample taxid=$taxid — skipping SPAdes"
            continue
        fi

        log "[INFO] SPAdes assembly for taxid=$taxid"

        set +e
        spades.py \
            -1 "$out_R1" \
            -2 "$out_R2" \
            -o "$sample_spades" \
            -t "$THREADS"
        spades_exit=$?
        set -e

        if [[ $spades_exit -ne 0 ]]; then
            log "[WARN] SPAdes failed for sample=$sample taxid=$taxid"
            continue
        fi

        log "[INFO] SPAdes completed for taxid=$taxid"
    done
    log "============================================"
    log "Sample $sample processing complete (all top taxids)"
    log "============================================"
done

log "============================================"
log "Pipeline complete for all samples!"
log "============================================"
