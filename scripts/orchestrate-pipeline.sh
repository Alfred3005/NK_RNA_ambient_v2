#!/bin/bash
# scripts/orchestrate-pipeline.sh
# PHOENIX_PROTOCOL: Phase 02 Automation

IN_DIR="data/processed/segments"
OUT_DIR="data/processed/corrected"
mkdir -p "$OUT_DIR"

echo "--- Starting PHOENIX Orchestrator: Phase 02 (Ambient Correction) ---"

FILES=$(ls $IN_DIR/*.h5ad)
COUNT=$(echo "$FILES" | wc -l)
CURRENT=0

for f in $FILES; do
    CURRENT=$((CURRENT + 1))
    BASENAME=$(basename "$f" .h5ad)
    OUT_FILE="$OUT_DIR/$BASENAME.rds"
    
    if [ -f "$OUT_FILE" ]; then
        echo "[$CURRENT/$COUNT] Skipping $BASENAME (Already exists)"
        continue
    fi
    
    echo "[$CURRENT/$COUNT] Correcting $BASENAME..."
    Rscript scripts/02-ambient-correction.R -i "$f" -o "$OUT_FILE"
done

echo "--- PHOENIX ORCHESTRATOR: Phase 02 Completed ---"
