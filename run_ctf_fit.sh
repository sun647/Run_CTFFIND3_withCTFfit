#!/bin/bash

SCRIPT=fit_ctf_resolution.py
PIXEL=2.5
KV=300
CS=2.7
AMP=0.10
THR=0.30

OUT_SUMMARY=ctf_fit_resolution_all.csv

# CSV header
echo "micrograph,ctf_fit_resolution_A" > "$OUT_SUMMARY"

for power_mrc in *patch_aligned_ctf.mrc; do
  base=${power_mrc%.mrc}
  log="${base%_ctf}_ctffind3.log"

  [[ ! -f "$log" ]] && { echo "Missing log for $power_mrc, skipping"; continue; }

  echo "Processing $power_mrc"

  result=$(python3 "$SCRIPT" \
    --power_mrc "$power_mrc" \
    --log "$log" \
    --pixel_A "$PIXEL" \
    --kv "$KV" \
    --cs_mm "$CS" \
    --amp_contrast "$AMP" \
    --thr "$THR" \
    --out_csv "${base}_fitcurve.csv" \
    | grep "CTF-fit resolution" \
    | awk '{print $(NF-1)}')

  [[ -z "$result" ]] && result="NA"

  echo "$power_mrc,$result" >> "$OUT_SUMMARY"
done

echo "Done. Wrote $OUT_SUMMARY"
