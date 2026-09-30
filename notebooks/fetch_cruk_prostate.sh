#!/usr/bin/env bash
# =============================================================================
# CRUK-ICGC Prostate Cancer WGS (EGA, DAC EGAC00001000010)
#   EGAD00001004125  normal prostatectomy project (30 PCa + 7 non-cancer men), 71 BAMs, ~5.5 TB
#   EGAD00001000689  3 men, multiple tumour regions, 18 BAMs, ~2.4 TB
#
#   metadata  : public, wget + jq            -> ./ega_meta/
#   download  : controlled, pyega3 only      -> ./ega/   (needs approved EGA account)
#
# Usage:
#   ./fetch_cruk_prostate.sh metadata
#   ./fetch_cruk_prostate.sh download [file_ids.txt]   # default: all files of both datasets
#   ./fetch_cruk_prostate.sh verify
# Needs: wget, jq, pyega3 (pip install pyega3), md5sum
# =============================================================================
set -euo pipefail

DATASETS=(EGAD00001004125 EGAD00001000689)
API=https://metadata.ega-archive.org
META=ega_meta
OUT=ega
CRED=ega_credentials.json        # {"username":"you@inst.org","password":"..."}  chmod 600
CONNECTIONS=8
PAGE=100

# ---- public metadata: paginated GET until an empty page ---------------------
get_all() {   # $1 endpoint path, $2 output json
  local path=$1 out=$2 offset=0 page
  echo "[]" > "$out"
  while :; do
    page=$(wget -q -O - "${API}${path}?limit=${PAGE}&offset=${offset}")
    [[ $(jq 'length' <<<"$page") -eq 0 ]] && break
    jq -s '.[0] + .[1]' "$out" <(echo "$page") > "$out.tmp" && mv "$out.tmp" "$out"
    offset=$((offset + PAGE))
  done
}

metadata() {
  mkdir -p "$META"
  wget -q -O "$META/dac_EGAC00001000010_datasets.json" "${API}/dacs/EGAC00001000010/datasets" || true
  for d in "${DATASETS[@]}"; do
    wget -q -O "$META/$d.json" "${API}/datasets/$d"
    get_all "/datasets/$d/files"   "$META/$d.files.json"
    get_all "/datasets/$d/samples" "$META/$d.samples.json"

    jq -r --arg d "$d" '.[] | [$d, .accession_id, .extension, .filesize,
                               .unencrypted_checksum] | @tsv' \
       "$META/$d.files.json"   > "$META/$d.files.tsv"
    jq -r --arg d "$d" '.[] | [$d, .accession_id, .subject_id, .title, .description,
                               .biological_sex, .phenotype] | @tsv' \
       "$META/$d.samples.json" > "$META/$d.samples.tsv"

    printf '%s\t%s files\t%s samples\t%.2f TB\n' "$d" \
      "$(wc -l < "$META/$d.files.tsv")" "$(wc -l < "$META/$d.samples.tsv")" \
      "$(jq '[.[].filesize] | add / 1e12' "$META/$d.files.json")"
  done
  cat "$META"/EGAD*.files.tsv   > "$META/all_files.tsv"     # dataset file_id ext bytes md5
  cat "$META"/EGAD*.samples.tsv > "$META/all_samples.tsv"   # dataset sample_id subject title desc sex phenotype
  echo "sample titles end in the tissue label (e.g. _N, _Blood, tumour codes): map them to T/B/H/N"
  echo "file<->sample links are not in the public API: use 'pyega3 files' after approval"
}

# ---- controlled data: pyega3 ------------------------------------------------
download() {
  [[ -f "$CRED" ]] || { echo "missing $CRED"; exit 1; }
  mkdir -p "$OUT"
  # file names (which carry sample IDs) are listed only for authorised users
  for d in "${DATASETS[@]}"; do
    pyega3 -cf "$CRED" files "$d" | tee "$META/$d.pyega3_files.txt"
  done
  local ids
  if [[ -n "${1:-}" ]]; then ids=$(cat "$1"); else ids=$(cut -f2 "$META/all_files.tsv"); fi
  for f in $ids; do
    pyega3 -cf "$CRED" -c "$CONNECTIONS" fetch "$f" --output-dir "$OUT" --max-retries 10
  done
}

verify() {   # pyega3 already checks MD5; this re-checks against the public metadata
  local bad=0
  while IFS=$'\t' read -r _ id _ _ md5; do
    local f
    f=$(find "$OUT/$id" -type f ! -name '*.md5' 2>/dev/null | head -1 || true)
    [[ -z "$f" ]] && { echo "MISSING $id"; bad=1; continue; }
    [[ $(md5sum "$f" | cut -d' ' -f1) == "$md5" ]] && echo "OK $id" || { echo "BAD $id"; bad=1; }
  done < "$META/all_files.tsv"
  return $bad
}

case "${1:-}" in
  metadata) metadata ;;
  download) download "${2:-}" ;;
  verify)   verify ;;
  *) sed -n '2,14p' "$0"; exit 1 ;;
esac
