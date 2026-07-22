#!/bin/bash
#PBS -N download_genomes_mrtA
#PBS -l select=1:ncpus=4:mem=16gb:scratch_ssd=50gb
#PBS -l walltime=2:00:00
#PBS -m ae
#PBS -M mail@mail.mail

# ============================================================================
# 02b_download_genomes.sh
# Download completed genomes and WGS assemblies that contain mrtA genes.
#   1. Search NCBI nuccore with the mrtA genome query
#   2. Fetch accession IDs
#   3. Download each accession with ncbi-acc-download
# ============================================================================

set -euo pipefail

# ---- CONFIGURATION ---------------------------------------------------------
STORAGE="/path/to/storage"
WORKDIR="${SCRATCHDIR}/mcrA_work"

QUERY_MRTA_GENOME='txid2157[ORGN] AND genome[Title] AND (mrtA[All fields] OR methyl coenzyme M reductase II subunit A[All fields] OR methyl coenzyme M reductase II subunit alpha[All fields] OR methyl coenzyme M reductase II alpha subunit[All fields] OR methyl-coenzyme M reductase II subunit A[All fields] OR methyl-coenzyme M reductase II alpha subunit[All fields] OR methyl-coenzyme M reductase II subunit alpha[All fields] OR methyl coenzyme M reductase II, subunit A[All fields] OR methyl coenzyme M reductase II, subunit alpha[All fields] OR methyl coenzyme M reductase II, alpha subunit[All fields] OR methyl-coenzyme M reductase II, subunit A[All fields] OR methyl-coenzyme M reductase II, alpha subunit[All fields] OR methyl-coenzyme M reductase II, subunit alpha[All fields]) NOT (uncultured archaeon[Title] OR uncultured methanogenic archaeon[Title] OR uncultured euryarchaeote[Title] OR uncultured methanogen[Title] OR euryarchaeota archaeon[Title] OR uncultured archeon[Title] OR uncultured rumen archaeon[Title] OR methanogenic archaeon enrichment culture[Title] OR uncultured marine archaeon[Title] OR uncultured soil archaeon[Title])'

# Number of parallel download jobs
PARALLEL_JOBS=4

# Optional: set your NCBI API key here
NCBI_API_KEY=""
# ----------------------------------------------------------------------------

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

# Activate ncbi_tools
module load mambaforge
mamba activate /path/to/.conda/envs/ncbi_tools

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Searching NCBI nuccore for mrtA genome entries ..."
esearch -db nuccore -query "${QUERY_MRTA_GENOME}" | \
    efetch -format acc > accessions.txt

TOTAL=$(wc -l < accessions.txt)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found ${TOTAL} accessions"

if [[ "${TOTAL}" -eq 0 ]]; then
    echo "No accessions found. Exiting."
    exit 0
fi

mkdir -p genbank_files
cd genbank_files

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading GenBank files with ncbi-acc-download (parallel=${PARALLEL_JOBS}) ..."

# ncbi-acc-download
# Sequential loop
if [[ "${PARALLEL_JOBS}" -eq 1 ]]; then
    while IFS= read -r acc; do
        ncbi-acc-download \
            --verbose \
            --format genbank \
            $(if [[ -n "${NCBI_API_KEY}" ]]; then echo "--api-key ${NCBI_API_KEY}"; fi) \
            "${acc}" || echo "FAIL ${acc}" >> ../failed_downloads.txt
    done < ../accessions.txt
else
    cat ../accessions.txt | xargs -P "${PARALLEL_JOBS}" -I {} \
        ncbi-acc-download \
            --verbose \
            --format genbank \
            $(if [[ -n "${NCBI_API_KEY}" ]]; then echo "--api-key ${NCBI_API_KEY}"; fi) \
            {} || echo "FAIL {}" >> ../failed_downloads.txt
fi

# Ensure extensions match GeneScoop (*.gbk) expectations
cd "${WORKDIR}"
for f in genbank_files/*.gbff genbank_files/*.gb; do
    [[ -f "$f" ]] && ln -s "$f" "${f%.*}.gbk" 2>/dev/null || true
done

DOWNLOADED=$(find genbank_files -type f \( -name "*.gbk" -o -name "*.gbff" -o -name "*.gb" \) | wc -l)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloaded ${DOWNLOADED} / ${TOTAL} files"

if [[ -f failed_downloads.txt ]]; then
    FAILED=$(wc -l < failed_downloads.txt)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: ${FAILED} failed downloads logged in failed_downloads.txt"
fi

# ---- Copy results back to storage -------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying results to storage ..."
OUT_STORAGE="${STORAGE}/02b_genomes"
mkdir -p "${OUT_STORAGE}"
rsync -a "${WORKDIR}/genbank_files/" "${OUT_STORAGE}/genbank_files/" || cp -r "${WORKDIR}/genbank_files" "${OUT_STORAGE}/"
cp "${WORKDIR}/accessions.txt" "${OUT_STORAGE}/"
[[ -f "${WORKDIR}/failed_downloads.txt" ]] && cp "${WORKDIR}/failed_downloads.txt" "${OUT_STORAGE}/"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 02b complete. Results in ${OUT_STORAGE}"
