#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PATH="${REPO_DIR}:$PATH"
unset SCAPTUREPATH

# 1. Verify shared helper functions (path resolution, strand validation, flag mapping)
source "${REPO_DIR}/scripts/scapture_common.sh"

scapture_resolve_path "${REPO_DIR}"
[ "${SCAPTUREPATH}" = "${REPO_DIR}/" ] || { echo "FAIL: SCAPTUREPATH=${SCAPTUREPATH}"; exit 1; }

[ "$(scapture_featurecounts_strand forward)" = "1" ] || exit 1
[ "$(scapture_featurecounts_strand reverse)" = "2" ] || exit 1
[ "$(scapture_featurecounts_strand unstranded)" = "0" ] || exit 1

[ "$(scapture_samtools_strand_flags + forward)" = "-F 0x10" ] || exit 1
[ "$(scapture_samtools_strand_flags - forward)" = "-f 0x10" ] || exit 1
[ "$(scapture_samtools_strand_flags + reverse)" = "-f 0x10" ] || exit 1
[ "$(scapture_samtools_strand_flags - reverse)" = "-F 0x10" ] || exit 1
[ "$(scapture_samtools_strand_flags + unstranded)" = "" ] || exit 1
[ "$(scapture_samtools_strand_flags - unstranded)" = "" ] || exit 1

if scapture_validate_strand invalid 2>/dev/null; then
    echo "FAIL: expected invalid strand to fail"
    exit 1
fi

# 2. Verify CLI help and invalid strand rejection on main entrypoints
for bin in scapture scapture_callpeak scapture_quant; do
    help_out=$("${REPO_DIR}/${bin}" --help 2>&1 || true)
    echo "$help_out" | grep -q -- "--strand"
done

if "${REPO_DIR}/scapture" -m PASquant --strand bogus -o /tmp/x -b /tmp/y --pas /tmp/z 2>/dev/null; then
    echo "FAIL: expected scapture --strand bogus to fail"
    exit 1
fi

# 3. Verify synthetic BAM strand filtering + featureCounts + DeepPASS + Normality.R
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

# Create synthetic tiny BAM with 1 '+' read (flag 0) and 1 '-' read (flag 16) and CB/UB tags
cat > "${WORKDIR}/tiny.sam" <<'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@RG	ID:rg1	SM:s1
read_fwd	0	chr1	101	255	20M	*	0	0	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	CB:Z:CELL1-1	UB:Z:UMI1	NH:i:1
read_rev	16	chr1	501	255	20M	*	0	0	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	CB:Z:CELL1-1	UB:Z:UMI2	NH:i:1
EOF
samtools view -bS "${WORKDIR}/tiny.sam" | samtools sort -o "${WORKDIR}/tiny.bam"
samtools index "${WORKDIR}/tiny.bam"

# Check samtools strand flags directly on tiny.bam
fwd_reads_plus_mode=$(samtools view -c $(scapture_samtools_strand_flags + forward) "${WORKDIR}/tiny.bam")
rev_reads_plus_mode=$(samtools view -c $(scapture_samtools_strand_flags + reverse) "${WORKDIR}/tiny.bam")
all_reads_unstranded=$(samtools view -c $(scapture_samtools_strand_flags + unstranded) "${WORKDIR}/tiny.bam")
[ "$fwd_reads_plus_mode" = "1" ] && samtools view $(scapture_samtools_strand_flags + forward) "${WORKDIR}/tiny.bam" | grep -q "read_fwd"
[ "$rev_reads_plus_mode" = "1" ] && samtools view $(scapture_samtools_strand_flags + reverse) "${WORKDIR}/tiny.bam" | grep -q "read_rev"
[ "$all_reads_unstranded" = "2" ]

# 4. Verify scapture_quant --strand forward vs reverse on synthetic PAS BED12
# PAS1 on '+' at [100,150) overlaps read_fwd; PAS2 on '-' at [500,550) overlaps read_rev.
# Wait: read_fwd is mapped to '+' (flag 0), so in 'forward' mode it matches PAS1 (+), NOT PAS2 (-).
# In 'reverse' mode (10x v1), '+' reads match '-' features (PAS2) and '-' reads match '+' features (PAS1).
cat > "${WORKDIR}/pas.bed" <<'EOF'
chr1	100	150	GenePlus|1|GenePlus|exon	100	+	100	150	0	1	50,	0,
chr1	500	550	GeneMinus|1|GeneMinus|exon	100	-	500	550	0	1	50,	0,
EOF

# Invert expected read alignment for PAS2 vs PAS1: both reads forward-mapped (flag 0)
cat > "${WORKDIR}/quant.sam" <<'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
read_p_fwd	0	chr1	101	255	20M	*	0	0	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	CB:Z:CELL1-1	UB:Z:UMI1	NH:i:1
read_m_fwd	0	chr1	501	255	20M	*	0	0	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	CB:Z:CELL1-1	UB:Z:UMI2	NH:i:1
EOF
samtools view -bS "${WORKDIR}/quant.sam" | samtools sort -o "${WORKDIR}/quant.bam"
samtools index "${WORKDIR}/quant.bam"

"${REPO_DIR}/scapture" -m PASquant -o "${WORKDIR}/out_fwd" -p 1 -b "${WORKDIR}/quant.bam" --pas "${WORKDIR}/pas.bed" --celllist FromBAM --strand forward
"${REPO_DIR}/scapture" -m PASquant -o "${WORKDIR}/out_rev" -p 1 -b "${WORKDIR}/quant.bam" --pas "${WORKDIR}/pas.bed" --celllist FromBAM --strand reverse

fwd_plus=$(zcat "${WORKDIR}/out_fwd.KeepCell.UMIs.tsv.gz" | awk '$1=="GenePlus-1"{print $2}')
fwd_minus=$(zcat "${WORKDIR}/out_fwd.KeepCell.UMIs.tsv.gz" | awk '$1=="GeneMinus-1"{print $2}')
rev_plus=$(zcat "${WORKDIR}/out_rev.KeepCell.UMIs.tsv.gz" | awk '$1=="GenePlus-1"{print $2}')
rev_minus=$(zcat "${WORKDIR}/out_rev.KeepCell.UMIs.tsv.gz" | awk '$1=="GeneMinus-1"{print $2}')

[ "$fwd_plus" = "1" ] && [ "$fwd_minus" = "0" ] || { echo "FAIL forward counts: +:${fwd_plus} -:${fwd_minus}"; exit 1; }
[ "$rev_plus" = "0" ] && [ "$rev_minus" = "1" ] || { echo "FAIL reverse counts: +:${rev_plus} -:${rev_minus}"; exit 1; }

# 5. Verify DeepPASS/Predict.py runs without ValueError on ragged/multi-seq input + unknown base N/R
python3 - <<'PY' "${REPO_DIR}/DeepPASS/Predict.py" "${WORKDIR}/seq.txt"
import sys
seq = "A" * 199 + "N"
with open(sys.argv[2], "w") as fh:
    fh.write(f"seq1\t{seq}\n")
    fh.write(f"seq2\t{'T'*199}R\n")
PY
python3 "${REPO_DIR}/DeepPASS/Predict.py" -m "${REPO_DIR}/DeepPASS/mm10_best_model.h5" -p "${WORKDIR}/seq.txt" -o "${WORKDIR}/deeppass_out"
[ "$(wc -l < "${WORKDIR}/deeppass_out/Predict_Result.txt")" = "2" ]

# 6. Verify Normality.R deterministic output with SCAPTURE_NORMALITY_SEED
python3 - <<'PY' "${WORKDIR}/cov.txt" "${WORKDIR}/peak.txt"
import sys
with open(sys.argv[1], "w") as fh:
    for i in range(1, 401):
        val = max(1, 50 - abs(i - 200) // 4)
        fh.write(f"chr1\t{i}\t{i+1}\t{i}\t{val}\t+\n")
with open(sys.argv[2], "w") as fh:
    fh.write("chr1\t50\t350\t200\t10\t+\n")
PY
cp "${WORKDIR}/peak.txt" "${WORKDIR}/peak1.txt"
cp "${WORKDIR}/peak.txt" "${WORKDIR}/peak2.txt"
SCAPTURE_NORMALITY_SEED=42 Rscript "${REPO_DIR}/scripts/Normality.R" "${WORKDIR}/peak1.txt" "${WORKDIR}/cov.txt" 400
SCAPTURE_NORMALITY_SEED=42 Rscript "${REPO_DIR}/scripts/Normality.R" "${WORKDIR}/peak2.txt" "${WORKDIR}/cov.txt" 400
# Different basename hashes should differ only if intentional; identical filename seed should match identically:
cp "${WORKDIR}/peak1.txt" "${WORKDIR}/peak1_copy.txt"
mkdir -p "${WORKDIR}/sub"
cp "${WORKDIR}/peak1.txt" "${WORKDIR}/sub/peak1.txt"
SCAPTURE_NORMALITY_SEED=42 Rscript "${REPO_DIR}/scripts/Normality.R" "${WORKDIR}/sub/peak1.txt" "${WORKDIR}/cov.txt" 400
diff -u "${WORKDIR}/peak1.Pvalue.txt" "${WORKDIR}/sub/peak1.Pvalue.txt"

echo "ALL_SCAPTURE_STRAND_CLI_TESTS_PASSED"
