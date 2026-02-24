set -euo pipefail

TAR="refdata-gex-GRCh38-2024-A.tar.gz"
WORKDIR="$PWD/GRCh38_2024A_MbpSplitWork"
mkdir -p "$WORKDIR"
tar -xzf "$TAR" -C "$WORKDIR"

# Find the unpacked ref directory
REFDIR="$(find "$WORKDIR" -maxdepth 1 -type d -name 'refdata-gex-GRCh38-2024-A*' | head -n 1)"
echo "REFDIR=$REFDIR"

REFDIR="/storage1/fs1/jin810/Active/References/GRCh38_10X/GRCh38_2024A_MbpSplitWork/refdata-gex-GRCh38-2024-A"

GTF_GZ="${REFDIR}/genes/genes.gtf.gz"
GTF_IN="${REFDIR}/genes/genes.gtf"   # we’ll create this
echo "GTF_GZ=$GTF_GZ"
ls -lh "$GTF_GZ"

# unzip to a plain .gtf for easier awk/grep
gunzip -c "$GTF_GZ" > "$GTF_IN"
ls -lh "$GTF_IN"

# Confirm expected files exist
ls -lh "$REFDIR/genes/genes.gtf" "$REFDIR/fasta/genome.fa"


REFDIR="/storage1/fs1/jin810/Active/References/GRCh38_10X/GRCh38_2024A_MbpSplitWork/refdata-gex-GRCh38-2024-A"
GTF_IN="${REFDIR}/genes/genes.gtf"

MBP_GENE_ID="$(awk -F'\t' '
  $0 !~ /^#/ && $3=="gene" && $9 ~ /gene_name "MBP"/ {
    if (match($9, /gene_id "([^"]+)"/, a)) { print a[1]; exit }
  }
' "$GTF_IN")"

echo "MBP_GENE_ID=$MBP_GENE_ID"

# Quick sanity check: show the gene line we matched
grep -m 1 'gene_name "MBP"' "$GTF_IN" | head

# Golli Exon
ENSE00003502747
ENSE00003479882
ENSE00001421865
ENSE00001530508

# Compact Exon
ENSE00003537820
ENSE00003476635
ENSE00003605230
ENSE00003666829
ENSE00003523858
ENSE00001530522
ENSE00001530539

cat > MBP_Golli.exon_ids.txt << 'EOF'
exon_id "ENSE00003502747"
exon_id "ENSE00003479882"
exon_id "ENSE00001421865"
exon_id "ENSE00001530508"
EOF

cat > MBP_Compact.exon_ids.txt << 'EOF'
exon_id "ENSE00003537820"
exon_id "ENSE00003476635"
exon_id "ENSE00003605230"
exon_id "ENSE00003666829"
exon_id "ENSE00003523858"
exon_id "ENSE00001530522"
exon_id "ENSE00001530539"
EOF

# Verify that they exist
REFDIR="/storage1/fs1/jin810/Active/References/GRCh38_10X/GRCh38_2024A_MbpSplitWork/refdata-gex-GRCh38-2024-A"
GTF_IN="${REFDIR}/genes/genes.gtf"

echo "Golli hits:"
grep -F -f MBP_Golli.exon_ids.txt "$GTF_IN" | head -n 5
echo "Count Golli hits:"
grep -F -f MBP_Golli.exon_ids.txt "$GTF_IN" | wc -l

echo "Compact hits:"
grep -F -f MBP_Compact.exon_ids.txt "$GTF_IN" | head -n 5
echo "Count Compact hits:"
grep -F -f MBP_Compact.exon_ids.txt "$GTF_IN" | wc -l


# Build minimal exon-only GTFs for your two pseudo-genes (dedupe by genomic coordinates)
REFDIR="/storage1/fs1/jin810/Active/References/GRCh38_10X/GRCh38_2024A_MbpSplitWork/refdata-gex-GRCh38-2024-A"
GTF_IN="${REFDIR}/genes/genes.gtf"

# Golli exon-only mini GTF
grep -F -f MBP_Golli.exon_ids.txt "$GTF_IN" \
  | awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="exon" {
      ex="NA";
      if (match($9, /exon_id "([^"]+)"/, a)) ex=a[1];
      key=$1 FS $4 FS $5 FS $7;   # dedupe by genomic interval (not exon_id)
      if (!(key in seen)) {
        seen[key]=1;
        print $1,"custom","exon",$4,$5,".",$7,".",
          "gene_id \"MBP_Golli\"; transcript_id \"MBP_Golli_tx\"; gene_name \"MBP_Golli\"; transcript_name \"MBP_Golli_tx\"; exon_id \""ex"\";"
      }
    }' > MBP_Golli.exons.min.gtf

# Compact exon-only mini GTF
grep -F -f MBP_Compact.exon_ids.txt "$GTF_IN" \
  | awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="exon" {
      ex="NA";
      if (match($9, /exon_id "([^"]+)"/, a)) ex=a[1];
      key=$1 FS $4 FS $5 FS $7;   # dedupe by genomic interval
      if (!(key in seen)) {
        seen[key]=1;
        print $1,"custom","exon",$4,$5,".",$7,".",
          "gene_id \"MBP_Compact\"; transcript_id \"MBP_Compact_tx\"; gene_name \"MBP_Compact\"; transcript_name \"MBP_Compact_tx\"; exon_id \""ex"\";"
      }
    }' > MBP_Compact.exons.min.gtf

echo "Unique genomic exon intervals retained:"
echo -n "  Golli:   "; wc -l MBP_Golli.exons.min.gtf
echo -n "  Compact: "; wc -l MBP_Compact.exons.min.gtf

echo "Peek Golli min GTF:"
head MBP_Golli.exons.min.gtf
echo "Peek Compact min GTF:"
head MBP_Compact.exons.min.gtf

# Create gene and transcript line per pseudogene
make_gene_tx_simple () {
  local exongtf="$1" gid="$2" txid="$3"
  local chrom strand start end
  chrom=$(awk -F'\t' 'NR==1{print $1}' "$exongtf")
  strand=$(awk -F'\t' 'NR==1{print $7}' "$exongtf")
  start=$(awk -F'\t' 'BEGIN{m=1e18} {if($4<m)m=$4} END{print m}' "$exongtf")
  end=$(awk -F'\t' 'BEGIN{m=0} {if($5>m)m=$5} END{print m}' "$exongtf")

  printf "%s\tcustom\tgene\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; gene_name \"%s\"; gene_type \"protein_coding\";\n" \
    "$chrom" "$start" "$end" "$strand" "$gid" "$gid"
  printf "%s\tcustom\ttranscript\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; transcript_name \"%s\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";\n" \
    "$chrom" "$start" "$end" "$strand" "$gid" "$txid" "$gid" "$txid"
}

make_gene_tx_simple MBP_Golli.exons.min.gtf   MBP_Golli   MBP_Golli_tx   > MBP_Golli.gene_tx.gtf
make_gene_tx_simple MBP_Compact.exons.min.gtf MBP_Compact MBP_Compact_tx > MBP_Compact.gene_tx.gtf

echo "Golli gene/tx:"
cat MBP_Golli.gene_tx.gtf
echo "Compact gene/tx:"
cat MBP_Compact.gene_tx.gtf



# Copy the ref, remove native MBP, append your two pseudo-genes, and write the new genes.gtf
set -euo pipefail

REFDIR="/storage1/fs1/jin810/Active/References/GRCh38_10X/GRCh38_2024A_MbpSplitWork/refdata-gex-GRCh38-2024-A"
GTF_IN="${REFDIR}/genes/genes.gtf"
MBP_GENE_ID="ENSG00000197971"

OUTREF="/storage1/fs1/jin810/Active/References/GRCh38_10X/refdata-gex-GRCh38-2024-A_MbpSplit"
rsync -a --info=progress2 "${REFDIR}/" "${OUTREF}/"

# Remove all records for native MBP gene_id
awk -F'\t' -v gid="$MBP_GENE_ID" '
  $0 ~ /^#/ {print; next}
  $9 !~ ("gene_id \""gid"\"") {print}
' "$GTF_IN" > genes.noMBP.gtf

# Append your custom genes
cat genes.noMBP.gtf \
  MBP_Golli.gene_tx.gtf MBP_Golli.exons.min.gtf \
  MBP_Compact.gene_tx.gtf MBP_Compact.exons.min.gtf \
  > genes.MbpSplit.gtf

# Validate structure
awk -F'\t' '
  $0 ~ /^#/ {next}
  NF!=9 {print "BAD", NR, NF; bad=1}
  END{ if(bad) exit 1; else print "OK: 9 columns for all non-comment lines" }
' genes.MbpSplit.gtf

# Install into OUTREF
cp -f genes.MbpSplit.gtf "${OUTREF}/genes/genes.gtf"
gzip -f "${OUTREF}/genes/genes.gtf"   # so the ref looks like a normal 10x refdata dir again

ls -lh "${OUTREF}/genes/genes.gtf.gz"


export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810:/storage1/fs1/jin810 /home/d.wendy:/home/d.wendy'

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub \
  -G compute-jin810 \
  -q general \
  -n 16 \
  -M 200GB \
  -R 'span[hosts=1] rusage[mem=150GB]' \
  -a 'docker(jinlab/velocytoxcellranger:vs0.17.17x9.0.1)' \
  cellranger mkref \
    --genome=GRCh38_MbpSplit \
    --fasta=/storage1/fs1/jin810/Active/References/GRCh38_10X/refdata-gex-GRCh38-2024-A_MbpSplit/fasta/genome.fa \
    --genes=/storage1/fs1/jin810/Active/References/GRCh38_10X/refdata-gex-GRCh38-2024-A_MbpSplit/genes/genes.gtf.gz


