set -euo pipefail

REFDIR="/storage1/fs1/jin810/Active/References/GRCm39_10X/refdata-gex-GRCm39-2024-A"
GTF_IN="${REFDIR}/genes/genes.gtf"
FA_IN="${REFDIR}/fasta/genome.fa"

MBP_GENE_ID="ENSMUSG00000041607"

# Optional: copy the 10x refdata directory (only needed if you want a refdata-like output)
OUTREF="refdata-gex-GRCm39-2024-A_MbpSplit"
rsync -a --info=progress2 "${REFDIR}/" "${OUTREF}/"

# Define exon IDs
cat > MbpGolli.patterns.txt << 'EOF'
exon_id "ENSMUSE00001256352"
exon_id "ENSMUSE00001241294"
EOF

cat > MbpCompact.patterns.txt << 'EOF'
exon_id "ENSMUSE00001293533"
exon_id "ENSMUSE00001227500"
EOF

# Build minimal exon GTFs directly from original genes.gtf
grep -F -f MbpGolli.patterns.txt "$GTF_IN" \
  | awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="exon" {
      ex="NA";
      if ($0 ~ /ENSMUSE00001256352/) ex="ENSMUSE00001256352";
      else if ($0 ~ /ENSMUSE00001241294/) ex="ENSMUSE00001241294";
      key=$1 FS $4 FS $5 FS $7 FS ex;
      if (!(key in seen)) {
        seen[key]=1;
        print $1,"custom","exon",$4,$5,".",$7,".",
              "gene_id \"Mbp_Golli\"; transcript_id \"Mbp_Golli_tx\"; gene_name \"Mbp_Golli\"; transcript_name \"Mbp_Golli_tx\"; exon_id \""ex"\";"
      }
    }' > Mbp_Golli.exons.min.gtf

grep -F -f MbpCompact.patterns.txt "$GTF_IN" \
  | awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="exon" {
      ex="NA";
      if ($0 ~ /ENSMUSE00001293533/) ex="ENSMUSE00001293533";
      else if ($0 ~ /ENSMUSE00001227500/) ex="ENSMUSE00001227500";
      key=$1 FS $4 FS $5 FS $7 FS ex;
      if (!(key in seen)) {
        seen[key]=1;
        print $1,"custom","exon",$4,$5,".",$7,".",
              "gene_id \"Mbp_Compact\"; transcript_id \"Mbp_Compact_tx\"; gene_name \"Mbp_Compact\"; transcript_name \"Mbp_Compact_tx\"; exon_id \""ex"\";"
      }
    }' > Mbp_Compact.exons.min.gtf

# Validate exon min files
awk -F'\t' '$0 !~ /^#/ && NF!=9{print "BAD",FNR,NF; exit 1}' Mbp_Golli.exons.min.gtf
awk -F'\t' '$0 !~ /^#/ && NF!=9{print "BAD",FNR,NF; exit 1}' Mbp_Compact.exons.min.gtf

# Build gene+transcript records
make_gene_tx_simple () {
  local exongtf="$1" gid="$2" gname="$3" txid="$4" txname="$5"
  local chrom strand start end
  chrom=$(awk -F'\t' 'NR==1{print $1}' "$exongtf")
  strand=$(awk -F'\t' 'NR==1{print $7}' "$exongtf")
  start=$(awk -F'\t' 'BEGIN{m=1e18} {if($4<m)m=$4} END{print m}' "$exongtf")
  end=$(awk -F'\t' 'BEGIN{m=0} {if($5>m)m=$5} END{print m}' "$exongtf")

  printf "%s\tcustom\tgene\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; gene_name \"%s\"; gene_type \"protein_coding\";\n" \
    "$chrom" "$start" "$end" "$strand" "$gid" "$gname"
  printf "%s\tcustom\ttranscript\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; transcript_name \"%s\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";\n" \
    "$chrom" "$start" "$end" "$strand" "$gid" "$txid" "$gname" "$txname"
}

make_gene_tx_simple Mbp_Golli.exons.min.gtf   Mbp_Golli   Mbp_Golli   Mbp_Golli_tx   Mbp_Golli_tx   > Mbp_Golli.gene_tx.gtf
make_gene_tx_simple Mbp_Compact.exons.min.gtf Mbp_Compact Mbp_Compact Mbp_Compact_tx Mbp_Compact_tx > Mbp_Compact.gene_tx.gtf

# Remove original Mbp and append split genes
awk -F'\t' -v gid="$MBP_GENE_ID" '
  $0 ~ /^#/ {print; next}
  $9 !~ ("gene_id \""gid"\"") {print}
' "$GTF_IN" > genes.noMbp.gtf

cat genes.noMbp.gtf \
    Mbp_Golli.gene_tx.gtf Mbp_Golli.exons.min.gtf \
    Mbp_Compact.gene_tx.gtf Mbp_Compact.exons.min.gtf \
  > genes.MbpSplit.gtf

# Structural validation (ignore comments)
awk -F'\t' '
  $0 ~ /^#/ {next}
  NF!=9 {print "BAD", NR, NF; bad=1}
  END {if (bad) exit 1; else print "OK: all non-comment lines have 9 columns"}
' genes.MbpSplit.gtf

# If you made a refdata copy, replace its GTF with the new one
cp -f genes.MbpSplit.gtf "${OUTREF}/genes/genes.gtf"



