# Start interactive session
bsub -Is -G compute-jin810-t3 -q subscription -sla jin810_t3 -n 8 \
    -R 'gpuhost rusage[mem=64GB]' \
    -gpu 'num=1:j_exclusive=yes' \
    -a 'docker(biocontainers/samtools:v1.9-4-deb_cv1)' \
/bin/bash


# 1) extract Jun region (+/- 5kb padding)
samtools view -b \
  C7_Injured_Sciatic_3_C7_Injured_Sciatic_3_98f996_7ac08f_ae73b0.bam \
  chr4:94932271-94945459 \
  > C7_Injured_Sciatic_3_Jun.bam

# 2) index the subset BAM
samtools index C7_Injured_Sciatic_3_Jun.bam


## Code to calculate # of full length reads at each truncated 3'UTR site
cd "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Jun_BAM_files"

OUT="Jun_site_counts_5primeFiltered.tsv"
echo -e "sample\tcondition\tsite\tcount" > "$OUT"


cd "/Users/wendydong/Documents/WUSM PhD/Long Read RNA Seq/Jun_BAM_files"

OUT="Jun_site_counts_5primeFiltered.tsv"

SITE=94937934
SITENAME="site1"
W=10

FIVE=94940459
W5=10

for bam in *Jun.bam; do
  cond="NA"
  [[ "$bam" == *"C0"* ]] && cond="C0"
  [[ "$bam" == *"C3"* ]] && cond="C3"
  [[ "$bam" == *"C7"* ]] && cond="C7"

  count=$(samtools view "$bam" | awk -v site="$SITE" -v w="$W" -v five="$FIVE" -v w5="$W5" '
    int($2/16)%2==1 {
      pos=$4
      c=$6
      ref=0
      # sum reference-consuming ops M,D,N,=,X
      while (match(c, /[0-9]+[MDN=X]/)) {
        tok=substr(c, RSTART, RLENGTH)
        len=substr(tok, 1, length(tok)-1)
        ref += len
        c = substr(c, RSTART + RLENGTH)
      }
      end = pos + ref - 1

      if (end >= five - w5 && end <= five + w5 && pos >= site - w && pos <= site + w) {
        cnt++
      }
    }
    END { print cnt+0 }
  ')

  echo -e "${bam}\t${cond}\t${SITENAME}\t${count}" >> "$OUT"
done

echo "Appended $SITENAME to $OUT"
