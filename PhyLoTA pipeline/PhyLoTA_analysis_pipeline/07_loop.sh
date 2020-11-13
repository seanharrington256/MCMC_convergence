for f in *.fasta; do
  Rscript --vanilla 08_trim_ends_fasta.R $f ${f%.fasta}.fas .25;
done

