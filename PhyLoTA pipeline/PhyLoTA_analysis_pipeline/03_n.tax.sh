#!/bin/bash

# count # > and append to filename
cd ../02_1TAXON_FASTA/
for f in *.output; do
  cnt=$(grep -c '>' $f)
  printf -v num "%07u" ${cnt% *}
  mv $f $num${f%_*}
done

