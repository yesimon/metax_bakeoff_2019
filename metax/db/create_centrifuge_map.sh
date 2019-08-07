#!/bin/bash
set -eo pipefail
export LC_ALL=C
NCBI_DIR="$1"

DOMAINS=(archaea bacteria viral)
for DOMAIN in ${DOMAINS[@]}; do
    echo $DOMAIN

    join -t$'\t' <(sort $NCBI_DIR/${DOMAIN}_selection.txt) <(sort $NCBI_DIR/assembly_summary_refseq.txt)  | cut -f1,6 > ${DOMAIN}_taxids.txt

    for i in $(cat $NCBI_DIR/${DOMAIN}_selection.txt) ; do
        for j in $(fasta_formatter -i $NCBI_DIR/refseq/$DOMAIN/$i/*_genomic.fna -t | cut -f1 | cut -f1 -d' ' ); do
            echo -e "$i\t$j";
        done ;
    done | sort > ${DOMAIN}_accessions.tsv

    join -t$'\t' ${DOMAIN}_taxids.txt ${DOMAIN}_accessions.tsv | awk 'BEGIN {FS=OFS="\t"} {print $3, $2}' > all-$DOMAIN.map
done


# centrifuge-compress.pl archaea taxonomy  -map all- archaea.map -o archaea-compressed.tmp -t 40 -maxG 50000000 2>&1 | tee centrifuge-compress-archaea.log
