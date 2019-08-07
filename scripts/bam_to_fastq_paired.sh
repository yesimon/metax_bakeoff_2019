#!/bin/bash

trap "trap - SIGTERM && pkill -P $$ 2> /dev/null || true" SIGINT SIGTERM EXIT
mkfifo one.pipe
mkfifo two.pipe

if [[ $# -gt 3 ]]; then
    samtools merge -O SAM -n - "${@:3}" |
        $PICARD -Djava.io.tmpdir=$TMPDIR SamToFastq INPUT=/dev/stdin CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X \
                FASTQ=one.pipe SECOND_END_FASTQ=two.pipe TMP_DIR=$TMPDIR &
else
    $PICARD -Djava.io.tmpdir=$TMPDIR SamToFastq INPUT="$3" CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X \
            FASTQ=one.pipe SECOND_END_FASTQ=two.pipe TMP_DIR=$TMPDIR &
fi

PSAM=$!

# sort -k1,1 -s -T 

cat one.pipe | lbzip2 > "$1" &
P1=$!
cat two.pipe | lbzip2 > "$2" &
P2=$!
if wait $PSAM && wait $P1 && wait $P2; then
    exit 0
else
    exit 1
fi
