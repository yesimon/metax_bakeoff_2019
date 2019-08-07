#!/bin/bash
# Get all refseq genomes GFF3 files.
HOST="ftp.ncbi.nlm.nih.gov"
# USER="username"
# PASS="password"
GFF_ROOT="$1"
echo "Downloading GFFs to $1"
RCD="/genomes/archive/old_refseq"

mkdir -p "$GFF_ROOT/Bacteria"
mkdir -p "$GFF_ROOT/Fungi"
mkdir -p "$GFF_ROOT/Protozoa"
mkdir -p "$GFF_ROOT/Plasmids/gff"
mkdir -p "$GFF_ROOT/Viruses"
mkdir -p "$GFF_ROOT/Human"
lftp -c "open '$HOST'; lcd $GFF_ROOT/Bacteria; cd $RCD/Bacteria; mirror --verbose --include-glob *.gff"
lftp -c "open '$HOST'; lcd $GFF_ROOT/Fungi; cd $RCD/Fungi; mirror --verbose --include-glob *.gff"
lftp -c "open '$HOST'; lcd $GFF_ROOT/Protozoa; cd $RCD/Protozoa; mirror --verbose --include-glob *.gff"
lftp -c "open '$HOST'; lcd $GFF_ROOT/Plasmids/gff; cd $RCD/Plasmids/gff; mirror --verbose --include-glob *.gff"

pushd "$GFF_ROOT/Viruses"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.gff.tar.gz
tar -zxf all.gff.tar.gz
rm all.gff.tar.gz
popd

HUMAN_GFF="ref_GRCh38.p2_top_level.gff3.gz"
pushd "$GFF_ROOT/Human"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/$HUMAN_GFF"
pigz -d "$HUMAN_GFF"
popd
