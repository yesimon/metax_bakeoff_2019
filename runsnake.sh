#!/bin/bash


MEMAVAIL=$(echo $(grep MemAvail /proc/meminfo | awk '{print $2}') '/' 1024 '/' 1024 | bc)
CPUS=$(cat /proc/cpuinfo | grep -i processor | wc -l)
snakemake --resources mem=$MEMAVAIL io=100 --cores $CPUS $@
