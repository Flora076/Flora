#!/bin/bash    paired with submit_script.sh
set -e

SRR=SRR26854895
BASE=/omics/groups/OE0436/data/linmq/Datasets/lung_adenocarcinoma

# Create output directory automatically
mkdir -p ${BASE}/${SRR}

# Run Telogator2
python /omics/groups/OE0436/data/linmq/Telogator2-software/telogator2.py \
    -i ${BASE}/fastq/${SRR}_1.fastq.gz \
    -o ${BASE}/${SRR} \
    -p 16 \
    -r ont \
    -n 10 \
    -d 6000 \
    -l 1000 \
    -c 5 \
    --minimap2 /software/minimap2/2.28-GCCcore-14.1.0/bin/minimap2
