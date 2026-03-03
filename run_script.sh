#!/bin/bash    paired with submit_script.sh
set -e

# Run Telogator2
python /omics/groups/OE0436/data/linmq/Telogator2-software/telogator2.py \
    -i /omics/groups/OE0436/data/linmq/Datasets/osteosarcoma/fastq/SRR26854882_1.fastq.gz\
    -o /omics/groups/OE0436/data/linmq/Datasets/osteosarcoma/SRR26854882\
    -p 16 \
    -r ont \
    -n 10 \
    -d 6000 \
    -l 1000 \
    -c 5 \
    --minimap2 /software/minimap2/2.28-GCCcore-14.1.0/bin/minimap2
