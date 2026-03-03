 submit_script.sh       for submitting jobs on the cluster, usually used to download SRR                                                                           
#!/bin/bash

#BSUB -J SRR26854882   
#BSUB -W 168:00 
#BSUB -n 10   
#BSUB -R "rusage[mem=128G]" 
#BSUB -q verylong 
#BSUB -o /omics/groups/OE0436/data/linmq/Datasets/osteosarcoma/this-log-file.log 
#BSUB -e /omics/groups/OE0436/data/linmq/Datasets/osteosarcoma/this-err-file.err 

module load Python/3.12.4-GCCcore-14.1.0
module load minimap2
pip install biopython matplotlib numpy pysam pywavelets scipy pandas
pip install --upgrade pip

cd /omics/groups/OE0436/data/linmq/Datasets/iPSC
 
./run_script.sh

