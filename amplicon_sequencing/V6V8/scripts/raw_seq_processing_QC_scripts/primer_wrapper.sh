#!/bin/bash

#wrapper script to run primer removal:
#run from just outside the directory of where the merged reads are stored
#Kristen R. Hunter-Cevera, Marine Biological Laboratory, January 2020

pathname=$(pwd)
inpath=$pathname/"merged_reads_primers_in/"
samplelist=$(ls $inpath | grep "MERGED")
outpath="/merged_reads_primers_out_with_python/"

if [ ! -d ${pathname}${outpath} ]
then
  echo "making output directory"
  mkdir ${pathname}${outpath}
fi

N=6 #refers to number of cores on your machine
(
for sample in $samplelist
do
  ((i=i%N)); ((i++==0)) && wait
  echo " ${sample}"
  python ~/Documents/V6V8_analysis/scripts/raw_seq_processing_QC_scripts/remove_forward_reverse_primers.py ${inpath}${sample} ${pathname}${outpath} ${sample} &

done
)
#
