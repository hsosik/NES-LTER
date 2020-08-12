#!/bin/bash

#script to run fastQC for quality assessment of reads for each sample
#because the number of reads in each file is not large, can just run this in serial
#Kristen R. Hunter-Cevera, Marine Biological Laboratory, January 2020

#RUN IN FiLE DIRECTORY....

samplelist=$(ls | grep -oP "[A-Z,0-9,\-,a-z]*(?=\_R1|\.R1)")
#
# #first check that directory exists to place fastQC output
if [ ! -d "fastqc_output" ]
then
  echo "making output directory"
  mkdir fastqc_output
fi

for sample in $samplelist
do
  echo "fastqc ${sample}"
  # fastqc ${sample}.R1 -extract -o fastqc_output
  # fastqc ${sample}.R2 -extract -o fastqc_output
  fastqc ${sample}_R1.fastq -extract -o fastqc_output
  fastqc ${sample}_R2.fastq -extract -o fastqc_output
done




#pull out snapshot of information....
#extract just the per base sequence quality and whether or not pass or failed...find first line when bp crosses 30 threshold...
#awk '/Per base sequence quality/{flag=1;next}/>>END_MODULE/{flag=0}flag' fastqc_data.txt

samplelist=$(ls | grep -oP "[A-Z,0-9,\-,\_,a-z]*(\_R1|\.R1)")

#for the forward reads:
echo -e "Sample ID \t BP for mean quality < 30 \t BP for median quality < 30 \t Pass or Fail" >> fastqc_output/out_R1fastqc_summary.txt

for sample in $samplelist
do

array[0]="${sample}_R1"
array[1]=$(awk '/Per base sequence quality/{flag=1;next}/>>END_MODULE/{flag=0}flag' fastqc_output/${sample}_fastqc/fastqc_data.txt | awk '$2 < 30' | cut -f1 | head -1) #when max quality decreases below 30
array[2]=$(awk '/Per base sequence quality/{flag=1;next}/>>END_MODULE/{flag=0}flag' fastqc_output/${sample}_fastqc/fastqc_data.txt | awk '$3 < 30' | cut -f1 | head -1) #when median quality decreases below 30
array[3]=$(grep ">>Per base sequence quality" fastqc_output/${sample}_fastqc/fastqc_data.txt | cut -f2)
echo ${array[@]} >> fastqc_output/out_R1fastqc_summary.txt

done
#awk is amazing!


#and now the reverse reads:
samplelist=$(ls | grep -oP "[A-Z,0-9,\-,\_,a-z]*(\_R2|\.R2)")
echo -e "Sample ID \t BP for mean quality < 30 \t BP for median quality < 30 \t Pass or Fail" >> fastqc_output/out_R2fastqc_summary.txt

for sample in $samplelist
do

array[0]="${sample}_R2"
array[1]=$(awk '/Per base sequence quality/{flag=1;next}/>>END_MODULE/{flag=0}flag' fastqc_output/${sample}_fastqc/fastqc_data.txt | awk '$2 < 30' | cut -f1 | head -1) #when max quality decreases below 30
array[2]=$(awk '/Per base sequence quality/{flag=1;next}/>>END_MODULE/{flag=0}flag' fastqc_output/${sample}_fastqc/fastqc_data.txt | awk '$3 < 30' | cut -f1 | head -1) #when median quality decreases below 30
array[3]=$(grep ">>Per base sequence quality" fastqc_output/${sample}_fastqc/fastqc_data.txt | cut -f2)
echo ${array[@]} >> fastqc_output/out_R2fastqc_summary.txt

done
