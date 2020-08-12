#!/bin/bash

#script to run merging of forward and reverse reads for each MVCO sample
#Kristen R. Hunter-Cevera, Marine Biological Laboratory, January 2020

#MUST CHANGE DATA PATHS INSIDE SCRIPT!!!!
#RUN FROM WITHIN DIRECTORY OF DEMULTIPLEXED READS
pathname="/home/kristen/Documents/V6V8_analysis/data/MVCO_1_31Oct2016_testrun/"

#IF WANT TO USE MEREN'S SCRIPT TO GENERATE .INI FILES:
# ls *_R1.fastq | tr '\t' '\n' > temp_file
# sed -e 's/_R1.fastq//g' temp_file > temp_file2
# ls *_R1.fastq | tr '\t' '\n' > read1_files
# ls *_R2.fastq | tr '\t' '\n' > read2_files
# paste temp_file2 read1_files read2_files > config_tab_file
# echo -e "sample\tr1\tr2" | cat - config_tab_file > /tmp/out && mv /tmp/out config_tab_file


# ls MV*.R1 | tr ’\t’ ’\n’ | grep -E -o "(MV\d{3}(A|B)?)|(MVM\d{1})" > sample_names.txt
# ls MV*.R1 | tr ’\t’ ’\n’ > read1.txt
# ls MV*.R2 | tr ’\t’ ’\n’ > read2.txt
# paste sample_names.txt read1.txt read2.txt > config_tab_file
# echo -e "sample\tr1\tr2" | cat - config_tab_file > /tmp/out && mv /tmp/out config_tab_file
# rm *.txt


source ~/virtual-envs/illumina-utils-v2.0.2/bin/activate

# mkdir merged_reads #FIX
# iu-gen-configs --e-mail bhamilton@mbl.edu --r1-prefix ˆ.....AAACT[CT]AAA[GT]GAATTGACGG --r2-prefix ˆACGGGCGGTGTGT[AG]C -o merged_reads config_tab_file

if [ ! -d ${pathname}merged_reads/primers_in ]
then
  echo "making output directory"
  mkdir ${pathname}merged_reads/primers_in
fi

filelist=$(ls | grep -oP "[A-Z,0-9,\-,a-z]*(?=\_R1|\.R1)")

for f in $filelist
do
  echo $f
  echo "[general]" > ${f}.ini
  echo "project_name = ${f}" >>  ${f}.ini
  echo "researcher_email = khunter-cevera@mbl.edu" >> ${f}.ini
  echo "input_directory = ${pathname}demultiplexed_reads/" >> ${f}.ini
  echo "output_directory = ${pathname}merged_reads/primers_in/" >> ${f}.ini

  echo "[files]" >> ${f}.ini

  # echo "pair_1 = ${f}_R1.fastq" >> ${f}.ini
  # echo "pair_2 = ${f}_R2.fastq" >> ${f}.ini

  echo "pair_1 = ${f}.R1" >> ${f}.ini
  echo "pair_2 = ${f}.R2" >> ${f}.ini

done

mv *.ini ${pathname}merged_reads/primers_in/ #the .ini files are made in the demultiplex folder  - move them to the priemrs_in folder
cd ${pathname}merged_reads/primers_in/

filelist=$(ls *.ini)

N=6 #refers to number of cores on your machine
(
for sample in $filelist
do
  ((i=i%N)); ((i++==0)) && wait
  echo $i " > " ${sample}
  iu-merge-pairs --max-num-mismatches 3 --enforce-Q30-check --ignore-Ns ${sample}  & # ignore N flag necessary for run 7; curious N in position ~20? for many reads
  #the background amperstand will allow them to go in parallel
done
)

deactivate
