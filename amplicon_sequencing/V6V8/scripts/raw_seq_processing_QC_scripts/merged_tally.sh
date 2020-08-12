#!/bin/bash

#script to examine number of sequences before and after merging and after primer trimming if doing this after merging:
#Kristen R. Hunter-Cevera, Marine Biological Laboratory, January 2020


#We would like:
#Sequence into merge
#Sequences recognized by merge program
#Sequences that passed prefix
#Sequences failed due to max num mismatch
#Sequences failed due to Q30 check
#Total merged sequences

#FOR PRIMERS OUT DURING MERGE:
#
# path="/home/kristen/Documents/V6V8_analysis/data/MVCO_7_13Jan2020_run/demultiplexed_reads/merged_reads/"
# outstats="MVCO_7_13Jan2020_run_merged_stats.txt"
#
# #filelist=$(ls MV*_MERGED-MAX-MISMATCH-3)
#
# filelist=$(ls ${path}*STATS) # | grep -E -o "MV(\d{3}|M2)")
# echo -e "Sample ID \t Total reads \t Passes prefix \t Failed max mismatch \t Failed Q30 \t Total merged" >> ${outstats}
#
# for sample in $filelist
#
# do
#
# echo ${sample}
# #P flag is for Perl regex, look ahead insertion: (?=..)
# array[0]=$(echo ${sample} | grep -oP "[A-Z, 0-9, \-]*(?=\_STATS)")
# array[1]=$(grep "Number of pairs analyzed" ${sample} | cut -f2)
# array[2]=$(grep "Passed prefix total" ${sample} | cut -f2)
# array[3]=$(grep "Merge discarded due to max num mismatches" ${sample} | cut -f2)
# array[4]=$(grep "Merge discarded due to Q30" ${sample} | cut -f2)
# array[5]=$(grep "Merged total" ${sample} | cut -f2)
#
# echo ${array[@]} >> $outstats
# done



#***************************************************************************************************************
#FOR PRIMERS OUT AFTER MERGE:

#run from top data folder:
folderlist=$(ls | grep -E "MVCO_[1-5,7]")

for folder in $folderlist
do

  echo "Processing ${folder} ..."

  num=$(echo ${folder} | grep -Eo "MVCO\_[0-9]" | grep -Eo "[0-9]{1}")
  pathname1=${folder}"/merged_reads/merged_reads_primers_out_with_python/"
  pathname2=${folder}"/merged_reads/merged_reads_primers_in/"

  outstats=${folder}"/merged_reads/primers_python_out_tally_stats.txt"

  filelist=$(ls ${pathname2}| grep -E "STATS" | grep -oP "((MV){0,1}[0-9]{3}[A-Z]{0,1}[Milli|Ster]{0,5})|(SynCommunity[A-D]{1})|(Dock[A,B]|((MV){0,1}M[1,2])|(MVBUCKET[P,V])|(C[0-9]{1,2}N[0-9]{0,2}-{0,1}[0-9]{0,2}[V,P]))")

  echo -e "Sample_ID \t Run_Number \t Total_reads \t Failed_max_mismatch \t Failed_Q30 \t Total_merged \t Primer_fail \t Primer_matched" > ${outstats}

  for sample in $filelist
  do

    echo ${sample}
    #P flag is for Perl regex, look ahead insertion: (?=..)
    array[0]=$(echo ${sample})
    array[1]=$(echo ${num})
    array[2]=$(grep "Number of pairs analyzed" ${pathname2}${sample}"_STATS" | cut -f2)
    array[3]=$(grep "Merge discarded due to max num mismatches" ${pathname2}${sample}"_STATS" | cut -f2)
    array[4]=$(grep "Merge discarded due to Q30" ${pathname2}${sample}"_STATS" | cut -f2)
    array[5]=$(grep "Merged total" ${pathname2}${sample}"_STATS" | cut -f2)

    array[6]=$(grep -c ">" ${pathname1}noprimer_match_${sample}_MERGED)
    array[7]=$(grep -c ">" ${pathname1}primer_match_${sample}_MERGED)

    echo ${array[@]} | tr ' ' '\t' >> $outstats

  done
done
