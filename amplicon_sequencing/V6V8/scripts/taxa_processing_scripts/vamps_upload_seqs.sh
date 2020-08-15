#!/bin/bash

#sort and organize reads from MVCO, OOI and mock communities, format sequence IDs with run included
#cat individual merged files into one file for vamps upload!

#Kristen R. Hunter-Cevera, Marine Biological Laboratory, Spring 2020

folderlist=$(ls | grep -E "MVCO_[1-5,7]")
echo $folderlist

for folder in $folderlist
do
echo $folder
num=$(echo ${folder} | grep -Eo "MVCO\_[0-9]" | grep -Eo "[0-9]{1}")
run="run"${num}
echo "Processing ${folder} ..."

pathname=${folder}"/merged_reads/merged_reads_primers_out_with_python/"
filelist=$(ls ${folder}"/merged_reads/merged_reads_primers_out_with_python/" | grep -E "^primer_match")

    for file in $filelist
    do
        mvco="[0-9]{3}[A-Z]{,1}"
        ooi="(C[0-9]{,2}N)|(BUCKET)"
         #(?=Dock|Syn)
        if [[ $file =~ $mvco ]]; then

            echo "found an MVCO file!"
            echo "File: " ${file} " from run " ${num}
            #echo ${pathname}${file}

            #special handling for each run: 1,2,4 have MV designations, 3,5,7 do not have it and run 5 have extra letters after
            if [[ $num == "1" || $num == "2" || $num == "4" ]]; then
                awk -v rr="$run" 'NR%2 {if ((getline tmp) > 0) {gsub(/\|/," "); sub(/\_/, "_"rr"."); print $0; print tmp}}' ${pathname}${file} >> MVCO_timeseries.fasta
            elif [[ $num == "3" || $num == "7" ]]; then
                awk -v rr="$run" 'NR%2 {if ((getline tmp) > 0) {sub(/>/, ">MV"); gsub(/\|/," "); sub(/\_/, "_"rr"."); print $0; print tmp}}' ${pathname}${file} >> MVCO_timeseries.fasta
            elif [[ $num == "5" ]]; then
                awk -v rr="$run" 'NR%2 {if ((getline tmp) > 0) {sub(/>/, ">MV"); sub(/[A,B]/, ""); gsub(/\|/," "); sub(/\_/, "_"rr"."); print $0; print tmp}}' ${pathname}${file} >> MVCO_timeseries.fasta
            fi

        elif [[ $file =~ $ooi ]]; then
            echo "found an OOI file!"
            echo ${file}
            awk -v rr="$run" 'NR%2 {if ((getline tmp) > 0) {gsub(/\|/," "); sub(/\_/, "_"rr"."); print $0; print tmp}}' ${pathname}${file} >> OOI_AR24B.fasta
        else
            echo "Mock community file?"
            echo ${file}
            awk -v rr="$run" 'NR%2 {if ((getline tmp) > 0) {gsub(/\|/," "); sub(/\_/, "_"rr"."); print $0; print tmp}}' ${pathname}${file} >> mock_syn_community.fasta
        fi
    done
#
done


# filelist=$(ls | grep -E "^primer_match") #not the most robust regex, ^ means match at the start of the line
#
# for file in $filelist
# do
#     echo "adding ${file}"
#     cat ${file} >> mvco_7_upload_to_vamps.fasta
# done
#
# gzip mvco_7_upload_to_vamps.fasta
