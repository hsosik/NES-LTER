#!/bin/bash

#Demultiplex files that still have different barcodes
#Rename files according to barcode and index mismatches
#Kristen R. Hunter-Cevera, Marine Biological Laboratory, January 2020


#If necessary, gunzip them:
#gunzip *.fastq.gz

#find the files:

filelist=$(ls *R?*.fastq)
regex="([A-Z]{6})[A-Z0-9\_]*.fastq"
declare -a barcode_array #store found barcodes in an array
echo -e "ID\tindex\tbarcode\tforward reads\treverse reads\tno barcode found" > demultiplex_logfile.txt #initialize log, the -e flag allows \t to be interpretted as a tab
tempfile=$(echo $filelist | cut -f 1 -d " ")
seqline_identifier=$(head -1 $tempfile |  cut -c 1-6) #just to get the start of each line, which has the sequencer id on it

i=0 #counter
for f in $filelist
do
    if [[ $f =~ $regex ]]
    then
        index="${BASH_REMATCH[1]}"
        echo $index
        index_array[$i]=$index
        #echo "${name}.jpg"    # concatenate strings
        #name="${name}.jpg"    # same thing stored in a variable
    else
        echo "$f doesn't seem to have an index in the filename? "
    fi
    i=$((i+1)) #increase counter
    #echo $i
done

index_array_unq=( $(echo "${index_array[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') )
#tr turns the spaces into lines, which the command sort is expecting, and then turns it back to spaces

#Demultiplex:

#First, add reverse complement onto searched indexes:
awk '{print $3}' sample_key.txt | rev | tr "ATCG" "TAGC" | paste sample_key.txt - > sample_key_wrev.txt
#reverse complement in bash courtesy of gist!

for ind in "${index_array_unq[@]}" #for each index file, go through and separate out the barcodes
do
  echo "demultiplexing ${ind} index"

  temparray=( $(echo ${filelist[@]} | tr ' ' '\n' | grep "$ind") ) #find the two files, R1 & R2

  #Combine the files to read on one line:

  cat ${temparray[0]} | pr -a -4 -J -t > temp_forward_reads
  cat ${temparray[1]} | pr -a -4 -J -t > temp_reverse_reads
  #the command \textbf{pr} stands for print and reshuffles the lines according to -a (across), -J (join) 4 lines and -t (exclude header info)).
  paste temp_forward_reads temp_reverse_reads > temp_total_reads
  #8 tab-delimited fields now: defline, read1\_forward, defline ('$+$'), qualscores, defline, read1\_reverse, defline ('$+$'), qualscores.

  #now, depending on how many indices are in here, split...
  barcode_array=( $(grep $ind sample_key_wrev.txt | awk '{print $2}') )
  echo "Found these barcodes: ${barcode_array[@]} for index ${ind}"

  for bc in "${barcode_array[@]}"
  do
      #incorporate variable outside of matching:
      regex="N:0:[0-9]{1,3}\t$bc" #check to make sure regex pattern matches!
      echo "searching for barcode: ${bc}"

      # commands contingent on finding at least 1 barcode match:
      #awk -v pat="$regex" -v bcmatch="${pathname}bc_match" -v nobcmatch="${pathname}no_bc_match" '{if ($0 ~ pat){print $0 >> bcmatch} else{print $0 >> nobcmatch}}' ${pathname}temp_total_reads
      awk -v pat="$regex" '{if ($0 ~ pat){print $0 >> "bc_match"} else{print $0 >> "no_bc_match"}}' temp_total_reads

      # #rename temp_total_reads to be the no_matches...
      cat no_bc_match > temp_total_reads
      #
      # #find the MV sample match:
      ID=$(grep $ind sample_key_wrev.txt | grep $bc | awk '{print $1}')
      # echo "matches ${ID}"
      #In tabular format, convert back with:
      cat bc_match | cut -f 1-4 | tr '\t' '\n' > ${ID}_R1.fastq
      cat bc_match | cut -f 5-8 | tr '\t' '\n' > ${ID}_R2.fastq

      length_r1=$(grep -c "$seqline_identifier" ${ID}_R1.fastq) # gives number of sequences
      length_r2=$(grep -c "$seqline_identifier" ${ID}_R2.fastq)
      # length_r1=$(wc -l ${ID}_R1.fastq | cut -f 1 -d " ") #log it!
      # length_r2=$(wc -l ${ID}_R2.fastq | cut -f 1 -d " ")
      echo -e "${ID}\t${ind}\t${bc}\t${length_r1}\t${length_r2}\t---" >> demultiplex_logfile.txt

      rm bc_match no_bc_match #clear these temporary files

  done

  mv temp_total_reads no_barcode_match_${ind}
  length_no_match=$(wc -l no_barcode_match_${ind} | cut -f 1 -d " ") #lines are okay here -not put back into fastq format
  echo -e "---\t${ind}\t---\t---\t---\t${length_no_match}" >> demultiplex_logfile.txt

  #This puts all files in current working directory, but probably should have them in a separate folder...

done


#reverse files are in the correct orientation!
