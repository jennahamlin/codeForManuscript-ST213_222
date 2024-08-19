#! /bin/bash 

## move files to same directory

for i in *_prokkaOut; do y=${i%_prokkaOut}; cp $i/$y.faa ~/mySandbox/prokkaWithLPGTags/allST213-222/faaFiles/ ; done

## text processing to generate a file with isolate specific locus tag 
## the associated and Toronto (lpt) locus tags. In some isoaltes the second
## column is gene, so wrangle those files, remove gene column and append the
## files together. 
for i in *.faa; do
 
  output_file="${i%.faa}_locusCleaned.txt"
  awk '/>/ { print $0 }' "$i" | \
  sort -k2 | \
  awk '$2 != "hypothetical" { print $0 }' | \
  awk '!/gene=/' | \
  awk '{ print $1 " " $2 }' | \
  sed 's/\[//g' | \
  awk '$2 ~/^locus_tag=lpt_*/' >> "$output_file"
 
  output_file="${i%.faa}_geneCleaned.txt"
  awk '/>/ { print $0 }' "$i" | \
  sort -k2 | \
  awk '$2 != "hypothetical" { print $0 }' | \
  awk '/gene=/' | \
  awk '{ print $1 " " $3 }' | \
  sed 's/\[//g' >> "$output_file"
done 

for gen in *_locusCleaned.txt; do 

    prefix=${gen%_locusCleaned.txt}
    #echo "$gen" "${prefix}_geneCleaned.txt"
    cat "$gen" "${prefix}_geneCleaned.txt" > "${prefix}_combinedCleaned.txt"
    #paste -d' ' <(cut -d'[' -f1 "${prefix}"_combinedCleaned.txt) <(cut -d'[' -f2 "${prefix}"_combinedCleaned.txt) > "${prefix}"_combinedColumns.txt
    cut -d '>' -f2 "${prefix}_combinedCleaned.txt"  >  "${prefix}_splitCleaned.txt"
    cut -d ']' -f1 "${prefix}_splitCleaned.txt" > "${prefix}_finalCleaned.txt"
    cat "${prefix}_finalCleaned.txt" | sort -k2 > "${prefix}_finalSortToHead.txt"
    cat "${prefix}_finalSortToHead.txt" | awk '{ print $2 ,  $1 }' >> "${prefix}_final.txt"
    sed -i "1i ${prefix}_TAG ${prefix}_ID" "${prefix}_final.txt"
    
done
