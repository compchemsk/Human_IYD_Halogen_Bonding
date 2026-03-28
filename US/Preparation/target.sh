#!/bin/bash

# Input and output files
input_file="dis.dat"  # Replace with your file name, this file will contain the trajectory frame vs COM distance
output_file="target_rows_1.txt"

# Define minima values
minima_values=$(seq 3.8 0.2 8.8)

# Initialize output file
> $output_file

# Loop through each minima value
for target in $minima_values; do
    # Find the row number with minimum distance to the target value
    row=$(awk -v target="$target" 'BEGIN {min_diff=999999; row_num=0} 
        {diff=sqrt(($2-target)^2); 
        if (diff < min_diff) {min_diff=diff; row_num=NR}} 
        END {print row_num}' "$input_file")
    
    # Print row number to output file
    echo "$target $row" >> $output_file
done

echo "Target rows written to $output_file"

# Input and output files
input_file="dis.dat"  # Replace with your file name
output_file="target_rows_2.txt"

# Define minima values
minima_values=$(seq 9.8 1.0 18.8)

# Initialize output file
> $output_file

# Loop through each minima value
for target in $minima_values; do
    # Find the row number with minimum distance to the target value
    row=$(awk -v target="$target" 'BEGIN {min_diff=999999; row_num=0} 
        {diff=sqrt(($2-target)^2); 
        if (diff < min_diff) {min_diff=diff; row_num=NR}} 
        END {print row_num}' "$input_file")
    
    # Print row number to output file
    echo "$target $row" >> $output_file
done

echo "Target rows written to $output_file"
awk '{printf "%s ", $2} END {print ""}' target_rows_1.txt > tgt-rst1.dat
awk '{printf "%s ", $2} END {print ""}' target_rows_2.txt > tgt-rst2.dat
rm target_rows_1.txt target_rows_2.txt
