#!/bin/bash

# Define the ranges for pt and cent
pt_ranges=( "(1 2)" "(2 3)" "(3 4)" "(4 5)" "(5 6)" "(6 8)" "(8 10)" "(10 15)" "(15 20)" )
cent_ranges=( "(0 10)" "(10 30)" "(30 50)" "(50 90)" )
#pt_ranges=( "(1 2)")
#cent_ranges=( "(0 10)")
#pt_ranges=( "(1 2)" "(2 3)" "(3 4)" "(4 5)" "(5 6)" "(6 8)" "(8 10)" "(10 15)" "(15 20)" "(20 30)" "(30 50)" "(50 100)" "(100 1000)" )
#cent_ranges=( "(0 10)" "(10 30)" "(30 50)" "(50 90)" )

# Path to the template file
template_file="TMVAClassification_template.C"

# Loop over each pt range
for pt_range in "${pt_ranges[@]}"; do
    # Remove parentheses from pt_range
    ptmin=$(echo "$pt_range" | tr -d '()' | awk '{print $1}')
    ptmax=$(echo "$pt_range" | tr -d '()' | awk '{print $2}')
    
    # Loop over each cent range
    for cent_range in "${cent_ranges[@]}"; do
        # Remove parentheses from cent_range
        centmin=$(echo "$cent_range" | tr -d '()' | awk '{print $1}')
        centmax=$(echo "$cent_range" | tr -d '()' | awk '{print $2}')
        
        # Generate the directory name  based on the ranges
        dir_name="pt_${ptmin}_${ptmax}_cent_${centmin}_${centmax}"

        # Create the new TMVAClassification.C file by modifying the template
        cp "$template_file" TMVAClassification.C

        # Replace placeholders in the copied template, without parentheses
        sed -i "s/pt_min/$ptmin/g" TMVAClassification.C
        sed -i "s/pt_max/$ptmax/g" TMVAClassification.C
        sed -i "s/cent_min/$centmin/g" TMVAClassification.C
        sed -i "s/cent_max/$centmax/g" TMVAClassification.C
        sed -i "s/directory_name/$dir_name/g" TMVAClassification.C
        
        # Run the script with ROOT
        root -b -q TMVAClassification.C

        rm -r BDT_training_output/$dir_name/./dataset
        mv dataset BDT_training_output/$dir_name/.
        echo "Finished processing $dir_name"
    done
done

