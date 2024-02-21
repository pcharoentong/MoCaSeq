#!/bin/bash

name=$1
tumor_name=$2
metrics_file_tumor=$3
metrics_file_normal=$4
unfiltered_delly_vcf_file=$5
variants_LOH_file=$6
hmmcopy_segments_file=$7
hmmcopy_log2RR_file=$8
repository_dir=$9
snpeff_dir=${10}
chromosomes=${11}
output_dir=${12}

rm -Rf ./$name

mkdir -p $output_dir
mkdir -p $name/results/Delly/

delly_tab_file=${name}/results/Delly/${name}.breakpoints.tab
filtered_delly_file=${name}/results/Delly/${name}.breakpoints.filtered.tab

format="pdf"
resolution=20000


coverage_tumor=$(ps -ef | awk '/MEAN_COVERAGE/ { getline; print $2 }' $metrics_file_tumor)
coverage_normal=$(ps -ef | awk '/MEAN_COVERAGE/ { getline; print $2 }' $metrics_file_normal)
coverage_sum=$(echo $coverage_normal + $coverage_tumor | bc)
coverage=$(echo "$coverage_sum / 2" | bc)

java -jar ${snpeff_dir}/SnpSift.jar extractFields ${unfiltered_delly_vcf_file} \
	CHROM POS CHR2 POS2 END PE SR PRECISE IMPRECISE \
	GEN[$tumor_name].DR GEN[$tumor_name].DV GEN[$tumor_name].RR GEN[$tumor_name].RV \
	GEN[normal].DR GEN[normal].DV GEN[normal].RR GEN[normal].RV \
	MAPQ CT \
	> ${delly_tab_file}



# writes $filtered_delly_file
Rscript $repository_dir/adapted_for_MoCaSnake/Chromothripsis_Delly_Annotate-Filter.R \
    -n $name -N normal -T $tumor_name -i ${delly_tab_file}  -o T -c $coverage -v 0.2 -d 6000

mkdir -p $name/results/Chromothripsis/

for chr in $chromosomes; do
    rearrangement_count=$(Rscript $repository_dir/Chromothripsis_RearrangementCounter.R -i $filtered_delly_file -c $chr)
    if [ ${rearrangement_count} -ge 4 ]; then
        echo 'Analysing Chromosome '$chr
	mkdir -p $name'/results/Chromothripsis/Chr'$chr

        echo '---- Hallmark: Clustering of breakpoints for Chr'$chr' ----'
        echo -e "$(date) \t timestamp: $(date +%s)"

	Rscript $repository_dir/Chromothripsis_DetectBreakpointClustering.R -i $filtered_delly_file -c $chr -n $name -f $format

        echo '---- Hallmark: Regularity of oscillating copy number states for Chr'$chr' ----'
        echo -e "$(date) \t timestamp: $(date +%s)"

	Rscript $repository_dir/Chromothripsis_SimulateCopyNumberStates.R -i $filtered_delly_file -o human -c $chr -n $name -s 1000 -a 1000 -f $format -v 0

	echo '---- Hallmark: Interspersed loss and retention of heterozygosity for Chr'$chr' ----'
	echo -e "$(date) \t timestamp: $(date +%s)"

        Rscript $repository_dir/Chromothripsis_PlotLOHPattern.R -s $hmmcopy_segments_file -d $hmmcopy_log2RR_file \
		-v $name/results/LOH/$name.VariantsForLOH.txt -o human -c $chr -n $name -f $format

        echo '---- Hallmark: Randomness of DNA fragment joins and segment order for Chr'$chr' ----'
	echo -e "$(date) \t timestamp: $(date +%s)" 

	Rscript $repository_dir/Chromothripsis_DetectRandomJoins.R -i $filtered_delly_file -c $chr -n $name -f $format
      
	echo '---- Hallmark: Ability to walk the derivative chromosome for Chr'$chr' ----' 
	echo -e "$(date) \t timestamp: $(date +%s)" 

	Rscript $repository_dir/Chromothripsis_WalkDerivativeChromosome.R -i $filtered_delly_file -c $chr -n $name -f $format

	echo '---- Visualisation: Copy number profile combined with complex rearrangements for Chr'$chr' ----' 
	echo -e "$(date) \t timestamp: $(date +%s)" 

	Rscript $repository_dir/Chromothripsis_PlotRearrangementGraph.R -i $filtered_delly_file -d $hmmcopy_log2RR_file -c $chr -n $name -f $format

    else
        # check if it did run successfully and print something to the output folder
	checkBreakpoints=$(cut -f 1 $name/results/Delly/$name.breakpoints.filtered.tab | sed -n 2p)
	if [[ $checkBreakspoints == "NA" ]]; then
		echo 'No breakpoints found for chromosome '$chr'.' | tee -a $name/results/Chromothripsis/$name.Chromothripsis.log
	else
		echo 'There are too few rearrangements in chromosome '$chr' ('$rearrangement_count').' | tee -a $name/results/Chromothripsis/$name.Chromothripsis.log
	fi
    fi
done

cp -R $name/results/Chromothripsis/ $output_dir
cp $name/results/Chromothripsis/$name.Chromothripsis.log $output_dir/$name.Chromothripsis.log
