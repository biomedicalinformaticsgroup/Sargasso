When separating mixed-species RNA-seq data according to true species of origin using Sargaaso, 
it is difficult to avoid making a small number of incorrect species assignments, 
and it is important to ensure that these reads do not lead to spurious conclusions.  
In order to avoid such misinterpretation, one approach was discussed in the https://www.nature.com/articles/s41596-018-0029-2, 
which can be used to estimate the amount of reads incorrectly assigned during the Sargasso separation process. 
This script is an implementation of the above approach. The details of the approach is describe in the above paper Procedure session step 30-34.


As an example, suppose we have samples contains cells from 3 species human, mouse and rat. In order to 
Here is a brief break-down of the steps:
-


Usage: 
Rscript calculate_misassignment.R  'mouse' '01_HMR,02_HMR,03_HMR,04_HMR' 'human,mouse,rat' '09_Pc,10_Pc,11_Pc,12_Pc,05_EC,06_EC,07_EC,08_EC' 'rat;rat;rat;rat;human;human;human;human' 0 '/home/xinhe/Projects/bbb_pbaxter/misassignment/gene_length.csv' '/home/xinhe/Projects/bbb_pbaxter/results/sargasso/filtered_reads/overall_filtering_summary.txt' '/home/xinhe/Projects/bbb_pbaxter/results/read_counts/' '/home/xinhe/Projects/bbb_pbaxter/results/misassignment'

Rscript calculate_misassignment.R <species_of_interest>
    <target_samples>
	<target_species>
	<reference_samples>
	<reference_species>
	<paired>
	<gene_lengths_file>
	<overall_filtering_summary_file>
	<read_count_dir>
	<result_dir>
    
Script parameters are described below:

- <species_of_interest>: the species of which misassignment % is calculated 

- <target_samples>:  comma separated mix-species sample names. The script will estimate the mis-assignment % of the species_of_interest in these samples.

- <target_species>:  semi-comma separated species name, matching the target_samples. If the target_samples contains more than one species, use comma to separate them. For example, 'rat;mouse' means the target_sample_1 is rat only sample, and target_sample_2 is mouse only sample. 'rat,human;mouse' means the the target_sample_1 is a co-culture sample contains rat and human cells while target_sample_2 is mono-culture contains mouse cell only.

- <reference_samples>: comma separated sample names. These samples usually cells from one or more species from the target_species, but do not contains the cells from species_of_interest. These samples are used to work out how likely a gene is to be assigned incorrectly from other species to the species_of_interest.

- <reference_species> :  similar to target_species. semi-comma separated species name, matching the reference_samples. If the reference_samples contains more than one species, use comma to separate them. 

- <paired>: whether the reference_samples are paired with the target_samples. Can be '0'(no pair) or '1'(paired)

- <gene_lengths_file>: csv file contains a gene column and a max_transcript_length column. 
  gene               max_transcript_length
  <chr>                              <dbl>
1 ENSMUSG00000000001                  3262
2 ENSMUSG00000000003                   902
...

- <overall_filtering_summary_file>: Sargasso output file, which can be found in filtered_reads folder under Sargasso result folder.

- <read_count_dir>: featureCount output folder. For each of the samples above, a read count file is expected under the name 'sample.special.counts'. 

- <result_dir>: output folder.


Output:
- gene_misassignment_percentage.csv
This files contains the pre gene misassignment percentage for the species_of_interest.
 
gene | misprec | misprec.human | misprec.rat
-- | -- | -- | --
ENSMUSG00000000056 | 0 | 0 | 0
ENSMUSG00000000058 | 0.0006917 | 0.0006917 | 0
ENSMUSG00000000078 | 0.02673777 | 0 | 0.02673777

- reference_sample_composition.csv
This file contains the Sargasso assigned reads percentage of each species in each sample

- target_sample_RNA_ratio.csv
This file contains, for each target_sample, the approximate ratio of RNA of target_species to species_of_interest in the sample. This corresponds to Procedure session step 32, or 'd' in step 33 in the paper.