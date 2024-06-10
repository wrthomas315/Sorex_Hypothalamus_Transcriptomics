while getopts f:s:g:o: flag
do
    case "${flag}" in
        f) filename=${OPTARG};;
        s) species_num=${OPTARG};;
	g) genelist=${OPTARG};;
	o) output=${OPTARG};;
    esac
done
var=1
header=$(($species_num+$var))

#how to convert the output
# take all of the spaces of LRT file and turn them into linebreaks, then save as a new file .LTR
tr ' ' '\n' < "/gpfs/scratch/withomas/eve/analysis/results_2022_10_03_DIVALV/"$filename > "/gpfs/scratch/withomas/eve/analysis/results_2022_10_03_DIVALV/"$filename".LRTs"
# use awk to convert it out of scientific notation
awk '{printf "%8.6f\n", $1}' "/gpfs/scratch/withomas/eve/analysis/results_2022_10_03_DIVALV/"$filename".LRTs" > "/gpfs/scratch/withomas/eve/analysis/results_2022_10_03_DIVALV/"$filename".LRTs.nosci"
#head -n$header ../14spec/14spec_500exp.txt | awk '{print $1}' > just14_7055
#paste it next to a file with just the ENST gene names
paste "/gpfs/scratch/withomas/eve/data/adult_hypo/"$genelist "/gpfs/scratch/withomas/eve/analysis/results_2022_10_03_DIVALV/"$filename".LRTs.nosci" > "/gpfs/scratch/withomas/eve/analysis/results_2022_10_03_DIVALV/"$output
