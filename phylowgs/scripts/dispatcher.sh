# This script dispatches a phylogeny reconstruction run.

usage() { echo -en "dispatcher. This script dispatches a phylogeny reconstruction run.\nExample: ./dispatcher.sh \n\t-s SAMPLE1\n\t-v SAMPLE.vcf\n\t-f SAMPLE1_facets_output.txt\n\t-o /output/\n\nFlags:\n" && grep " .)\ #" $0; exit 0; } 
[ $# -eq 0 ] && usage
while getopts ":hs:v:f:o:" arg; do
    case $arg in
        s) #Sample Name
            sample_name=${OPTARG}
            ;;
        v) #Sample VCF
            vcf=${OPTARG}
            ;;
        f) #Sample Facets output file
            facets=${OPTARG}
            ;;
        o) #Output Directory
            output_path=${OPTARG}
            ;;
        h | *) # Display help.
        usage
            exit 0
            ;;
    esac
done

# use full paths for the vcf, facets, and output_path
#sample_name="MMRF_2692_1_BM"
#vcf="$(pwd)""/../data/""$sample_name"/"$sample_name"".vcf"
#facets="$(pwd)""/../data/""$sample_name"/"$sample_name""_facets_output.txt"
#output_path="$(pwd)""/../results/"

/scripts/run_phylogeny.sh \
    "$sample_name" \
    "$vcf" \
    "$facets" \
    "$output_path"
