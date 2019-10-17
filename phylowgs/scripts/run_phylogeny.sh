# This script runs tree reconstruction on a single sample.

# input arguments
sample_name=$1
vcf=$2
facets=$3
output_path=$4

########## PRE-DEFINE PATHS ##########

# export gsl lib/ directory to LD_LIBRARY_PATH
LD_LIBRARY_PATH="/usr/include/gsl/lib/"
export LD_LIBRARY_PATH

# full paths for certain scripts
phylowgs_master_path="/scripts/phylowgs-master/"
src_path="/scripts/src/"

######################################

# create output path
output_path=$output_path/$sample_name/
mkdir -p $output_path

# enter output directory and process
cd $output_path

# copy vcf file
vcf_file="$sample_name"".vcf"
cp $vcf $vcf_file

# copy cnv file
cnv_file="$sample_name""_facets.txt"
cp $facets $cnv_file

# parse input files
python "$phylowgs_master_path""/parser/create_phylowgs_inputs.py" $vcf_file -v meta --cnvs $cnv_file

# make phylowgs output directory
mkdir "outdir"/

# move parsed files to phylowgs output directory
mv "cnv_data.txt" "outdir/""$sample_name""_cnv.txt"
mv "ssm_data.txt" "outdir/""$sample_name""_ssm.txt"

# correct phylowgs parser output
cd "outdir/"
python "$src_path""/mainPhyloWGS.py" 2 ../ .

# run PhyloWGS
ssm="$sample_name""_ssm.txt"
cnv="$sample_name""_cnv.txt"

python "$phylowgs_master_path""/evolve.py" $ssm $cnv

# generate JSON results
python "$phylowgs_master_path""/write_results.py" \
    $sample_name \
    "trees.zip" \
    "summ_""$sample_name"".json.gz" \
    "muts_""$sample_name"".json.gz" \
    "mutass_""$sample_name"".zip"

rm "trees.zip"*
