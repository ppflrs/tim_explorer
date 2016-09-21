file=$1
gff_file=$2
output_dir=$3

dataset=$(basename $file | cut -d'_' -f1)
#pre_bam_file="$bam_dir"/"$dataset"*.bam
bam_file=$file
tsv_file="$output_dir"/"$dataset"".tsv"
error_log_file="$output_dir"/"$dataset"".error.log"
if [ ! -f "$dataset".tsv ]; then
	echo "$(date)"": ""Parsing ""$dataset"
	samtools view "$bam_file" | python count_reads.py --gff "$gff_file" -i 80 --order n -o "$tsv_file" - 2> "$error_log_file"
else
	echo "$dataset"" already processed."
fi
	
#find . -type l -exec basename {} \; | cut -d'_' -f1
