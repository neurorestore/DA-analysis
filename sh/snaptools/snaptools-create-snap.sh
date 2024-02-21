base_dir=/scratch/st-bkkwon-1/DA-analysis/
dev_dir=/scratch/st-bkkwon-1/dev/

mkdir -p ${processed_data_dir}
cd ${processed_data_dir}

input_bam_file=${bam_dir}/${name}.bam
output_snap_file=${snap_dir}/${name}.snap

mkdir -p ${bam_dir}
cd ${bam_dir}

if [ -f $output_snap_file ]; then
	echo $output_snap_file already exists
else
	echo making $output_snap_file
	snaptools snap-pre  \
		--input-file=$input_bam_file  \
		--output-snap=$output_snap_file  \
		--genome-name=$ref_genome  \
		--genome-size=${base_dir}/ref_genome/${ref_genome}/${ref_genome}.chrom.sizes  \
		--min-mapq=30  \
		--min-flen=0  \
		--max-flen=1000  \
		--keep-chrm=TRUE  \
		--keep-single=TRUE  \
		--keep-secondary=False  \
		--overwrite=True  \
		--max-num=1000000  \
		--min-cov=100  \
		--verbose=True
fi