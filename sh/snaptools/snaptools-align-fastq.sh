base_dir=/scratch/st-bkkwon-1/DA-analysis/
dev_dir=/scratch/st-bkkwon-1/dev/

mkdir -p ${processed_data_dir}
cd ${processed_data_dir}

output_bam_file=${bam_dir}/${name}.bam
output_snap_file=${snap_dir}/${name}.snap

r1_file=${name}_1.dex.trimmed.fq
r2_file=${name}_2.dex.trimmed.fq
if [  -f ${r1_file} ] && [ -f ${r2_file} ]; then
	echo trimmed files already exists
elif [ -f ${r1_filename} ] && [ -f ${r3_filename} ]; then
	echo paired files exist
	trim_galore --paired ${r1_filename} ${r3_filename} -o ${processed_data_dir} --basename ${name}
	mv ${name}_val_1.fq.gz ${r1_file}
	mv ${name}_val_2.fq.gz ${r2_file}
else
	echo no file exists
fi

r1_file=${processed_data_dir}/${r1_file}
r2_file=${processed_data_dir}/${r2_file}

mkdir -p ${bam_dir}
cd ${bam_dir}

if [ -f $output_bam_file ]; then
	echo $output_bam_file already exists
else
	echo Aligning and making $output_bam_file
	if [ -f "$r1_file" ] && [ -f "$r2_file" ]; then
		echo paired-reads
		snaptools align-paired-end	\
			--input-reference=${base_dir}/ref_genome/${ref_genome}/${ref_genome}.fa	\
			--input-fastq1=$r1_file	\
			--input-fastq2=$r2_file	\
			--output-bam=$output_bam_file	\
			--aligner=bwa	\
			--path-to-aligner=${dev_dir}/bwa-0.7.17/	\
			--read-fastq-command=cat	\
			--min-cov=0	\
			--num-threads=5	\
			--if-sort=True	\
			--tmp-folder=./	\
			--overwrite=TRUE
	else
		echo No reads
	fi
fi