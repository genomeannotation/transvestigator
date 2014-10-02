CPU=32
module load trinity


mkdir ../RSEM_Abundance
cd ../RSEM_Abundance
cat ../RawData/*R1*fastq > all.R1.fastq
cat ../RawData/*R2*fastq > all.R2.fastq
ln -s ../Trinity_Assembly/trinity_out_dir/Trinity.fasta . 
align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left all.R1.fastq --right all.R2.fastq --est_method RSEM --aln_method bowtie --thread_count $CPU --trinity_mode --prep_reference
rm all*R*fastq 
cd ../Scripts
