#Need Trinity installed and in path

#Some Variables:

ASSEM_MEMORY=120G #MAKE SURE TO PUT G
ASSEM_THREADS=32
ASSEM_MINKMER=1
SS_LIB_TYPE=RF
JACCARD_CLIP=--jaccard_clip  #FILL IN WITH "--jaccard_clip" if you want to use, blank will skip
NORMALIZE_READS=--normalize_reads
##ALSO CAN MODIFY TO ADD BUTTERFLY OPTIONS
mkdir ../Trinity_Assembly
cd ../Trinity_Assembly

module load trinity
ulimit

Trinity --seqType fq --JM $ASSEM_MEMORY --left ../RawData/*.R1*fastq --right ../RawData/*.R2*fastq --CPU $ASSEM_THREADS $JACCARD_CLIP $NORMALIZE_READS --normalize_max_read_cov 50 --min_kmer_cov $ASSEM_MINKMER --bflyCalculateCPU --bflyHeapSpaceInit 4G

ln -s trinity_out_dir/Trinity.fasta . 

if [ ! -f trinity_out_dir/Trinity.fasta ] ; then
   echo "Problem with Assembly!";
   cd ..
   exit
fi
cd ..
echo "Assembly ran successfully!"
echo "proceed with Assembly Filtering"

