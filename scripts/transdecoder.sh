CPU=32
mkdir ../Transdecoder_noPfam
cd ../Transdecoder_noPfam
module load trinity
module load hmmer
module load cd-hit

ln -s ../Trinity_Assembly/trinity_out_dir/Trinity.fasta

TransDecoder -t Trinity.fasta  --CPU $CPU
 
cd ../Scripts
