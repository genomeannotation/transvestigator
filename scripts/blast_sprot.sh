CPU=32
DB=/data0/opt/BlastSoftware/DB/uniprot_sprot.fasta
mkdir ../BlastSprot
cd ../BlastSprot
ln -s ../Trinity_Assembly/trinity_out_dir/Trinity.fasta . 

module load blast 
module load trinity

blastx -query Trinity.fasta -db $DB -out blastx.outfmt6 -evalue 1e-20 -num_threads $CPU -max_target_seqs 1 -outfmt 6
analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 Trinity.fasta $DB


cd ../Scripts
