#!/bin/bash

#$ -q rcc-30d

# Moger-Reischer_Lennon_Bif
#this script uses Slurm. TORQUE/Moab are disabled on Carbonate.
date
cd /N/project/Bifidobacterium

AR=( $(seq 1 2 ) )

for i in "${AR[@]}"
do
	module load breseq
	cd Sample_${i}
	
	#do breseq
	echo "#!/bin/bash" > breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH -J Sample_${i}_breseq" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH -p general" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH -o Sample_${i}_breseq_stdout.txt" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH -e Sample_${i}_breseq_error.txt" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH --mail-type=ALL" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH --mail-user=rzmogerr@indiana.edu" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH --nodes=1" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH --ntasks-per-node=1" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH --cpus-per-task=12" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH --time=20:00:00" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "#SBATCH --mem=100G" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "module load breseq" >> breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	echo "cd /N/project/Bifidobacterium" > breseq_${i}_.sh
	echo "" >> breseq_${i}_.sh
	#echo "breseq -j 8 -p -o 3B_${i}_breseq -r /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.gb Sample_${i}/Sample_${i}_R1_trimmed.fastq" >> 3B_breseq_${i}
	echo "breseq -j 8 -p -o breseq_${i} -r /N/project/Bifidobacterium/BB12_Jensen_2021.gb Sample_${i}/Sample_${i}_R1_trimmed.fastq Sample_${i}/Sample_${i}_R2_trimmed.fastq" >> breseq_${i}_.sh
	
	chmod u+x breseq_${i}_.sh
	
	sbatch breseq_${i}_.sh
	done