#!/bin/bash

#$ -q rcc-30d

# Moger-Reischer_Lennon_Bif
#this script uses Slurm. TORQUE/Moab are disabled on Carbonate.
date
cd /N/project/Bifidobacterium
AR=( $(seq 1 25 ) )

for i in "${AR[@]}"
do
    module load gatk/3.8; module load samtools; module unload python; module load python/2.7.16; module load cutadapt; module load perl/5.30.1; module load fastqc; module load bwa; module load picard
	cd Sample_${i}
	cat *R1_001.fastq > Sample_${i}_R1.fastq
	cat *R2_001.fastq > Sample_${i}_R2.fastq
	
	# run qc
	echo "#!/bin/bash" > Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH -J Sample_${i}_qc" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH -p general" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH -o Sample_${i}_qc_stdout.txt" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH -e Sample_${i}_qc_error.txt" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH --mail-type=ALL" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH --mail-user=rzmogerr@indiana.edu" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH --nodes=1" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH --ntasks-per-node=1" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH --cpus-per-task=12" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH --time=2:00:00" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "#SBATCH --mem=20G" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "module load samtools; module unload python; module load python/2.7.16; module load cutadapt; module load perl/5.30.1; module load fastqc; module load bwa; module load java; module load picard" >> Bif_${i}_qc.sh
	echo ""
	echo "cd /N/project/Bifidobacterium/Sample_${i}" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "mkdir fastqc_untrimmed" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "fastqc -o fastqc_untrimmed/ Sample_${i}_R1.fastq" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "fastqc -o fastqc_untrimmed/ Sample_${i}_R2.fastq" >> Bif_${i}_qc.sh
	echo "" >> Bif_${i}_qc.sh
	echo "exit" >> Bif_${i}_qc.sh
	
		
	chmod u+x Bif_${i}_qc.sh

		
	sbatch Bif_${i}_qc.sh
	cd ..
	done