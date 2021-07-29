#!/bin/bash

#$ -q rcc-30d

# Moger-Reischer_Lennon_Bif
#this script uses Slurm. TORQUE/Moab are disabled on Carbonate.
date
cd /N/project/Bifidobacterium/breseq
AR=( $(seq 1 25 ) )

for i in "${AR[@]}"
do
    module load breseq
	cd .
	
	# clean and trim
	echo "#!/bin/bash" > Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH -J Sample_${i}_breseq" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH -p general" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH -o Sample_${i}_breseq_stdout.txt" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH -e Sample_${i}_breseq_error.txt" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH --mail-type=ALL" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH --mail-user=rzmogerr@indiana.edu" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH --nodes=1" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH --ntasks-per-node=1" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH --cpus-per-task=12" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH --time=28:00:00" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "#SBATCH --mem=120G" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "module load breseq" >> Bif_${i}_breseq.sh
	echo ""
	echo "cd /N/project/Bifidobacterium/breseq" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "breseq -j 8 -p -o breseq_${i} -r /N/project/Bifidobacterium/BB12_Jensen_2021.gb ../Sample_${i}/Sample_${i}_R1_trimmed.fastq ../Sample_${i}/Sample_${i}_R2_trimmed.fastq" >> Bif_${i}_breseq.sh
	echo "" >> Bif_${i}_breseq.sh
	echo "exit" >> Bif_${i}_breseq.sh
	
		
	chmod u+x Bif_${i}_breseq.sh
		
	sbatch Bif_${i}_breseq.sh
	#cd ..
	done