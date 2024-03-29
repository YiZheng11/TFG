#!/bin/bash
# Javier Ortiz de Artiñano	2023	Trabajo de Fin de Grado
# This script analyzes a dataset of Bradyrhizobium assemblies.
# Requires a dataset (found in data/), generated by DBGenerator.bash:
#	- index.txt -> index of RefSeq Bradyrhizobium assemblies to be analyzed. Every row contains the following fields separated by tabs:
#	AssemblyAccn	AssemblyName	Organism	AssemblyStatus	BioprojectAccn	BioSampleAccn	ContigN50	ScaffoldN50	FtpPath_RefSeq
#	- assemblies/	-> directory containing the indexed and downloaded assemblies. The .fna files are also decompressed.
# Analysis (found in results/):
#	- BlastP of a query file in FASTA format containing one or more protein sequences. 
#		· Its results are filtered and a human readable tsv file is produced, with the following fields:
#		AssemblyAccn	hitsGeneA	hitsGeneB	hitsGeneC	...


#################################################
# FUNCTIONS
#################################################

#ayuda() prints a help message
ayuda() {
   echo 'BLASTP Analyzer'
   echo 'javiortizdea 2023'
   echo -e '\nusage: BLASTPAnalyzer.bash -o <output_folder> [-q <query.fa>] [-I <min_identity>] [-C <min_coverage>] [-E <min_evalue>] [-h]\n'
   echo -e '-o output_folder: proyect folder results will be stored. It must contain a data/ folder inside.'
   echo -e "-q query.fa	: fasta/multifasta file containing protein query sequences for BlastP."
   echo -e "-I min_identity	: minimum identity allowed in BlastP results. The default value is 30%."
   echo -e "-C min_coverage	: minimum query coverage allowed in BlastP results. The default value is 50%."
   echo -e "-E min_evalue	: maximum expected value (E-value) allowed in BlastP results. The default value is 1e-5."
   echo -e "-h		: print this help and exit."
}

#join_rec() and recurrentJoin() are two recursive join functions.
#recurrentJoin() uses join_rec().
#When recurrentJoin() is called, it joins as many (sorted) files as it is provided.
join_rec() {
    if [ $# -eq 1 ]; then
        join - "$1"
    else
        f=$1; shift
        join - "$f" | join_rec "$@"
    fi
}

recurrentJoin() {
	if [ $# -le 2 ]; then
		join "$@"
	else
		f1=$1; f2=$2; shift 2
		join "$f1" "$f2" | join_rec "$@"
	fi
}

#BlastPAnalysis() performs the BlastP analysis of the query sequences over the BLAST database.
#Also, the corresponding subject AssemblyAccn is printed next to every hit. This is done to facilitate
#the following data filtering and formating.
#It also filters the BLAST result with the set thresholds -> blastResultFiltered.tsv
BlastPAnalysis() {
	if [ ! -d $projectName/results ]; then
		mkdir $projectName/results
	fi
	rm -f $projectName/results/*
	blastp -query $query -db $projectName/data/brady_T6SS_context_DB -outfmt '6 qseqid sseq pident qcovs evalue sseqid' -max_target_seqs '2000' > $projectName/results/blastResult.tsv 2>> $projectName/log.txt
	#blastp -query $query -db $projectName/data/brady_DB.fa -outfmt '6 qseqid sseq pident qcovs evalue sseqid' -max_target_seqs '2000' > $projectName/results/blastResult.tsv 2>> $projectName/log.txt
	echo -e "\tDone!"
	echo -e "\tFiltering results..."
	gawk -v min_iden=$minimumIdentity -v min_cov=$minimumCoverage -v max_evalue=$maximumEValue '$3 >= min_iden && $4 >= min_cov && $5 <= max_evalue {print}' $projectName/results/blastResult.tsv > $projectName/results/blastResultFiltered.tsv 2>> $projectName/log.txt
	echo -e "\tDone!" 
}

#hitsCounter() filters the results of the BLAST analysis and gives them a human-readable format.
#To do so:
#	1. Counts the number of hits for every combination of query sequence and assembly -> totalHits.tsv
#	2. Generates a hitsX.tsv file with the number of hits for every assembly (including the ones without significative hits) for every X query.
#	3. With recurrentJoin() the hitsX.tsv files are joined in a human-readable table -> readableTotalHits.tsv
#	4. Organism names are added -> readableTotalHits.csv
#	5. BioSample classifications are added -> bioSamplesAndHits.csv (IT REQUIRES A classes.tsv FILE IN THE DATA FOLDER)
hitsCounter() {
	echo -e "\tGenerating human-readable table..."
	touch $projectName/results/totalHits.tsv
	cut -f1,6 $projectName/results/blastResultFiltered.tsv | sed 's/!.*$//g' | sort | uniq | while read match; do
		nHits=`cat $projectName/results/blastResultFiltered.tsv | cut -f1,6 | grep -c "$match"` 2>> $projectName/log.txt
		echo "$match" | sed -e "s/$/\t$nHits/" | tr " " "\t" >> $projectName/results/totalHits.tsv 2>> $projectName/log.txt
	done
	coolHeader='#Accn'
	cat $projectName/results/totalHits.tsv | cut -f1 | sort | uniq | ( while read geneName; do
		join -1 1 -2 1 -a1 -a2 -o '2.1,1.2' -e 0 <(grep "^$geneName" $projectName/results/totalHits.tsv | cut -f2,3 | sort) <(tail +2 $projectName/data/index.txt | cut -f1 | sort) > $projectName/results/hits$geneName.tsv 2>> $projectName/log.txt
		coolHeader=`echo -e "$coolHeader\t$geneName"`
	done; echo -e "$coolHeader\t" > $projectName/results/readableTotalHits.tsv )
	recurrentJoin $projectName/results/hits*.tsv | tr " " "\t" >> $projectName/results/readableTotalHits.tsv 2>> $projectName/log.txt
	echo -e "\tDone! Results can be found in $projectName/results/readableTotalHits.tsv"
}

#################################################
# DEFAULT PARAMETERS
#################################################

doBlastPAnalysis=false
minimumIdentity=30
minimumCoverage=50
maximumEValue='1e-5'

while getopts ho:q:I:C:E: option; do
	case $option in
		h) ayuda; exit 0;;
		o) projectName=$OPTARG;;
		q) doBlastPAnalysis=true; query=$OPTARG;;
		I) minimumIdentity=$OPTARG;;
		C) minimumCoverage=$OPTARG;;
		E) maximumEValue=$OPTARG;;
		*) echo "error: unknown option" 1>&2; ayuda 1>&2; exit 1;;
	esac
done

#To make sure that at least the output folder name has been provided.
if [ "$#" -lt 1 -o "$projectName" = false ]; then
	echo -e "error: too few arguments" 1>&2
	ayuda 1>&2
	exit 1 
fi

#################################################
# BLASTP ANALYSIS
#################################################

if [ -d $projectName/results ]; then
	echo -e "\tWARNING: the /results folder will be overwritten. We recommend changing its name \n\ty/n"
	read respuesta
	case $respuesta in
		y|Y) echo -e "\tThe specified directory will be used and its files will be overwritten.";;
		n|N) echo -e "\tStopping the process..."; exit 0;;
		*) echo -e "\tInvalid answer.\n\tIf you want to continue with the indicated directory, enter \"y\" o \"Y\"." 1>&2 ; echo -e "\tStopping the process.." 1>&2; exit 1;;
	esac
fi

if [ $doBlastPAnalysis = true ]; then
	BlastPAnalysis
	hitsCounter
fi
