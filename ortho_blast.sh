#!/bin/sh

echo "Starting $0"
echo "ID File=$1"
echo "Blast task type=$2"
echo "Output directory=$3"
echo "Annotation=$4"
echo "Genome=$5\n"

echo "\nCreating query directory"
mkdir query
echo "Creating output directory"
mkdir $3
line_total=$(wc -l $1 | cut -d" " -f1)
minus_1=$(($line_total - 1))
#first_line=$(tail -n$minus_1 $1 | cut -d"," -f1 | head -n1)
col=$(tail -n$minus_1 $1 | cut -d"," -f1)
echo "$col"
for entry in $col
do
	grep "$entry" $4 | head -n1 | tee "./query/$entry.bed"
	echo "\nEcho query directory contents"
	ls ./query
	echo ""
	bedpath="./query/$entry.bed" 
	fastapath="./query/$entry.fasta"
	outpath="./$3/$entry-blast.csv"
	bedargs="-fi $5 -bed $bedpath -fo $fastapath"
	echo "Creating query fasta:"
	bedtools getfasta $bedargs
	cat $fastapath
	start_date=$(date)
	echo "\nStarting blast for: $entry at $start_date"
	blastn -db core_nt -query $fastapath -task $2 -dust no -num_threads 4 -outfmt "7 delim=, sacc qstart qend evalue score length staxid ssciname qcovs"| tee $outpath
	end_date=$(date)
	echo "\nFinished blast for: $entry at $end_date\n"
done

