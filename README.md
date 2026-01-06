**ncbi-datasets-batch.sh** : loop for ncbi datasers program to run on batch genomes. Runs as ./ncbi-datasets-batch.sh <file-containing-GCF-ids>

**mash-to-clusters.py**: Converts a mash dist matrix to a file with genome ids and a corresponding cluster number, for dereplication of a dataset.
Runs as python mash-to-clusters.py mash-dist-table.tab <cluster-mash-dist> output.file.txt
*cluster-mash-dist*: Choose the mash dist you want to use to cluster your genomes for dereplication of a dataset.


**koFamScan.sh** : Processes your protein fasta files (GCF_xxx.faa) and then runs https://github.com/takaram/kofam_scan as a loop on each single file. https://github.com/takaram/kofam_scan needs to be pre-installed

**process-metabolic-scorecard.py**: Python script to process KEGG-decoder analysis scorecard and drop all columns that have "0", retaining only metabolic pathways that have a value >0  in at least one row of the scorecard

**CorrelDefenceMetabol.py**: run as python CorrelDefenceMetabol.py <genome-defence-file> <processed-scorecard-from-process-metabolic-scorecard>

**PADLOC-DefFinder-combine.py** :Combines PADLOC and Defensefinder outputs. 
Run as python PADLOC-DefFinder-combine.py PADLOC.tsv DEFENSEFINDER.TSV OUTPUTFILE.tsv

**nuccore-to-assembly.py**: Runs as python nuccore-to-assembly.py padloc-defense-combined.xlsx OutputfileName.xlsx

In the plots folder:
MetabolvsDefence.py: takes the tables produced by CorrelDefenceMetabol.py to produce plots.
