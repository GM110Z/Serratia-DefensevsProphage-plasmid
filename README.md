**ncbi-datasets-batch.sh** : loop for ncbi datasers program to run on batch genomes. Runs as ./ncbi-datasets-batch.sh <file-containing-GCF-ids>

**mash-to-clusters.py**: Converts a mash dist matrix to a file with genome ids and a corresponding cluster number, for dereplication of a dataset.
Runs as python mash-to-clusters.py mash-dist-table.tab <cluster-mash-dist> output.file.txt
*cluster-mash-dist*: Choose the mash dist you want to use to cluster your genomes for dereplication of a dataset.


**koFamScan.sh** : Processes your protein fasta files (GCF_xxx.faa) and then runs https://github.com/takaram/kofam_scan as a loop on each single file. https://github.com/takaram/kofam_scan needs to be pre-installed
