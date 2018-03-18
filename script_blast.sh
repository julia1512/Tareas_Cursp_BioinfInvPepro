#!/bin/bash
# Este script baja la secuencia de una proteína prión de humano y corre un blast para ver si se parece a alguna proteína en el pez zebra. 


#Bajar la secuencia de la proteína
wget http://www.uniprot.org/uniprot/P04156.fasta 

#Bajar y descomprimir la base de datos de pez zebra de NCBI
curl -O ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.1.protein.faa.gz
gunzip zebrafish.1.protein.faa.gz

#Declarar variable con la ruta absoluta a donde se va a utilizar para montar un volumen.
mi_ruta=/home/julia/Documents/LIBB/8o_semestre/Bioinf_Repro/Tareas/Unidad5
#Montar el volumen, preparar la base de datos y correr el blast. Los resultados se guardan al archivo results.txt
docker run -v $mi_ruta:/data/ biocontainers/blast makeblastdb -in zebrafish.1.protein.faa -dbtype prot
docker run -v $mi_ruta:/data/ biocontainers/blast blastp -query P04156.fasta -db zebrafish.1.protein.faa -out results.txt
