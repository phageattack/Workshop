**Workshop de Metagenómica: Análisis y Procesamiento de Datos**

**Dr. Sergio Guajardo-Leiva**

Lo primero sera descargar los datos desde onedrive

[WORKSHOP_Microbial_Frontiers](https://udetalca-my.sharepoint.com/:f:/g/personal/sergio_guajardo_utalca_cl/EihvlvVzs-ZBj4NAqAYOAQMBqkR8XWB9awSlAGniU4kKdQ?e=ukW450)

**PARTE 1: ANÁLISIS CENTRADO EN EL GENOMA**

**1. Revisar calidad de la secuencias**

Para esto usaremos FastQC que pueden descargar aquí: [Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```
mkdir workshop
cd workshop
mkdir fastqc # primero creamos el directorio de output

fastqc -o ./fastqc -t 2 *.fq

# -t numero de hilos

ls -lh
```

Dentro del directorio fastqc encontraremos una serie de archivos. Los mas relevantes son los de extensión HTML . Desde los reportes se desprende que no existen adaptadores en nuestras secuencias y que la calidad promedio por secuencia es bastante buena.

" Q 30" : Probabilidad de una determinación de base incorrecta de 1 en 1000 y una precisión de la determinación de base del 99.9%.

En consecuencia, solo iremos eliminando bases con puntajes Q < 30 y no realizaremos ningún "hard-clipping".

Para esto utilizaremos el software cutadapt [Installation &mdash; Cutadapt 0.1 documentation](https://cutadapt.readthedocs.io/en/stable/installation.html)

**2. Recortar las secuencias**

```
mkdir trimed 
ls -lh
cutadapt -j 2 -q 30,30 -Q 30,30 -m 50 --max-n 0 -o ./trimed/Sample1_trimed1.fq -p ./trimed/Sample1_trimed2.fq Sample1_1.fq Sample1_2.fq 
# -j número de procesadores 
# -q calidad por base en par 1 
# -Q calidad por base en par 2
# -m largo minímo luego del triming
# --max-n máximo número de bases sin identificar 
```

Ahora revisaremos el proceso de limpieza de nuestras secuencias

```
#Contar el numero de secuencias antes y despues del recorte
grep -c "^@" Sample1_1.fq
wc -l Sample1_1.fq | awk '{print $1/4}' 
grep -c "^@" trimed/Sample1_trimmed1.fq 
#Que diferencias hay?
#Contar el numero de bases totales antes y despues del recorte
grep -A 1 "^@" Sample1_1.fq| grep -v "^@" | grep -v "^+" | tr -d '\n' | wc -c
grep -A 1 "^@" Sample1_trimmed1.fq| grep -v "^@" | grep -v "^+" | tr -d '\n' | wc -c
#Que direncias hay?
```

**3. Ensamblar las secuencias**

Una vez que las reads han sido recortadas por calidad realizamos en el ensamble con megahit [GitHub - voutcn/megahit: Ultra-fast and memory-efficient (meta-)genome assembler](https://github.com/voutcn/megahit)

```
megahit -m 4e9 -t 2 --min-contig-len 1500 -1 ./trimed/Sample1_trimed1.fq -2 ./trimed/Sample1_trimed2.fq -o ./Sample1_assembly
# -m memoria para el ensamble 2Gb
# -t numero de procesadores
```

**4. Binning metagenomico**

Antes del binning debemos alinear las reads contra los contigs para conocer su profundidad de secuenciación. Para esto ultilizamos Bowtie2 [Bowtie 2: fast and sensitive read alignment](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

```
# Primero hacemos una base de datos con los contigs
bowtie2-build  ./Sample1_assembly/final.contigs.fa  ./Sample1_assembly/final.contigs

#ahora mapeamos las reads contra la base de datos

bowtie2 --verbose -p 2 --sensitive --end-to-end --no-unal -x ./Sample1_assembly/final.contigs -1 ./trimed/Sample1_trimed1.fq -2 ./trimed/Sample1_trimed2.fq -S ./profundidad/final.contigs.sam
sed -i '1,5d' ./profundidad/final.contigs.sam #removemos las primeras 5 lineas que escribe bowtie
```

Ahora utilizaremos Metabat2 (https://bitbucket.org/berkeleylab/metabat/src/master/) como herramienta para hacer el binning. Lo primero es convertir nuestro archivo SAM a un formato que Metabat2 pueda leer. Para esto utilizaremos samtools (http://www.htslib.org)

```
samtools view -b ./profundidad/final.contigs.sam | samtools sort > ./profundidad/final.contigs.bam
```

Ahora utilizaremos un script incluido en Metabat2 para general un archivo de texto con la profundidad de cada contig.

```
 jgi_summarize_bam_contig_depths --outputDepth ./profundidad/depth.txt ./profundidad/final.contigs.bam 
```

Ahora procedemos al binning

```
metabat2 -i ./Sample1_assembly/final.contigs.fa -a ./profundidad/depth.txt -m 1500  -o bins/bin
```

Ahora que tenemos nuestros bins debemos revisar su grado de completitud y contaminación. Para esto utilizaremos CheckM (https://github.com/chklovski/CheckM2)

```
checkm lineage_wf -t 2 -x fa ./bins/ ./bins/checkm
```

Finalmente cuantificaremos la covertura de cada bien en nuestra muestra, para esto utilizaremos CoverM [GitHub - wwood/CoverM: Read coverage calculator for metagenomics](https://github.com/wwood/CoverM)

```
nano rename.sh
```

Luego pegaremos el siguiente codigo

```
#!/bin/bash

# Define the directory containing your .fa files
input_dir="./"

# Loop through each .fa file in the directory
for file in "$input_dir"/*.fa; do
    # Check if the file exists
    if [ -e "$file" ]; then
        # Get the filename without the extension
        filename=$(basename -- "$file")
        filename_no_ext="${filename%.*}"

        # Create a temporary file to store the modified sequences
        tmp_file="${filename_no_ext}_tmp.fa"

        # Process each line in the input file
        while IFS= read -r line; do
            # Check if the line starts with '>'
            if [[ $line == '>'* ]]; then
                # Remove '>' and replace spaces with underscores
                modified_line=">${filename_no_ext}_${line:1}"
            else
                modified_line="$line"
            fi
            # Append the modified line to the temporary file
            echo "$modified_line" >> "$tmp_file"
        done < "$file"

        # Replace the original file with the temporary file
        mv "$tmp_file" "$file"

        echo "Modified headers in $file"
    else
        echo "File not found: $file"
    fi
done

echo "Done"
```

Ahora con los contigs ya renombrados, utilizamos CoverM

```
coverm genome -t 2 -m rpkm -d ./bins -x fa -1 ./trimed/Sample1_trimed1.fq  -2 ./trimed/Sample1_trimed2.fq

-m relative_abundance
-p bwa-mem 
```

Finalmente asignaremos una clasificaciones taxonómica "objetiva" a los bins basándonos en la Taxonomía de la Base de Datos de Genomas (GTDB). Para esto utilizaremos GTDB-Tk [GitHub - Ecogenomics/GTDBTk: GTDB-Tk: a toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes.](https://github.com/Ecogenomics/GTDBTk)

Lamentablemente, este paso es solo teórico, ya que GTDB-Tk consume muchos recursos, especialmente requiere mucha memoria RAM. Este software debe utilizarse en un servidor o estación de trabajo con al menos 128-256 GB de RAM.

```
gtdbtk classify_wf --genome_dir ./bins --out_dir ./gtdb-tk --skip_ani_screen --extension fa
```

En el directorio gtdb-tk, estan los archivos de salida, revisar el gtdbtk.bac120.summary.tsv

**PARTE 2: ANÁLISIS CENTRADO EN GENES**

En esta segunda parte utilizaremos varios archivos generados en la primera parte del presente curso. En específico nos serán útiles las secuencias recortadas y también los
contigs ensamblados.

Lo primero que haremos es crear una base de datos a partir de secuencias aminoacídicas de genes claves en el metabolismo microbiano (que sacamos de aquí https://bridges.monash.edu/collections/_/5230745 ) contra las cuales alinearemos nuestras reads. Esta base de datos corresponde a la siguiente publicación https://doi.org/10.1038/s41564-023-01322-0.

En este proceso usaremos Diamond [GitHub - bbuchfink/diamond: Accelerated BLAST compatible local sequence aligner.](https://github.com/bbuchfink/diamond)

```
diamond makedb --in metabolic_concatenated_genes.fa -d metabolic_concatenated_genes
```

Con la base de datos ya creada en el formato que Diamond requiere, procedemos a mapear las reads recortadas.

```

mkdir blastx

diamond blastx --fast --threads 2 --outfmt 6 qtitle stitle pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp --max-target-seqs 1 --max-hsps 1 --id 50 --query-gencode 11 -d ./metabolic_concatenated_genes.dmnd -q ./trimed/Sample1_trimmed1.fq -o ./blastx/Sample1.blastx.txt

# --ultra-sensitive si tuviesemos mucho tiempo
```

Ahora haremos lo mismo pero utilizando una aproximacion diferente. En vez de utilizar directamente las secuencias recortadas utilizaremos las secuencias aminoacídicas codificadas en nuestros contigs. Para esto utilizaremos primero prodigal [GitHub - hyattpd/Prodigal: Prodigal Gene Prediction Software](https://github.com/hyattpd/Prodigal) que nos permitira predecir secuencias codificantes de proteinas.

```
prodigal -i Sample1_assembly/final.contigs.fa -a Sample1_proteins.faa -d Sample1_genes.fna -p meta -n
```

Luego alinearemos estas secuencias contra nuestra base de datos utilizando diamond nuevamente.

```
mkdir blastp 
diamond blastp --fast --outfmt 6 qtitle stitle pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp --max-target-seqs 1 --max-hsps 1 --id 50 -d ./metabolic_concatenated_genes.dmnd -q ./Sample1_proteins.faa -o ./blastp/Sample1_proteins.blastp.txt
# Como este proceso es mas rapido, podemos probar la opción --ultra-sensitive
```

Estos dos analisis nos entregaran una archivo de texto en formato tabulado el cual deberemos inspeccionar y resumir para convertirlo en un formato que nos haga sentido.

Para esto utilizaremos una serie de scripts cortos en bash que nos haran la tarea mas facil.

```
cd blastx
#primero filtraremos todos los alineamientos de < 40 aa
awk -F '\t' '{if ($4 >= 40) print $0}' Sample1.blastx.txtt  > Sample1.blastx_40aa.txt

# nos quedaremos con las filas 2,3 y agregaremos una columna llena de 1
cut -f2,3 Sample1.blastx_40aa.txt | awk -F'\t' 'BEGIN {OFS=FS} {print $0, 1}' > Sample1.blastx_cut.txt

# Ahora filtraremos por cada gen deacuerdo a un porcentaje de identidad distinto
awk -F"\t" '$2 >= 60 || !($1 ~ /(_CoxL_|_MmoA_|_AmoA_|_NxrA_|_RbcL_|_NuoF_|_FeFe_|_NiFe_)/)' Sample1.blastx_cut.txt > Sample1.blastx_cut_f1-60%.txt

awk -F'\t' '$1 !~ /.*_HbsT_.*|^\*/ || $2 >= 75' Sample1.blastx_cut_f1-60%.txt > Sample1.blastx_cut_f2-75%.txt

awk -F'\t' '$1 !~ /.*_PsaA_.*|^\*/ || $2 >= 80' Sample1.blastx_cut_f2-75%.txt > Sample1.blastx_cut_f3-80%.txt

awk -F"\t" '$2 >= 70 || !($1 ~ /(_PsbA_|_IsoA_|_AtpA_|_YgfK_|_ARO_)/)' Sample1.blastx_cut_f3-80%.txt > Sample1.blastx_cut_f4-70%.txt

awk -F '\t' '{if ($2 >= 50) print $0}' Sample1.blastx_cut_f4-70%.txt > Sample1.blastx_cut_f5-50%.txt

#vamos a ver como fuimos filtrando secuencias
wc-l *.txt 

#Sumaremos todas las secuencias que corresponde a la misma proteina

awk 'BEGIN{FS="\t"; OFS="\t"} {sum[$1] += $3} END{for (key in sum) print key, sum[key]}' Sample1.blastx_cut_f5-50%.txt > Sample1.blastx_cut_filtered_sum.txt

```

Ahora aplicaremos una serie de filtros a los resultados de blastp para luego cuantificar los genes que pasen los filtros.

```
cd blastp
#primero filtraremos todos los alineamientos de con una cobertura < 80%
awk -F '\t' '{if ($11 >= 80) print $0}' Sample1_proteins.blastp.txt  > Sample1_proteins.blastp_80qcov.txt

# nos quedaremos con las filas 1-3 
cut -f1-3 Sample1_proteins.blastp_80qcov.txt > Sample1_proteins.blastp_80qcov_cut.txt

# Ahora filtraremos por cada gen deacuerdo a un porcentaje de identidad distinto
awk -F"\t" '$3 >= 60 || !($2 ~ /(_CoxL_|_MmoA_|_AmoA_|_NxrA_|_RbcL_|_NuoF_|_FeFe_|_NiFe_|_AtpA)_|_PsbA_/)' Sample1_proteins.blastp_80qcov_cut.txt > Sample1_proteins.blastp_cut_f1-60%.txt

awk -F'\t' '$2 !~ /.*_HbsT_.*|^\*/ || $3 >= 75' Sample1_proteins.blastp_cut_f1-60%.txt > Sample1_proteins.blastp_cut_f2-75%.txt

awk -F'\t' '$2 !~ /.*_PsaA_.*|^\*/ || $3 >= 80' Sample1_proteins.blastp_cut_f2-75%.txt > Sample1_proteins.blastp_cut_f3-80%.txt

awk -F"\t" '$3 >= 70 || !($2 ~ /(_IsoA_|_YgfK_|_ARO_)/)' Sample1_proteins.blastp_cut_f3-80%.txt > Sample1_proteins.blastp_cut_f4-70%.txt

awk -F'\t' '$2 !~ /.*_RdhA_.*|^\*/ || $3 >= 45' Sample1_proteins.blastp_cut_f4-70%.txt > Sample1_proteins.blastp_cut_f5-45%.txt

awk -F'\t' '$2 !~ /.*_Cyc2_.*|^\*/ || $3 >= 35' Sample1_proteins.blastp_cut_f5-45%.txt > Sample1_proteins.blastp_cut_f6-35%.txt

awk -F'\t' '$2 !~ /.*_RHO_.*|^\*/ || $3 >= 30' Sample1_proteins.blastp_cut_f6-35%.txt > Sample1_proteins.blastp_cut_f7-30%.txt

#vamos a ver como fuimos filtrando secuencias
wc-l *.txt 

#finalmente resumimos la tabla y extraemos la lista de genes a cuantificar

cut -f1,2 Sample1_proteins.blastp_cut_f7-30%.txt | awk -F'\t' '{split($1, a, " "); $1 = a[1]}1' OFS='\t' > final_blastp.txt
cut -f1 final_blastp.txt > headers.txt
awk 'BEGIN{FS=" "}{if(FNR==NR){headers[$1]; next} } /^>/{header=substr($1,2)} header in headers' headers.txt ../Sample1_genes.fna > final_blastp_genes.fna
```

Ahora para cuantificar utilizaremos bowtie2

```
# Primero hacemos una base de datos con los genes
bowtie2-build  ./blastp/final_blastp_genes.fna  ./blastp/final_blastp_genes

#ahora mapeamos las reads contra la base de datos

bowtie2 --verbose -p 2 --sensitive --end-to-end --no-unal -x ./blastp/final_blastp_genes -1 ./trimed/Sample1_trimed1.fq -2 ./trimed/Sample1_trimed2.fq -S ./blastp/final_blastp_genes.sam
sed -i '1,5d' ./blastp/final_blastp_genes.sam #removemos las primeras 5 lineas que escribe bowtie

```

Ahora usaremos bbmap para generar un texto plano que nos permita cuantificar cada gen

```
pileup.sh in=./blastp/final_blastp_genes.sam  out=./blastp/final_blastp_genes_coverage.txt rpkm=./blastp/final_blastp_genes_rpkm.txt
```

Copiaremos estos archivos a windows para poder verlos en excel

final_blastp_genes_coverage.txt

final_blastp_genes_RPKM.txt

final_blastp.txt

Sample1.blastx_cut_filtered_sum.txt

**Ahora en R**

Ahora importaremos este el ultimo archivo de blastx en R para hacer la suma de las diferentes variantes del la misma proteina y graficar.

```

setwd("~/workshop_SGL/blastx_R/")


length <-read.table(file = "./metabolic_concatenated_genes.length.tsv", header = T)
Sample_1 <-read.table(file = "Sample1.blastx_cut_filtered_sum.txt", header = F)

join1 <- merge(length,Sample_1, by = 1, all.y =TRUE)
write.table(join1, file = "genes_per_sample.csv", sep = ",", row.names = F, col.names = T)
# en bash 
awk -F, '{split($1, arr, "_"); $1 = arr[1]; for (i=2; i<length(arr); i++) $1 = $1 "_" arr[i]; print $1 "," $2 "," $3}' genes_per_sample.csv > genes_per_sample_to_TPM.csv
# y calculamos los TPM en excel
 # volvemos a R
 TPM <- read.table("Sample1_TPM.csv", sep = ",", header = T)

library(dplyr)
library(tidyverse)
df_summed <- TPM %>%
  group_by(gene) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("gene")
write.csv(df_summed, file = "TPM_sumada.csv")



log_TPM <- log1p(df_summed)

library(ComplexHeatmap)


matrix1 <- as.matrix(data.frame(log_TPM))
ComplexHeatmap::pheatmap(matrix1, name = "Log(TPM)")


colours = colorRampPalette(c("lightyellow","orange1","red","purple4","black"),space="rgb")(100)
ComplexHeatmap::pheatmap(matrix1, name = "Log(TPM)", color = colours)



setwd("~/workshop_SGL/blastp_R/")

en bash vamos a preparar el archivo de RPKM para utilizarlo
sed -i '1,4d' final_blastp_genes_rpkm.txt
cut -f 1,2,5 final_blastp_genes_rpkm.txt > final_blastp_genes_rpkm_cut.txt

rpkm <-read.table(file = "final_blastp_genes_rpkm_cut.txt", header = T)
Sample_1 <-read.table(file = "final_blastp.txt", header = F)

join1 <- merge(rpkm, Sample_1, by = 1, all.y =TRUE)
write.table(join1, file = "genes_per_sample.csv", sep = ",", row.names = F, col.names = T)
# en bash 
awk -F, 'BEGIN {OFS=","} {sub(/_[^_]*$/, "", $2); print}' genes_per_sample.csv > genes_per_sample_to_TPM.csv
# y calculamos los TPM en excel
 # volvemos a R
 TPM <- read.table("Sample1_TPM.csv", sep = ",", header = T)

library(dplyr)
library(tidyverse)
df_summed <- TPM %>%
  group_by(gene) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("gene")
write.csv(df_summed, file = "TPM_sumada.csv")



log_TPM <- log1p(df_summed)

library(ComplexHeatmap)


matrix1 <- as.matrix(data.frame(log_TPM))
ComplexHeatmap::pheatmap(matrix1, name = "Log(TPM)")


colours = colorRampPalette(c("lightyellow","orange1","red","purple4","black"),space="rgb")(100)
ComplexHeatmap::pheatmap(matrix1, name = "Log(TPM)", color = colours)

```
