**Workshop de Metagenómica: Análisis y Procesamiento de Datos**

Profesor encargado: **Dr. Sergio Guajardo-Leiva**

Profesor Ayudante: **Dr. (c) Valentín Berríos-Farías** 

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

Para esto utilizaremos el software cutadapt [Installation — Cutadapt 0.1 documentation](https://cutadapt.readthedocs.io/en/stable/installation.html)

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

Ahora utilizaremos Metabat2 (https://bitbucket.org/berkeleylab/metabat/src/master/) como herramienta para hacer el binning. Lo primero es convertir nuestro archivo SAM a un formato que Metabat2 pueda leer. Para esto utilizaremos samtools ([http://www.htslib.org](http://www.htslib.org))

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

Ahora que tenemos nuestros bins debemos revisar su grado de completitud y contaminación. Para esto utilizaremos CheckM ([GitHub - chklovski/CheckM2: Assessing the quality of metagenome-derived genome bins using machine learning](https://github.com/chklovski/CheckM2))

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

En el directorio gtdb-tk, estan los archivos de salida, revisar el gtdbtk.bac120.summary.tsv**Workshop de Metagenómica: Análisis y Procesamiento de Datos**

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

Para esto utilizaremos el software cutadapt [Installation — Cutadapt 0.1 documentation](https://cutadapt.readthedocs.io/en/stable/installation.html)

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

Ahora utilizaremos Metabat2 ([Bitbucket](https://bitbucket.org/berkeleylab/metabat/src/master/)) como herramienta para hacer el binning. Lo primero es convertir nuestro archivo SAM a un formato que Metabat2 pueda leer. Para esto utilizaremos samtools ([http://www.htslib.org](http://www.htslib.org))

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

Ahora que tenemos nuestros bins debemos revisar su grado de completitud y contaminación. Para esto utilizaremos CheckM ([GitHub - chklovski/CheckM2: Assessing the quality of metagenome-derived genome bins using machine learning](https://github.com/chklovski/CheckM2))

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
