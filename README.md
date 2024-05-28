# TFG
_English explanation below_

Código desarrollado para el trabajo de fin de grado "Identificación de efectores dependientes de sistemas de secreción tipo VI en rizobios mediante análisis genómico de adaptadores"

Grado de Biotecnología en la Universidad Politécnica de Madrid

Orden del proceso del trabajo:
1. db_generator.py, para descargar los genomas anotados en NCBI de _Rhizobium_, _Bradyrhizobium_, _Mesorhizobium_ y _Sinorhizobium_, y crear una base de datos para BLAST.
2. blastanalysis.py, para cada una de los secuencias referencia (queries) que se quiera buscar.
3. rpstblastn.py, para cada una de las secuencias de proteínas adaptadoras que se quiera buscar.
4. count_organisms.py, para analizar los resultados obtenidos de blastanalysis.py y rpstblastn.py.
5. find_adjacent_genes.py
6. find_cluster.py

Code developed for the Bachelor's thesis titled "Bioinformatic strategy for searching Type VI Secretion System (T6SS) dependent effectors in rizobia by identifying adaptor proteins"

BSc in Biotechnology in Polytechnic University of Madrid.

Workflow:
1. db_generator.py
