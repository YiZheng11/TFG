# TFG
Código desarrollado para el trabajo de fin de grado "Identificación de efectores dependientes de sistemas de secreción tipo VI en rizobios mediante análisis genómico de adaptadores"

- main
1. dbgenerator.py: Same as DBGenerator.bash but in Python

- Linux_old
1. DBGenerator.bash es el código de Javier, sin embargo, le he hecho unos retoques ya que el código original de Javier no me iba correctamente
2. DBGenerator_Yi.bash: Mismo código que DBGenerator.bash pero para más géneros de la familia de Rhizobiaceae además de Bradyrhizobium
3. decompress.bash: Usar antes de DB1_MultiFastaGenerator.ipynb para descomprimir los archivos dentro del directorio de cada assembly.
4. DB1_MultiFastaGenerator.ipynb
5. T6SSgenes_MultiFastaGenerator.ipynb: Igual que el código de DB1_MultiFastaGenerator.ipynb solo que para generar un MultiFasta de genes concretos de un genoma, por ejemplo, para hacer un MultiFasta de los genes de tssB y tssC dado el genoma de Bradyrhizobium valentinum LmjM3.
