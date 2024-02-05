#!/bin/bash

help(){
	echo "Usage: $0 <folder>"
}

# Verifying if the user has provided the folder as an argument
if [ $# -eq 0 ]; then
	echo "Error: no argument provided" 1>&2
	help 1>&2
	exit 1
fi

## VARIABLES
folder=$1

## DECOMPRESS PROCESS

# Check if the folder provided exists
if [ ! -d $folder ]; then
    echo "Folder '$folder' does not exist."
    exit 1
fi

# Go to the directory of the folder given by the user
cd $folder/data/assemblies || exit 1

# Descompress all the files
for directory in *; do
	for compressed_file in $directory/*; do
		if [ -f $compressed_file ]; then
			extension="${compressed_file##*.}"
			case $extension in
				gz) gzip -d $compressed_file ;;
				*)
			esac
		fi
	done
done

echo "Process completed"
