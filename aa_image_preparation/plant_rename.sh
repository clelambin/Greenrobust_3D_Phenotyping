#!/bin/bash
# Script which rename group folder based on input file list

# Input variable
RENAME_FILE="$1"
BACKGROUND_FOLDER="Background"
BACKGROUND_THRESH=5
PROCESSING_FOLDER="2026*"

# Set bash to stop when facing an error
set -euo pipefail

# Compare number of folder to content of input file
nb_folder=$(ls -d -c1 $PROCESSING_FOLDER/ | wc -l)
nb_input=$(cat "$RENAME_FILE" | wc -l)
processed_count=$(ls --ignore="2026*" --ignore="*.*" -c1 | wc -l)
if [ $nb_folder -ne $((nb_input-processed_count)) ]; then
  echo "Error, different file count in input file than working directory"
  echo "Number of folder in working directory: $nb_folder"
  echo "Number of lines in input file: $((nb_input-processed_count))"
  exit 1
fi

# Initialise count of processed folder and start from next folder
folder_count=$((processed_count+1))
# Look through folder and rename based on input folder name
for folder in $(ls -d $PROCESSING_FOLDER/); do
  # Count number of file in folder, compare with expected folder name
  folder_content=$(ls -c1 "$folder" | wc -l)

  # Get next name in input file
  rename_folder=$(cat "$RENAME_FILE" | head -n $folder_count | tail -n 1)

  # Check folder to input matching based on folder content
  # If less file than specified in input background threshold, check if name is input background name
  if [ $folder_content -le $BACKGROUND_THRESH ] && [ "$rename_folder" != "$BACKGROUND_FOLDER" ]; then
    # Expecting Background file alert the user and cancel the script
    echo "Error, $folder only include $folder_content files but expected name ($rename_folder) is not a background"
    exit 1
  # If more than specified in input background threshold, check if name is input background name
  elif [ $folder_content -gt $BACKGROUND_THRESH ] && [ "$rename_folder" == "$BACKGROUND_FOLDER" ]; then
    echo "Error, $folder is a background folder but contain more files than input for background"
    exit 1
  fi

  # If background folder, append original file name to folder
  if [ "$rename_folder" == "$BACKGROUND_FOLDER" ]; then
    rename_folder="${BACKGROUND_FOLDER}_${folder}"
  fi

  # Check if folder already exist
  if [ -e "$rename_folder" ]; then
    echo "Error, attempt to rename $folder into $rename_folder, but folder already exists"
    exit 1
  fi

  # Rename folder to match input file and increment folder count
  mv "$folder" $rename_folder
  folder_count=$((folder_count+1))
done
