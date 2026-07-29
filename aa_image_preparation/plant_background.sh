#!/bin/bash
# Script which move backgound files into corresponding folder

# Input variable
BACKGROUND_RENAME="$1"
BACKGROUND_FOLDER="Background"

# Set bash to stop when facing an error
set -euo pipefail

# Backup previous background timestamp file
if [ -f timestamp_info.txt ]; then
  mv timestamp_info.txt timestamp_info_backup.txt
fi

# Loop through background folder and map folder name to modified date
for folder in $(ls -d ${BACKGROUND_FOLDER}*); do
  # Check that only one file in folder
  background_file=$(ls -c1 $folder)
  if [ $(echo "$background_file" | wc -l) -gt 1 ]; then
    echo "Error, more than one file in $folder"
    exit 1
  fi
  # If input argument specify, rename files in background folders
  if [ ! -z "$BACKGROUND_RENAME" ]; then
    mv "$folder/$background_file" "$folder/$BACKGROUND_RENAME"
    background_file="$BACKGROUND_RENAME"
  fi
  # Append timestamp into separate file
  background_timestamp=$(stat --format="%Y" "$folder/$background_file")
  echo "$background_timestamp,$folder/$background_file" >> timestamp_info.txt
done

# Sort timestamp_info.txt by timestamp
sort -n timestamp_info.txt -o timestamp_info.txt

# Loop through non background folders and add matching background
for folder in $(ls --ignore="${BACKGROUND_FOLDER}*" --ignore="*.*"); do
  # Check timestamp of one picture within processing folder
  sample_picture=$(ls -c1 $folder | head -n 1)
  sample_timestamp=$(stat --format="%Y" $folder/$sample_picture)
  # Find matching background file: last background file with timestamp smaller than sample picture
  matching_background=$(awk -v target="$sample_timestamp" -F , '$1 < target { var = $2 } END {print var}' timestamp_info.txt)
  # If matching background file exist, copy it into current folder
  if [ -f $matching_background ]; then
    #echo "Moving $matching_background to $folder"
    cp "$matching_background" "$folder"
  else
    echo "Warning, $matching_background does not exists, skipping $folder"
  fi
done
