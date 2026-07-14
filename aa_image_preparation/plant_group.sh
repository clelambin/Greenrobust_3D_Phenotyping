#!/bin/bash
# Script which group files in separate folder based on relative time between each file

# Input variable
FILE_PATTERN="*.JPG"
WORKING_DIR="."
TIME_DELTA_SEC=10

# Set bash to stop when facing an error
set -euo pipefail

# Initialise timestamp and destination_folder
prev_timestamp=0
destination_folder=""

# Script function
move_to_destination () {
  # Script which check if file exist at destination, if not move it
  # Input: - file
  #        - destination
  if [ ! -f "$2/$1" ]; then
    mv "$1" "$2/$1"
  else
    # If destination file already exist, return code error 1
    echo "Error, file $1 already exists in destination"
    return 1 
  fi
}

# Loop through file and check stat between them
for file in $(ls -t --time=birth -r $WORKING_DIR/$FILE_PATTERN); do
  # Extract modification time from current file
  # timestamp: number of second since EPOCH
  # datetime:  YYYY-MM-DD_HH:MM:SS
  file_timestamp=$(stat --format="%Y" "$file")
  file_datetime=$(stat --format="%y" "$file" | sed "s/ /_/" | cut -d "." -f 1 | sed "s/:/-/g")

  if [ $((file_timestamp-prev_timestamp)) -gt $TIME_DELTA_SEC ]; then
    # Time difference bigger than input, start new group of file
    destination_folder="$WORKING_DIR/$file_datetime"
    mkdir -p "$destination_folder"
  fi

  # Move file to destination
  move_to_destination "$file" "$destination_folder"

  # Update previous timestamp
  prev_timestamp="$file_timestamp"
done
