#!/bin/bash

#for d in $(ls -d Camera_*/); do
for d in $(ls -d */); do
  echo "- $d -"
  for e in $(ls $d/); do
    echo "--- $e ---"
    mkdir -p $e/$d
    mv $d/$e/*JPG $e/$d/
  done
done
