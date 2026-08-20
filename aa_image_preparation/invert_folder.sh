#!/bin/bash

for d in $(ls -d Camera_*/); do
  echo "- $d -"
  for e in $(ls $d/); do
    echo "--- $e ---"
    mkdir -p $e/$d
    mv $d/$e/*JPG $e/$d/
  done
done

# Move background images one folder down (exclude the Camera folder)
for background in $(find -name "Camera_*_background.JPG"); do
  target=$(echo $background | sed "s#/Camera_0./#/#")
  mv $background $target
done
