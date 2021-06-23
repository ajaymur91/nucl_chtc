#!/bin/bash
RETRY=$(($1+1))
mkdir -p work
mkdir -p log
mkdir -p err
mkdir -p out

rm -rf work/*
rm -rf Pr.txt
cp -r N_$RETRY/*.gro ./work/
