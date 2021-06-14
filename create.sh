#!/bin/bash
RETRY=$(($1+1))
mkdir -p work
rm -rf work/*
cp -r N_$RETRY/*.gro ./work/
