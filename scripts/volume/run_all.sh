#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo $DIR
echo $(pwd)
bash $DIR/sasa.sh $1
bash $DIR/insert.sh $1
