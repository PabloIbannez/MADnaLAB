#!/bin/bash

BIN_path="/home/pablo/bin"
UAMMD_path="/home/pablo/UAMMD"

python ./setEnvironment.py $BIN_path $UAMMD_path

python -m pip install -r requirements.txt

make -C ./src
ln -s $MADNALABPATH/bin/MADnaLAB $BIN_path/MADnaLAB

echo Type: source ~/.bashrc
