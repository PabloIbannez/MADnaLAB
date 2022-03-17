#!/bin/bash

BIN_path="/home/pablo/bin"
UAMMD_path="/home/pablo/Desktop/UAMMD"

python ./src/setENV.py $BIN_path $UAMMD_path
source ~/.bashrc

python -m pip install -r requirements.txt

make
cp ./bin/MADnaLAB $BIN_path
