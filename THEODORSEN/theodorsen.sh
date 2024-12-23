#!/bin/bash

# Directory dello Script
cd "$(dirname "$0")"

# Script di Comandi da Eseguire su XFoil
echo "LOAD BL430.dat
PPAR
N 200


OPER
ALFA 1.5
CPWR th15.dat
ALFA 1.7
CPWR th17.dat
ALFA 1.9
CPWR th19.dat
ALFA 3.0
CPWR th30.dat

QUIT" > theodorsen.txt

# Esecuzione Sequenza di Comandi su XFoil
xfoil < theodorsen.txt

# Rimozione Script
rm theodorsen.txt
