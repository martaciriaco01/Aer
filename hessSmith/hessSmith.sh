#!/bin/bash

# Directory dello Script
cd "$(dirname "$0")"

# Script di Comandi da Eseguire su XFoil
echo "NACA 0012
GDES
TGAP 0 0
EXEC

PPAR
N 200


SAVE NACA_0012.dat
y
OPER
PACC
polar.dat

ASEQ -4 8 2
PACC

NACA 0012
OPER
ALFA 2
CPWR cp.dat

QUIT" > hessSmith.txt

# Esecuzione Sequenza di Comandi su XFoil
xfoil < hessSmith.txt

# Rimozione Script
rm hessSmith.txt
