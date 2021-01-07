#!/bin/bash

#programa para ejecutar los casos sucesivos de Mercury.
#para ejecutar multiples casos, usar el formato

###################
#    echo "Running case: $CASE_NUMBER" 
##### Ejecutar el generador para cada caso
#    ./generator -$CASE_NUMBER > generator.out
##### Copiar el resultado a Results
#    cp big.in ./Results/big_$CASE_NUMBER.in
##### Ejecutar el simulador, guardar resultados en output.txt
#    ./mercury > output.txt
##### Copiar los resultados dentro de Results
#    cp xv.out ./Results/xv_$CASE_NUMBER.out
#    cp info.out ./Results/info_$CASE_NUMBER.out
#    cp output.txt ./Results/output_$CASE_NUMBER.txt
##### Borrar los temporales del integrador para volver a ejecutar 
#    rm *.out *.txt *.dmp *.tmp
###################

for CASE_NUMBER in $(seq 1 50)
do
    echo "Running case: $CASE_NUMBER" 
    ./generator -$CASE_NUMBER > generator.out
    cp big.in ./Results/big_$CASE_NUMBER.in

    ./mercury > output.txt

    cp xv.out ./Results/xv_$CASE_NUMBER.out
    cp info.out ./Results/info_$CASE_NUMBER.out
    cp output.txt ./Results/output_$CASE_NUMBER.txt 

    rm *.out *.txt *.dmp *.tmp

    ./mercuryMod > output.txt

    cp xv.out ./Results/xv_mod_$CASE_NUMBER.out
    cp info.out ./Results/info_mod_$CASE_NUMBER.out
    cp output.txt ./Results/output_mod_$CASE_NUMBER.txt 

    rm *.out *.txt *.dmp *.tmp embriones big.in
done
