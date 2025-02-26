#!/bin/bash
set -ex

quick_test=0

if [[ ($# -eq 1) && ($1 == '-q') ]]; then
    quick_test=1
fi

make clean
make tests

rm -rf test_output/*.dat
./tests


make

if [ "$quick_test" == 0 ]; then
    rm -rf exampleRun
    rm -rf exampleRun_nml
    
    omv test -V neuromlLocal/.test.w2d.nrn.omt
    omv test -V neuromlLocal/.test.w2d.omt
    
    omv test -V .test.example.omt
    omv test -V .test.nmlNS.omt
    
    #omv test -V neuromlLocal/.test.w2d.omt
    #omv test -V .test.example.omt

    #time python run_main.py -R 1233 -p 96 --doEvol --folderName exampleRun

    #time omv all -V 
fi

make tests2

rm -rf test_output_2/*.dat
./tests2