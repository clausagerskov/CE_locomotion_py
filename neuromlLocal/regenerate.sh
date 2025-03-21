#!/bin/bash
set -ex

./clean.sh 

ruff format *py

python regenerate.py exampleRunCEW2D

omv test -V .test.21w2d.nrn.omt

python regenerate.py exampleRun21W2D

omv test -V .test.21w2d.nrn.omt

cd testc302SigSim
./testc302NervousSystem --popString "AS DA DB DD VD VB VA" --popSize 7 --datString "21W2D"
cd ..

python regenerate.py exampleRun

#omv all -V 
omv test -V .test.w2d.nrn.omt
omv test -V .test.w2d.omt

cd testc302SigSim
make clean all
./testc302NervousSystem --popString "DA DB DD VD VA VB" --popSize 10 --datString "CEOrig"
cd ..

