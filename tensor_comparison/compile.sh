#!/bin/bash

g++ 2plane.cpp -lgsl -lgslcblas -lm
LD_LIBRARY_PATH=/usr/local/lib
export LD_LIBRARY_PATH
./a.out > interp.dat
graph -T ps < interp.dat > interp.ps
evince interp.ps &
