#!/bin/bash

POLY="4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[3 1 0 0] (2  0 1)[1 0 1 2] (2  0 1)[0 1 0 3]"
N=1
OUT=gmconnection-03

./gmconnection-p "$POLY" $N > $OUT

