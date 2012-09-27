#!/bin/bash

POLY="3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[1 1 1]"
N=1000
OUT=gmconnection-01

./gmconnection-p "$POLY" $N > $OUT

