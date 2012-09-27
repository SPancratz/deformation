#!/bin/bash

POLY="4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[1 1 1 1]"
N=100
OUT=gmconnection-02

./gmconnection-p "$POLY" $N > $OUT

