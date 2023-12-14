#!/bin/bash
(
    CMD="time ./lamine -t 1 -o 2 -r 2 -no-po -m square_2mat_per.msh"
    echo $CMD >&2; $CMD
    CMD="time ./lamine -t 1 -o 2 -r 2 -no-po -m square_2mat_per.msh -pa"
    echo $CMD >&2; $CMD
    CMD="time ./lamine -t 1 -o 2 -r 1 -no-po -m square_2mat_per.msh -pa -d cuda"
    echo $CMD >&2; $CMD
    CMD="time ./lamine -t 3 -o 2 -r 1 -no-po -m cube_2mat_per.msh"
    echo $CMD >&2; $CMD
    CMD="time ./lamine -t 3 -o 2 -r 1 -no-po -m cube_2mat_per.msh -pa"
    echo $CMD >&2; $CMD
    CMD="time ./lamine -t 3 -o 2 -r 1 -no-po -m cube_2mat_per.msh -d cuda -pa"
    echo $CMD >&2; $CMD
) 1>out1 2>out2 
