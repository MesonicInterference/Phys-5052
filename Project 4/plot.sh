#!/bin/bash

if [[ $1 == 1 ]]
  then
    /usr/bin/gnuplot phase_3d.gnu
fi

if [[ $1 == 2 ]]
  then
    /usr/bin/gnuplot phase_2d.gnu
fi

if [[ $1 == 3 ]]
  then
    /usr/bin/gnuplot pos_versus_time.gnu
fi

if [[ $1 == 4 ]]
  then
    /usr/bin/gnuplot vel_versus_time.gnu
fi