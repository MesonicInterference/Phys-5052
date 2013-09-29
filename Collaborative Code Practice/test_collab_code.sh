#!/bin/bash

# This script is written for use with bash in Linux.  It compiles the main
# program and each of our .cpp files into object files so they can all be included
# in the final executable.  That part is accomplished in the line before the
# program's execution,  where g++ takes the object files as arguments instead of
# .cpp files.
#
# I've made this script executable on my laptop, but I'm not sure if that will
# transfer once I commit this to GitHub.  If you can't run this, check to make
# sure this file is set to be executable.

# compile the main program
g++ altered_number.cpp -c

# compile kyle's code
g++ kyle.cpp -c

# compile arjun's code
g++ arjun.cpp -c

# compile jeremy's code
g++ jeremy.cpp -c

# combine the object files
g++ altered_number.o kyle.o arjun.o jeremy.o -o altered_number.out

# run the program
./altered_number.out