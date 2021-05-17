#!/bin/bash

# Exit the script when pressing Ctrl+C
trap "exit" SIGINT

# Run the program
make
./run.x

mv data2.ecsv ./data/data2.ecsv
rm ./data/prueba*