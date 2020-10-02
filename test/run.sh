#!/bin/bash

../bin/gen_input 1 1 1 inputs ../data/bus_data ../data/line_data
../bin/gen_fail inputs inputs_fail ../data/line_data 1
../bin/alg2 ../data/bus_data ../data/line_data 0 0 0