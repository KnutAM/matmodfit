#!/bin/bash
# Script to run all examples

examples=( example01 example02 example04 example05 external_script user_error user_mpar_scaling user_optimization user_simulation )

for example in "${examples[@]}"
do
	cd $example
	matmodfit *.inp
	cd ..
done
