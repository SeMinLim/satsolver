#!/bin/bash

for var in {333..400}
do
	echo "benchmark$var"
	./obj/main ../../../benchmark/benchmark$var.cnf
done
