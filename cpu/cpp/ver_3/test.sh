#!/bin/bash

for var in {1..50}
do
	echo "benchmark$var"
	./obj/main ../../../benchmark/benchmark$var.cnf
done
