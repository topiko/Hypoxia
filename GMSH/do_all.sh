#!/bin/bash
NP=100
for layer in {1..14}
do
	python gmsh_from_file.py $layer $NP
        python fenics_msh_from_gmsh.py $layer $NP
done
