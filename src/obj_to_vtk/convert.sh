#!/bin/bash
 
path=../../data/mesh_test
path_new=../../data/mesh_test
files=$(ls $path)

i=0
for filename in ${files}
do
    python obj_to_vtk.py ${path}/${filename} ${path_new}/DB_${i}.vtk
    mv ${path}/${filename} ${path_new}/DB_${i}.obj
    let i+=1
    echo $filename >> filename.txt
done
