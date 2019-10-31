#! /bin/bash

project_dir=${HOME}/project/mytest/autohair_test/
mesh_obj_tag=${project_dir}/data/mesh_new/DB
mesh_vtk_tag=${project_dir}/data/mesh_new/DB
mesh_vtk_new_tag=${project_dir}/data/mesh_data_aug/DB
out_dir_tag=${project_dir}/data/mesh_data_aug/DB
python_path=${project_dir}/src/data_aug/


########################################
#define the cluster method 
cluster_method=agglomerativeclustering
#########################################
mesh_label_tag=${project_dir}/data/mesh_data_aug/DB_label_${cluster_method}

for ((i=0;i<20;i+=1))
do
    mesh_obj=${mesh_obj_tag}_${i}.obj
    mesh_vtk=${mesh_vtk_tag}_${i}.vtk
    mesh_vtk_new=${mesh_vtk_new_tag}_${i}.vtk
    out_dir=${out_dir_tag}_${i}
    distance_path=${out_dir}/distance.txt
    mesh_label=${mesh_label_tag}_${i}.vtk
    mkdir -p ${out_dir}
    ${project_dir}/src/data_aug/build/bin/data_aug cluster_method=${cluster_method} prog=decomposition mesh_vtk=${mesh_vtk} mesh_vtk_new=${mesh_vtk_new} out_dir=${out_dir} mesh_obj=${mesh_obj} python_path=${python_path} distance_path=${distance_path} mesh_label=${mesh_label}
done
