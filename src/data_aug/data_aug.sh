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

for ((i=20;i<20;i+=1))
do
    mesh_obj_0=${mesh_obj_tag}_${i}.obj
    mesh_vtk_0=${mesh_vtk_tag}_${i}.vtk
    mesh_vtk_new_0=${mesh_vtk_new_tag}_${i}.vtk
    out_dir_0=${out_dir_tag}_${i}
    distance_path_0=${out_dir_0}/distance.txt
    mesh_label_0=${mesh_label_tag}_${i}.vtk
    mkdir -p ${out_dir_0}
    ${project_dir}/src/data_aug/build/bin/data_aug cluster_method=${cluster_method} prog=decomposition mesh_vtk_0=${mesh_vtk_0} mesh_vtk_new_0=${mesh_vtk_new_0} out_dir_0=${out_dir_0} mesh_obj_0=${mesh_obj_0} python_path=${python_path} distance_path_0=${distance_path_0} mesh_label_0=${mesh_label_0}
done

begin_num=100
max_num=15
for((i=${begin_num};i<${max_num};i+=1))
do
    for((j=$[$i+1];j<${max_num};j+=1))
    do
	#echo ${i} ${j}
	mesh_obj_0=${mesh_obj_tag}_${i}.obj
	mesh_vtk_0=${mesh_vtk_tag}_${i}.vtk
	mesh_vtk_new_0=${mesh_vtk_new_tag}_${i}.vtk
	out_dir_0=${out_dir_tag}_${i}
	distance_path_0=${out_dir_0}/distance.txt
	mesh_label_0=${mesh_label_tag}_${i}.vtk
	mkdir -p ${out_dir_0}

	mesh_obj_1=${mesh_obj_tag}_${j}.obj
	mesh_vtk_1=${mesh_vtk_tag}_${j}.vtk
	mesh_vtk_new_1=${mesh_vtk_new_tag}_${j}.vtk
	out_dir_1=${out_dir_tag}_${j}
	distance_path_1=${out_dir_1}/distance.txt
	mesh_label_1=${mesh_label_tag}_${j}.vtk
	mkdir -p ${out_dir_1}

	combine_result_path=${out_dir_tag}_combine_${i}_${j}/
	mkdir -p ${combine_result_path}

	${project_dir}/src/data_aug/build/bin/data_aug   prog=recombination mesh_vtk_0=${mesh_vtk_0} mesh_vtk_new_0=${mesh_vtk_new_0} out_dir_0=${out_dir_0} mesh_obj_0=${mesh_obj_0}  distance_path_0=${distance_path_0} mesh_label_0=${mesh_label_0}       mesh_vtk_1=${mesh_vtk_1} mesh_vtk_new_1=${mesh_vtk_new_1} out_dir_1=${out_dir_1} mesh_obj_1=${mesh_obj_1} distance_path_1=${distance_path_1} mesh_label_1=${mesh_label_1}  cluster_method=${cluster_method} python_path=${python_path} combine_result_path=${combine_result_path}
    done
done


i=13
j=14
#echo ${i} ${j}
mesh_obj_0=${mesh_obj_tag}_${i}.obj
mesh_vtk_0=${mesh_vtk_tag}_${i}.vtk
mesh_vtk_new_0=${mesh_vtk_new_tag}_${i}.vtk
out_dir_0=${out_dir_tag}_${i}
distance_path_0=${out_dir_0}/distance.txt
mesh_label_0=${mesh_label_tag}_${i}.vtk
mkdir -p ${out_dir_0}

mesh_obj_1=${mesh_obj_tag}_${j}.obj
mesh_vtk_1=${mesh_vtk_tag}_${j}.vtk
mesh_vtk_new_1=${mesh_vtk_new_tag}_${j}.vtk
out_dir_1=${out_dir_tag}_${j}
distance_path_1=${out_dir_1}/distance.txt
mesh_label_1=${mesh_label_tag}_${j}.vtk
mkdir -p ${out_dir_1}

combine_result_path=${out_dir_tag}_combine_${i}_${j}/
mkdir -p ${combine_result_path}

${project_dir}/src/data_aug/build/bin/data_aug   prog=recombination mesh_vtk_0=${mesh_vtk_0} mesh_vtk_new_0=${mesh_vtk_new_0} out_dir_0=${out_dir_0} mesh_obj_0=${mesh_obj_0}  distance_path_0=${distance_path_0} mesh_label_0=${mesh_label_0}       mesh_vtk_1=${mesh_vtk_1} mesh_vtk_new_1=${mesh_vtk_new_1} out_dir_1=${out_dir_1} mesh_obj_1=${mesh_obj_1} distance_path_1=${distance_path_1} mesh_label_1=${mesh_label_1}  cluster_method=${cluster_method} python_path=${python_path} combine_result_path=${combine_result_path}
