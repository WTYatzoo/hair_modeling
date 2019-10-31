#! /bin/bash

height=${height}
width=${width}

project_dir=${HOME}/project/mytest/autohair_test/
sub_dir=${sub_dir}
sparse_theta_map_path=${project_dir}/data/${sub_dir}/sparse_theta_map_used.vtk

dense_theta_map_path=${project_dir}/data/${sub_dir}/dense_theta_map_visual_10.vtk

dense_theta_map_py_path=${project_dir}/data/${sub_dir}/dense_theta_map_py.txt

${project_dir}/src/diffuse/build/bin/diffuse height=${height} width=${width} sparse_theta_map_path=${sparse_theta_map_path} dense_theta_map_path=${dense_theta_map_path} dense_theta_map_py_path=${dense_theta_map_py_path}
