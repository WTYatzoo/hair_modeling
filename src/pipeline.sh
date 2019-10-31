#! /bin/bash
height=725
width=500
sub_dir=test3
project_dir=${HOME}/project/mytest/autohair_test/
iteration_num=8
threshold=0.4

#python ./cut_image.py img_data_path ../data/${sub_dir} threshold ${threshold} iteration_num ${iteration_num}
#python ./prepro_image.py img_data_path ../data/${sub_dir} threshold ${threshold} iteration_num ${iteration_num}
python3 ./annotate/annotate_dir.py img_data_path ${project_dir}/data/${sub_dir} width ${width} height ${height}
python ./solve_amb.py img_data_path ../data/${sub_dir} iteration_num ${iteration_num}
source ./diffuse/diffuse.sh 
python ./visualize_dense.py img_data_path ../data/${sub_dir}
