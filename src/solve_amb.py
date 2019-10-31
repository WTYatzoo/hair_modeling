from scipy import misc
import cv2
import numpy as np
import math
import paramhelpers as ph
import matplotlib.image as image
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from PIL import Image

img_data_path=ph.getParam("img_data_path","../data/test1")
iteration_num=ph.getParam("iteration_num",6)
def get_theta_map():
    sparse_theta_map=np.loadtxt(img_data_path+"/"+"sparse_theta_map_"+str(iteration_num)+".txt")
    dir_map_x=np.loadtxt(img_data_path+"/"+"dir_map_x.txt")
    dir_map_y=np.loadtxt(img_data_path+"/"+"dir_map_y.txt")
    whether_hair=np.loadtxt(img_data_path+"/"+"whether_hair.txt")
    #print(dir_map_x.shape,dir_map_x[280,367],dir_map_y[280,367])
    #print(sparse_theta_map.shape,sparse_theta_map[280,367])
    #print(whether_hair.shape,whether_hair[280,367])

    mat_shape=dir_map_x.shape

    ####################################
    #write the sparse theta map into vtk for vector visualization
    f = open(img_data_path+"/sparse_theta_map_origin.vtk",'w')
    f.write('# vtk DataFile Version 2.0\n')
    f.write('point\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS '+str(mat_shape[0]*mat_shape[1])+' double\n')
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write(str(i)+" "+str(j)+" 0\n")

    f.write("CELLS "+str(mat_shape[0]*mat_shape[1])+" "+str(mat_shape[0]*mat_shape[1]*2)+"\n")

    index=0
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write("1 "+str(index)+"\n")
            index=index+1

    f.write("CELL_TYPES "+str(mat_shape[0]*mat_shape[1])+"\n")
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write("1\n")
    
    f.write("POINT_DATA "+str(mat_shape[0]*mat_shape[1])+"\n")
    f.write("SCALARS theta double \n")
    f.write("LOOKUP_TABLE default\n")
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            if(sparse_theta_map[i,j]>-0.00001):
                #f.write(str(dir_map_x[i,j])+" "+str(dir_map_y[i,j])+" 0\n")
                f.write(str(sparse_theta_map[i,j])+"\n")
            else:
                f.write("0\n")

    f.write("SCALARS whether_hair double \n")
    f.write("LOOKUP_TABLE default\n")
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write(str(whether_hair[i,j])+"\n")          

    f.close()
    ###############################################################
    index_trace=0
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):

            #if(index_trace==170858):
                #print(np.cos(sparse_theta_map[i,j])*dir_map_x[i,j]+np.sin(sparse_theta_map[i,j])*(-1)*dir_map_y[i,j])
                #print(sparse_theta_map[i,j],"trace_here")
                #print(whether_hair[i][j])
                #print(i,j)
                #print(dir_map_x[i,j],dir_map_y[i,j])
                
            index_trace=index_trace+1
            if(whether_hair[i,j]>0 and sparse_theta_map[i,j]>-0.00001):
                
                if(np.cos(sparse_theta_map[i,j])*dir_map_x[i,j]+np.sin(sparse_theta_map[i,j])*(-1)*dir_map_y[i,j]>0):                    
                    zt=1
                else:
                    sparse_theta_map[i,j]=sparse_theta_map[i,j]+np.pi
            

    ####################################
    #write the sparse theta map into vtk for vector visualization
    f = open(img_data_path+"/sparse_theta_map.vtk",'w')
    f.write('# vtk DataFile Version 2.0\n')
    f.write('point\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS '+str(mat_shape[0]*mat_shape[1])+' double\n')
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write(str(i)+" "+str(j)+" 0\n")

    f.write("CELLS "+str(mat_shape[0]*mat_shape[1])+" "+str(mat_shape[0]*mat_shape[1]*2)+"\n")

    index=0
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write("1 "+str(index)+"\n")
            index=index+1

    f.write("CELL_TYPES "+str(mat_shape[0]*mat_shape[1])+"\n")
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write("1\n")
    
    f.write("POINT_DATA "+str(mat_shape[0]*mat_shape[1])+"\n")
    f.write("SCALARS theta double \n")
    f.write("LOOKUP_TABLE default\n")
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            if(sparse_theta_map[i,j]>-0.00001):
                #f.write(str(dir_map_x[i,j])+" "+str(dir_map_y[i,j])+" 0\n")
                f.write(str(sparse_theta_map[i,j])+"\n")
            else:
                f.write("0\n")
                
    f.write("SCALARS whether_hair double \n")
    f.write("LOOKUP_TABLE default\n")
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write(str(whether_hair[i,j])+"\n")
            
    '''
    f.write("VECTORS v double\n")
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            if(sparse_theta_map[i,j]>-0.00001):
                #f.write(str(dir_map_x[i,j])+" "+str(dir_map_y[i,j])+" 0\n")
                f.write(str(np.cos(sparse_theta_map[i,j]))+" "+str(np.sin(sparse_theta_map[i,j]))+" 0\n")
            else:
                f.write("0 0 0\n")
                

    '''
    
    f.close()
    ###############################################################
    

    ####################################
    #write the sparse theta map into vtk for vector visualization
    f = open(img_data_path+"/sparse_theta_map_used.vtk",'w')
    f.write('# vtk DataFile Version 2.0\n')
    f.write('point\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS '+str(mat_shape[0]*mat_shape[1])+' double\n')
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write(str(i)+" "+str(j)+" 0\n")

    f.write("CELLS "+str(mat_shape[0]*mat_shape[1])+" "+str(mat_shape[0]*mat_shape[1]*2)+"\n")

    index=0
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write("1 "+str(index)+"\n")
            index=index+1

    f.write("CELL_TYPES "+str(mat_shape[0]*mat_shape[1])+"\n")
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write("1\n")
    
    f.write("POINT_DATA "+str(mat_shape[0]*mat_shape[1])+"\n")
    f.write("SCALARS theta double \n")
    f.write("LOOKUP_TABLE default\n")
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            if(sparse_theta_map[i,j]>-0.00001):
                #f.write(str(dir_map_x[i,j])+" "+str(dir_map_y[i,j])+" 0\n")
                f.write(str(sparse_theta_map[i,j])+"\n")
            else:
                f.write("-1\n")
                
    f.write("SCALARS whether_hair double \n")
    f.write("LOOKUP_TABLE default\n")
    for i in range(mat_shape[0]):
        for j in range(mat_shape[1]):
            f.write(str(whether_hair[i,j])+"\n")
    
    f.close()
    ###############################################################
    

get_theta_map()
    
