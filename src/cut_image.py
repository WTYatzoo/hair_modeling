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
threshold=float(ph.getParam("threshold",0.4))
iteration_num=int(ph.getParam("iteration_num",6))

def prepro_img():

    '''
    #test for numpy's function
    test = np.full((2, 2), np.nan)
    test[0,0]=0.12
    test[1,1]=1.222

    A = np.array( [[1,1],
               [0,1]] )
    B = np.array( [[2,0],
               [3,4]] )
    print ("A=","\n",A)
    print ("B=","\n",B)
    print ("A*B=","\n",A*B) # Aij * Bij
    print ("A.dot(B)=","\n",np.dot(A,B)) # matrix multiply
    print(np.min([1,2,3]))
    print(np.nanmin(test))
    print(np.nanmax(test))

    '''
    

    ###################################################
    
    hair_img_path=img_data_path+"/"+"hair.png"
    conf_img_path=img_data_path+"/"+"conf.png"
    mask_img_path=img_data_path+"/"+"mask.png"
    ori_dir_img_path=img_data_path+"/"+"ori_smooth_2d.png"
    hair_img_new_path=img_data_path+"/"+"hair_new.png"
    hair_img_gray_new_path=img_data_path+"/"+"hair_gray_new.png"
    conf_img_new_path=img_data_path+"/"+"conf_new.png"
    mask_img_new_path=img_data_path+"/"+"mask_new.png"
    ori_dir_img_new_path=img_data_path+"/"+"ori_smooth_2d_new.png"


    theta_img_path=img_data_path+"/"+"theta_img.png"
    my_conf_img_path=img_data_path+"/"+"conf_img.png"


    hair_img_rgb = Image.open(hair_img_path)
    #   hair_img_rgb.show()
    hair_img_gray = np.array(hair_img_rgb.convert('L'))
    hair_img_rgb=np.array(hair_img_rgb)
    #   hair_img_gray.show()
    print("hair_img_rgb shape",hair_img_rgb.shape)
    print("hair_img_gray shape",hair_img_gray.shape)

    hair_img=plt.imread(hair_img_path)
    hair_img_shape=hair_img.shape
    print("hair_img_shape",hair_img_shape)

    conf_img=plt.imread(conf_img_path)
    conf_img_shape=conf_img.shape
    print("conf_img_shape",conf_img_shape)
    #print(conf_img[280,367],conf_img[0,0])
    
    mask_img=plt.imread(mask_img_path)
    mask_img_shape=mask_img.shape
    print("mask_img dtype",mask_img.dtype)
    print("mask_img_shape",mask_img_shape)
    #print(mask_img[280,367],conf_img[0,0])
    
    ori_dir_img=plt.imread(ori_dir_img_path)
    ori_dir_img_shape=ori_dir_img.shape
    print("ori_dir_img_shape",ori_dir_img_shape)
    #print(ori_dir_img[280,367,0],ori_dir_img[280,367,1],ori_dir_img[280,367,2])
    #print(ori_dir_img[281,367,0],ori_dir_img[281,367,1],ori_dir_img[280,367,2])
    #print(ori_dir_img[282,367,0],ori_dir_img[282,367,1],ori_dir_img[280,367,2])
    
    corner_point=np.zeros(shape=(4,2),dtype=np.int16)

    find_corner=0
    for i in range(0,hair_img_shape[0]):
        for j in range(0,hair_img_shape[1]):
            if (hair_img[i,j,0]>0.00001):
                find_corner=1
                corner_point[0,0]=i
                corner_point[0,1]=j
                break
        if (find_corner==1):
            break
    #print(corner_point[0,0],corner_point[0,1])

    find_corner=0
    for i in range(0,hair_img_shape[0]):
        for j in range(hair_img_shape[1]-1,-1,-1):
            if (hair_img[i,j,0]>0.00001):
                find_corner=1
                corner_point[1,0]=i
                corner_point[1,1]=j
                break
        if (find_corner==1):
            break
    #print(corner_point[1,0],corner_point[1,1])

    find_corner=0
    for i in range(hair_img_shape[0]-1,-1,-1):
        for j in range(0,hair_img_shape[1]):
            if (hair_img[i,j,0]>0.00001):
                find_corner=1
                corner_point[2,0]=i
                corner_point[2,1]=j
                break
        if (find_corner==1):
            break
    #print(corner_point[2,0],corner_point[2,1])

    '''
    find_corner=0
    for i in range(hair_img_shape[0]-1,-1,-1):
        for j in range(hair_img_shape[1]-1,-1,-1):
            if (hair_img[i,j,0]>0.00001):
                find_corner=1
                corner_point[3,0]=i
                corner_point[3,1]=j
                break
        if (find_corner==1):
            break
    print(corner_point[3,0],corner_point[3,1])
    '''
    
    hair_img_new=hair_img[corner_point[0,0]:corner_point[2,0],corner_point[0,1]:corner_point[1,1],:]
    image.imsave(hair_img_new_path,hair_img_new)
    hair_img_gray_new=hair_img_gray[corner_point[0,0]:corner_point[2,0],corner_point[0,1]:corner_point[1,1]]
    image.imsave(hair_img_gray_new_path,hair_img_gray_new,cmap='gray')
    conf_img_new=conf_img[corner_point[0,0]:corner_point[2,0],corner_point[0,1]:corner_point[1,1]]
    image.imsave(conf_img_new_path,conf_img_new,cmap='gray')
    mask_img_new=mask_img[corner_point[0,0]:corner_point[2,0],corner_point[0,1]:corner_point[1,1]]
    image.imsave(mask_img_new_path,mask_img_new,cmap='gray')
    ori_dir_img_new=ori_dir_img[corner_point[0,0]:corner_point[2,0],corner_point[0,1]:corner_point[1,1],:]
    image.imsave(ori_dir_img_new_path,ori_dir_img_new)
        
prepro_img()
 
