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

def get_theta_map():
    dense_theta_map=np.loadtxt(img_data_path+"/"+"dense_theta_map_py.txt",dtype=np.float32)
    print(dense_theta_map.shape)
    shape_here=dense_theta_map.shape
    visualize_img=np.zeros(shape=(dense_theta_map.shape[0],dense_theta_map.shape[1],3),dtype=np.float32)
    for i in range (shape_here[0]):
        for j in range  (shape_here[1]):
            if(dense_theta_map[i,j]>-0.0001):
                visualize_img[i,j,0]=(np.cos(dense_theta_map[i,j])+1)*0.5
                visualize_img[i,j,1]=(np.sin(dense_theta_map[i,j])+1)*0.5
                visualize_img[i,j,2]=0
    image.imsave(img_data_path+"/"+"dense_theta_map_visalize.png",visualize_img)

get_theta_map()
    
