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

def get_theta_conf(half_kernel_size,gabor_kernel,input_content,mask_img,iteration_num):
    input_content_shape=input_content.shape

    #theta_img=np.full((input_content_shape[0],input_content_shape[1]),np.nan)
    theta_img=np.zeros(shape=(input_content_shape[0],input_content_shape[1]),dtype=np.float32)
    
    theta_img[:,:]=-1 # invalid value
    #conf_img=np.full((input_content_shape[0],input_content_shape[1]),np.nan)
    conf_img=np.zeros(shape=(input_content_shape[0],input_content_shape[1]),dtype=np.float32) 
    
    F_sub=np.zeros(shape=(32),dtype=np.float32)
    
    for i in range(0,input_content_shape[0]):
        for j in range(0,input_content_shape[1]):
            if(np.fabs(mask_img[i,j]-1.0)<=0.00001):
                input_content_now=input_content[i-half_kernel_size:i+half_kernel_size+1,j-half_kernel_size:j+half_kernel_size+1]              
                score_now=-1
                theta_final=0                
                for theta_num in range(0,32):
                    theta_now=np.pi/31*theta_num
                    score=np.sum(input_content_now[:,:]*gabor_kernel[theta_num,:,:])                    
                    F_sub[theta_num]=score
                    if(theta_num==0 or score>score_now):
                        score_now=score
                        theta_final=theta_now

                conf_now=0
                for theta_num in range(0,32):
                    theta_now=np.pi/31*theta_num
                    mini_angle=np.min([np.abs(theta_now-theta_final),np.abs(theta_now-theta_final+np.pi),np.abs(theta_now-theta_final-np.pi)])
                    conf_now=conf_now+np.sqrt(mini_angle*((F_sub[theta_num]-score_now)**2))                    
                theta_img[i,j]=theta_final
                conf_img[i,j]=conf_now                                  
                #print(theta_img[i,j],conf_img[i,j])

    conf_img_max=np.nanmax(conf_img)
    conf_img_min=np.nanmin(conf_img)
    print("conf_img_min,conf_img_max",conf_img_min,conf_img_max)

    max_min=conf_img_max-conf_img_min
    thre=conf_img_min+threshold*(conf_img_max-conf_img_min)
    for i in range(0,input_content_shape[0]):
        for j in range(0,input_content_shape[1]):
            if(conf_img[i,j]<thre and conf_img[i,j]!=np.nan):
                #conf_img[i,j]=np.nan
                theta_img[i,j]=-1
                x=1
            elif(conf_img[i,j]>thre and conf_img[i,j]!=np.nan):                
                conf_img[i,j]=thre

    # visualize the theta_img & conf_img
    conf_img=(((conf_img-np.nanmin(conf_img))/(np.nanmax(conf_img)-np.nanmin(conf_img))))

    np.savetxt(img_data_path+"/"+"sparse_theta_map_"+str(iteration_num)+".txt",theta_img)
    
    theta_img_max=np.nanmax(theta_img)
    theta_img_min=np.nanmin(theta_img)
    theta_img=(theta_img-theta_img_min)/(theta_img_max-theta_img_min)
    image.imsave(img_data_path+"/theta_"+str(iteration_num)+".png",theta_img,cmap="jet")
    image.imsave(img_data_path+"/conf_"+str(iteration_num)+".png",conf_img,cmap="gray")
    return theta_img,conf_img


def refine_theta_conf(half_kernel_size,dis_kernel,conf_img_old,theta_img_old,hair_img,mask_img,iteration_num):
    theta_img_old_shape=theta_img_old.shape
    #theta_img=np.full((theta_img_old_shape[0],theta_img_old_shape[1]),np.nan)
    theta_img=np.zeros(shape=(theta_img_old_shape[0],theta_img_old_shape[1]),dtype=np.float32)
    conf_img=np.zeros(shape=(theta_img_old_shape[0],theta_img_old_shape[1]),dtype=np.float32) 
    #conf_img=np.full((input_content_shape[0],input_content_shape[1]),np.nan)
    
    for i in range(0,theta_img_old_shape[0]):
        for j in range(0,theta_img_old_shape[1]):
            if(np.fabs(mask_img[i,j]-1.0)<=0.00001):
                theta_content_old_now=theta_img_old[i-half_kernel_size:i+half_kernel_size+1,j-half_kernel_size:j+half_kernel_size+1]
                conf_content_old_now=conf_img_old[i-half_kernel_size:i+half_kernel_size+1,j-half_kernel_size:j+half_kernel_size+1]
                hair_content_now=hair_img[i-half_kernel_size:i+half_kernel_size+1,j-half_kernel_size:j+half_kernel_size+1]
                bilateral_kernel=dis_kernel*np.exp(-1*conf_content_old_now/conf_img_old[i,j])*np.exp(-1*(hair_content_now-hair_img[i,j])**2/0.01)
                theta_img[i,j]=np.sum(theta_content_old_now*bilateral_kernel)                
                conf_img[i,j]=np.sum(conf_content_old_now*bilateral_kernel)

    conf_img_max=np.nanmax(conf_img)
    conf_img_min=np.nanmin(conf_img)
    print("conf_img_min,conf_img_max",conf_img_min,conf_img_max)

    thre=conf_img_min+0.02*(conf_img_max-conf_img_min)
    for i in range(0,theta_img_old_shape[0]):
        for j in range(0,theta_img_old_shape[1]):
            if(conf_img[i,j]<thre and conf_img[i,j]!=np.nan):
                #conf_img[i,j]=np.nan
                #theta_img[i,j]=np.nan
                x=1

                
    # visualize the theta_img & conf_img
    conf_img=(conf_img-np.nanmin(conf_img))/(np.nanmax(conf_img)-np.nanmin(conf_img))
    theta_img_max=np.nanmax(theta_img)
    theta_img_min=np.nanmin(theta_img)
    theta_img=(theta_img-theta_img_min)/(theta_img_max-theta_img_min)
    image.imsave(img_data_path+"/theta_"+str(iteration_num)+".png",theta_img,cmap="jet")
    image.imsave(img_data_path+"/conf_"+str(iteration_num)+".png",conf_img,cmap="gray")
    return theta_img,conf_img


def prepro_img():

    ###################################################
    hair_img_new_path=img_data_path+"/"+"hair_new.png"
    hair_img_gray_new_path=img_data_path+"/"+"hair_gray_new.png"
    conf_img_new_path=img_data_path+"/"+"conf_new.png"
    mask_img_new_path=img_data_path+"/"+"mask_new.png"
    ori_dir_img_new_path=img_data_path+"/"+"ori_smooth_2d_new.png"

    theta_img_path=img_data_path+"/"+"theta_img.png"
    my_conf_img_path=img_data_path+"/"+"conf_img.png"


    hair_img_rgb_new = Image.open(hair_img_new_path)
    #hair_img_rgb_new.show()
    hair_img_gray_new = np.array(hair_img_rgb_new.convert('L'))
    hair_img_rgb_new=np.array(hair_img_rgb_new)
    

    hair_img_new=plt.imread(hair_img_new_path)
    hair_img_new_shape=hair_img_new.shape
    print("hair_img_new_shape",hair_img_new_shape)

    #conf_img_new=plt.imread(conf_img_new_path)
    #conf_img_new=conf_img_new[:,:,1]
    #conf_img_new_shape=conf_img_new.shape
    #print("conf_img_new_shape",conf_img_new_shape)

    #print(conf_img_new[280,367],conf_img_new[0,0])
    
    mask_img_new=plt.imread(mask_img_new_path)
    mask_img_new=mask_img_new[:,:,1]
    mask_img_new_shape=mask_img_new.shape
    print("mask_img_new dtype",mask_img_new.dtype)
    print("mask_img_new_shape",mask_img_new_shape)
    #print(mask_img_new[280,367],conf_img_new[0,0])
    
    #ori_dir_img_new=plt.imread(ori_dir_img_new_path)
    #ori_dir_img_new_shape=ori_dir_img_new.shape
    #print("ori_dir_img_new_shape",ori_dir_img_new_shape)
    
    ##########################################################

    kernel_size=11 #9
    hair_img_new_shape=hair_img_new.shape
    
    gabor_kernel=np.zeros(shape=(32,kernel_size,kernel_size),dtype=np.float32)
    #theta_img=np.zeros(shape=(hair_img_new_shape[0],hair_img_new_shape[1]),dtype=np.float32)
    #conf_img=np.zeros(shape=(hair_img_new_shape[0],hair_img_new_shape[1]),dtype=np.float32) 
    
    sigma_u=1.8 #1.4
    sigma_v=2.4 #1.6
    Lambda=4 #5

   
    ##############################
    
    for i in range(0,hair_img_new_shape[0]):
        for j in range(0,hair_img_new_shape[1]):
            if(np.fabs(mask_img_new[i,j]-0.0)<=0.00001):
                hair_img_gray_new[i,j]=0.0
                #   hair_img_new[i,j,:]=0.0
                #print(hair_img_new[i,j])

    mask_img_max=np.nanmax(mask_img_new)
    mask_img_min=np.nanmin(mask_img_new)
    print("mask_img_min,mask_img_max",mask_img_min,mask_img_max)
    mask_img_new=(mask_img_new-mask_img_min)/(mask_img_max-mask_img_min)
    mask_img_max=np.nanmax(mask_img_new)
    mask_img_min=np.nanmin(mask_img_new)
    print("mask_img_min,mask_img_max",mask_img_min,mask_img_max)
    
    hair_img_gray_new_max=float(np.nanmax(hair_img_gray_new))
    hair_img_gray_new_min=float(np.nanmin(hair_img_gray_new))
    print("hair_img_gray_new_min,hair_img_gray_new_max",hair_img_gray_new_min,hair_img_gray_new_max)
    hair_img_gray_new=(hair_img_gray_new-hair_img_gray_new_min)/(hair_img_gray_new_max-hair_img_gray_new_min)

    half_kernel_size=(kernel_size-1)/2
    for theta_num in range(0,32):
        theta_now=np.pi/31*theta_num
        #print("theta_now",theta_now)
        (y,x)=np.meshgrid(np.arange(-1*half_kernel_size,half_kernel_size+1),np.arange(-1*half_kernel_size,half_kernel_size+1))
        u=x*np.cos(theta_now)+y*np.sin(theta_now)
        v=-x*np.sin(theta_now)+y*np.cos(theta_now)
        gabor_kernel[theta_num,:,:]=np.exp(-0.5*(u**2/(sigma_u**2)+v**2/(sigma_v**2)))*np.cos(2.0*np.pi*u/Lambda)

    (theta_img,conf_img)=get_theta_conf(half_kernel_size,gabor_kernel,hair_img_gray_new,mask_img_new,0)
    for i in range(1,iteration_num+1,1):
        (theta_img,conf_img)=get_theta_conf(half_kernel_size,gabor_kernel,conf_img,mask_img_new,i)
    
    
    ############################################################

    '''
    kernel_size=11
    hair_img_new_shape=hair_img_new.shape
    bilateral_kernel=np.zeros(shape=(kernel_size,kernel_size),dtype=np.float32)
    sigma_d=3 
    sigma_p=1

    half_kernel_size=(kernel_size-1)/2
    (y,x)=np.meshgrid(np.arange(-1*half_kernel_size,half_kernel_size+1),np.arange(-1*half_kernel_size,half_kernel_size+1))
    dis_kernel=np.exp(-1*((y**2+x**2)/(sigma_d**2)))
    (theta_img,_)=refine_theta_conf(half_kernel_size,dis_kernel,conf_img,theta_img,hair_img_gray_new,mask_img_new,9)
    (theta_img,_)=refine_theta_conf(half_kernel_size,dis_kernel,conf_img,theta_img,hair_img_gray_new,mask_img_new,10)
    (theta_img,_)=refine_theta_conf(half_kernel_size,dis_kernel,conf_img,theta_img,hair_img_gray_new,mask_img_new,11)
    (theta_img,_)=refine_theta_conf(half_kernel_size,dis_kernel,conf_img,theta_img,hair_img_gray_new,mask_img_new,12)
    (theta_img,_)=refine_theta_conf(half_kernel_size,dis_kernel,conf_img,theta_img,hair_img_gray_new,mask_img_new,13)
    '''
        
prepro_img()
 
