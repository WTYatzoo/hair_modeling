import cv2
import matplotlib.image as image
import matplotlib.pyplot as plt
import numpy as np
import math
from multiprocessing.dummy import Pool as ThreadPool
from numpy.linalg import solve

lambda_normal=10
lambda_warping=1000
incre_s=20
sample_hair=100  #100#200
sample_image=1080 #1080#1000#530#2000+#530

path_dic_result="./data_0/"
scale=4
sigma=scale/2.0/math.sqrt(math.pi)
kernal_size=(int(6*sigma+1),int(6*sigma+1))
thres=0.4
iter_num=25

def gf(input):    
    input=cv2.GaussianBlur(input,kernal_size,sigma)
    input_here=input.copy()
    input_here[input_here<thres]=0
    return input_here

input_model=plt.imread(path_dic_result+"hair_without_face.png")
print(input_model.shape)
input_image_mask= plt.imread(path_dic_result+"mask_new.png")
print(input_image_mask.shape)

input_model_here=input_model[:,:,1]
input_image_mask=input_image_mask[:,:,1]
input_model_blur=input_model_here.reshape((input_model.shape[0],input_model.shape[1]))

for i in range(iter_num):
    input_model_blur=gf(input_model_blur)

input_model_blur[input_model_blur<thres]=0
input_model_blur[input_model_blur>=thres]=1
input_model_blur[input_model_blur==1]=-1
input_model_blur[input_model_blur==0]=1
input_model_blur[input_model_blur==-1]=0

input_model[:,:,0]=input_model_blur
input_model[:,:,1]=input_model_blur
input_model[:,:,2]=input_model_blur
        
image.imsave(path_dic_result+'hair_without_face_blur.png',input_model)

####################################################################################

def getContour(input):
    gray=cv2.cvtColor(input,cv2.COLOR_BGR2GRAY)
    print(gray.shape)
    ret,binary=cv2.threshold(gray,100,255,cv2.THRESH_BINARY)
    binary,contours,hierarchy = cv2.findContours(binary,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    #cv2.drawContours(input,contours,-1,(0,0,255),5)
    #cv2.imshow("img",input)
    #cv2.waitKey(0)
    cnts_size=len(contours)
    print(cnts_size)
    max_len=-1
    index_max=-1
    for i in range(cnts_size):
        print(contours[i].shape)
        if(contours[i].shape[0]>max_len):
            max_len=contours[i].shape[0]
            index_max=i
    return contours[index_max]

#input_model_blur=cv2.imread(path_dic_result+'mask_new.png')
input_model_blur=cv2.imread(path_dic_result+'hair_without_face_blur.png')
cont_hair=getContour(input_model_blur)
input_mask=cv2.imread(path_dic_result+'mask_new.png')
cont_image=getContour(input_mask)

cont_hair=cont_hair.reshape(-1,2)

for i in range(cont_hair.shape[0]):
    input_model_blur[cont_hair[i,1],cont_hair[i,0]]=[0,0,255]
cv2.imshow("img",input_model_blur)
cv2.waitKey(0)
cont_image=cont_image.reshape(-1,2)

for i in range(cont_image.shape[0]):
    input_mask[cont_image[i,1],cont_image[i,0]]=[0,0,255]
cv2.imshow("img",input_mask)
cv2.waitKey(0)

print(cont_hair.shape)
print(cont_image.shape)
########################################################################
def getPoint(input):
    point_set=[]
    point_set.append(input[0])
    for i in range(1,input.shape[0]):
        point_before=point_set[-1]
        point_now=input[i]
        if(point_now[0]!=point_before[0] and point_now[1]!=point_before[1]):
            if(point_before[0]>point_now[0]):
                dx=-1
            elif(point_before[0]<point_now[0]):
                dx=1
            if(point_before[1]>point_now[1]):
                dy=-1
            elif(point_before[1]<point_now[1]):
                dy=1

            max_dis=np.abs(point_now[0]-point_before[0])
            for j in range(1,max_dis+1):
                point_set.append(np.array([point_before[0]+j*dx,point_before[1]+j*dy]))
        elif(point_now[0]==point_before[0]):
            if(point_before[1]>point_now[1]):
                dx=-1
            else:
                dx=1
            for j in range(point_before[1]+dx,point_now[1]+dx,dx):
                point_set.append(np.array([point_now[0],j]))           
           
        elif(point_now[1]==point_before[1]):            
            if(point_before[0]>point_now[0]):
                dx=-1
            else:
                dx=1
            for j in range(point_before[0]+dx,point_now[0]+dx,dx):
                point_set.append(np.array([j,point_now[1]]))
                
    point_before=point_set[-1]
    point_now=input[0]
    if(point_now[0]!=point_before[0] and point_now[1]!=point_before[1]):
        if(point_before[0]>point_now[0]):
            dx=-1
        elif(point_before[0]<point_now[0]):
            dx=1
        if(point_before[1]>point_now[1]):
            dy=-1
        elif(point_before[1]<point_now[1]):
            dy=1
        max_dis=np.abs(point_now[0]-point_before[0])
        for j in range(1,max_dis):
            point_set.append(np.array([point_before[0]+j*dx,point_before[1]+j*dy]))
    elif(point_now[0]==point_before[0]):
        if(point_before[1]>point_now[1]):
            dx=-1
        else:
            dx=1
        for j in range(point_before[1]+dx,point_now[1],dx):
            point_set.append(np.array([point_now[0],j]))
    elif(point_now[1]==point_before[1]):
        if(point_before[0]>point_now[0]):
            dx=-1
        else:
            dx=1
        for j in range(point_before[0]+dx,point_now[0],dx):
            point_set.append(np.array([j,point_now[1]]))     
    return point_set


point_hair=getPoint(cont_hair)
#for i in range(len(point_hair)):
    #print(point_hair[i])
print(len(point_hair),cont_hair.shape[0])

point_image=getPoint(cont_image)
#for i in range(20):
    #print(point_image[i])
print(len(point_image),cont_image.shape[0])

#############################################################

def getNormal(point_set):
    normal_set=[]
    len_list=len(point_set)
    for i in range(len_list):
        idx_a=(i-1)%len_list
        idx_b=(i+1)%len_list
        a=point_set[i][:]-point_set[idx_a][:]
        #print("a",a)
        b=point_set[idx_b][:]-point_set[i][:]
        #print("b",b)
        x=float(-a[1]-b[1])
        y=float(a[0]+b[0])
        tt=np.sqrt(x*x+y*y)
        if(tt==0.0):
            help=normal_set[-1].copy()
            normal_set.append(help)
        else:
            x=x/tt
            y=y/tt
            #print(x,y)
            normal_set.append(np.array([x,y]))
    return normal_set


normal_hair=getNormal(point_hair)
print(len(normal_hair))

normal_image=getNormal(point_image)
print(len(normal_image))

############################################################################
'''
#get boundary correspondence
test_correspondence=input_model_blur.copy()
test_correspondence[:,:,:]=255

for i in range(len(point_hair)):
    test_correspondence[point_hair[i][1],point_hair[i][0]]=[255,0,0]
for i in range(len(point_image)):
    test_correspondence[point_image[i][1],point_image[i][0]]=[0,0,255]
cv2.imshow("img",test_correspondence)
cv2.waitKey(0)

incre_h=np.abs(len(normal_hair)/sample_hair)
index_hair=np.zeros(shape=(sample_hair),dtype=np.int)
for i in range(sample_hair):
    index_hair[i]=i*incre_h

incre_i=np.abs(len(normal_image)/sample_image)
index_image=np.zeros(shape=(sample_image),dtype=np.int)
for i in range(sample_image):
    index_image[i]=i*incre_i

#K*T
T1=np.zeros(shape=(sample_image,sample_hair),dtype=np.float)
#K*T
T2=np.zeros(shape=(sample_image,sample_hair),dtype=np.int)
#Ep or B #K*T
B=np.zeros(shape=(sample_image,sample_hair),dtype=np.float)
#Ee or A #K*K
A=np.zeros(shape=(sample_hair,sample_image,sample_image),dtype=np.float)

#point_image point_hair normal_image normal_hair
#index_hair index_image
for i in range(sample_hair):
    for j in range(sample_image):
        pp=point_hair[index_hair[i]]-point_image[index_image[j]]
        nn=(normal_hair[index_hair[i]]*normal_image[index_image[j]]).sum()
        B[j,i]=(pp*pp).sum()+lambda_normal*(1-nn)*(1-nn)

for i in range(sample_hair):
    print("yaya")
    for j in range(sample_image):
        for k in range(sample_image):
            pp_1=point_hair[index_hair[i]]-point_hair[index_hair[(i-1)%sample_hair]]
            pp_2=point_image[index_image[k]]-point_image[index_image[j]]
            tt=np.sqrt((pp_1*pp_1).sum())-np.sqrt((pp_2*pp_2).sum())
            A[i,j,k]=tt*tt

for i in range(sample_image):
    T1[i,0]=B[i,0]
    T2[i,0]=0

for j in range(1,sample_hair):
    for i in range(sample_image):
        min_now=-1
        index_k=-1
        for k in range(sample_image):
            help=T1[k,j-1]+A[j,k,i]+B[i,j]
            if(index_k==-1 or help<min_now):
                min_now=help
                index_k=k
        T1[i,j]=min_now
        T2[i,j]=index_k


ZT=np.argmin(T1[:,sample_hair-1])
print(ZT)
correspondence=[]
#hair image
correspondence.append(np.array([sample_hair-1,ZT]))

Zj=ZT
for j in range(sample_hair-1,1,-1):
    Zjb=T2[Zj,j]
    correspondence.append(np.array([j-1,Zjb]))
    Zj=Zjb

green = (0, 255, 0) 
for i in range(len(correspondence)):
    index_h=correspondence[i][0]
    index_i=correspondence[i][1]
    p_h=point_hair[index_hair[index_h]]
    p_i=point_image[index_image[index_i]]
    cv2.line(test_correspondence, (p_h[0],p_h[1]), (p_i[0],p_i[1]), green)

cv2.imshow("img",test_correspondence)
cv2.waitKey(0)

'''

#get boundary correspondence
test_correspondence=input_model_blur.copy()
test_correspondence[:,:,:]=255

for i in range(len(point_hair)):
    test_correspondence[point_hair[i][1],point_hair[i][0]]=[255,0,0]
for i in range(len(point_image)):
    test_correspondence[point_image[i][1],point_image[i][0]]=[0,0,255]
cv2.imshow("img",test_correspondence)
cv2.waitKey(0)

incre_h=np.abs(len(normal_hair)/sample_hair)
index_hair=np.zeros(shape=(sample_hair),dtype=np.int)
for i in range(sample_hair):
    index_hair[i]=i*incre_h

incre_i=np.abs(len(normal_image)/sample_image)
index_image=np.zeros(shape=(sample_image),dtype=np.int)
for i in range(sample_image):
    index_image[i]=i*incre_i

#K*T
T1=np.zeros(shape=(sample_image,sample_hair),dtype=np.float)
#K*T
T2=np.zeros(shape=(sample_image,sample_hair),dtype=np.int)
#Ep or B #K*T
B=np.zeros(shape=(sample_image,sample_hair),dtype=np.float)
#Ee or A #K*K

A1=np.zeros(shape=(sample_hair),dtype=np.float)
A2=np.zeros(shape=(sample_image,sample_image),dtype=np.float)  

#point_image point_hair normal_image normal_hair
#index_hair index_image
for i in range(sample_hair):
    for j in range(sample_image):
        pp=point_hair[index_hair[i]]-point_image[index_image[j]]
        nn=(normal_hair[index_hair[i]]*normal_image[index_image[j]]).sum()
        B[j,i]=(pp*pp).sum()+lambda_normal*(1-nn)*(1-nn)

for i in range(sample_hair):
    pp_1=point_hair[index_hair[i]]-point_hair[index_hair[(i-1)%sample_hair]]
    A1[i]=np.sqrt((pp_1*pp_1).sum())
    
for j in range(sample_image):
    for k in range(sample_image):
        pp_2=point_image[index_image[k]]-point_image[index_image[j]]
        A2[j,k]=np.sqrt((pp_2*pp_2).sum())

A3=np.zeros(shape=(sample_image,sample_image),dtype=np.float)  
A3_help=np.zeros(shape=(sample_image),dtype=np.float)

sum_len=0
for i in range(sample_image):
    pp_1=point_image[index_image[i]]-point_image[index_image[(i+1)%sample_image]]
    A3_help[i]=np.sqrt((pp_1*pp_1).sum())
    sum_len=sum_len+A3_help[i]

A3[0,0]=0
for i in range(1,sample_image):
    A3[0,i]=A3[0,i-1]+A3_help[i-1]
for j in range(1,sample_image):
    for k in range(j,sample_image):
        A3[j,k]=A3[j-1,k]-A3_help[j-1]

for j in range(sample_image):
    print(A3[j][j])

for j in range(sample_image):
    for k in range(j,sample_image):
        A3[j,k]=min(A3[j,k],sum_len-A3[j,k])
        A3[k,j]=A3[j,k]
        
#final A = (A1-A2)*(A1-A2)

for i in range(sample_image):
    T1[i,0]=B[i,0]
    T2[i,0]=0

for j in range(1,sample_hair):
    for i in range(sample_image):
        min_now=-1
        index_k=-1
        for k in range(sample_image):
            A_here=np.square(A1[j]-A2[k,i])
            help=T1[k,j-1]+A_here+B[i,j]
            if(index_k==-1 or help<min_now):
                min_now=help
                index_k=k
        T1[i,j]=min_now
        T2[i,j]=index_k


ZT=np.argmin(T1[:,sample_hair-1])
print(ZT)
correspondence=[]
#hair image
correspondence.append(np.array([sample_hair-1,ZT]))

Zj=ZT
for j in range(sample_hair-1,0,-1):
    Zjb=T2[Zj,j]
    correspondence.append(np.array([j-1,Zjb]))
    Zj=Zjb

green = (0, 255, 0)
len_cor=len(correspondence)
print("len_cor",len_cor)
for i in range(len_cor):
    index_h=correspondence[i][0]
    index_i=correspondence[i][1]
    p_h=point_hair[index_hair[index_h]]
    p_i=point_image[index_image[index_i]]
    cv2.line(test_correspondence, (p_h[0],p_h[1]), (p_i[0],p_i[1]), green)

cv2.imshow("img",test_correspondence)
cv2.waitKey(0)

#cv2.imshow("img",test_correspondence)
#cv2.waitKey(0)

#################################################################
#get warping function
K_P=np.zeros(shape=(len_cor+3,len_cor+3),dtype=np.float)
y_0=np.zeros(shape=(len_cor+3),dtype=np.float)
y_1=np.zeros(shape=(len_cor+3),dtype=np.float)
for i in range(len_cor):
    for j in range(i,len_cor):
        index_h_0=correspondence[i][0] # source
        index_h_1=correspondence[j][0] # source
        p_h_0=point_hair[index_hair[index_h_0]]
        p_h_1=point_hair[index_hair[index_h_1]]
        p=p_h_1-p_h_0
        len_here=np.sqrt(np.sum(p*p))

        if(i==j):
            len_here=len_here+1e-5
        K_P[i,j]=np.square(len_here)*np.log(len_here)
        K_P[j,i]=K_P[i,j]
#print(K_P)
for i in range(len_cor):
    K_P[i,i]=K_P[i,i]+1000

for i in range(len_cor):
    index_h_0=correspondence[i][0] # source
    p_h_0=point_hair[index_hair[index_h_0]]
    K_P[i,len_cor+0]=K_P[len_cor+0,i]=1
    K_P[i,len_cor+1]=K_P[len_cor+1,i]=p_h_0[0]
    K_P[i,len_cor+2]=K_P[len_cor+2,i]=p_h_0[1]

    index_i=correspondence[i][1]
    p_i=point_image[index_image[index_i]]

    y_0[i]=p_i[0]
    y_1[i]=p_i[1]


w_0=solve(K_P,y_0)
w_1=solve(K_P,y_1)


print("w_0",w_0)
print("w_1",w_1)

test_correspondence[:,:,:]=255
sample_i=test_correspondence.shape[0]/incre_s
sample_j=test_correspondence.shape[1]/incre_s
for i in range(sample_i):
    for j in range(sample_j):
        cv2.circle(test_correspondence, (j*incre_s,i*incre_s), 2, [0,0,0],-1)  
        #test_correspondence[i*incre_s,j*incre_s,:]=[0,0,0]

for i in range(len_cor):
    index_i=correspondence[i][1]
    p_i=point_image[index_image[index_i]]
    cv2.circle(test_correspondence, (p_i[0], p_i[1]), 5, [0,0,255],-1)  #26
    #test_correspondence[p_i[1],p_i[0],:]=[255,0,0]
cv2.imshow("img",test_correspondence)
cv2.waitKey(0)

test_correspondence[:,:,:]=255
sample_i=test_correspondence.shape[0]/incre_s
sample_j=test_correspondence.shape[1]/incre_s
for i in range(sample_i):
    for j in range(sample_j):
        cv2.circle(test_correspondence, (j*incre_s,i*incre_s), 2, [0,0,0],-1)  
        #test_correspondence[i*incre_s,j*incre_s,:]=[0,0,0]

for i in range(len_cor):
    index_h=correspondence[i][0]
    p_h=point_hair[index_hair[index_h]]
    cv2.circle(test_correspondence, (p_h[0], p_h[1]), 5, [0,255,0],-1)  #26
    #test_correspondence[p_i[1],p_i[0],:]=[255,0,0]
cv2.imshow("img",test_correspondence)
cv2.waitKey(0)


p_h_ori_sample=np.zeros(shape=(2,len_cor),dtype=np.float)
for i in range(len_cor):
    index_h=correspondence[i][0]
    p_h=point_hair[index_hair[index_h]]
    p_h_ori_sample[:,i]=[p_h[0],p_h[1]]


def get_warping_point(p_ori,p_h_ori_sample,w_0,w_1):
    print(p_h_ori_sample.shape)
    p_ori.resize(2,1)
    p_h_ori_sample_new=p_h_ori_sample-p_ori
    p_h_ori_sample_new=p_h_ori_sample_new*p_h_ori_sample_new
    p_h_ori_sample_new=np.sqrt(np.sum(p_h_ori_sample_new,0))

    #help
    p_h_ori_sample_new=p_h_ori_sample_new+1e-5
    print(p_h_ori_sample_new.shape)
    p_h_ori_sample_new_1=np.log(p_h_ori_sample_new)
    p_h_ori_sample_new_2=p_h_ori_sample_new*p_h_ori_sample_new
    p_h_ori_sample_new_3=p_h_ori_sample_new_1*p_h_ori_sample_new_2
    
    p_h_ori_sample_new_3=np.append(p_h_ori_sample_new_3,1.0)
    p_h_ori_sample_new_3=np.append(p_h_ori_sample_new_3,p_ori[0])
    p_h_ori_sample_new_3=np.append(p_h_ori_sample_new_3,p_ori[1])

    a=np.sum(p_h_ori_sample_new_3*w_0)
    b=np.sum(p_h_ori_sample_new_3*w_1)
    return a,b

test_correspondence[:,:,:]=255
#get_warping_point(np.array([0,0]),p_h_ori_sample,w_0,w_1)
#get_warping_point(p_ori,p_h_ori_sample,w_0,w_1)
for i in range(sample_i):
    for j in range(sample_j):
        p_ori=np.array([j*incre_s,i*incre_s])
        a,b=get_warping_point(p_ori,p_h_ori_sample,w_0,w_1)
        cv2.circle(test_correspondence, (int(a),int(b)), 2, [0,0,0],-1)  
        #test_correspondence[i*incre_s,j*incre_s,:]=[0,0,0]

for i in range(len_cor):
    index_h=correspondence[i][0]
    p_h=point_hair[index_hair[index_h]]
    a,b=get_warping_point(p_h,p_h_ori_sample,w_0,w_1)
    cv2.circle(test_correspondence, (int(a),int(b)), 5, [0,255,0],-1)  #26
    #test_correspondence[p_i[1],p_i[0],:]=[255,0,0]
cv2.imshow("img",test_correspondence)
cv2.waitKey(0)

cv2.imshow("img",test_correspondence)
cv2.waitKey(0)


