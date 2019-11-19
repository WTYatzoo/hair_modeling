import matplotlib.image as image
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import pyplot
from PIL import Image
from PIL import ImageOps
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize
from scipy.optimize import least_squares
from scipy import misc
import scipy.io as io
import cv2
import numpy as np
import math
import random
import vtk

from OpenGL.GL import *
from OpenGL.arrays import vbo
from OpenGL.GLU import *
from OpenGL.GLUT import *
#from skimage import io

path_dic_result="./data_0/"
data=io.loadmat("./BFM_2009.mat")
len_data=len(data)
print(len(data))

for key in data:
    print(key)

topology=data["tl"]
print("topology")
print(topology.shape)

shape_avg=data["shapeMU"]
print("shape_avg")
print(shape_avg.shape)
shape_avg=shape_avg.reshape((-1,3))
print(shape_avg.shape)
x = shape_avg[:,0]
y = shape_avg[:,1]
z = shape_avg[:,2]

####################################
f = open(path_dic_result+"face_avg.vtk",'w')
f.write('# vtk DataFile Version 2.0\n')
f.write('point\n')
f.write('ASCII\n')
f.write('DATASET UNSTRUCTURED_GRID\n')
f.write('POINTS '+str(x.shape[0])+' double\n')
for i in range(x.shape[0]):
    f.write(str(x[i])+" "+str(y[i])+" "+str(z[i])+"\n")
f.write("CELLS "+str(topology.shape[0])+" "+str(topology.shape[0]*4)+"\n")
for i in range(topology.shape[0]):
    f.write("3 "+str(topology[i,0]-1)+" "+str(topology[i,1]-1)+" "+str(topology[i,2]-1)+"\n")
f.write("CELL_TYPES "+str(topology.shape[0])+"\n")
for i in range(topology.shape[0]):
    f.write("5\n")
f.close()
###############################################################

shape_EV=data["shapeEV"] # standard deviation 
print("shape_EV")
print(shape_EV.shape)
#print(shape_EV)
shape_EV=shape_EV.reshape((-1))

shape_pc=data["shapePC"]
print("shape_pc")
print(shape_pc.shape)
shape_pc=shape_pc.reshape((-1,3,199))
shape_pc_0=shape_pc[:,:,1]
print(shape_pc_0.shape)

#sd=np.sqrt(shape_EV[0])
sd=shape_EV[1]
x=x+5*sd*shape_pc_0[:,0]
y=y+5*sd*shape_pc_0[:,1]
z=z+5*sd*shape_pc_0[:,2]

print(x[0],y[0],z[0])

####################################
f = open(path_dic_result+"face_test_d.vtk",'w')
f.write('# vtk DataFile Version 2.0\n')
f.write('point\n')
f.write('ASCII\n')
f.write('DATASET UNSTRUCTURED_GRID\n')
f.write('POINTS '+str(shape_pc_0.shape[0])+' double\n')
for i in range(shape_pc_0.shape[0]):
    f.write(str(x[i])+" "+str(y[i])+" "+str(z[i])+"\n")
f.write("CELLS "+str(topology.shape[0])+" "+str(topology.shape[0]*4)+"\n")
for i in range(topology.shape[0]):
    f.write("3 "+str(topology[i,0]-1)+" "+str(topology[i,1]-1)+" "+str(topology[i,2]-1)+"\n")
f.write("CELL_TYPES "+str(topology.shape[0])+"\n")
for i in range(topology.shape[0]):
    f.write("5\n")
f.close()
###############################################################

num_basis=100
num_var=6+num_basis
begin_landmark_index=0

shape_EV_now=np.zeros(shape=(num_basis),dtype=np.float32)
shape_EV_now=shape_EV[0:num_basis]

#get 3d landmarks
landmark_3d=np.loadtxt(path_dic_result+"Landmarks68_BFM.txt")
print(landmark_3d.shape)

uu=np.zeros(shape=(landmark_3d.shape[0]-begin_landmark_index,3,num_basis),dtype=np.float32)
print(uu.shape)
#print(int(landmark_3d[0]))
#print(shape_pc[int(landmark_3d[0]),0,0],shape_pc[int(landmark_3d[0]),1,0],shape_pc[int(landmark_3d[0]),2,0])

for i in range(begin_landmark_index,landmark_3d.shape[0]):
    for j in range(num_basis):
        uu[i-begin_landmark_index,:,j]=shape_pc[int(landmark_3d[i]),:,j]

vv=np.zeros(shape=(landmark_3d.shape[0]-begin_landmark_index,3),dtype=np.float32)
print(vv.shape)
for i in range(begin_landmark_index,landmark_3d.shape[0]):
    vv[i-begin_landmark_index,:]=shape_avg[int(landmark_3d[i]),:]

#get 2d landmarks
landmark_2d=np.loadtxt(path_dic_result+"2d_landmark.txt")
print(landmark_2d.shape)
landmark_2d=landmark_2d[begin_landmark_index:,:]

xx_0=np.zeros(shape=(num_var),dtype=np.float32) #1+3+num_basis+2
print("xx_0.shape",xx_0.shape)
lower_bound=np.zeros(shape=(num_var),dtype=np.float32)
for i in range(num_var):
    if(i==0):
        lower_bound[0]=0
    else:
        lower_bound[i]=-np.inf

###############################################################################
#fig = pyplot.figure()
#ax = Axes3D(fig)
#ax.scatter(x,y,z)
#pyplot.show()
##################################################################################

def proj_step(xx,uu,vv,shape_EV_now):
    proj=np.array([[1,0,0],[0,1,0]])
    R_0=np.array([[1,0,0],[0,np.cos(xx[1]),np.sin(xx[1])],[0,-np.sin(xx[1]),np.cos(xx[1])]])
    R_1=np.array([[np.cos(xx[2]),0,-np.sin(xx[2])],[0,1,0],[np.sin(xx[2]),0,np.cos(xx[2])]])
    R_2=np.array([[np.cos(xx[3]),np.sin(xx[3]),0],[-np.sin(xx[3]),np.cos(xx[3]),0],[0,0,1]])

    R=np.matmul(R_0,R_1)
    R=np.matmul(R,R_2)
    #print("R\n",R)
    proj_all=xx[0]*np.matmul(proj,R)
    #print("proj_all\n",proj_all)

    here_0=xx[4:4+num_basis]
    #here_0=xx[4:4+num_basis]*shape_EV_now
    #print(here_0.shape)
    here_1=np.matmul(uu,here_0)
    #print(here_1.shape)
    here_2=here_1+vv
    #print(here_2.shape)
    
    #for i in range(0,u.shape[1]):
    #    print(i)
    #    print(np.matmul(proj_all,u[i,:].transpose()))

    
    here_2_t=here_2.transpose()
    ans=np.matmul(proj_all,here_2_t[None,:])
    #print(ans.shape)
    ans=ans.reshape((2,-1))
    x_translate=xx[num_basis+4:num_basis+6]
    x_translate.resize(2,1)
    ans=ans+x_translate
    return ans

#proj_step(xx_0,uu,vv)

def fun(xx,uu,vv,landmark_2d,shape_EV_now):
    print("fun")
    landmark_2d_t=landmark_2d.transpose()
    print(landmark_2d_t.shape)
    ans_fun=proj_step(xx,uu,vv,shape_EV_now)-landmark_2d_t
    print(ans_fun.shape)
    return ans_fun

#fun(xx_0,uu,vv,landmark_2d)

def fun_reg(xx,uu,vv,landmark_2d,shape_EV_now):
    print("fun_reg")
    reg=20 #0.001
    f=fun(xx,uu,vv,landmark_2d,shape_EV_now)
    print(f.shape)
    f=f.reshape((-1))
    print(f.shape)
    reg_help=np.sqrt(0.5*reg*np.square(xx[4:4+num_basis]/shape_EV_now))
    #print(shape_EV_now.shape)
    #print(reg_help.shape)
    f=np.append(f,reg_help)
    return f
###########################################################################
#fun_reg(xx_0,uu,vv,landmark_2d,shape_EV_now)
ans_proj_before=proj_step(xx_0,uu,vv,shape_EV_now)
input = plt.imread(path_dic_result+"hair_new_landmarks.png")
print(input.shape)
font = cv2.FONT_HERSHEY_SIMPLEX
for i in range(ans_proj_before.shape[1]):
    a=int(ans_proj_before[0,i])
    b=int(ans_proj_before[1,i])
    #print(a,b,"biubiubiu")
    input[b,a,:]=[0,0,1]
    cv2.putText(input, str(i), (a,b), font, 0.3, (187, 0, 255), 1, cv2.LINE_AA)
image.imsave(path_dic_result+'hair_new_landmarks_before_opt.png',input)

###########################################################################
#f_before=fun_reg(xx_0,uu,vv,landmark_2d,shape_EV_now)
#print(0.5*np.sum(np.square(f_before)))

res=least_squares(fun_reg,xx_0,method="trf",args=(uu,vv,landmark_2d,shape_EV_now),bounds=(lower_bound,np.inf),verbose=1)

print(res.x)
print(res.cost)
print(res.optimality)

x_result=res.x
for i in range(num_basis):
    #shape_avg=shape_avg+x_result[4+i]*shape_EV_now[i]*shape_pc[:,:,i]
    shape_avg=shape_avg+x_result[4+i]*shape_pc[:,:,i]
    zt_here=1

####################################
#ok
f = open(path_dic_result+"face_result.vtk",'w')
f.write('# vtk DataFile Version 2.0\n')
f.write('point\n')
f.write('ASCII\n')
f.write('DATASET UNSTRUCTURED_GRID\n')
f.write('POINTS '+str(shape_pc_0.shape[0])+' double\n')
for i in range(shape_pc_0.shape[0]):
    f.write(str(shape_avg[i,0])+" "+str(shape_avg[i,1])+" "+str(shape_avg[i,2])+"\n")
f.write("CELLS "+str(topology.shape[0])+" "+str(topology.shape[0]*4)+"\n")
for i in range(topology.shape[0]):
    f.write("3 "+str(topology[i,0]-1)+" "+str(topology[i,1]-1)+" "+str(topology[i,2]-1)+"\n")
f.write("CELL_TYPES "+str(topology.shape[0])+"\n")
for i in range(topology.shape[0]):
    f.write("5\n")
f.close()

##############################################################
#ok
proj=np.array([[1,0,0],[0,1,0]])
R_0=np.array([[1,0,0],[0,np.cos(x_result[1]),np.sin(x_result[1])],[0,-np.sin(x_result[1]),np.cos(x_result[1])]])
R_1=np.array([[np.cos(x_result[2]),0,-np.sin(x_result[2])],[0,1,0],[np.sin(x_result[2]),0,np.cos(x_result[2])]])
R_2=np.array([[np.cos(x_result[3]),np.sin(x_result[3]),0],[-np.sin(x_result[3]),np.cos(x_result[3]),0],[0,0,1]])
R=np.matmul(R_0,R_1)
R=np.matmul(R,R_2)

ss_1=x_result[0].copy()
RR_1=R.copy()

R=x_result[0]*R
shape_avg_t=shape_avg.transpose()
ans_3d=np.matmul(R,shape_avg_t[None,:])
ans_3d=ans_3d.reshape((3,-1))
ans=np.matmul(proj,ans_3d) #2D result after projection
x_translate=x_result[num_basis+4:num_basis+6]
x_translate.resize(2,1)

############
x_translate_3d=np.matrix([x_translate[0,0],x_translate[1,0],0])
x_translate_3d=x_translate_3d.transpose()
print("x_translate_3d.shape",x_translate_3d.shape)
global ans_3d_tra
ans_3d_tra=ans_3d+x_translate_3d

print("ans_3d_tra.shape",ans_3d_tra.shape)
print("max x",np.max(ans_3d_tra[0,:]))
print("min x",np.min(ans_3d_tra[0,:]))
print("max y",np.max(ans_3d_tra[1,:]))
print("min y",np.min(ans_3d_tra[1,:]))
print("max z",np.max(ans_3d_tra[2,:]))
print("min z",np.min(ans_3d_tra[2,:]))
##############

tt_1=x_translate.copy()

ans=ans+x_translate
print("ans",ans.shape)

ans_color=np.zeros(shape=(ans.shape[1],3),dtype=np.float)
input_tex=plt.imread(path_dic_result+"hair_new.png")
zbuffer=np.full((input_tex.shape[0],input_tex.shape[1]),np.inf)

for i in range(ans.shape[1]):
    a=int(ans[0,i])
    b=int(ans[1,i])
    ans_color[i,0:3]=input_tex[b,a,0:3]
print(input_tex[0,0,:])
#######################################################
#ok
f = open(path_dic_result+"face_result_color.vtk",'w')
f.write('# vtk DataFile Version 2.0\n')
f.write('point\n')
f.write('ASCII\n')
f.write('DATASET UNSTRUCTURED_GRID\n')
f.write('POINTS '+str(shape_pc_0.shape[0])+' double\n')
for i in range(shape_pc_0.shape[0]):
    f.write(str(shape_avg[i,0])+" "+str(shape_avg[i,1])+" "+str(shape_avg[i,2])+"\n")
f.write("CELLS "+str(topology.shape[0])+" "+str(topology.shape[0]*4)+"\n")
for i in range(topology.shape[0]):
    f.write("3 "+str(topology[i,0]-1)+" "+str(topology[i,1]-1)+" "+str(topology[i,2]-1)+"\n")
f.write("CELL_TYPES "+str(topology.shape[0])+"\n")
for i in range(topology.shape[0]):
    f.write("5\n")

f.write("CELL_DATA "+str(topology.shape[0])+"\n")
f.write("SCALARS sample_cell double 1\n")
f.write("LOOKUP_TABLE rgbtable_cell\n")
for i in range(topology.shape[0]):
    f.write(str(float(i)/topology.shape[0])+"\n")
    
f.write("LOOKUP_TABLE rgbtable_cell "+str(topology.shape[0])+"\n")
for i in range(topology.shape[0]):
    index_v_0=topology[i,0]-1
    index_v_1=topology[i,1]-1
    index_v_2=topology[i,2]-1
    if(i==101619 or i==101620):
        print(index_v_0,index_v_1,index_v_2)
        print(ans_color[index_v_0,:])
    ans_color_now=(ans_color[index_v_0,:]+ans_color[index_v_1,:]+ans_color[index_v_2,:])/3.0

    #f.write(str(ans_color[index_v_0,0])+" "+str(ans_color[index_v_0,1])+" "+str(ans_color[index_v_0,2])+" 1.0\n")
    f.write(str(ans_color_now[0])+" "+str(ans_color_now[1])+" "+str(ans_color_now[2])+" 1\n")

f.close()

###############################################################
#ok
ans_proj_result=proj_step(x_result,uu,vv,shape_EV_now)
input = plt.imread(path_dic_result+"hair_new_landmarks.png")
print(input.shape)
font = cv2.FONT_HERSHEY_SIMPLEX
for i in range(ans_proj_result.shape[1]):
    a=int(ans_proj_result[0,i])
    b=int(ans_proj_result[1,i])

    input[b,a,:]=[0,0,1]
    cv2.putText(input, str(i+begin_landmark_index), (a,b), font, 0.3, (187, 0, 255), 1, cv2.LINE_AA)

image.imsave(path_dic_result+'hair_new_landmarks_predict.png',input)

##############################################################
#ok
input = plt.imread(path_dic_result+"hair_new_landmarks.png")
print(input.shape)
font = cv2.FONT_HERSHEY_SIMPLEX
for i in range(ans.shape[1]):
    a=int(ans[0,i])
    b=int(ans[1,i])

    input[b,a,:]=[0,1,0]

image.imsave(path_dic_result+'hair_new_seg.png',input)

###################################################
#cal hair seg
reader = vtk.vtkUnstructuredGridReader()


#reader=vtk.vtkPolyDataReader()
reader.SetFileName(path_dic_result+'hair.vtk')
reader.Update()
output = reader.GetOutput()

hair_vertex=np.zeros(shape=(output.GetNumberOfPoints(),3),dtype=np.float)
print("hair_vertex.shape",hair_vertex.shape)
for i in range(output.GetNumberOfPoints()):
    pt = output.GetPoint(i)
    hair_vertex[i,:]=pt[:]
print(hair_vertex[0:22,:])

hair_cell=np.zeros(shape=(output.GetNumberOfCells(),3),dtype=np.int)
print("hair_cell.shape",hair_cell.shape)
for i in range(output.GetNumberOfCells()):
    cell=output.GetCell(i)
    hair_cell[i,0:3]=[cell.GetPointId(0),cell.GetPointId(1),cell.GetPointId(2)]
print(hair_cell[0:22,:])


np.savez(path_dic_result+'hair',a=hair_vertex,b=hair_cell)

#code the transform
ss_2=1.4e-06
tt_2=np.array([-0.000477652762006945,1.69798898873505,-0.038696890169479])
tt_2.resize(3,1)

proj=np.array([[1,0,0],[0,1,0]])
ss_3=ss_1/ss_2
RR_3=RR_1
help_1=ss_3*np.matmul(proj,RR_3)
help_2=np.matmul(help_1,tt_2)
tt_3=tt_1-help_2
print(tt_3)

hair_vertex_t=hair_vertex.transpose()

hair_3d=ss_3*np.matmul(RR_3,hair_vertex_t[None,:])
hair_3d=hair_3d.reshape((3,-1))

ssRR_3=ss_3*RR_3
np.savetxt(path_dic_result+'hair_sR.txt',ssRR_3)

hair_2d=np.matmul(help_1,hair_vertex_t[None,:])
hair_2d=hair_2d.reshape((2,-1))
hair_2d=hair_2d+tt_3

############
x_translate_3d_hair=x_translate_3d-ss_3*np.matmul(RR_3,tt_2)
np.savetxt(path_dic_result+'hair_translate.txt',x_translate_3d_hair)

#x_translate_3d_hair=x_translate_3d_hair.transpose()
print("x_translate_3d_hair.shape",x_translate_3d_hair.shape)
global hair_3d_tra
hair_3d_tra=hair_3d+x_translate_3d_hair

print("hair_3d_tra.shape",hair_3d_tra.shape)
print("max x",np.max(hair_3d_tra[0,:]))
print("min x",np.min(hair_3d_tra[0,:]))
print("max y",np.max(hair_3d_tra[1,:]))
print("min y",np.min(hair_3d_tra[1,:]))
print("max z",np.max(hair_3d_tra[2,:]))
print("min z",np.min(hair_3d_tra[2,:]))
##############

max_depth=max(np.max(hair_3d_tra[2,:]),np.max(ans_3d_tra[2,:]))
min_depth=min(np.min(hair_3d_tra[2,:]),np.min(ans_3d_tra[2,:]))

print("max min",max_depth,min_depth)
#############################################
#ok
input = plt.imread(path_dic_result+"hair_new_landmarks.png")
print(input.shape)
font = cv2.FONT_HERSHEY_SIMPLEX
for i in range(ans.shape[1]):
    a=int(ans[0,i])
    b=int(ans[1,i])

    input[b,a,:]=[0,1,0]
for i in range(hair_2d.shape[1]):
    a=int(hair_2d[0,i])
    b=int(hair_2d[1,i])

    if(input[b,a,0]!=0 or input[b,a,0]!=1 or input[b,a,0] !=0):
        input[b,a,:]=[1,0,0]

image.imsave(path_dic_result+'hair_new_hair_region.png',input)

###########################################
#opengl
global R_opengl, W, H
R_opengl=input.shape[1]/2
W=input.shape[1]
H=input.shape[0]

def init():
    glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA)
    glutInitWindowPosition(20, 20)
    glutInitWindowSize(W,H)
    glutCreateWindow("render")

    #for creating image file
    #glutHideWindow()

    glPointSize(4.0)
    glEnable(GL_NORMALIZE)
    glEnable(GL_DEPTH_TEST)
    glDepthFunc(GL_LESS)
    glClearDepth(1.0)
    glClearColor(1.0, 1.0, 1.0, 1.0)

    glEnable(GL_LINE_SMOOTH)
    glEnable(GL_POINT_SMOOTH)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
    glHint(GL_LINE_SMOOTH_HINT,GL_DONT_CARE)
    glLineWidth(1.0)		

    glEnable(GL_POLYGON_SMOOTH)

def reshape(w, h):
    if h <= 0:
        h = 1
    glViewport(0, 0, w, h)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    if w <= h:
        glOrtho(-R_opengl,R_opengl,-R_opengl* h / w, R_opengl * h / w,min_depth-5,max_depth+5)
    else:
        glOrtho(-R_opengl * w / h, R_opengl * w / h, -R_opengl, R_opengl,min_depth-5,max_depth+5)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    gluLookAt(W/2,H/2,-1, W/2,H/2,0,  0,-1,0)

    
def keyboard(key, x, y):
    if key == chr(27) or key == "q":  # Esc is 27
        sys.exit()
    
def drawfunc():
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()
    gluLookAt(W/2,H/2,-1, W/2,H/2,0,  0,-1,0)
    glColor3f(1.0, 1.0, 1.0) 

    '''
    glBegin(GL_POINTS)
    for i in range(ans_3d_tra.shape[1]):
        glVertex3fv(ans_3d_tra[:,i])
    glEnd()

    '''
    #print(topology.shape[0])
    for i in range(topology.shape[0]):
        glBegin(GL_POLYGON)
        index_0=topology[i,0]-1
        glVertex3fv(ans_3d_tra[:,index_0])
        index_1=topology[i,1]-1
        glVertex3fv(ans_3d_tra[:,index_1])
        index_2=topology[i,2]-1
        glVertex3fv(ans_3d_tra[:,index_2])
        glEnd()
    
    '''
    glColor3f(0.0, 1.0, 0.0) 
    glBegin(GL_POINTS)
    for i in range(hair_3d_tra.shape[1]):
        glVertex3fv(hair_3d_tra[:,i])
    glEnd()
    '''
    glColor3f(0.0, 0, 1.0) 
    
    #print(hair_cell.shape[0])
    for i in range(hair_cell.shape[0]):
        glBegin(GL_POLYGON)
        index_0=hair_cell[i,0]
        glVertex3fv(hair_3d_tra[:,index_0])
        index_1=hair_cell[i,1]
        glVertex3fv(hair_3d_tra[:,index_1])
        index_2=hair_cell[i,2]
        glVertex3fv(hair_3d_tra[:,index_2])
        glEnd()
       
    glutSwapBuffers()
    glFlush()

    x, y, width, height = glGetIntegerv(GL_VIEWPORT)
    print(x,y,width,height)
    data=glReadPixels(x,y,width,height, GL_RGB, GL_UNSIGNED_BYTE)
    image = Image.frombytes("RGB", (W, H), data)
    image = ImageOps.flip(image)
    image.save(path_dic_result+'hair_without_face.png', 'PNG')
    sys.exit()
        
    glutPostRedisplay()
    
def main():
    glutInit(sys.argv)
    init()
    glutReshapeFunc(reshape)
    glutDisplayFunc(drawfunc)
    glutKeyboardFunc(keyboard)
    glutMainLoop()
    
############################################
#ok
test_0=np.array([[1,2],[0,1]])
test_1=np.array([4,5])
test_3=np.matmul(test_0,test_1)
print(test_3)
test_1.resize(2,1)
test_2=test_0+test_1.transpose()
print(test_2)
#############################################

main()


