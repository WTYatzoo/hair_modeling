import face_alignment
from skimage import io
import cv2
import numpy as np

fa_3d = face_alignment.FaceAlignment(face_alignment.LandmarksType._3D, flip_input=False)
fa_2d = face_alignment.FaceAlignment(face_alignment.LandmarksType._2D, flip_input=False)
input = io.imread('./test/hair_new.jpg')

input_for_3d = input.copy()
input_for_2d = input.copy()
print (input.shape)

preds_3d = fa_3d.get_landmarks(input)
preds_2d = fa_2d.get_landmarks(input)
print(preds_3d)
print(preds_2d)
print(len(preds_2d[0]))

landmark_2d=np.zeros(shape=(len(preds_2d[0]),2),dtype=np.int32)
landmark_3d=np.zeros(shape=(len(preds_3d[0]),3),dtype=np.int32)
landmark_3d_proj=np.zeros(shape=(len(preds_3d[0]),2),dtype=np.int32)

font = cv2.FONT_HERSHEY_SIMPLEX

for i in range(len(preds_3d[0])):
    a=int(preds_3d[0][i][0])
    b=int(preds_3d[0][i][1])
    c=int(preds_3d[0][i][2])
    landmark_3d[i,:]=[a,b,c]
    landmark_3d_proj[i,:]=[a,b]
    print(a,b,c)
    input_for_3d[b,a,0]=255
    input_for_3d[b,a,1]=0
    input_for_3d[b,a,2]=0
    #if(i==64):
    #    cv2.putText(input, str(i), (a,b), font, 0.3, (187, 255, 255), 1, cv2.LINE_AA)

    cv2.putText(input_for_3d, str(i), (a,b), font, 0.3, (187, 255, 255), 1, cv2.LINE_AA)
#    print(input[a,b])

io.imsave('./test/hair_new_landmarks_3d.png',input_for_3d)
np.savetxt("./test/3d_landmark.txt",landmark_3d)
np.savetxt("./test/3d_landmark_proj.txt",landmark_3d_proj)

for i in range(len(preds_2d[0])):
    a=int(preds_2d[0][i][0])
    b=int(preds_2d[0][i][1])
    landmark_2d[i,:]=[a,b]
    print(a,b)
    input_for_2d[b,a,0]=255
    input_for_2d[b,a,1]=0
    input_for_2d[b,a,2]=0

    cv2.putText(input_for_2d, str(i), (a,b), font, 0.3, (187, 255, 255), 1, cv2.LINE_AA)

io.imsave('./test/hair_new_landmarks_2d.png',input_for_2d)
np.savetxt("./test/2d_landmark.txt",landmark_2d)




