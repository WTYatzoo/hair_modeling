'''
import sys
from PyQt5.QtWidgets import (QWidget, QToolTip, 
    QPushButton, QApplication)
from PyQt5.QtGui import QFont     
class Example(QWidget):
    
    def __init__(self):
        super().__init__()        
        self.initUI()               
    def initUI(self):
        QToolTip.setFont(QFont('SansSerif', 10))
        self.setToolTip('This is a <b>QWidget</b> widget')
        btn = QPushButton('Button', self)
        btn.setToolTip('This is a <b>QPushButton</b> widget')
        btn.resize(btn.sizeHint())
        btn.move(50, 50)               
        self.setGeometry(300, 300, 300, 200)
        self.setWindowTitle('Tooltips')    
        self.show()
                
if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
'''

'''

import sys
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

class picture(QWidget):
    def __init__(self):
        super(picture, self).__init__()
        self.resize(820, 1000)
        self.setWindowTitle("test")
        self.label = QLabel(self)
        #self.label.setText(" show image")
        self.label.setFixedSize(500, 625)
        self.label.move(160, 160)
        self.label.setStyleSheet("QLabel{background:white;}" "QLabel{color:rgb(300,300,300,120);font-size:10px;font-weight:bold;font-family:;}" )
        btn = QPushButton(self)
        btn.setText("open image")
        btn.move(10, 30)
        btn.clicked.connect(self.openimage)
    def openimage(self):
        imgName, imgType = QFileDialog.getOpenFileName(self, "open image", "", "*.jpg;;*.png;;All Files(*)")
        jpg = QtGui.QPixmap(imgName).scaled(self.label.width(), self.label.height())
        self.label.setPixmap(jpg)
        print(jpg.shape)

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    my_pic = picture()
    my_pic.show()
    sys.exit(app.exec_())
'''

import sys
import copy
import numpy as np
import paramhelpers as ph
from PyQt5.QtWidgets import QApplication, QLabel, QSlider
from PyQt5 import QtWidgets, QtCore, QtGui
import PyQt5.QtCore 
from PyQt5.QtCore import Qt, QPoint
from PyQt5.QtGui import QPixmap,  QImage, QPainter
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
import cv2

img_data_path=ph.getParam("img_data_path","../../data/test1")
width=int(ph.getParam("width",0))
height=int(ph.getParam("height",0))

class BtnLabel(QLabel):
    global my_dialog
    def __init__(self,parent=None):  
        super(BtnLabel,self).__init__(parent)

        self.if_mouse_press = False
        self.segment_status=False
        self.dir_status=False
        self.pix = QPixmap() #empty painterdevice
        self.tempix=QPixmap() #empty temporary painterdevice
        self.lastPoint = QPoint() 
        self.endPoint = QPoint()

        self.pos_xy=[] #
        self.dir_set=[] #

    def set_segment(self,segment_status):
        self.segment_status=segment_status
        self.dir_status=False

    def set_dir(self,dir_status):
        self.dir_status=dir_status
        print(self.dir_status)
        self.segment_status=False
        
    def paintEvent(self, event):
        if(self.segment_status==True):
            print("segment")
            painter = QPainter(self)
            painter.begin(self)

            pen = QPen(Qt.red, 2, Qt.SolidLine)
            pp=QPainter(self.pix)
            pp.setPen(pen)
            pos_size=len(self.pos_xy)
            if(pos_size>=3):
                for i in range(0,pos_size-1):
                    #print(i,"i")
                    point_start = self.pos_xy[i]
                    point_end = self.pos_xy[i+1]
                    if point_end == (-1, -1) and i==pos_size-2:
                        break
                    elif point_end==(-1,-1) and i!=pos_size-2:
                        i+=1 # no use
                    elif point_start==(-1,-1):
                        continue
                    else:
                        print(point_start,"start")
                        print(point_end,"end")
                        pp.drawLine(point_start[0], point_start[1], point_end[0], point_end[1])

            painter.drawPixmap(0,0,self.pix)
            painter.end()   
                                                                  
            
        elif(self.dir_status==True):
            print("dir")
            painter = QPainter(self)
            painter.begin(self)

            # draw the direction with arrow
            pp=QPainter(self.pix)
            pen = QPen(Qt.green, 2, Qt.SolidLine)
            pp.setPen(pen)
            brush=QBrush(Qt.green,Qt.SolidPattern)
            pp.setBrush(brush)

            line=QLineF(QPointF(self.lastPoint.x(),self.lastPoint.y()),QPointF(self.endPoint.x(),self.endPoint.y()))
            v=line.unitVector()
            v.setLength(8)
            v.translate(QPointF(line.dx(),line.dy()))
            n=v.normalVector()
            n.setLength(n.length() * 0.5)
            n2 = n.normalVector().normalVector()
            p1 = v.p2()
            p2 = n.p2()
            p3 = n2.p2()        
            pp.drawLine(self.lastPoint.x(),self.lastPoint.y(), self.endPoint.x(), self.endPoint.y())
            pp.drawPolygon(p1, p2, p3)
            
            painter.drawPixmap(0,0,self.pix)
            painter.end()
        else:
            painter = QPainter(self)
            painter.begin(self)
                        
            painter.drawPixmap(0,0,self.pix)
            painter.end()

        '''
        painter = QPainter(self)
        painter.begin(self)
        x=self.lastPoint.x()
        y=self.lastPoint.y()
        w=self.endPoint.x()-x
        h=self.endPoint.y()-y

        if self.if_mouse_press:
            pp=QPainter(self.tempix)
            pp.drawRect(x,y,w,h)
            painter.begin(self)
            painter.drawPixmap(0,0,self.tempix)
        else:
            pp=QPainter(self.pix)
            pp.drawRect(x,y,w,h)
            painter.begin(self)
            painter.drawPixmap(0,0,self.pix)           
        '''                
        
    def getPixmap(self,ori_pix):
        self.pix=ori_pix
        self.tempix=ori_pix
        self.segment_status=False
        self.dir_status=False
        self.pos_xy=[]
        self.dir_set=[]

    def getFinalPixmap(self,ori_pix):
        self.pix=ori_pix
        self.tempix=ori_pix
        self.segment_status=False
        self.dir_status=False
        
    def mouseMoveEvent(self,event):
        if (event.buttons()):
            if(self.dir_status==True and self.segment_status==False):
                self.endPoint = event.pos()
            elif(self.dir_status==False and self.segment_status==True):
                pos_tmp = (event.pos().x(), event.pos().y())
                print("append",pos_tmp)
                self.pos_xy.append(pos_tmp)
                
            #self.update()            
        #print ('mouse move:(%d,%d)\n'%(event.pos().x(),event.pos().y()))
        #if self.if_mouse_press:
        #    my_dialog.move_point(event.pos().x(),event.pos().y())
    def mousePressEvent(self,event):
        if(event.button() == Qt.LeftButton):
            self.if_mouse_press = True
            self.update()
            
            if(self.dir_status==True and self.segment_status==False):
                self.lastPoint = event.pos()
                self.endPoint = self.lastPoint

                pos_start=(event.pos().x(),event.pos().y())
                self.dir_set.append(pos_start)
                
            elif(self.dir_status==False and self.segment_status==True):
                pos_tmp = (event.pos().x(), event.pos().y())
                print("append",pos_tmp)
                self.pos_xy.append(pos_tmp)                           
                
        #print ('mousePressEvent(%d,%d)\n'%(event.pos().x(),event.pos().y()))
        
        #my_dialog.move_point(event.pos().x(),event.pos().y())
 
    def mouseReleaseEvent(self,event):
        print("release")
        if event.button() == Qt.LeftButton:
            
            self.update()
            self.if_mouse_press = False

            if(self.dir_status==True and self.segment_status==False):
                self.endPoint = event.pos()
                pos_end=(event.pos().x(),event.pos().y())
                self.dir_set.append(pos_end)

            elif(self.dir_status==False and self.segment_status==True):
                pos_test = (-1, -1)
                print("append",pos_test)
                self.pos_xy.append(pos_test)
                
        #print ('mouseReleaseEvent(%d,%d)\n'%(event.pos().x(),event.pos().y()))
        

class MainDialog(QDialog):  
    def __init__(self,parent=None):  
        super(MainDialog,self).__init__(parent)

        #########################################
        self.hair_path=0 #
        self.mask_path=0 #

        self.cnts=0
        self.max_index=0
        #########################################
        self.hair_label = BtnLabel(self)  
        self.hair_label.setGeometry(160, 40,width,height)  
        self.mask_label = BtnLabel(self)  
        self.mask_label.setGeometry(800, 40, width,height)  
        #set open hair image button
        self.hair_btn = QtWidgets.QPushButton(self)  
        self.hair_btn.setObjectName("open_hair_btn")  
        self.hair_btn.setGeometry(1, 0, 100, 40)
        self.hair_btn.setText("open hair ")  
        self.hair_btn.clicked.connect(self.open_hair_image)
        
 
        #set open mask image button
        self.mask_btn = QtWidgets.QPushButton(self)  
        self.mask_btn.setObjectName("open_mask_btn")  
        self.mask_btn.setGeometry(1, 60, 100, 40)
        self.mask_btn.setText("open mask")  
        self.mask_btn.clicked.connect(self.open_mask_image)  
 
        #set apply mask button
        self.apply_btn = QtWidgets.QPushButton(self)  
        self.apply_btn.setObjectName("apply_mask_btn")  
        self.apply_btn.setGeometry(1, 120, 100, 40)
        self.apply_btn.setText("apply mask")  
        self.apply_btn.clicked.connect(self.apply_mask)  
 
        #set segment button
        self.segment_btn = QtWidgets.QPushButton(self)  
        self.segment_btn.setObjectName("segment_btn")  
        self.segment_btn.setGeometry(1, 180, 100, 40)
        self.segment_btn.setText("segment")  
        self.segment_btn.clicked.connect(self.segment)
        
        #set direct button
        self.direct_btn = QtWidgets.QPushButton(self)  
        self.direct_btn.setObjectName("direct_btn")  
        self.direct_btn.setGeometry(1, 240, 100, 40)
        self.direct_btn.setText("direct")  
        self.direct_btn.clicked.connect(self.direct)  
 
        #set save dir button
        self.save_btn = QtWidgets.QPushButton(self)  
        self.save_btn.setObjectName("save_btn")  
        self.save_btn.setGeometry(1, 300, 100, 40)
        self.save_btn.setText("save direction")  
        self.save_btn.clicked.connect(self.save_dir)  
 
        #set result text
        self.text_label = QLabel(self)
        self.text_label.setAlignment(Qt.AlignCenter)  
        self.text_label.setGeometry(1, 360, 140, 100)
        
    def open_hair_image(self):
        imgName, imgType = QFileDialog.getOpenFileName(self, "open image", "", "*.png;;*.jpg;;All Files(*)")
        self.hair_path=imgName
        pix = QtGui.QPixmap(imgName).scaled(self.hair_label.width(), self.hair_label.height())
        self.hair_label.getPixmap(pix)
        
    def open_mask_image(self):
        imgName, imgType = QFileDialog.getOpenFileName(self, "open image", "", "*.png;;*.jpg;;All Files(*)")
        self.mask_path=imgName
        pix = QtGui.QPixmap(imgName).scaled(self.mask_label.width(), self.mask_label.height())
        self.mask_label.getPixmap(pix)
        
    def apply_mask(self):
        hair_img=cv2.imread(self.hair_path)
        mask_img=cv2.imread(self.mask_path,0)
        ret,thresh=cv2.threshold(mask_img,150,255,0)
        self.cnts,hierarchy=cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
        #cv2.imshow('oringin image',thresh)
        #cv2.waitKey()
        cnts_size=np.size(self.cnts)
        print(np.size(self.cnts))
        max_size=-1
        max_index=-1
        for i in range(cnts_size):
            if(np.size(self.cnts[i])>max_size):
                max_index=i
                max_size=np.size(self.cnts[i])
        print(max_size,max_index)
        print(self.cnts[0])
        #print(hierarchy)
        self.max_index=max_index
        img=cv2.drawContours(hair_img,self.cnts,self.max_index,(0,0,255),3)

        #convert opencv image to qt image
        height, width, channel = img.shape
        bytesPerLine = 3 * width
        qImg = QImage(img.data, width, height, bytesPerLine, QImage.Format_RGB888)
        #convert qt image to qt pixmap
        qPix=QtGui.QPixmap.fromImage(qImg)
        self.hair_label.getPixmap(qPix)
        #print(img.shape)
        #cv2.imshow('image with contours',img)
        #cv2.waitKey()
        
    def segment(self):
        self.hair_label.set_segment(True)
        
    def direct(self):
        self.hair_label.set_dir(True)
        
    def save_dir(self):
        hair_img=cv2.imread(self.hair_path)
        hair_img[:,:,:]=0
        img=cv2.drawContours(hair_img,self.cnts,self.max_index,(0,0,255),1)
        #convert opencv image to qt image
        height, width, channel = img.shape
        bytesPerLine = 3 * width
        qImg = QImage(img.data, width, height, bytesPerLine, QImage.Format_RGB888)
        #convert qt image to qt pixmap
        qPix=QtGui.QPixmap.fromImage(qImg)

        pen = QPen(Qt.blue, 1, Qt.SolidLine)
        pp=QPainter(qPix)
        pp.setPen(pen)
        pos_size=len(self.hair_label.pos_xy)
        if(pos_size>=3):
            for i in range(0,pos_size-1):
                #print(i,"i")
                point_start = self.hair_label.pos_xy[i]
                point_end = self.hair_label.pos_xy[i+1]
                if point_end == (-1, -1) and i==pos_size-2:
                    break
                elif point_end==(-1,-1) and i!=pos_size-2:
                    i+=1 # no use
                elif point_start==(-1,-1):
                    continue
                else:
                    print(point_start,"start")
                    print(point_end,"end")
                    pp.drawLine(point_start[0], point_start[1], point_end[0], point_end[1])

        self.hair_label.getFinalPixmap(qPix)
        
        image = qPix.toImage()
        image = image.rgbSwapped()
        ptr = image.bits()
        ptr.setsize(image.byteCount())
        arr = np.array(ptr).reshape(height, width, 4)  #  Copies the data
        gray = cv2.cvtColor (arr,cv2.COLOR_BGR2GRAY)
        #cv2.imshow('gray',gray)
        #cv2.waitKey()
        ret,thresh=cv2.threshold(gray,10,255,0)
        #cv2.imshow('thresh',thresh)
        #cv2.waitKey()
        cnts_here,hierarchy=cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
        cnts_size=np.size(cnts_here)
        print(cnts_size)
        img=cv2.drawContours(hair_img,cnts_here,0,(0,255,0),2)
        #cv2.imshow('image with contours',img)
        #cv2.waitKey()

        dir=np.zeros((cnts_size,2),dtype=np.float32)
        sum=np.zeros((cnts_size))

        #data need to be saved
        print(height,width)
        dir_map=np.zeros((height,width,2),dtype=np.float32)
        whether_hair=np.zeros((height,width),dtype=np.float32)
        
        len_dir=len(self.hair_label.dir_set)
        for i in range(0,len_dir-1,2):
            for j in range(cnts_size):
                dist=cv2.pointPolygonTest(cnts_here[j],self.hair_label.dir_set[i],True)
                if(dist>0):
                    sum[j]=sum[j]+1
                    print(i,j)
        for j in range(cnts_size):
            for i in range(0,len_dir-1,2):
                if(sum[j]==1):
                    dist=cv2.pointPolygonTest(cnts_here[j],self.hair_label.dir_set[i],True)
                    if(dist>0):
                        print(i,j)
                        dir[j,0]=self.hair_label.dir_set[i+1][0]-self.hair_label.dir_set[i][0]
                        dir[j,1]=self.hair_label.dir_set[i+1][1]-self.hair_label.dir_set[i][1]
                        print(dir[j,0],dir[j,1])
                        
        for i in range(height): #625
            for j in range(width): #500
                for k in range(cnts_size):
                    if (sum[k]==1 and dir_map[i,j,0]==0 and dir_map[i,j,1]==0):
                        dist=cv2.pointPolygonTest(cnts_here[k],(j,i),True)
                        if(dist>0):
                            dir_map[i,j]=[dir[k,0],dir[k,1]]
                    elif (sum[k]>=2):
                        dist=cv2.pointPolygonTest(cnts_here[k],(j,i),True)
                        if(dist>0):
                            whether_hair[i,j]=200
        #cv2.imshow('whether_hair',whether_hair)
        #cv2.waitKey()
        #cv2.imshow('dir_map',dir_map)
        #cv2.waitKey()

        #for fix isolated hair problem
        whether_hair[553,461]=0
        # dir map is measure in image space
        np.savetxt(img_data_path+"/"+"dir_map_x.txt",dir_map[:,:,0])
        np.savetxt(img_data_path+"/"+"dir_map_y.txt",dir_map[:,:,1])
        np.savetxt(img_data_path+"/"+"whether_hair.txt",whether_hair)
        print("save annotate dir_map into file")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    my_dialog = MainDialog()
    my_dialog.show()
    sys.exit(app.exec_())
