#include "myvector.h"
#include "object.h"
#include "io.h"

using namespace std;
using namespace Eigen;

const static double help[8][3]= {
  {
    -1,-1,-1
  },
  {
    -1,-1,1
  },
  {
    -1,1,-1
  },
  {
    -1,1,1
  },
  {
    1,-1,-1
  },
  {
    1,-1,1
  },
  {
    1,1,-1
  },
  {
    1,1,1
  }
};

const static double quadrature[8][3]={
  {
    -0.57735, -0.57735, -0.57735
  },
  {
    -0.57735, -0.57735, 0.57735
  },
  {
    -0.57735, 0.57735, -0.57735
  },
  {
    -0.57735, 0.57735, 0.57735
  },
  {
    0.57735, -0.57735, -0.57735
  },
  {
    0.57735, -0.57735, 0.57735
  },
  {
    0.57735, 0.57735, -0.57735
  },
  {
     0.57735, 0.57735, 0.57735
  }
};

const static size_t for_face[6][4]={
  {
    0,4,1,5
  },
  {
    0,4,6,2
  },
  {
    6,2,3,7
  },
  {
    3,7,5,1
  },
  {
    5,4,6,7
  },
  {
    0,2,3,1
  }
};

object::object()
{
  while(!myvertexs.empty())
    {
      myvertexs.pop_back();
    }
  while(!myhexahedrons.empty())
    {
      myhexahedrons.pop_back();
    }
  while(!myfaces.empty())
    {
      myfaces.pop_back();
    }

  while(!mytetrahedrons.empty())
    {
      mytetrahedrons.pop_back();
    }
  while(!mytri_faces.empty())
    {
      mytri_faces.pop_back();
    }

while(!time_norm_pair.empty())
    {
      time_norm_pair.pop_back();
    }
  time_all=norm_Jacobian_cal=0;
  this->index_for_vertex=NULL;
  this->mapIndexToLocInMartix=NULL;
  this->mapLocInMatrixToIndex=NULL;
  center_loc=myvector(0,0,0);
  calShapeFuncGrad();
}

object::~object()
{
  // because the harmonic deformation do not need these two array, they can be null
  if(mapIndexToLocInMartix!=NULL)
    {
      delete[] mapIndexToLocInMartix;
    }
  if(mapLocInMatrixToIndex!=NULL)
    {
      delete[] mapLocInMatrixToIndex;
    }
  if(index_for_vertex!=NULL)
    {
      size_t i,j;
      for(i=0;i<lc;++i)
	{
	  for(j=0;j<wc;++j)
	    {
	      delete[] index_for_vertex[i][j];
	    }
	  delete[] index_for_vertex[i];
	}
      delete[] index_for_vertex;
    }
}
object::object(double dmetric,std::string input_dir,double dt,double density,int line_search,double weight_line_search,std::string constitutive_model)
{
  while(!myvertexs.empty())
    {
      myvertexs.pop_back();
    }
  while(!myhexahedrons.empty())
    {
      myhexahedrons.pop_back();
    }
  while(!myfaces.empty())
    {
      myfaces.pop_back();
    }

  while(!mytetrahedrons.empty())
    {
      mytetrahedrons.pop_back();
    }
  while(!mytri_faces.empty())
    {
      mytri_faces.pop_back();
    }
  
while(!time_norm_pair.empty())
    {
      time_norm_pair.pop_back();
    }
  time_all=norm_Jacobian_cal=0;

  mpFromIndexVertexToIndexFace.clear();
  converge=0;
  center_loc=myvector(0,0,0);
  this->dt=dt;
  this->density=density;
  this->line_search=line_search;
  this->weight_line_search=weight_line_search;
  this->constitutive_model=constitutive_model;
  this->dmetric=dmetric;
  
  this->index_for_vertex=NULL;
  this->mapIndexToLocInMartix=NULL;
  this->mapLocInMatrixToIndex=NULL;
  size_t i,j,k;
  size_t x,y,z;

  io myio=io();
  myio.getVertexAndHex(myvertexs,myhexahedrons,input_dir);
  num_vertex=myvertexs.size();
  printf("num_vertex ::%u \n",num_vertex);
  num_hexahedrons=myhexahedrons.size();
  printf("num_hex::%u\n",num_hexahedrons);
  
  num_all_dof=3*num_vertex;
  for(i=0;i<num_vertex;++i)
    {
      center_loc+=myvertexs[i].location_original;
    }
  center_loc/=num_vertex;
  size_t index_vertex_now[2][2][2];
  size_t index_vertex_for_face_now[4];  //for  each face, there are four vertexs
  size_t index_face[6]; //当前六个面的index
  //traverse all the grid
  size_t a,b,c;
  for(i=0;i<num_hexahedrons;++i)
    {
      for(a=0;a<6;++a)
	{
	  for(b=0;b<4;++b)
	    {
	      x=for_face[a][b]/4;
	      y=(for_face[a][b]%4)/2;
	      z=(for_face[a][b]%4)%2;
	      index_vertex_for_face_now[b]=myhexahedrons[i].index_vertex[x][y][z];
	    }
	  sort(index_vertex_for_face_now,index_vertex_for_face_now+4);
	  if(mpFromIndexVertexToIndexFace.find(make_pair(index_vertex_for_face_now[0],make_pair(index_vertex_for_face_now[1],index_vertex_for_face_now[2])))==mpFromIndexVertexToIndexFace.end())
	    {
	      index_face[a]=myfaces.size();
	      mpFromIndexVertexToIndexFace[make_pair(index_vertex_for_face_now[0],make_pair(index_vertex_for_face_now[1],index_vertex_for_face_now[2]))]=index_face[a];
	      myfaces.push_back(face(index_vertex_for_face_now));
	      myfaces[index_face[a]].num_hex=1;
	      // cal normal
	      myvector v1=myvertexs[index_vertex_for_face_now[0]].location_original-myvertexs[index_vertex_for_face_now[1]].location_original;
	      myvector v2=myvertexs[index_vertex_for_face_now[0]].location_original-myvertexs[index_vertex_for_face_now[2]].location_original;
	      myvector v3=v1.cross(v2); v3.normalize();
	      myfaces[index_face[a]].normal_ori=v3;
	      myvector v4=myvertexs[index_vertex_for_face_now[0]].location_original-center_loc;
	      if(v4.dot(v3)<0.0)
		{
		  myfaces[index_face[a]].normal_ori*=-1;
		}
	    }
	  else
	    {
	      int index_thisface=mpFromIndexVertexToIndexFace[make_pair(index_vertex_for_face_now[0],make_pair(index_vertex_for_face_now[1],index_vertex_for_face_now[2]))];
	      index_face[a]=index_thisface;
	      myfaces[index_face[a]].num_hex=2;
	    }
	}
    }	      
  num_faces=myfaces.size();
  printf("num_face::%u\n",num_faces);

  prepare();
  
}
object::object( double dmetric, size_t l_size, size_t w_size, size_t h_size, double dt, double density,int line_search, double weight_line_search, std::string constitutive_model)
{
  while(!myvertexs.empty())
    {
      myvertexs.pop_back();
    }
  while(!myhexahedrons.empty())
    {
      myhexahedrons.pop_back();
    }
  while(!myfaces.empty())
    {
      myfaces.pop_back();
    }

  while(!mytetrahedrons.empty())
    {
      mytetrahedrons.pop_back();
    }
  while(!mytri_faces.empty())
    {
      mytri_faces.pop_back();
    }
  
while(!time_norm_pair.empty())
    {
      time_norm_pair.pop_back();
    }
  time_all=norm_Jacobian_cal=0;

  mpFromIndexVertexToIndexFace.clear();
  converge=0;
  center_loc=myvector(0,0,0);
  
  this->dt=dt;
  this->density=density;
  this->line_search=line_search;
  this->weight_line_search=weight_line_search;
  this->constitutive_model=constitutive_model;
  length=dmetric*l_size; width=dmetric*w_size; height=dmetric*h_size; 
  this->dmetric=dmetric;
  
  this->index_for_vertex=NULL;
  this->mapIndexToLocInMartix=NULL;
  this->mapLocInMatrixToIndex=NULL;
  double length_now,width_now,height_now;
  size_t i,j,k;
  size_t x,y,z;

  lc=l_size+1;
  wc=w_size+1;
  hc=h_size+1;
  
  // 分配内存空间 保存顶点索引
  index_for_vertex=(size_t***)new size_t**[lc];
  for(i=0;i<lc;++i)
    {
      index_for_vertex[i]=(size_t**)new size_t*[wc];
      for(j=0;j<wc;++j)
	{
	  index_for_vertex[i][j]=new size_t[hc];
	}
    }

  size_t index_now_vertex=0;
  //set index for all vertexs and EPS is to make double`error correction and set it as local constant
  double EPS=1e-6;
  for(i=0,length_now=-length/2.0;length_now<length/2.0+EPS;length_now+=dmetric,++i)
    {
      for(j=0,width_now=-width/2.0;width_now<width/2.0+EPS;width_now+=dmetric,++j)
	{
	  for(k=0,height_now=-height/2.0;height_now<height/2.0+EPS;height_now+=dmetric,++k)
	    {
	      index_for_vertex[i][j][k]=index_now_vertex;
	      myvertexs.push_back(vertex(myvector(length_now,width_now,height_now)));
	      index_now_vertex++;
	    }
	}
    }

  num_vertex=myvertexs.size();
  num_all_dof=3*num_vertex;
  
  printf("num_vertex ::%u \n",num_vertex);

  for(i=0;i<num_vertex;++i)
    {
      center_loc+=myvertexs[i].location_original;
    }
  center_loc/=num_vertex;
  size_t index_vertex_now[2][2][2];
  size_t index_vertex_for_face_now[4];  //for  each face, there are four vertexs
  size_t index_face[6]; //当前六个面的index
  //traverse all the grid
  size_t a,b,c;
   for(i=0;i<lc-1;++i)
    {   
      for(j=0;j<wc-1;++j)
	{
	  for(k=0;k<hc-1;++k)
	    {  	 
	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  index_vertex_now[a][b][c]=index_for_vertex[i+a][j+b][k+c];
			}
		    }
		}
	      myhexahedrons.push_back(hexahedron(index_vertex_now));
	      for(a=0;a<6;++a)
		{
		  for(b=0;b<4;++b)
		    {
		      x=for_face[a][b]/4;
		      y=(for_face[a][b]%4)/2;
		      z=(for_face[a][b]%4)%2;
		      index_vertex_for_face_now[b]=index_vertex_now[x][y][z];
		    }
		  sort(index_vertex_for_face_now,index_vertex_for_face_now+4);
		  if(mpFromIndexVertexToIndexFace.find(make_pair(index_vertex_for_face_now[0],make_pair(index_vertex_for_face_now[1],index_vertex_for_face_now[2])))==mpFromIndexVertexToIndexFace.end())
		    {
		      index_face[a]=myfaces.size();
		      mpFromIndexVertexToIndexFace[make_pair(index_vertex_for_face_now[0],make_pair(index_vertex_for_face_now[1],index_vertex_for_face_now[2]))]=index_face[a];
		      myfaces.push_back(face(index_vertex_for_face_now));
		      myfaces[index_face[a]].num_hex=1;
		      // cal normal
		      myvector v1=myvertexs[index_vertex_for_face_now[0]].location_original-myvertexs[index_vertex_for_face_now[1]].location_original;
		      myvector v2=myvertexs[index_vertex_for_face_now[0]].location_original-myvertexs[index_vertex_for_face_now[2]].location_original;
		      myvector v3=v1.cross(v2); v3.normalize();
		      myfaces[index_face[a]].normal_ori=v3;
		      myvector v4=myvertexs[index_vertex_for_face_now[0]].location_original-center_loc;
		      if(v4.dot(v3)<0.0)
			{
			  myfaces[index_face[a]].normal_ori*=-1;
			}
       		    }
		  else
		    {
		      int index_thisface=mpFromIndexVertexToIndexFace[make_pair(index_vertex_for_face_now[0],make_pair(index_vertex_for_face_now[1],index_vertex_for_face_now[2]))];
		      index_face[a]=index_thisface;
		      myfaces[index_face[a]].num_hex=2;
       		    }
		}
	    }    
	}
    }
  num_hexahedrons=myhexahedrons.size();
  num_faces=myfaces.size();
  printf("num_hex::%u\n",num_hexahedrons);
  printf("num_face::%u\n",num_faces);
  
  prepare();
}

int object::prepare()
{
  calShapeFuncGrad();
  init_Energy_now_ForHex();
  calMassForVertex();

  if(constitutive_model=="linear"||constitutive_model=="linear_with_stiffness_tensor")
    {
      line_search=0;
      max_iteration_num=1;
    }
  else if(constitutive_model=="co_rotated_linear"||constitutive_model=="co_rotated_linear_with_stiffness_tensor")
    {
      max_iteration_num=40;
    }
  else if(constitutive_model=="stvk"||constitutive_model=="neo_hookean"||constitutive_model=="stvk_with_stiffness_tensor")
    {
      // do not change line_search
    }
  return 0;
}
int object::checkFixedOrFree()
{
  int i,j,k,a,b,c;  // 这里使用int因为如果是size_t的话使用abs()函数有问题
  num_cal_dof=0;
  mapIndexToLocInMartix=new size_t[num_all_dof];
  mapLocInMatrixToIndex=new size_t[num_all_dof];
  for(i=0;i<num_vertex;++i)
    {
      if(myvertexs[i].isFixed==0)
	{
	  for(j=0;j<3;++j)
	    {
	      mapIndexToLocInMartix[i*3+j]=num_cal_dof;
	      mapLocInMatrixToIndex[num_cal_dof]=i*3+j;
	      ++num_cal_dof;
	    }
	}
    }
  return 0;
}

int object::calShapeFuncGrad()
{
  size_t i,j,k;
  size_t whichShapeFunc;
  size_t a,b,c;
  size_t whichQuadrature;
  for(i=0;i<2;++i)
    {
      for(j=0;j<2;++j)
	{
	  for(k=0;k<2;++k)
	    {
	      whichShapeFunc=4*i+2*j+k;
	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  whichQuadrature=4*a+2*b+c;
			  
			  shapeFuncGrad[i][j][k][a][b][c][0]=help[whichShapeFunc][0]*(1+help[whichShapeFunc][1]*quadrature[whichQuadrature][1])*(1+help[whichShapeFunc][2]*quadrature[whichQuadrature][2])*0.125;
			  shapeFuncGrad[i][j][k][a][b][c][1]=help[whichShapeFunc][1]*(1+help[whichShapeFunc][0]*quadrature[whichQuadrature][0])*(1+help[whichShapeFunc][2]*quadrature[whichQuadrature][2])*0.125;
			  shapeFuncGrad[i][j][k][a][b][c][2]=help[whichShapeFunc][2]*(1+help[whichShapeFunc][1]*quadrature[whichQuadrature][1])*(1+help[whichShapeFunc][0]*quadrature[whichQuadrature][0])*0.125;
			      
			}
		    }
		}
	    }
	}
    }
  return 0;
}

int object::init_Energy_now_ForHex()
{
  size_t i;
  for(i=0;i<num_hexahedrons;++i)
    {
      myhexahedrons[i].energy_now=0;
    }
  return 0;
}

int object::calMassForVertex()
{
  size_t i;
  size_t a,b,c;
  for(i=0;i<num_vertex;++i)
    {
      myvertexs[i].mass=0;
    }
  double vol=pow(dmetric,3.0);
  for(i=0;i<num_hexahedrons;++i)
    {
      for(a=0;a<2;++a)
	{
	  for(b=0;b<2;++b)
	    {
	      for(c=0;c<2;++c)
		{
		  myvertexs[myhexahedrons[i].index_vertex[a][b][c]].mass+=(0.125*vol*density);
		}
	    }
	}
    }
  return 0;
}

double object::calElasticEnergy()
{
  size_t i;
  double elasticE=0;
  for(i=0;i<num_hexahedrons;++i)
    {
      elasticE+=myhexahedrons[i].energy_now;
    }
  return elasticE;
}

int object::dynamicSimulator()
{
  size_t i;
  //printf("num_fixed: %u\n",num_fixed);
  
  MatrixXd Hessian=MatrixXd::Random(num_all_dof,num_all_dof);
  VectorXd Jacobian(num_all_dof);
  iteration_num=0;
  converge=0;
  for(i=0;i<num_vertex;++i)
    {
      myvertexs[i].location_lastFrame=myvertexs[i].location;
      myvertexs[i].velocity_lastFrame=myvertexs[i].velocity;
    }
  clock_t start,finish;
  double totaltime;
  while(!converge)
    { 
      Hessian.fill(0);
      Jacobian.fill(0);
      
      printf("--[Inf]:calJacobianAndHessianForHex\n");
      start=clock();
      calJacobianAndHessianForHex(Hessian,Jacobian);
      finish=clock();
      totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
      time_all+=totaltime; 
      printf("Assemble Time Cost: %lf\n",totaltime);

      
      printf("--[Inf]:solve Matrix\n");
      start=clock();
      //  printf("here converge before solve: %d\n",converge);
      solve(Hessian,Jacobian);
      finish=clock();
      totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
      time_all+=totaltime;
      printf("Solve Time Cost: %lf\n",totaltime);
      //  printf("here converge after solve: %d\n",converge);
    }
  return 0;
}

int object::harmonic_def_static()
{
  if(constitutive_model!="linear")
    {
      return 1; //nonlinear strain can not be used  
    }
  MatrixXd Hessian=MatrixXd::Random(num_all_dof,num_all_dof); // 3n dof 
  VectorXd Jacobian(num_all_dof);

  iteration_num=0;
  converge=0;
  //解static 只有一帧
  while(!converge)
    {
      Hessian.fill(0);
      Jacobian.fill(0);
       printf("--[Inf]:calJacobianAndHessianForHex\n");
      calJacobianAndHessianForHex(Hessian,Jacobian);
       printf("--[Inf]:solve Matrix for harmonic static deformation\n");
      solve_harmonic_def_static(Hessian,Jacobian);
    }
  return 0;
}
int object::calJacobianAndHessianForHex(Eigen::MatrixXd &Hessian,Eigen::VectorXd &Jacobian)
{
  size_t i,j;   
  for(i=0;i<num_hexahedrons;++i)
    {
      myhexahedrons[i].calJacobianAndHessian(shapeFuncGrad,myvertexs,dmetric,Hessian,Jacobian,constitutive_model);
    }

  // diff by hand 
  {
     double factor=dmetric*dmetric*0.25;
     Hessian*=factor;
     Jacobian*=(factor*-1);
   }
  
  //auto diff
  // Jacobian*=-1;
  
  for(i=0;i<num_vertex;++i)
    {
      for(j=0;j<3;++j)
	{
	  Jacobian(i*3+j)+=myvertexs[i].force_external(j); 
	}
    }

  double EPS=1e-5; // local constant to judge zero for dt-100
  double d1dt;
  if(fabs(dt-100)<EPS) 
    {
      printf("--[INF]:: No mass matrix\n");
      d1dt=0;
    }
  else
    {
      d1dt=1.0/dt;
    }
  double mdt,mdtdt;
  myvector vmdt;
  for(i=0;i<num_vertex;++i)
    {
      mdt=myvertexs[i].mass*d1dt;
      mdtdt=mdt*d1dt;
      vmdt=(myvertexs[i].velocity_lastFrame-(myvertexs[i].location-myvertexs[i].location_lastFrame)*d1dt)*mdt;
      for(j=0;j<3;++j)
	{
	  Hessian(i*3+j,i*3+j)+=mdtdt;
	  Jacobian(i*3+j)+=vmdt(j);
	}
    }
   
  return 0;
}
 
bool object::checkInversion(Eigen::VectorXd &dx,Eigen::VectorXd &Jacobian)
{
  size_t i;
  size_t x,y;
  Jacobian.fill(0);
  for(i=0;i<num_cal_dof;++i)
    {
      y=mapLocInMatrixToIndex[i]%3;
      x=mapLocInMatrixToIndex[i]/3;
      myvertexs[x].location_maybe(y)=dx(i)+myvertexs[x].location(y);
    }
  for(i=0;i<num_hexahedrons;++i)
    {         
      if(myhexahedrons[i].checkInversion(shapeFuncGrad,myvertexs,dmetric,Jacobian,constitutive_model))
	{
	  return true;
	}      
    }  
  return false;
}

bool object::checkInversion(Eigen::VectorXd &dx)
{
  size_t i;
  size_t x,y;
  for(i=0;i<num_cal_dof;++i)
    {
      y=mapLocInMatrixToIndex[i]%3;
      x=mapLocInMatrixToIndex[i]/3;
      myvertexs[x].location_maybe(y)=dx(i)+myvertexs[x].location(y);
    }
  for(i=0;i<num_hexahedrons;++i)
    {         
      if(myhexahedrons[i].checkInversion(shapeFuncGrad,myvertexs,dmetric,constitutive_model))
	{
	  return true;
	}      
    }  
  return false;
}

bool object::checkInversion_for_harmonic_def(Eigen::VectorXd &dx)
{
  size_t i;
  size_t x,y;
  for(i=0;i<num_all_dof;++i)
    {
      y=i%3;
      x=i/3;
      myvertexs[x].location_maybe(y)=dx(i)+myvertexs[x].location(y);
    }
  for(i=0;i<num_hexahedrons;++i)
    {
      if(myhexahedrons[i].checkInversion(shapeFuncGrad,myvertexs,dmetric,constitutive_model))
	{
	  return true;
	}
    }  
  return false;
}

double object::calEnergyDif()
{
  size_t i;
  double energy_old,energy_new;
  energy_new=energy_old=0;
  for(i=0;i<num_hexahedrons;++i)
    {
      energy_new+=myhexahedrons[i].energy_maybe;
      energy_old+=myhexahedrons[i].energy_now;
    }

  printf("the first part old energy:%lf new energy:%lf  \n",energy_old,energy_new);
  double EPS=1e-5; // local constant to judge zero for dt-100
  double d1dt;
  if(fabs(dt-100)<EPS) 
    {  
      d1dt=0;
    }
  else
    {
      d1dt=1.0/dt;
    }
  
  myvector help1,help2;
  for(i=0;i<num_vertex;++i)
    {
      help1=(myvertexs[i].location_maybe-myvertexs[i].location_lastFrame)*d1dt-myvertexs[i].velocity_lastFrame;
      help2=(myvertexs[i].location-myvertexs[i].location_lastFrame)*d1dt-myvertexs[i].velocity_lastFrame;
      
      energy_new+=help1.len_sq()*myvertexs[i].mass*0.5;
      energy_old+=help2.len_sq()*myvertexs[i].mass*0.5;
      
      energy_new+=-1*myvertexs[i].location_maybe.dot(myvertexs[i].force_external);
      energy_old+=-1*+myvertexs[i].location.dot(myvertexs[i].force_external);
    }
   printf("old energy:%lf new energy:%lf \n",energy_old,energy_new);
  return energy_new-energy_old;
}

int object::solve_harmonic_def_static(Eigen::MatrixXd &Hessian,Eigen::VectorXd &Jacobian)
{
  size_t row=num_all_dof+6;
  size_t col=row;
  SparseMatrix<double > Hessianspa(row,col);
  VectorXd Jacobian_cal(row),dx(row);
  vector< Triplet<double > > tripletsForHessianspa;

  SparseLU<SparseMatrix<double>,COLAMDOrdering<int>> linearSolver;//拉格朗日乘子系统是不定的，有负特征根，llt ldlt分解都不可以用，ldlt要求矩阵正定　llt要求矩阵半正定　lu适用于 ambitrary　matrix

  //共轭梯度法求解线性方程组从原先的适用于正定matrix拓展到适用于ambitrary matrix
  ConjugateGradient<SparseMatrix<double> > cg;
  cg.setMaxIterations(50);
  
  size_t i,j;

  double EPS=1e-10; // local constant to judge zero for K(i,j)

  for(i=0;i<num_all_dof;++i)
    {
      for(j=0;j<num_all_dof;++j)
	{
	  if(fabs(Hessian(i,j))>= EPS)
	    {
	      tripletsForHessianspa.emplace_back(i,j,Hessian(i,j));
	    }
	}
    }
  for(i=0;i<num_all_dof;++i)
    {
      Jacobian_cal(i)=Jacobian(i);
    }

  for(i=0;i<6;++i)
    {
      Jacobian_cal(num_all_dof+i)=0; // because the constraint is satisfied at the beginning 
    }
    
  MatrixXd Jacobian_for_constraint=MatrixXd::Random(6,num_all_dof);
  Jacobian_for_constraint.fill(0);
  // build KKT system
  {
    for(i=0;i<3;++i)
      {
	for(j=0;j<num_vertex;++j)
	  {
	    Jacobian_for_constraint(i,j*3+i)=1;
	  }
      }
    size_t up,down;
    for(i=0;i<3;++i)
      {
	for(j=0;j<num_vertex;++j)
	  {
	    up=(i+2)%3; down=(i+1)%3;
	    Jacobian_for_constraint(i+3,j*3+down)=myvertexs[j].location_original(up)-center_loc(up);
	    Jacobian_for_constraint(i+3,j*3+up)=center_loc(down)-myvertexs[j].location_original(down);
	  }
      }
  }
  for(i=0;i<6;++i)
    {
      for(j=0;j<num_all_dof;++j)
	{
	   if(fabs(Jacobian_for_constraint(i,j))>= EPS)
	    {
	      tripletsForHessianspa.emplace_back(i+num_all_dof,j,Jacobian_for_constraint(i,j));
	      tripletsForHessianspa.emplace_back(j,i+num_all_dof,Jacobian_for_constraint(i,j));
	    }
	}
    }
  
  Hessianspa.setFromTriplets(tripletsForHessianspa.begin(),tripletsForHessianspa.end());
  Hessianspa.makeCompressed();

  //共轭梯度法
  /*
  {
    cg.compute(Hessianspa);
    dx=cg.solve(Jacobian_cal);
  }
  */
  
  //LU分解法
  {
    linearSolver.compute(Hessianspa);
    dx=linearSolver.solve(Jacobian_cal);
  }
  
  checkInversion_for_harmonic_def(dx);

  for(i=0;i<num_vertex;++i)
    {
      myvertexs[i].location=myvertexs[i].location_maybe;
    }
  for(i=0;i<num_hexahedrons;++i)
    {
      myhexahedrons[i].energy_now=myhexahedrons[i].energy_maybe;
    }
  converge=1;
  return 0;
}

int object::solve(Eigen::MatrixXd &Hessian,Eigen::VectorXd &Jacobian)
{
  size_t row=num_cal_dof;
  size_t col=row;
  SparseMatrix<double > Hessianspa(row,col);
  VectorXd Jacobian_cal(row),dx(row),dx_now(row);
  vector< Triplet<double > > tripletsForHessianspa;

  SimplicialLLT<SparseMatrix<double>> linearSolver;
   
  size_t i,j;
  size_t ii,jj;

  double EPS=1e-10; // local constant to judge zero for K(i,j)

  for(i=0;i<row;++i)
    {
      for(j=0;j<col;++j)
	{
	  if(fabs(Hessian(mapLocInMatrixToIndex[i],mapLocInMatrixToIndex[j]))>= EPS)
	    {
	      tripletsForHessianspa.emplace_back(i,j,Hessian(mapLocInMatrixToIndex[i],mapLocInMatrixToIndex[j]));
	    }
	}
    }
  for(i=0;i<row;++i)
    {
      Jacobian_cal(i)=Jacobian(mapLocInMatrixToIndex[i]);
    }

  double norm_Jacobian_cal=0;
  for(i=0;i<row;++i)
    {
      norm_Jacobian_cal+=(Jacobian_cal(i)*Jacobian_cal(i));
    }
  norm_Jacobian_cal=sqrt(norm_Jacobian_cal);
  printf("norm_Jacobian_cal: %lf\n",norm_Jacobian_cal);  
  time_norm_pair.push_back(make_pair(time_all,norm_Jacobian_cal));

  //正定化 当hessian负定时
  SparseMatrix<double > diagMatrix_spa(row,col);
  vector< Triplet<double > > tripletsFordiagMatrix_spa;

  for(i=0;i<row;i++)
    {
      tripletsFordiagMatrix_spa.emplace_back(i,i,1e-8);
    }
  diagMatrix_spa.setFromTriplets(tripletsFordiagMatrix_spa.begin(),tripletsFordiagMatrix_spa.end());
  diagMatrix_spa.makeCompressed();
  Hessianspa.setFromTriplets(tripletsForHessianspa.begin(),tripletsForHessianspa.end());
  Hessianspa.makeCompressed();

  while(1)
    {
      linearSolver.compute(Hessianspa);
      int info=(int)linearSolver.info();
      // cout<<info<<"info"<<endl;
      if(info==0)
	{
	  break;
	}
      else if(info==1)
	{
	  Hessianspa=Hessianspa+diagMatrix_spa;
	  diagMatrix_spa=diagMatrix_spa*2;
	}
    }
// cout<<Eigen::Success<<endl;
  dx=linearSolver.solve(Jacobian_cal);
  
  double h=2;
  size_t max_lineSearch=30;
  double energy_dif,threshold;
  bool find=0;

  VectorXd Jacobian_loc_maybe(num_all_dof);
  VectorXd Jacobian_cal_loc_maybe(row);
  
  if(line_search==1)
    {
      //printf("[control info]:: line_search open\n");
      for(i=0;i<max_lineSearch;++i)
	{
	  h*=0.5;
	  dx_now=dx*h;
	  if(checkInversion(dx_now,Jacobian_loc_maybe)==1)
	    {
	      continue;
	    }
	  else
	    {
	      // diff by hand
	      {
		double factor=dmetric*dmetric*0.25;
		Jacobian_loc_maybe*=(factor*-1);
	      }  
	      for(ii=0;ii<num_vertex;++ii)
		{
		  for(jj=0;jj<3;++jj)
		    {
		      Jacobian_loc_maybe(ii*3+jj)+=myvertexs[ii].force_external(jj); 
		    }
		}
	      for(ii=0;ii<row;++ii)
		{
		  Jacobian_cal_loc_maybe(ii)=Jacobian_loc_maybe(mapLocInMatrixToIndex[ii]);
		}

	      double Jdx=dx.dot(Jacobian_cal_loc_maybe);
	      double Jdx_now=dx.dot(Jacobian_cal);
	      energy_dif=calEnergyDif();
	      threshold=weight_line_search*dx_now.dot(Jacobian_cal)*-1;
	      if((energy_dif<=threshold)/*&&(fabs(Jdx)<=fabs(Jdx_now))*/)
		{
		  find=1;
		  break;
		}
	    }
	}
    }
  else if(line_search==0)
    {
      printf("[control info]:: line_search closed\n");
      h=1;
      dx_now=dx*h;
// 即使这里没有使用line search 但是为了保证鲁棒性：程序不是因为nan 等原因终止，因此仍插入可能的line search
      if(checkInversion(dx_now)==1)
	{
	  for(i=0;i<max_lineSearch;++i)
	    {
	      h*=0.5;
	      dx_now=dx*h;
	      if(checkInversion(dx_now)==1)
		{
		  continue;
		}
	      else
		{
		  energy_dif=calEnergyDif();
		  printf("energy_dif %lf\n",energy_dif);
		  threshold=weight_line_search*dx_now.dot(Jacobian_cal)*-1;
		  if(energy_dif<=threshold)
		    {
		      find=1;
		      break;
		    }
		}
	    }
	}
      find=1;  
      energy_dif=calEnergyDif();
    }  
      
  EPS=1e-4; // local constant to judge zero for energy_dif
  if(find==1)
    {
// fabs(energy_dif)<=EPS 即走不下去时为终止条件  ； fabs(norm_Jacobian_cal)<=EPS 即梯度很小时为终止条件
      if(fabs(norm_Jacobian_cal)<=EPS&&(constitutive_model=="stvk"||constitutive_model=="neo_hookean"||constitutive_model=="stvk_with_stiffness_tensor"||constitutive_model=="co_rotated_linear"||constitutive_model=="co_rotated_linear_with_stiffness_tensor"))
	{
          // printf("bingo\n");
	  converge=1;
	}
      for(i=0;i<num_vertex;++i)
	{
	  myvertexs[i].location=myvertexs[i].location_maybe;
	}
      for(i=0;i<num_hexahedrons;++i)
	{
	  myhexahedrons[i].energy_now=myhexahedrons[i].energy_maybe;
	}
    }
  else if(find==0)
    {
      converge=1;
    }

  ++iteration_num;
  if(iteration_num==max_iteration_num&&(constitutive_model=="linear"||constitutive_model=="co_rotated_linear"||constitutive_model=="linear_with_stiffness_tensor"||constitutive_model=="co_rotated_linear_with_stiffness_tensor"))
    {
      converge=1;
    }
  //co_rotated 既可以用 line_search 也可以不用 line_search; 不用line_search到达最大迭代次数仍会被认为是收敛　因此会造成翻转以后能量值巨大未收敛却被认为是收敛的情况
  if(converge==1)
    {
      double EPS=1e-5; // local constant to judge zero for dt-100
      double d1dt;
      if(fabs(dt-100)<EPS) 
	{
	  d1dt=0;
	}
      else
	{
	  d1dt=1.0/dt;
	}
      for(i=0;i<num_vertex;++i)
	{
	  //注意这里是更新速度,而不是累加速度,原因在于我们是用当前计算的这一帧的location减去前一帧的location,这里是经过了dt时间location的变化,因此获得的是一个位置对时间的变化率,因此是速度
	  //而需要注意的是在SG 2012 course中它在计算当前这一帧的收敛时,X(k+1)-x(k)=dt(v(k+1)-v(k))=> dx=dtdv 这里的dv是在这一帧的迭代过程的速度随时间的矫正率,而非速度
	  //因此这里不是+=而是=
	  myvertexs[i].velocity=(myvertexs[i].location-myvertexs[i].location_lastFrame)*d1dt;
	}
    }
  
  return 0;
}
