#include "head.h"
#include "io.h"
#include "autodiff.h"
using namespace std;
using namespace Eigen;
using namespace cv;
DECLARE_DIFFSCALAR_BASE();
void readcmdline(int argc, char* argv[],boost::property_tree::ptree &para_tree)
{
  size_t i;
  for(i=1;i<argc;++i)
    {
      string para_here=argv[i];
      size_t pos=para_here.find("=");
      if(pos!= string::npos)
	{
	  string key=para_here.substr(0,pos);
	  string value=para_here.substr(pos+1);
	  para_tree.put(key+".value",value);
	  printf("--[cmdline para] %s %s \n",key.c_str(),value.c_str());
	}
    }
  return;
}
int main(int argc, char *argv[])
{
  // test access point
  boost::property_tree::ptree para_tree;
  readcmdline(argc,argv,para_tree);
  string sparse_theta_map_path=para_tree.get<string>("sparse_theta_map_path.value");
  string dense_theta_map_path=para_tree.get<string>("dense_theta_map_path.value");
  string dense_theta_map_py_path=para_tree.get<string>("dense_theta_map_py_path.value");
 
  int height=para_tree.get<int>("height.value");
  int width=para_tree.get<int>("width.value");

  MatrixXd sparse_theta_map=MatrixXd::Random(height,width);
  MatrixXd dense_theta_map=MatrixXd::Random(height,width);
  
  MatrixXd whether_hair=MatrixXd::Random(height,width);
  
  io* my_io=new io();
  my_io->getData(sparse_theta_map,whether_hair,height,width,sparse_theta_map_path);

  int i,j,a,b;
  int aa,bb;
  int row,col;
  int ct,ct2;
  int num_hair=0;
  

  int index[height][width];
  double dir[height][width][2];
  for(i=0;i<height;++i)
    {
      for(j=0;j<width;++j)
	{
	  if(whether_hair(i,j)>10)
	    {
	      // trace here
	      if(i==341&&j==358)
		{
		  printf("theta %lf\n",sparse_theta_map(i,j));
		}
	      index[i][j]=num_hair;
	      num_hair++;
	      if(sparse_theta_map(i,j)>-0.0001)
		{
		  dir[i][j][0]=cos(sparse_theta_map(i,j));
		  dir[i][j][1]=sin(sparse_theta_map(i,j));

		}
	      else
		{
		  dir[i][j][0]=cos(0);
		  dir[i][j][1]=sin(0);
		}
	    }
	}
    }

  printf("num_hair:%d\n",num_hair);
  int num_hair_2=num_hair*2;

  vector<Triplet<double> > tripletsForHessian;
  while(!tripletsForHessian.empty())
    {
      tripletsForHessian.pop_back();
    }
  VectorXd Jacobian(num_hair_2);
  Jacobian.fill(0);

  typedef DScalar2<double,VectorXd, MatrixXd> DScalar; // use DScalar2 for calculating gradient and hessian and use DScalar1 for calculating gradient
  for(i=0;i<height;++i)
    {
      for(j=0;j<width;++j)
	{
	  if(whether_hair(i,j)>10)
	    {
	      int num_hair_now=0;

	      int index_center;
	      for(a=-1;a<2;a++)
		{
		  for(b=-1;b<2;b++)
		    {
		      if(whether_hair(i+a,j+b)>10)
			{
			  if(a==0&&b==0)
			    {
			      index_center=num_hair_now;
			    }
			  num_hair_now++;
			}		      
		    }
		}
	      VectorXd x(num_hair_now*2);

	      int index_local=0;
	      for(a=-1;a<2;a++)
		{
		  for(b=-1;b<2;b++)
		    {
		      if(whether_hair(i+a,j+b)>10)
			{
			  x(index_local*2+0)=dir[i+a][j+b][0];
			  x(index_local*2+1)=dir[i+a][j+b][1];
			  index_local++;
			}
		    }
		}

	      DiffScalarBase::setVariableCount(num_hair_now*2);
	      DScalar x_d[num_hair_now*2];
	      DScalar energy_d=DScalar(0);
	      for(a=0;a<num_hair_now*2;a++)
		{
		  x_d[a]=DScalar(a,x(a));
		}
	      DScalar dif_d[2]; dif_d[0]=DScalar(0); dif_d[1]=DScalar(0);
	      for(a=0;a<num_hair_now;a++)
		{
		  if(a!=index_center)
		    {
		      for(b=0;b<2;b++)
			{
			  dif_d[b]+=(x_d[a*2+b]-x_d[index_center*2+b]);
			}
		    }
		}
	      if(num_hair_now==1)
		{
		  printf("isolated hair %d %d\n",i,j);
		}
	      DScalar dif_1_d[2]; dif_1_d[0]=dif_d[0]/(num_hair_now-1); dif_1_d[1]=dif_d[1]/(num_hair_now-1);
	      energy_d= dif_1_d[0]*dif_1_d[0]+dif_1_d[1]*dif_1_d[1];
	      
	      if(sparse_theta_map(i,j)>-0.0001)
		{
		  DScalar dif_dir_d[2];
		  dif_dir_d[0]=x_d[index_center*2+0]-dir[i][j][0];
		  dif_dir_d[1]=x_d[index_center*2+1]-dir[i][j][1];
		  energy_d+=0.1*(dif_dir_d[0]*dif_dir_d[0]+dif_dir_d[1]*dif_dir_d[1]);
		}
	      MatrixXd grad(num_hair_now*2,1);
	      MatrixXd hes(num_hair_now*2,num_hair_now*2);
	      grad=energy_d.getGradient();
	      hes=energy_d.getHessian();

	      ct=0;

	      for(a=-1;a<2;a++)
		{
		  for(b=-1;b<2;b++)
		    {
		      if(whether_hair(i+a,j+b)>10)
			{
			  for(row=0;row<2;++row)
			    {
			      Jacobian(index[i+a][j+b]*2+row)+=grad(ct*2+row);
			      ct2=0;
			      for(aa=-1;aa<2;aa++)
				{
				  for(bb=-1;bb<2;bb++)
				    {
				      if(whether_hair(i+aa,j+bb)>10)
					{
					  for(col=0;col<2;col++)
					    {
					      tripletsForHessian.emplace_back(index[i+a][j+b]*2+row,index[i+aa][j+bb]*2+col,hes(ct*2+row,ct2*2+col));
					    }
					  ct2++;
					   
					}
				    }
				}
			    }			 
			  ct++;
			  
			}
		    }
		}
	    }
	}
    }

  Jacobian*=-1;

  printf("assemble complete\n");
  SparseMatrix<double,Eigen::RowMajor > Hessianspa(num_hair_2,num_hair_2);
  Hessianspa.setFromTriplets(tripletsForHessian.begin(),tripletsForHessian.end());
  Hessianspa.makeCompressed();
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > linearSolver;

  linearSolver.compute(Hessianspa);

  VectorXd dx(num_hair_2);
  dx=linearSolver.solve(Jacobian);

  double dot,det,angle;
  for(i=0;i<height;++i)
    {
      for(j=0;j<width;++j)
	{
	  if(whether_hair(i,j)>10)
	    {
	      if(i==341&&j==358)
		{
		  cout<<"direction before optimization: "<<dir[i][j][0]<<" "<<dir[i][j][1]<<endl;
		}
	     
	      for(a=0;a<2;a++)
		{
		  dir[i][j][a]=dir[i][j][a]+dx(index[i][j]*2+a);		   
		}
	      if(i==341&&j==358)
		{
		  cout<<"dx: "<<dx(index[i][j]*2+0)<<" "<<dx(index[i][j]*2+1)<<endl;
		  cout<<"direction after optimization: "<<dir[i][j][0]<<" "<<dir[i][j][1]<<endl;

		}

	      double x=dir[i][j][0]; double y=dir[i][j][1];
	      double norm=sqrt(x*x+y*y);
	      x=x/norm; y=y/norm;
	      angle=acos(x);

	      if(fabs(sin(angle)-y)>0.0001)
		{
		  angle=2*pi-angle;
		}
	      dense_theta_map(i,j)=angle;
	    }
	  else
	    {
	      dense_theta_map(i,j)=0;
	    }
	  if(i==341&&j==358)
	    {
	      printf("theta after optimization %lf\n",dense_theta_map(i,j));
	    }
	}
    }

  FILE* fp;
  fp=fopen(dense_theta_map_py_path.c_str(),"w");
  for(i=0;i<height;i++)
    {
      for(j=0;j<width;j++)
	{
	  if(whether_hair(i,j)>10)
	    {
	      fprintf(fp,"%lf ",dense_theta_map(i,j));
	    }
	  else
	    {
	      fprintf(fp,"-1 ");
	    }
	  
	}
      fprintf(fp,"\n");
    }
  fclose(fp);
  /*
  // for better visualization 
  for(i=0;i<height;i++)
    {
      for(j=0;j<width;j++)
	{
	  if(sparse_theta_map(i,j)<0&&whether_hair(i,j)>10)
	    {
	      sparse_theta_map(i,j)+=pi*2;
	    }
	  else if(whether_hair(i,j)<0)
	    {
	      sparse_theta_map(i,j)=0;
	    }
	}
    }
  my_io->saveAsVTK(sparse_theta_map,whether_hair,height,width,dense_theta_map_path);
  */
  ////////////////////////////////////////////////
  my_io->saveAsVTK(dense_theta_map,whether_hair,height,width,dense_theta_map_path);
  delete my_io;
  return 0;
}
