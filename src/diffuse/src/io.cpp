#include "io.h"
#include "head.h"
using namespace std;
using namespace Eigen;

int io::saveAsVTK(Eigen::MatrixXd &dense_theta_map,Eigen::MatrixXd &whether_hair,int height,int width,const std::string name)
{
  FILE *fp;
  int i,j;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"point\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n",height*width);
  for(i=0;i<height;++i)
    {
      for(j=0;j<width;++j)
	{
	  fprintf(fp,"%d %d 0\n",i,j);
	}
    }
  fprintf(fp,"CELLS %d %d\n",height*width,height*width*2);
  int index=0;
  for(i=0;i<height;++i)
    {
      for(j=0;j<width;++j)
	{
	  fprintf(fp,"1 %d\n",index);
	  index++;
	}
    }
  
  fprintf(fp,"CELL_TYPES %d\n",height*width);
  for(i=0;i<height;++i)
    {
      for(j=0;j<width;++j)
	{
	  fprintf(fp,"1\n");
	}
    }

  fprintf(fp,"POINT_DATA %d\n",height*width);
  fprintf(fp,"SCALARS theta double\n");
  fprintf(fp,"LOOKUP_TABLE default\n");

  for(i=0;i<height;++i)
    {
      for(j=0;j<width;++j)
	{
	  // -1 or real theta
	  fprintf(fp,"%lf\n",dense_theta_map(i,j));
	}
    }
  fprintf(fp,"SCALARS whether_hair double\n");
  fprintf(fp,"LOOKUP_TABLE default\n");

  for(i=0;i<height;++i)
    {
      for(j=0;j<width;++j)
	{
	  fprintf(fp,"%lf\n",whether_hair(i,j));
	}
    }
  
  fclose(fp);
  return 0;
}

int io::getData(Eigen::MatrixXd &sparse_theta_map,Eigen::MatrixXd &whether_hair,int height,int width,const std::string name)
{
  printf("%s \n",name.c_str());
  FILE* fp;
  fp=fopen(name.c_str(),"r");
  char filter[10];
  do
   {
     fscanf(fp,"%s",filter);
   }while(filter[0]!='P');
  int num_vertex_file;
  int index_vertex_file[4]; size_t index_vertex[4];
  int filter_1,filter_2;
  // 从文件中读取整数用int接受%d ,即使使用时是用size_t ,如果直接用size_t接受%u使用是会出错，此处可能是本身的bug

  int x,y,z;
  int num_simplexs_file,num_simplexs_dim_2_file;
  fscanf(fp,"%d%s",&num_vertex_file,filter);
  int i,j;
  for(i=0;i<num_vertex_file;++i)
    {
      fscanf(fp,"%d%d%d",&x,&y,&z);
    }
  fscanf(fp,"%s%d%d",filter,&num_simplexs_file,&num_simplexs_dim_2_file);

  for(i=0;i<num_simplexs_file;++i)
    {
      fscanf(fp,"%d %d",&filter_1,&filter_2);
    }

  fscanf(fp,"%s%d",filter,&num_vertex_file);
  for(i=0;i<num_vertex_file;++i)
    {
      fscanf(fp,"%d\n",&x);
    }
  
  fscanf(fp,"%s%d",filter,&num_vertex_file);
  fscanf(fp,"%s",filter);fscanf(fp,"%s",filter);fscanf(fp,"%s",filter);
  fscanf(fp,"%s",filter);fscanf(fp,"%s",filter);

  double theta_now,whether_hair_now;
  for(i=0;i<height;++i)
    {
      for(j=0;j<width;++j)
	{
	  fscanf(fp,"%lf",&theta_now);
	  sparse_theta_map(i,j)=theta_now;
	}
    }

   fscanf(fp,"%s",filter);fscanf(fp,"%s",filter);fscanf(fp,"%s",filter);
  fscanf(fp,"%s",filter);fscanf(fp,"%s",filter);


  for(i=0;i<height;++i)
    {
      for(j=0;j<width;++j)
	{
	  fscanf(fp,"%lf",&whether_hair_now);
	  whether_hair(i,j)=whether_hair_now;
	}
    }
  
  fclose(fp);
  return 0;
}


