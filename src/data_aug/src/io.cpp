#include "io.h"
using namespace std;

int io::saveArrayAsVTK(std::vector<vertex> &myvertexs,std::vector<tri_face> &mytris,const std::string name)
{
  FILE *fp;
  int i,j;
  //printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"vtk output\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET POLYDATA\n");
  fprintf(fp,"POINTS %d float\n",myvertexs.size());
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%f %f %f\n",myvertexs[i].location.x,myvertexs[i].location.y,myvertexs[i].location.z);
    }
  fprintf(fp,"POLYGONS %d %d\n",mytris.size(),mytris.size()*4);
  for(i=0;i<mytris.size();++i)
    {
      fprintf(fp,"3");
      for(j=0;j<3;++j)
	{
	  fprintf(fp," %d",mytris[i].index_vertex[j]);
	}
      fprintf(fp,"\n");

    }
  fprintf(fp,"POINT_DATA %d\n",myvertexs.size());
  fprintf(fp,"TEXTURE_COORDINATES TCoords 2 float\n");
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%f %f\n",myvertexs[i].u_x,myvertexs[i].u_y);
    }
  fclose(fp);
  return 0;
}

int io::savePolylineAsVTK(std::vector<vertex> &polyline,const std::string name)
{
  FILE *fp;
  int i,j;
  //printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"line\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n",polyline.size());
  
  for(i=0;i<polyline.size();++i)
    {
      fprintf(fp,"%f %f %f\n",polyline[i].location.x,polyline[i].location.y,polyline[i].location.z);
    }
  fprintf(fp,"CELLS %d %d\n",polyline.size()-1,(polyline.size()-1)*3);
  for(i=0;i<polyline.size()-1;++i)
    {
      fprintf(fp,"2 %d %d\n",i,i+1);
    }
  fprintf(fp,"CELL_TYPES %d\n",polyline.size()-1);
  for(i=0;i<polyline.size()-1;++i)
    {
      fprintf(fp,"3\n");
    }

   fclose(fp);
   return 0;
}


int io::getVertexAndTri(std::vector<vertex> &myvertexs,std::vector<tri_face > &mytris,const std::string name)
{
  printf("%s \n",name.c_str());
  FILE* fp;
  fp=fopen(name.c_str(),"r");
  char filter[15],filter_1[15],filter_2[15];
  
  do
   {
     fscanf(fp,"%s",filter);
   }while(filter[0]!='P'||filter[1]!='O'||filter[2]!='I');

  int num_vertex_file;
  int index_vertex_file[3]; int index_vertex[3];
  int num_tri_file,num_tri4_file;
  int filter3;
  float x,y,z;

  fscanf(fp,"%d%s",&num_vertex_file,filter);
  int i,j;
  for(i=0;i<num_vertex_file;++i)
    {
      fscanf(fp,"%f%f%f",&x,&y,&z);
      myvertexs.push_back(vertex(myvector(x,y,z)));
    }
  fscanf(fp,"%s%d%d",filter,&num_tri_file,&num_tri4_file);
 
  for(i=0;i<num_tri_file;++i)
    {
      fscanf(fp,"%d",&filter3);
      for(j=0;j<3;++j)
	{
	  fscanf(fp,"%d",&index_vertex_file[j]);
	  index_vertex[j]=index_vertex_file[j];
	}
      mytris.push_back(tri_face(index_vertex));
    }

  fscanf(fp,"%s%d",filter,&num_vertex_file);
  fscanf(fp,"%s",filter);
  if(filter[0]=='N')
    {
      fscanf(fp,"%s",filter); fscanf(fp,"%s",filter);
      for(i=0;i<num_vertex_file;i++)
	{
	  fscanf(fp,"%f%f%f",&x,&y,&z);
	}
      fscanf(fp,"%s%s%d%s",filter,filter_1,&x,filter_2);
      printf("%s %s %s\n",filter,filter_1,filter_2);
      for(i=0;i<num_vertex_file;i++)
	{
	  fscanf(fp,"%f%f",&x,&y);
	  myvertexs[i].u_x=x;
	  myvertexs[i].u_y=y;
	}
    }
  else
    {
      fscanf(fp,"%s%d%s",filter_1,&x,filter_2);
      printf("%s %s\n",filter_1,filter_2);
      for(i=0;i<num_vertex_file;i++)
	{
	  fscanf(fp,"%f%f",&x,&y);
	  myvertexs[i].u_x=x;
	  myvertexs[i].u_y=y;
	}
    }
  fclose(fp);
  return 0;
}

/*

// assimp can not be used because of the lack of robustness
int io::getVertexAndTriFromObj(std::vector<vertex> &myvertexs,std::vector<tri_face > &mytris,const std::string name)
{
  Assimp::Importer import;
  const aiScene *scene = import.ReadFile(name, aiProcess_Triangulate | aiProcess_FlipUVs);	
  if(!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) 
    {
      cout << "ERROR::ASSIMP::" << import.GetErrorString() << endl;
      return 1;
    }
  int i,j,k;
  float x,y,z;
  // size_t index_vertex[3];
  printf("have mesh: %d\n",scene->HasMeshes());
  printf("num meshes: %d\n",scene->mNumMeshes);
  for(i=0;i<scene->mNumMeshes;i++)
    {
      aiMesh *mesh = scene->mMeshes[i];
      aiVector3D* mVertices=mesh->mVertices;
      aiFace* mFaces=mesh->mFaces;
      printf("num vertex: %d\n",mesh->mNumVertices);
      printf("num face: %d\n",mesh->mNumFaces);
      for(j=0;j<10;j++)
	{
	  x=mVertices[j].x; y=mVertices[j].y; z=mVertices[j].z;
	  printf("x: %lf y: %lf z: %lf\n",x,y,z);
	  myvertexs.push_back(vertex(myvector(x,y,z)));
	}

      for(j=0;j<0;j++)
	{	  
	  x=mesh->mTextureCoords[j][0].x; y=mesh->mTextureCoords[j][0].y;
	  printf("x: %lf y: %lf\n",x,y);
	  myvertexs[j].u_x=x; myvertexs[j].u_y=y;
	}
      for(j=0;j<10;j++)
	{
	  for(k=0;k<3;++k)
	    {
	      index_vertex[k]=mFaces[j].mIndices[k];
	    }
	  printf("index %d %d %d\n",index_vertex[0],index_vertex[1],index_vertex[2]);
	  mytris.push_back(tri_face(index_vertex));
	}
    }
  
  return 0;
  }*/

int io::getVertexAndTriFromObj(std::vector<vertex> &myvertexs,std::vector<tri_face > &mytris,const std::string name)
{
  typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
  MyMesh mesh;
  mesh.request_vertex_texcoords2D();
  mesh.request_face_texture_index();
  OpenMesh::IO::Options opt = OpenMesh::IO::Options::VertexTexCoord;
  OpenMesh::IO::read_mesh(mesh,name,opt);

  int now=0;
  for (MyMesh::VertexIter v_it=mesh.vertices_begin();v_it!=mesh.vertices_end();++v_it)
    {
      OpenMesh::Vec3f &point_now=mesh.point(*v_it);
      vertex vertex_now=vertex(myvector(point_now[0],point_now[1],point_now[2]));
      if (mesh.has_vertex_texcoords2D()) 
	{
	  vertex_now.u_x=mesh.texcoord2D(*v_it).data()[0];
	  vertex_now.u_y=mesh.texcoord2D(*v_it).data()[1];
	  /*
	  if(now<10)
	    {
	      printf("%f %f\n",vertex_now.u_x,vertex_now.u_y);
	      }*/
	}
      myvertexs.push_back(vertex_now);    
      now++;
    }

  int index_vertex[3];

  int num_face=0;
  for (MyMesh::FaceIter f_it=mesh.faces_begin();f_it!=mesh.faces_end();++f_it)
    {      
      int index=0;
      for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it) 
	{ 
	  int id = fv_it->idx(); 
	  index_vertex[index] = id; 
	  index++; 
	}
      mytris.push_back(tri_face(index_vertex));            
      num_face++;
    }
  printf("num vertex: %d num_face: %d\n",now,num_face);  
  return 0;
}


int io::saveArrayAsVTKwithLabel(int num_strip,int num_tri,int which_vertex[],int label_strip[],std::vector<tri_face > mp[],std::vector<vertex > &myvertexs, const std::string name)
{
  FILE *fp;
  int i,j,k;
  printf("%s \n",name.c_str());
  
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"vtk output\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET POLYDATA\n");
  fprintf(fp,"POINTS %d float\n",myvertexs.size());
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%f %f %f\n",myvertexs[i].location.x,myvertexs[i].location.y,myvertexs[i].location.z);
    }
  fprintf(fp,"POLYGONS %d %d\n",num_tri,num_tri*4);
  for(i=0;i<num_strip;++i)
    {
      int which_vertex_now=which_vertex[i];
      int num_tri_now=mp[which_vertex_now].size();
      for(j=0;j<num_tri_now;j++)
	{
	  fprintf(fp,"3");
	  for(k=0;k<3;++k)
	    {
	      fprintf(fp," %d",mp[which_vertex_now][j].index_vertex[k]);
	    }
	  fprintf(fp,"\n");
	}
    }
  fprintf(fp,"CELL_DATA %d\n",num_tri);
  fprintf(fp,"SCALARS labels int 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");

  for(i=0;i<num_strip;++i)
    {
      int which_vertex_now=which_vertex[i];
      int label_now=label_strip[i];
      int num_tri_now=mp[which_vertex_now].size();
      for(j=0;j<num_tri_now;j++)
	{
	  fprintf(fp,"%d\n",label_now);
	}
    }
  fclose(fp);
  
  return 0;
}

int io::saveClusterPolylineAsVTK(std::vector<int > &strip_cluster,std::vector<vertex> polyline[],std::vector<vertex > &myvertexs,const std::string name)
{
  FILE *fp;
  int i,j;

  int num_vertex=0;
  int num_line=0;
  for(i=0;i<strip_cluster.size();i++)
    {
      int index_polyline=strip_cluster[i];
      num_vertex+=polyline[index_polyline].size();
      num_line+=polyline[index_polyline].size()-1;
    }
  //printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"line\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n",num_vertex);

  for(i=0;i<strip_cluster.size();i++)
    {
      int index_polyline=strip_cluster[i];
      for(j=0;j<polyline[index_polyline].size();++j)
	{
	  fprintf(fp,"%f %f %f\n",polyline[index_polyline][j].location.x,polyline[index_polyline][j].location.y,polyline[index_polyline][j].location.z);
	}
    }
  
  fprintf(fp,"CELLS %d %d\n",num_line,num_line*3);

  int num_vertex_now=0;
  for(i=0;i<strip_cluster.size();i++)
    {
      int index_polyline=strip_cluster[i];
      for(j=0;j<polyline[index_polyline].size()-1;++j)
	{
	  fprintf(fp,"2 %d %d\n",num_vertex_now+j,num_vertex_now+j+1);
	}
      num_vertex_now+=polyline[index_polyline].size();
    }
  fprintf(fp,"CELL_TYPES %d\n",num_line);
  for(i=0;i<num_line;++i)
    {
      fprintf(fp,"3\n");
    }

  fclose(fp);

  return 0;
}


int io::saveClusterAsVTK(std::vector<int > &strip_cluster,std::vector<vertex > &myvertexs,int which_vertex[],std::vector<tri_face > mp[],const std::string name)
{
  FILE *fp;
  int i,j,k;

  int num_tri=0;
  for(i=0;i<strip_cluster.size();i++)
    {
      int strip_now=strip_cluster[i]; int vertex_now=which_vertex[strip_now];
      num_tri+=mp[vertex_now].size();
    }
  //printf("%s \n",name.c_str());
  
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"vtk output\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET POLYDATA\n");
  fprintf(fp,"POINTS %d float\n",myvertexs.size());
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%f %f %f\n",myvertexs[i].location.x,myvertexs[i].location.y,myvertexs[i].location.z);
    }
  fprintf(fp,"POLYGONS %d %d\n",num_tri,num_tri*4);
  for(i=0;i<strip_cluster.size();++i)
    {
      int strip_now=strip_cluster[i]; int vertex_now=which_vertex[strip_now];
      int num_tri_now=mp[vertex_now].size();
      for(j=0;j<num_tri_now;j++)
	{
	  fprintf(fp,"3");
	  for(k=0;k<3;++k)
	    {
	      fprintf(fp," %d",mp[vertex_now][j].index_vertex[k]);
	    }
	  fprintf(fp,"\n");
	}
    }
  fclose(fp);
  return 0;
}

int io::saveCombine(object &myobject_target,object &myobject_source,std::vector<int > &cluster_now_used,const std::string name)
{
  FILE *fp;
  int i,j,k,l;

  int num_tri=0;
  num_tri+=myobject_target.num_tri;
  for(i=0;i<cluster_now_used.size();i++)
    {
      int cluster_now=cluster_now_used[i];
      for(j=0;j<myobject_source.strip_cluster[cluster_now].size();j++)
	{
	  int strip_now=myobject_source.strip_cluster[cluster_now][j]; int vertex_now=myobject_source.which_vertex[strip_now];
	  num_tri+=myobject_source.mp[vertex_now].size();
	}      
    }
  //printf("%s \n",name.c_str());
  
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"vtk output\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET POLYDATA\n");
  fprintf(fp,"POINTS %d float\n",myobject_target.num_vertex+myobject_source.num_vertex);
  for(i=0;i<myobject_target.myvertexs.size();++i)
    {
      fprintf(fp,"%f %f %f\n",myobject_target.myvertexs[i].location.x,myobject_target.myvertexs[i].location.y,myobject_target.myvertexs[i].location.z);
    }

  for(i=0;i<myobject_source.myvertexs.size();++i)
    {
      fprintf(fp,"%f %f %f\n",myobject_source.myvertexs[i].location.x,myobject_source.myvertexs[i].location.y,myobject_source.myvertexs[i].location.z);
    }
  
  fprintf(fp,"POLYGONS %d %d\n",num_tri,num_tri*4);

  for(i=0;i<myobject_target.num_tri;++i)
    {
      fprintf(fp,"3");
      for(j=0;j<3;++j)
	{
	  fprintf(fp," %d",myobject_target.mytris[i].index_vertex[j]);
	}
      fprintf(fp,"\n");
    }

  for(i=0;i<cluster_now_used.size();i++)
    {
      int cluster_now=cluster_now_used[i];
      for(j=0;j<myobject_source.strip_cluster[cluster_now].size();j++)
	{
	  int strip_now=myobject_source.strip_cluster[cluster_now][j]; int vertex_now=myobject_source.which_vertex[strip_now];	  
	  int num_tri_now=myobject_source.mp[vertex_now].size();
	  for(k=0;k<num_tri_now;k++)
	    {
	      fprintf(fp,"3");
	      for(l=0;l<3;++l)
		{
		  fprintf(fp," %d",myobject_source.mp[vertex_now][k].index_vertex[l]+myobject_target.num_vertex);
		}
	      fprintf(fp,"\n");
	    }
	}      
    }

  fprintf(fp,"CELL_DATA %d\n",num_tri);
  fprintf(fp,"SCALARS combine int 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");


  for(i=0;i<myobject_target.num_tri;++i)
    {
      fprintf(fp,"0\n");      
    }

  for(i=0;i<cluster_now_used.size();i++)
    {
      int cluster_now=cluster_now_used[i];
      for(j=0;j<myobject_source.strip_cluster[cluster_now].size();j++)
	{
	  int strip_now=myobject_source.strip_cluster[cluster_now][j]; int vertex_now=myobject_source.which_vertex[strip_now];	  
	  int num_tri_now=myobject_source.mp[vertex_now].size();
	  for(k=0;k<num_tri_now;k++)
	    {
	      fprintf(fp,"1\n");
	    }
	}      
    }

  fclose(fp);
  return 0;
}

int io::saveCombine(object &myobject_target,object &myobject_source,std::vector<int > &cluster_now_used_target,std::vector<int > &cluster_now_used,const std::string name)
{
  FILE *fp;
  int i,j,k,l;

  int num_tri=0;

  for(i=0;i<cluster_now_used_target.size();i++)
    {
      int cluster_now=cluster_now_used_target[i];
      for(j=0;j<myobject_target.strip_cluster[cluster_now].size();j++)
	{
	  int strip_now=myobject_target.strip_cluster[cluster_now][j]; int vertex_now=myobject_target.which_vertex[strip_now];
	  num_tri+=myobject_target.mp[vertex_now].size();
	}      
    }
 
  for(i=0;i<cluster_now_used.size();i++)
    {
      int cluster_now=cluster_now_used[i];
      for(j=0;j<myobject_source.strip_cluster[cluster_now].size();j++)
	{
	  int strip_now=myobject_source.strip_cluster[cluster_now][j]; int vertex_now=myobject_source.which_vertex[strip_now];
	  num_tri+=myobject_source.mp[vertex_now].size();
	}      
    }
  //printf("%s \n",name.c_str());
  
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"vtk output\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET POLYDATA\n");
  fprintf(fp,"POINTS %d float\n",myobject_target.num_vertex+myobject_source.num_vertex);
  for(i=0;i<myobject_target.myvertexs.size();++i)
    {
      fprintf(fp,"%f %f %f\n",myobject_target.myvertexs[i].location.x,myobject_target.myvertexs[i].location.y,myobject_target.myvertexs[i].location.z);
    }

  for(i=0;i<myobject_source.myvertexs.size();++i)
    {
      fprintf(fp,"%f %f %f\n",myobject_source.myvertexs[i].location.x,myobject_source.myvertexs[i].location.y,myobject_source.myvertexs[i].location.z);
    }
  
  fprintf(fp,"POLYGONS %d %d\n",num_tri,num_tri*4);

  for(i=0;i<cluster_now_used_target.size();i++)
    {
      int cluster_now=cluster_now_used_target[i];
      for(j=0;j<myobject_target.strip_cluster[cluster_now].size();j++)
	{
	  int strip_now=myobject_target.strip_cluster[cluster_now][j]; int vertex_now=myobject_target.which_vertex[strip_now];	  
	  int num_tri_now=myobject_target.mp[vertex_now].size();
	  for(k=0;k<num_tri_now;k++)
	    {
	      fprintf(fp,"3");
	      for(l=0;l<3;++l)
		{
		  fprintf(fp," %d",myobject_target.mp[vertex_now][k].index_vertex[l]);
		}
	      fprintf(fp,"\n");
	    }
	}      
    }

  for(i=0;i<cluster_now_used.size();i++)
    {
      int cluster_now=cluster_now_used[i];
      for(j=0;j<myobject_source.strip_cluster[cluster_now].size();j++)
	{
	  int strip_now=myobject_source.strip_cluster[cluster_now][j]; int vertex_now=myobject_source.which_vertex[strip_now];	  
	  int num_tri_now=myobject_source.mp[vertex_now].size();
	  for(k=0;k<num_tri_now;k++)
	    {
	      fprintf(fp,"3");
	      for(l=0;l<3;++l)
		{
		  fprintf(fp," %d",myobject_source.mp[vertex_now][k].index_vertex[l]+myobject_target.num_vertex);
		}
	      fprintf(fp,"\n");
	    }
	}      
    }

  fprintf(fp,"CELL_DATA %d\n",num_tri);
  fprintf(fp,"SCALARS combine int 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");


  for(i=0;i<cluster_now_used_target.size();i++)
    {
      int cluster_now=cluster_now_used_target[i];
      for(j=0;j<myobject_target.strip_cluster[cluster_now].size();j++)
	{
	  int strip_now=myobject_target.strip_cluster[cluster_now][j]; int vertex_now=myobject_target.which_vertex[strip_now];	  
	  int num_tri_now=myobject_target.mp[vertex_now].size();
	  for(k=0;k<num_tri_now;k++)
	    {
	      fprintf(fp,"0\n");
	    }
	}      
    }

  for(i=0;i<cluster_now_used.size();i++)
    {
      int cluster_now=cluster_now_used[i];
      for(j=0;j<myobject_source.strip_cluster[cluster_now].size();j++)
	{
	  int strip_now=myobject_source.strip_cluster[cluster_now][j]; int vertex_now=myobject_source.which_vertex[strip_now];	  
	  int num_tri_now=myobject_source.mp[vertex_now].size();
	  for(k=0;k<num_tri_now;k++)
	    {
	      fprintf(fp,"1\n");
	    }
	}      
    }

  fclose(fp);
  return 0;
}
