#include "head.h"
#include "vertex.h"
#include "tri_face.h"
#include "io.h"
#include "object.h"
using namespace std;

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

int Find(int a,int* father)
{
  if(father[a]!=a)
    {
      father[a]=Find(father[a],father);
    }
  return father[a];
}
void Union(int a,int b,int* father)
{
  a=Find(a,father);
  b=Find(b,father);
  if(a==b)
    {
      return;
    }
  else
    {
      if(a<b)
	{
	  father[b]=a;
	}
      else
	{
	  father[a]=b;
	}
    }
  return;
}

bool cmp(const vertex&a,const vertex&b)
{
  if(fabs(a.u_y-b.u_y)<0.007)
    {
      return a.u_x<b.u_x;
    }
  else
    {
      return a.u_y>b.u_y;
    }
}


int decompose(boost::property_tree::ptree &para_tree,object &myobject,int which_object,PyObject* pFunc)
{
  while(!myobject.mytris.empty())
    {
      myobject.mytris.pop_back();
    }
  while(!myobject.myvertexs.empty())
    {
      myobject.myvertexs.pop_back();
    }

  stringstream ss; string now_string; ss<<which_object; ss>>now_string;
  string help="mesh_obj_"+now_string+".value";
  string mesh_obj_name=para_tree.get<string>(help,"");
  help="mesh_vtk_"+now_string+".value";
  string mesh_vtk_name=para_tree.get<string>(help,"");
  help="mesh_vtk_new_"+now_string+".value";
  string mesh_vtk_new_name=para_tree.get<string>(help,"");
  help="out_dir_"+now_string+".value";
  string out_dir=para_tree.get<string>(help,"");

  help="distance_path_"+now_string+".value";
  string distance_path=para_tree.get<string>(help,"");
  help="mesh_label_"+now_string+".value";
  string mesh_label_name=para_tree.get<string>(help,"");

  io* myio=new io();
  
  //myio->getVertexAndTri(myvertexs,mytris,mesh_vtk_name);
  myio->getVertexAndTriFromObj(myobject.myvertexs,myobject.mytris,mesh_obj_name);
  myobject.num_vertex=myobject.myvertexs.size(); myobject.num_tri=myobject.mytris.size();

  int i,j,k;

  ////////////////////////////////////////////////////////////
  //Union-Find set for finding the connected component
  int num_vertex=myobject.num_vertex; int num_tri=myobject.num_tri;
  int father[num_vertex];
  int ct[num_vertex];

  for(i=0;i<num_vertex;i++)
    {
      ct[i]=-1;
      father[i]=i;
    }
  printf("num_vertex: %d num_tri: %d\n",num_vertex,num_tri);
  for(i=0;i<num_tri;i++)
    {
      Union(myobject.mytris[i].index_vertex[0],myobject.mytris[i].index_vertex[1],father);
      Union(myobject.mytris[i].index_vertex[0],myobject.mytris[i].index_vertex[2],father);
      Union(myobject.mytris[i].index_vertex[2],myobject.mytris[i].index_vertex[1],father);
    }

  
  for(i=0;i<num_vertex;++i)
    {
      while(!myobject.mp[i].empty())
	{
	  myobject.mp[i].pop_back();
	}
    }
  for(i=0;i<num_tri;i++)
    {
      int kk=Find(myobject.mytris[i].index_vertex[0],father);
      //printf("kk: %d\n",kk);
      myobject.mp[kk].push_back(myobject.mytris[i]);
    }

  int now=0;
  for(i=0;i<num_vertex;i++)
    {
      if(!myobject.mp[i].empty())
	{
	  stringstream ss; string now_string; ss<<now; ss>>now_string;
	  string vtk_now=out_dir+"/test_"+now_string+".vtk";
	  myio->saveArrayAsVTK(myobject.myvertexs,myobject.mp[i],vtk_now);
	  now++;
	}
    }

  myobject.num_strip=now;
  int num_strip=now;
  
  for(i=0;i<num_strip;i++)
    {
      while(!myobject.polyline[i].empty())
	{
	  myobject.polyline[i].pop_back();
	}
    }
  
  
  printf("the total num of connected component: %d\n",now);
  /////////////////////////////////////////////////////
  myio->saveArrayAsVTK(myobject.myvertexs,myobject.mytris,mesh_vtk_new_name);

  ///////////////////////////////////////////////////////////
  //  get simplified representation of strip

  int ct_now=0;
  for(i=0;i<num_vertex;i++)
    {
      if(!myobject.mp[i].empty())
	{
	  myobject.which_vertex[ct_now]=i;
	  vector<vertex > vertex_set;
	  while(!vertex_set.empty())
	    {
	      vertex_set.pop_back();
	    }
	  int num_face_now=myobject.mp[i].size();
	  for(j=0;j<num_face_now;j++)
	    {
	      tri_face tri_face_now=myobject.mp[i][j];
	      for(k=0;k<3;k++)
		{
		  int index_vertex_now=tri_face_now.index_vertex[k];
		  if(myobject.myvertexs[index_vertex_now].strip==-1)
		    {
		      myobject.myvertexs[index_vertex_now].strip=ct_now;
		      vertex_set.push_back(myobject.myvertexs[index_vertex_now]);
		    }
		}	      
	    }
	  int size_vertex_set=vertex_set.size();
	  
	  sort(vertex_set.begin(),vertex_set.end(),cmp);

	  /*
	    if(ct_now==128)
	    {
	    printf("size_vertex_set: %d\n",size_vertex_set);
	      
	    for(j=0;j<size_vertex_set;j++)
	    {
	    printf("x: %f y: %f z: %f\n",vertex_set[j].location.x,vertex_set[j].location.y,vertex_set[j].location.z);
	    printf("u_y: %f u_x: %f\n",vertex_set[j].u_y,vertex_set[j].u_x);
	    }
	    }*/
	  double radius_all=0;
	  int num_level=0;

	  vertex v_begin,v_end;
	  for(j=0;j<size_vertex_set;j++)
	    {
	      v_begin=vertex_set[j];
	      v_end=vertex_set[j];
	      for(k=j+1;k<size_vertex_set;k++)
		{
		  if(fabs(vertex_set[k].u_y-v_begin.u_y)<0.007)
		    {
		      v_end=vertex_set[k];
		    }
		  else
		    {
		      break;
		    }
		}
	      j=k-1;

	      if(ct_now==128)
		{
		  //printf("v_begin: %f %f %f v_end: %f %f %f\n",v_begin.location.x,v_begin.location.y,v_begin.location.z,v_end.location.x,v_end.location.y,v_end.location.z);
		}
	      myvector loc_now=0.5*(v_begin.location+v_end.location);
	      myobject.polyline[ct_now].push_back(vertex(loc_now));
	      
	      myvector radius_v=(v_begin.location-v_end.location)*0.5;
	      
	      radius_all+=radius_v.len();
	      num_level+=1;
	    }

	  
	  if(ct_now==128)
	    {
	      //    printf("num_level: %d\n",num_level);    
	    }	  
	  myobject.radius_strip[ct_now]=radius_all/num_level;
	  if(ct_now==128)
	    {
	      //  printf("radius_strip: %f\n",radius_strip[ct_now]);
	    }
	  ct_now++;
	}
    }
  for(i=0;i<num_strip;i++)
    {
      for(j=0;j<num_strip;j++)
	{
	  if(i==j)
	    {
	      myobject.dis[i][j]=0;
	    }
	  else if(j<i)
	    {
	      myobject.dis[i][j]=myobject.dis[j][i];
	    }
	  else if(j>i)
	    {
	      int min_siz=myobject.polyline[i].size()>myobject.polyline[j].size()?myobject.polyline[j].size():myobject.polyline[i].size();
	      float xx=0; 
	      for(k=0;k<min_siz;k++)
		{
		  myvector loc=myobject.polyline[i][k].location-myobject.polyline[j][k].location;
		  //float here=fabs(loc.len()); // ignore the radius  		  
		  float here=fabs(loc.len())-myobject.radius_strip[i]-myobject.radius_strip[j];
		  float xx_now=here>0.0?here:0.0;
		  xx+=xx_now;
		}
	      myobject.dis[i][j]=xx/min_siz;
	    }
	}      
    }

  for(i=0;i<num_strip;i++)
    {
      for(j=0;j<num_strip;j++)
	{
	  //printf("%f ",myobject.dis[i][j]);
	}
      //printf("\n");
    }

  
  FILE* fp;
  fp=fopen(distance_path.c_str(),"w");
  for(i=0;i<num_strip;i++)
    {
      for(j=0;j<num_strip;j++)
	{
	  fprintf(fp,"%f ",myobject.dis[i][j]);	  
	}
      fprintf(fp,"\n");
    }
  fclose(fp);


  //  printf("d 1\n");
  ////////////////////////////////////////////////////////////////////////
  //cluster test    

  //set paras
  PyObject *pArgsT = PyTuple_New(1);
  PyObject* pArgsD = PyDict_New();
  PyDict_SetItemString(pArgsD, "distance_path", Py_BuildValue("s", distance_path.c_str()));
  PyTuple_SetItem(pArgsT, 0, pArgsD);
  //调用py方法

  PyObject *pReturn = PyObject_CallObject(pFunc,pArgsT); //PyEval_CallObject(pFunc, pArgsT);
  //  pReturn = PyObject_CallObject(pFunc,pArgsT); 

  // get the returned value from Python
  int nTupleSize = PyTuple_Size(pReturn);
  int nValue = -2;

  int num_cluster=0;
  for (int l = 0; l < nTupleSize; l++)
    {
      PyObject *pTuple_now = PyTuple_GET_ITEM(pReturn, l);
      PyArg_Parse(pTuple_now, "i", &nValue);
      if(nValue>num_cluster)
	{
	  num_cluster=nValue;
	}
      myobject.label_strip[l]=nValue;
      //      std::cout << nValue << std::endl;

      /*std::cout<< (*pTuple_now)<<std::end;
	int nTupleListSize = PyList_Size(pTuple_now);
	printf("dddd %d\n",nTupleListSize);
	for (int m = 0; m < nTupleListSize; m++)
        {
	PyObject* pTupleListValue = PyList_GetItem(pTuple_now, m);
	  
	}*/
    }
  
  num_cluster++;
  myobject.num_cluster=num_cluster;
  //PyRun_SimpleString("import cluster");
  //PyRun_SimpleString("cluster.test()");
  

  
  for(i=0;i<num_cluster;i++)
    {
      while(!myobject.strip_cluster[i].empty())
	{
	  myobject.strip_cluster[i].pop_back();
	}
    }
  for(i=0;i<num_strip;i++)
    {
      int cluster_now=myobject.label_strip[i];
      myobject.strip_cluster[cluster_now].push_back(i);
    }

  myio->saveArrayAsVTKwithLabel(num_strip,num_tri,myobject.which_vertex,myobject.label_strip,myobject.mp,myobject.myvertexs,mesh_label_name);

  //printf("d 4\n");
  /////////////////////////////////////////////////////////////////////////

  for(i=0;i<num_cluster;i++)
    {
      stringstream ss; string now_string; ss<<i; ss>>now_string;
      string polyline_vtk_now=out_dir+"/cluster_polyline_"+now_string+".vtk";
      string strip_vtk_now=out_dir+"/cluster_strip_"+now_string+".vtk";
      myio->saveClusterPolylineAsVTK(myobject.strip_cluster[i],myobject.polyline,myobject.myvertexs,polyline_vtk_now);
      myio->saveClusterAsVTK(myobject.strip_cluster[i],myobject.myvertexs,myobject.which_vertex,myobject.mp,strip_vtk_now);
    }
  
  for(i=0;i<num_strip;i++)
    {
      stringstream ss; string now_string; ss<<i; ss>>now_string;
      string vtk_now=out_dir+"/polyline_"+now_string+".vtk";
      myio->savePolylineAsVTK(myobject.polyline[i],vtk_now);
    }

  //printf("d 5\n");
  delete myio;
  return 0;
}


void combine( int a[], int n, int m,  int b[], const int M, vector<int > &help )
{ 
  for(int i=n; i>=m; i--)   // 注意这里的循环范围
    {
      b[m-1] = i - 1;
      if (m > 1)
	combine(a,i-1,m-1,b,M,help);
      else                     // m == 1, 输出一个组合
	{   
	  for(int j=M-1; j>=0; j--)
	    {
	      help.push_back(a[b[j]]);
	    }
	}
    }
}


float hausdorff_dis(vector<myvector > &location_source, vector<myvector > &location_target)
{
  int i,j;
  int size_source=location_source.size();
  int size_target=location_target.size();
  //printf("%d %d \n",size_source,size_target);

  float dis_max=-1;
  for(i=0;i<size_source;i++)
    {
      float dis_min=100000.0;
      for(j=0;j<size_target;j++)
	{
	  float dis_now=(location_source[i]-location_target[j]).len();
	  if(dis_now<dis_min)
	    {
	      dis_min=dis_now;
	    }
	}
      if(dis_min>dis_max)
	{
	  dis_max=dis_min;
	}
    }
  //printf("dis_max: %f\n",dis_max);

  
  float dis_max_1=-1;
  for(i=0;i<size_target;i++)
    {
      float dis_min_1=100000.0;
      for(j=0;j<size_source;j++)
	{
	  float dis_now=(location_target[i]-location_source[j]).len();
	  if(dis_now<dis_min_1)
	    {
	      dis_min_1=dis_now;
	    }
	}
      if(dis_min_1>dis_max_1)
	{
	  dis_max_1=dis_min_1;
	}
    }

  printf("hausdorff %lf %lf \n",dis_max,dis_max_1);
  if(dis_max_1>dis_max)
    {
      return dis_max_1;
    }
  else
    {
      return dis_max;
    }
  
  return dis_max;
}

int recombine(boost::property_tree::ptree &para_tree)
{
  Py_Initialize();
  string cluster_method=para_tree.get<string>("cluster_method.value","");
  string python_path=para_tree.get<string>("python_path.value","");

  string chdir_cmd = string("sys.path.append(\"") + python_path + "\")";
  const char* cstr_cmd = chdir_cmd.c_str();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString(cstr_cmd);

  PyObject* pModule = PyImport_ImportModule("cluster");
  PyObject* pFunc = PyObject_GetAttrString(pModule, cluster_method.c_str());

  /////////////////////////////////////////////////////////////////
  object myobject_target;
  object myobject_source;
  decompose(para_tree,myobject_target,0,pFunc);
  decompose(para_tree,myobject_source,1,pFunc);
  /////////////////////////////////////////////////////////////////
  printf("result num_strip: %d %d\n",myobject_target.num_strip,myobject_source.num_strip);
   
  Py_Finalize();

  //////////////////////////////////////////////////////////////////
  //recombination from here

  int i,j,k;
  int use[myobject_source.num_cluster];
  int ct=0;

  vector<myvector > location_source[myobject_source.num_cluster];
  vector<myvector > location_target[myobject_target.num_cluster];

  for(i=0;i<myobject_source.num_vertex;i++)
    {
      int strip_now=myobject_source.myvertexs[i].strip;
      location_source[myobject_source.label_strip[strip_now]].push_back(myobject_source.myvertexs[i].location);
    }
  for(i=0;i<myobject_target.num_vertex;i++)
    {
      int strip_now=myobject_target.myvertexs[i].strip;
      location_target[myobject_target.label_strip[strip_now]].push_back(myobject_target.myvertexs[i].location);
    }
  
  for(i=0;i<myobject_source.num_cluster;i++)
    {
      printf("strip num:%d ",myobject_source.strip_cluster[i].size());
      if(myobject_source.strip_cluster[i].size()<5||location_source[i].size()<400)
	{
	  printf("no use here",i);
	}
      else
	{
	  use[ct]=i;
	  ct++;
	}
      printf("\n");
    }
  int num_cluster_use=ct;
  
  int conflict[myobject_source.num_cluster]; //mark which cluster of target is conflicting seriously
  float hausdorff_now[myobject_source.num_cluster];

  for(i=0;i<num_cluster_use;i++)
    {
      for(j=0;j<myobject_target.num_cluster;j++)
	{
	  float result_dis=hausdorff_dis(location_source[use[i]],location_target[j]);
	  if(j==0||(j>0&&result_dis<hausdorff_now[use[i]]))
	    {
	      hausdorff_now[use[i]]=result_dis;
	      conflict[use[i]]=j;
	      //printf("conflict %d haus %f\n",conflict[use[i]],hausdorff_now[use[i]]);
	    }	  
	}
    }

  for(i=0;i<num_cluster_use;i++)
    {
      printf("conflict %d: %d %lf\n",use[i],conflict[use[i]],hausdorff_now[use[i]]);
    }

  /*
  //https://www.jb51.net/article/54443.htm
  const int N = 4;
  const int M = 2;
  int a[N],b[N];
  for(int i=0;i<N;i++)
  a[i] = i+1;
  // 回溯方法
  combine(a,N,2,b,2);
  */
  vector<int > combine_help[num_cluster_use];  
  io* myio=new io();
  string combine_result_path=para_tree.get<string>("combine_result_path.value","");
  int b_help[num_cluster_use];
  int ct_obj=0;
  for(i=3;i<4;i++)// C(num_cluster_use,1),C(num_cluster_use,2),C(num_cluster_use,3)
    {
      combine(use,num_cluster_use,i,b_help,i,combine_help[i]);
      for(j=0;j<combine_help[i].size();j+=i)
	{
	  vector<int > cluster_now_used;
	  vector<int > cluster_now_used_target;
	  while(!cluster_now_used.empty())
	    {
	      cluster_now_used.pop_back();
	    }
	  while(!cluster_now_used_target.empty())
	    {
	      cluster_now_used_target.pop_back();
	    }
	  int conflict_tag[myobject_target.num_cluster];
	  for(k=0;k<myobject_target.num_cluster;k++)
	    {
	      conflict_tag[k]=-1;
	    }
	  for(k=j;k<j+i;k++)
	    {
	      cluster_now_used.push_back(combine_help[i][k]);
	      conflict_tag[conflict[combine_help[i][k]]]=1;
	      //printf("%d ",combine_help[i][k]);
	    }

	  for(k=0;k<myobject_target.num_cluster;k++)
	    {
	      if(conflict_tag[k]==-1)
		{
		  cluster_now_used_target.push_back(k);
		}
	    }
	  //printf("\n");

	  stringstream ss; string now_string; ss<<ct_obj; ss>>now_string;
	  string name=combine_result_path+"/result_"+now_string+".vtk";
	  myio->saveCombine(myobject_target,myobject_source,cluster_now_used_target,cluster_now_used,name);
	  ct_obj++;
	}
    }
  delete myio;
  
  return 0;
}



int main(int argc, char *argv[])
{
  // test access point
  boost::property_tree::ptree para_tree;
  readcmdline(argc,argv,para_tree);

  string prog=para_tree.get<string>("prog.value");
  if(prog=="decomposition")
    {
      Py_Initialize();
      string cluster_method=para_tree.get<string>("cluster_method.value","");
      string python_path=para_tree.get<string>("python_path.value","");
  
      string chdir_cmd = string("sys.path.append(\"") + python_path + "\")";
      const char* cstr_cmd = chdir_cmd.c_str();
      PyRun_SimpleString("import sys");
      PyRun_SimpleString(cstr_cmd);

      PyObject* pModule = PyImport_ImportModule("cluster");
      PyObject* pFunc = PyObject_GetAttrString(pModule, cluster_method.c_str());

      
      ///////////////////////////////////////////
      object myobject_now;
      decompose(para_tree,myobject_now,0,pFunc);      
      ///////////////////////////////////////////

      
      Py_Finalize();

    }
  else if(prog=="recombination")
    {     
      recombine(para_tree);
    }
  return 0;
}
