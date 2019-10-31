#ifndef _IO_
#define _IO_
#include "head.h"
#include "vertex.h"
#include "tri_face.h"
#include "object.h"
class io
{
 public:
  io(){}
  ~io(){}
  //name 不用引用，因为可能实参直接是"xxx" 
  int saveArrayAsVTK(std::vector<vertex> &myvertexs,std::vector<tri_face> &mytris,const std::string name);
  
  int getVertexAndTri(std::vector<vertex> &myvertexs,std::vector<tri_face > &mytris,const std::string name);
  int getVertexAndTriFromObj(std::vector<vertex> &myvertexs,std::vector<tri_face > &mytris,const std::string name);

  int savePolylineAsVTK(std::vector<vertex> &polyline,const std::string name);

  int saveArrayAsVTKwithLabel(int num_strip,int num_tri,int which_vertex[],int label_strip[],std::vector<tri_face > mp[],std::vector<vertex > &myvertexs, const std::string name);
  //int saveArrayAsVTKwithLabel(int num_strip,int num_tri,int (&which_vertex)[],int (&label_strip)[],std::vector<tri_face > (&mp)[],std::vector<vertex > &myvertexs, const std::string name);

  int saveClusterPolylineAsVTK(std::vector<int > &strip_cluster,std::vector<vertex> polyline[],std::vector<vertex > &myvertexs,const std::string name);
  int saveClusterAsVTK(std::vector<int > &strip_cluster,std::vector<vertex > &myvertexs,int which_vertex[],std::vector<tri_face > mp[],const std::string name);

  int saveCombine(object &myobject_target,object &myobject_source,std::vector<int > &cluster_now_used,const std::string name);

  int saveCombine(object &myobject_target,object &myobject_source,std::vector<int > &cluster_now_used_target,std::vector<int > &cluster_now_used,const std::string name);
};
#endif
