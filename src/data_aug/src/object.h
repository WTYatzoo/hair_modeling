#ifndef _OBJECT_
#define _OBJECT_

#include "head.h"
#include "myvector.h"
#include "vertex.h"
#include "tri_face.h"

class object
{
 public:
  std::vector<vertex > myvertexs; //always
  std::vector<tri_face > mytris; //always
  int num_vertex,num_tri; //always
  std::vector<tri_face > mp[15000]; // always 15000 is max_tolerance of num_vertex
  int num_strip; //always
  std::vector<vertex > polyline[200];  //always 200 is max_tolerance of num_strip
  float radius_strip[200]; //always
  int which_vertex[200]; //always
  float dis[200][200]; // save the distance measure between strip a & strip b //always

  int label_strip[200]; //always
  int num_cluster; //always
  std::vector<int > strip_cluster[20]; //always 20 is max_tolerance of num_cluster
  object();
  ~object();
};
#endif
