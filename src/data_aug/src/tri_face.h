#ifndef _TRI_FACE_
#define _TRI_FACE_
#include "myvector.h"

class tri_face
{
 public:
  int index_vertex[3];
  tri_face(){}
  ~tri_face(){}
  tri_face(const int (&index_vertex)[3]);
};
#endif
