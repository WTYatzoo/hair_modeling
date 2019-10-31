#include "tri_face.h"

tri_face::tri_face(const int (&index_vertex)[3])
{
  int i;
  for(i=0;i<3;++i)
    {
      this->index_vertex[i]=index_vertex[i];
    }
  return;
}
