#ifndef _VERTEX_
#define _VERTEX_

#include "head.h"
#include "myvector.h"

class vertex
{
 public:
  myvector location;
  float u_x,u_y;
  int strip;
  vertex();
  vertex(const myvector &location);
  ~vertex();
};

#endif
