#include "vertex.h"
using namespace std;

vertex::vertex()
{
  ;
}

vertex::vertex(const myvector &location)
{
  this->location=location;
  this->u_x=-1;
  this->u_y=-1;
  this->strip=-1;
}

vertex::~vertex()
{
  ;
}
