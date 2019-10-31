#ifndef _MYVECTOR_
#define _MYVECTOR_

#include "head.h"

class myvector
{
 public :
 float x,y,z;
 myvector():x(0),y(0),z(0){}
 myvector(const myvector&a):x(a.x),y(a.y),z(a.z){}
 myvector(float x,float y,float z):x(x),y(y),z(z){}
  ~myvector() {}
  void set(float x,float y,float z)
  {
    this->x=x; this->y=y; this->z=z;
    return ;
  }

  float& operator()(int key)
  {
    switch(key)
      {
      case 0:
	return x;
	break;
      case 1:
	return y;
	break;
      case 2:
	return z;
	break;
      }
  }
	
  myvector& operator=(const myvector& a)
    {
      this->x=a.x; this->y=a.y; this->z=a.z; 
      return *this;
    }
	
  myvector& operator+=(float a)
    {
      this->x+=a; this->y+=a; this->z+=a;
      return *this;
    }
	
  myvector& operator-=(float a)
    {
      this->x-=a; this->y-=a; this->z-=a;
      return *this;
    }
	
  myvector& operator*=(float a)
    {
      this->x*=a; this->y*=a; this->z*=a;
      return *this;
    }
	
  myvector& operator/=(float a)
    {
      this->x/=a; this->y/=a; this->z/=a;
      return *this;
    }
	
  myvector& operator+=(const myvector& a)
    {
      this->x+=a.x; this->y+=a.y; this->z+=a.z;
      return  *this;
    }
	
  myvector& operator-=(const myvector& a)
    {
      this->x-=a.x; this->y-=a.y; this->z-=a.z;
      return  *this;
    }
	
  float dot(const myvector& a) 
  {
    return this->x*a.x+this->y*a.y+this->z*a.z;
  }
	
  myvector cross(const myvector& a)
  {
    return myvector(this->y*a.z-this->z*a.y,this->z*a.x-this->x*a.z,this->x*a.y-this->y*a.x);
  }

  float len_sq(void)
  {
    return x*x+y*y+z*z;
  }
	
  float len(void)
  {
    return sqrt(x*x+y*y+z*z);
  }

  void normalize(void)
  {
    float len_sq=x*x+y*y+z*z;
    if(len_sq!=0.0) 
      {
	float len=sqrt(len_sq);
	x/=len; y/=len; z/=len;
      }
  }
	
};

inline myvector operator + (const myvector& a, const myvector& b)
{
  return myvector(a.x+b.x, a.y+b.y, a.z+b.z);
}

inline myvector operator - (const myvector& a, const myvector& b)
{
  return myvector(a.x-b.x, a.y-b.y, a.z-b.z);
}

inline myvector operator * (float s, const myvector& a)
{
  return myvector(s*a.x, s*a.y, s*a.z);
}

inline myvector operator * (const myvector& a,float s)
{
  return myvector(s*a.x, s*a.y, s*a.z);
}

inline myvector operator + (float s, const myvector& a)
{
  return myvector(s+a.x, s+a.y, s+a.z);
}

inline myvector operator + (const myvector& a,float s)
{
  return myvector(s+a.x, s+a.y, s+a.z);
}

inline myvector operator / (const myvector& a,float s)
{
  return myvector(a.x/s, a.y/s, a.z/s);
}

#endif
