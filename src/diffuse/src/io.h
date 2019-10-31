#ifndef _IO_
#define _IO_
#include "head.h"
class io
{
 public:
  io(){}
  ~io(){}
  int saveAsVTK(Eigen::MatrixXd &dense_theta_map,Eigen::MatrixXd &whether_hair,int height,int width,const std::string name);
  int getData(Eigen::MatrixXd &sparse_theta_map,Eigen::MatrixXd &whether_hair,int height,int width,const std::string name);
};
#endif
