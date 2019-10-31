#ifndef _HEAD_
#define _HEAD_

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <sstream>
#include <random>
#include <time.h>

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/CholmodSupport>

#include <SymGEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/info_parser.hpp>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh> 


#include "Python.h"
#include <thread>

#define pi 3.1415926535898
#define min(a,b) (((a) < (b)) ? (a) : (b))

#endif
