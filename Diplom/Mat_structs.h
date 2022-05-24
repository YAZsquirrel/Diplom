#pragma once
#include "GridMaker.h"
#include <iostream>
typedef double real;
using namespace mesh_comps;
namespace mats{

   struct Matrix {
      real* l, * u, * di;
      int* ig, * jg;
      size_t dim;
   };

   Matrix* MakeSparseFormat(int localsize, int elemsize, int FEsize, bool isknots, mesh_comps::Mesh* mesh);
   void copy(std::vector<real>& v, std::vector<real>& u);
   real scalar(std::vector<real> &v, std::vector<real> &u);
   
   void AddElement(Matrix* M, int i, int j, real elem);
   void MatxVec(std::vector<real> &v, Matrix* A, std::vector<real>& b);
   void SolveSLAE(Matrix* M, std::vector<real>& q, std::vector<real> &b);
}