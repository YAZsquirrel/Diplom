#pragma once
#include "GridMaker.h"
#include <iostream>
#include <fstream>

typedef double real;
using namespace mesh_comps;
namespace mats{

   struct Matrix {
      std::vector<real> l, u, di;
      std::vector<int> ig, jg;
      size_t dim;
   };

   Matrix* MakeSparseFormat(int localsize, int elemsize, int FEsize, bool isknots, mesh_comps::Mesh* mesh);
   void copy(std::vector<real>& v, std::vector<real>& u);
   real scalar(std::vector<real> &v, std::vector<real> &u);
   void AddElement(Matrix* M, int i, int j, real elem);
   void MatxVec(std::vector<real> &v, Matrix* A, std::vector<real>& b);
   void SolveSLAE(Matrix* M, std::vector<real>& q, std::vector<real> &b);
   void WriteMatrix(Matrix* M);
   void MatSymmetrisation(Matrix* M, std::vector<real>& b, int i);
}