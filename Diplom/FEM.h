#pragma once
//#define DEBUG
//#define DEBUG1
#include <iomanip>
#include <cmath>
#include <iostream>
#include <functional>
#include <list>
#include "Math_structs.h"
using namespace maths;
//using namespace filtr;
using namespace mesh_comps;
namespace FEMns
{
class FEM
{
   private:
   real f(knot &knot_);
   real ug(knot &knot_);
   int num_of_knots, num_of_FE, un;

   real localM2d[4][4];
   real localM[8][8]; // 8*8
   real localG[8][8];
   real localA[8][8];
   real reversed_J[3][3];
   real J[3][3];
   real J2D[2][2];
   real Jgrad_i[3];
   real gradi[3];
   real Jgrad_j[3];
   real gradj[3];
   inline real det_J();
   real prime_by_var(int what, int varOnFE, int knot_num[8], real ksi, real etta, real tetha);
   inline int mu(int index);
   inline int v(int index);
   inline int nu(int index);
   real W(int index, real alpha);
   real d_phi(int index, int what, real ksi, real etta, real tetha);
   inline real phi(int index, real ksi, real etta, real tetha);
   void calc_grad(int ij, int index, real ksi, real etta, real tetha);

   void check_test();
   void AddFirstBounds();
   void AddSecondBounds();
   void AddToA(element &hexa);
   void CreateSLAE();
   void CreateM(element &hexa);
   void CreateG(element &hexa); 
   void Createb(element &hexa);

   std::shared_ptr<Matrix> A;
   std::vector<real>b;
   std::vector<real> q; 
   std::shared_ptr<Mesh> mesh;
   
   real Integrate(const std::function<real(real, real, real, int, int, int[8])> f, int i, int j, int knot_num[8]);

   //real Integrate2D(const std::function<real(real, real, int, int, int[4])> f, int i, int j, int knot_num[4]);

   std::function<real(real, real, real, int, int, int[8])> Gij;
   std::function<real(real, real, real, int, int, int[8])> Mij;

   public:
   int GetKnotsNum () {return num_of_knots;}
   int GetHexasNum () {return num_of_FE;}

   std::vector<real>& GetKnots() { return q; };
   FEM(std::shared_ptr<Mesh> _mesh);
   void SolveElliptic();
   void GetSolutionOnPlane(real z);
   void Output(std::ofstream &out);
};
}