#pragma once
#include "GridMaker.h"
#include "FEM.h"
#include "Streams.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
typedef double real;
//using namespace mesh_comps;
using namespace streams;

namespace filtr
{
   //struct component
   //{
   //   std::vector<std::vector<real>> proportion; // M * L  - кол-во
   //};

   struct phase
   {
      real h;
      real n;      // n (etta) viscosity
      real k;      // k (kappa) penetrability
      //static int L;                   // size of comp array  
   };

   struct solid
   {
      real K;
      real Fi;
      knot corner0;
      knot corner1;
   };

   class Filtration
   {
      public:
      //Filtration();
      std::vector<solid> pors;
      Mesh* mesh;
      FEMns::FEM* fem;
      Streams* str;

      //real* flow_in_face;  // size face_size
      std::vector<phase> phases;
      std::vector<std::vector<phase>> elem_phases;       // size M * Ne
      //component comps;    // size 1..2..3....
      std::vector<std::vector<real>> comps_in_phases;    // size M * 1..2..3....

      void SetDiffKoef();
      void Start();
   };
}
