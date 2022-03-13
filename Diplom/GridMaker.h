#pragma once
#include <vector>
#include <set>
//#include <>

typedef double real;

namespace mesh_comps
{

   struct knot
   {
      unsigned int knot_num;
      knot(real _x, real _y, real _z) : x(_x), y(_y), z(_z) {}
      knot(real xyz) : x(xyz), y(xyz), z(xyz) {}
      knot() : x(0.0), y(0.0), z(0.0) {}
      real x, y, z;
   };

   struct face
   {
      unsigned int face_num;
      unsigned int hexas_num[2]{};
      unsigned int knots_num[4]{};
      real face_area = 0.0;
      knot normal = knot(0.0);
      //real Vin, Vout;
   };

   struct hexahedron {
      unsigned int hexa_num;
      unsigned int knots_num[8]{};
      unsigned int faces_num[6]{};
      real dif = 0.0;
   };

   class Mesh
   {
      public: 
	    Mesh();

       std::vector<knot> knots;
       std::vector<face> faces;
       hexahedron** hexas;

       void GenerateMesh(); //
       private: 
       
       void GenerateKnots();
       void FindAllFaces();
       void FindAllHexas();
       void 
   };

}