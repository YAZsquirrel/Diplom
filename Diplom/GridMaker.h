#pragma once
#include <vector>
#include <set>
#include <map>
#include <list>
typedef double real;

namespace mesh_comps
{
   struct bound {
      int face_num;
      int knots_num[4];
      real value;
   };
   struct knot
   {
      //unsigned int knot_num;
      knot(real _x, real _y, real _z) : x(_x), y(_y), z(_z) {}
      knot() : x(0.0), y(0.0), z(0.0) {}
      ~knot(){}
      real x, y, z;
   };

   struct well
   {
      real x, y;
      real radius;
      real h1, h2; 
      real intake;
      std::vector<int*> faces_knots_num;
      std::vector<int> bound_num;
      well(real _x, real _y, real rad, real _h1, real _h2, real th) :
         x(_x), y(_y), radius(rad), h1(_h1), h2(_h2), intake(th) {}
   };

   struct face 
   {
      std::vector<int> hexa_nums;
      int knots_num[4];
      knot normal;
   };

   struct hexahedron
   {
      int knots_num[8];
      int faces_sign[6];
      std::vector<int> faces_num;
      std::vector<int> phases_num;
      std::vector<real> Sm; // 
      real lam = 0.;
      // n (etta) viscosity
      // k (kappa) penetrability

      std::set<int> heighbors;
      int containsFace(int n)
      {
         bool found = false;
         int i = 0;
         for (;i < 6 && !found; i++)
            found = faces_num[i] == n;
         return found ? i-- : -1;

      }
      bool containsKnot(int n)
      {
         bool found = false;
         for (int i = 0; i < 8 && !found; i++)
            found = knots_num[i] == n;
         return found;
      }
   };


   class Mesh
   {
      public: 
	    Mesh();

       std::vector<knot*> knots;
       std::vector<hexahedron*> hexas;
       std::vector<face*> faces;
       std::vector<std::set<int>> neighbors;
       std::vector<bound*> bounds1;
       std::vector<bound*> bounds2;

       void FindNeighborsAndFaces();

       void set_layers()
       {
            for (int i = 0; i < w_info.wells.size(); i++)
            {
               bool found1 = false;
               bool found2 = false;

               for (auto l : layers)
               {
                  if (!found1) 
                     found1 = abs(l - w_info.wells[i].h1) < 1e-10;
                  if (!found2) 
                     found2 = abs(l - w_info.wells[i].h2) < 1e-10;
               }
               if (!found1)
                  layers.push_back(w_info.wells[i].h1);
               if (!found2)
                  layers.push_back(w_info.wells[i].h2);
            }
       }

       
       void FindFaceNormals();
       void GenerateMesh();
       void SetSignsForHexaNormals();
       struct well_info
       {
          std::vector<well> wells;
          real conc_rad, conc;
          int rad_knots;
       } w_info;

       private: 

       short sign(real x) { return -(x < 0.) + (x > 0); }
       int xn = 0, yn = 0, zn = 0;
       std::vector<real> layers;
       knot env_corner1, env_corner2, step;

       void FindAllHexasAndBounds(int plain_size, int* well_inds, std::vector<int>** inwell_indecies, real* zs);
       //void 
   };

}