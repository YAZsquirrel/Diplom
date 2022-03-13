#pragma once
typedef double real;

namespace filtration
{
   struct component
   {
      int num;          // num of a component
      real proportion;  // in a phase
   };

   struct phase
   {
      real density;        // ro
      real viscosity;      // etta
      real penetrability;  // k (kappa)
      real saturation;     // Sm

      component* compsInPhase; 
      int L;                   // size of comp array  
   };


   real* flow_in_face;  // size face_size

   phase* phases;       // size M
   component* comps;    // size 1..3...

   real* alpha;   // size 
   real* beta;    // size knot_num
}
