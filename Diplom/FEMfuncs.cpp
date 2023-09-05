#include "FEM.h"
#include <math.h>

//real FEM::lambda(real knot[2], real t)
//{
//	return real();
//}
namespace FEMns{
real FEM::f(knot &knot)
{
	real x = knot.x;
	real y = knot.y;
	real z = knot.z;
	//return 0;
	switch (un)
	{
	case 1: return 0;
	case 2: return 0;
	case 3: return 0;
	case 4: return -6 * x;
	case 5: return sin(x); //

	default:
		return 0;
	}
}

//real FEM::ug(knot &knot)
//{
//	real x = knot->x;
//	real y = knot->y;
//	real z = knot->z;
//
//	switch (un)
//	{
//	case 1: return 1;
//	case 2: return x + y + z;
//	case 3: return x;
//	case 4: return x * x * x;
//	case 5: return sin(x); // x^3 + y^3 + xy + 1
//
//	default:
//		return 0;
//	}
//}

}
