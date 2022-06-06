//#define DEBUG

#ifdef DEBUG
#include "FEM.h"
using namespace FEMns;
#else
#include "time.h"
#include "Filtration.h"
using namespace filtr;
#endif // DEBUG

int main()
{
#ifdef DEBUG
	FEM* fem = new FEM(nullptr);
	fem->SolveElliptic();
#else 
	std::clock_t start = clock();
	Filtration* filtrn = new Filtration();
	filtrn->Start();
	std::clock_t end = clock();
	std::cout << std::defaultfloat;
	std::cout << "Time taken: " << end - start << " ms/ " << (end - start) / 1000. << "secs/ " << (end - start) / 1000. / 60. << "min \nKnots num: " << filtrn->fem->GetKnotNum() << "\n";
	filtrn->fem->GetSolutionOnPlane(0);
#endif // DEBUG

}
