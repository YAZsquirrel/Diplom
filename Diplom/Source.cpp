#include "FEM.h"
using namespace filtration;
#include "time.h"
int main()
{
	std::clock_t start = clock();
	FEM *fem = new FEM();
	fem->SolveElliptic();
	std::clock_t end = clock();
	std::cout << "Time taken: " << end - start << " ms/ " << (end - start)/1000. << "secs/ " << (end - start) / 1000. / 60. << "min \nKnots num: " << fem->GetKnotNum() << "\n";
	fem->GetSolutionOnPlane(51);


}
