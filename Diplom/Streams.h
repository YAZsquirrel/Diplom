#pragma once
#include <vector>
typedef double real;

#include <functional>
#include "Mat_structs.h"
using namespace mats;

typedef double real;

namespace streams
{

	class Streams
	{
		public:
			std::vector<real> Q;
			Streams(mesh_comps::Mesh* mesh, std::vector<real> &P);
			void FindStreams();
		private:
			std::function<real(real, real, int, int, int, int)> G2Dij;
			real Integrate2D(const std::function<real(real, real, int, int, int, int)> f, int face_num, int e, int xyz, int opposite);
			real prime_by_var(int what, int varOnFE, int nFE, real ksi, real etta, real theta);
			inline int mu(int index);
			inline int v(int index);
			inline int nu(int index);
			real W(int index, real alpha);
			real d_phi(int index, int what, real ksi, real etta, real tetha);
			void calc_grad(int index, int Pindex, real ksi, real etta, int xyz, real opposite);
			void AssembleMatrix();
			void AssembleRightPart();
			bool CheckAccuracy(real eps);

			void FindAverageStreams();
			void BalanceStreams();
			mesh_comps::Mesh* mesh;
			real reversedJ2D[2][2];
			real J2D[2][2];
			real Jgrad[2];
			real grad[2];
			std::vector<real> P;
			//std::vector<std::vector<real>> B;
			Matrix *B;
			std::vector<real> d;
			std::vector<real> Q_av;
			std::vector<real> dQ;
			std::vector<real> alpha;
			std::vector<real> beta;

	};
}