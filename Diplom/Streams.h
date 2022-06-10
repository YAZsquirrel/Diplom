#pragma once
//#define DEBUG0

#include <vector>
#include <cmath>
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
			void OutputStreams(std::ofstream& out, std::ofstream& outone, int face_num_to_out);
		private:
			short sign(real x) {
				return -(x < 0.) + (x > 0);
			}
			real clamp(real max, real min, real val) {
				return (val > max) * max + (val < min) * min + (val < max && val > min) * val;
			}
			void SetAlpha();

			std::function<real(real, real, int, int, int, int)> G2Dij;
			real Integrate2D(const std::function<real(real, real, int, int, int, int)> f, int face_num, int e, int xyz, int opposite);
			real prime_by_var(int what, int varOnFE, int nFE, real ksi, real etta, real theta);
			real prime_by_var2D(int what, int varOnFace, int nFace, real u, real v, real x[4], real y[4]);
			inline int mu(int index);
			inline int v(int index);
			inline int nu(int index);
			real W(int index, real alpha);
			real d_phi(int index, int what, real ksi, real etta, real tetha);
			real d_phi2D(int index, int what, real u, real v);
			real det_J3D();
			void calc_grad(int index, int Pindex, real ksi, real etta, real tetha);

			void AssembleMatrix();
			void AssembleRightPart();
			bool AdjustBeta(real eps);					

			void FindAverageStreams();
			void BalanceStreams();
			mesh_comps::Mesh* mesh;
			real reversedJ2D[2][2];
			real J2D[2][2];
			real reversedJ3D[3][3];
			real J3D[3][3];
			real Jgrad[3];
			real grad[3];
			std::vector<real> P;
			//std::vector<std::vector<real>> B;
			Matrix *B;
			std::vector<real> d;
			std::vector<real> Q_av;
			std::vector<real> dQ;
			std::vector<real> alpha;
			std::vector<real> beta;
			//std::vector<real> disbalance;

	};
}