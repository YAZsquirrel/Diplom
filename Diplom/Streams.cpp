#include "Streams.h"

namespace streams
{
   Streams::Streams(mesh_comps::Mesh* _mesh, std::vector<real>& _P)
   {
      mesh = _mesh; 
      
      Q.resize(mesh->faces.size());
      d.resize(mesh->faces.size());
      Q_av.resize(mesh->faces.size());
      dQ.resize(mesh->faces.size());
      alpha.resize(mesh->faces.size(), 1e+10);
      beta.resize(mesh->hexas.size(), 1e+10);
      //B.reserve(mesh->faces.size());
      //for (int i = 0; i < mesh->faces.size(); i++)
      //   B[i].reserve(mesh->faces.size());
      P = _P;

      // TODO: проверить
      // 
      // 
      G2Dij = [this](real ksi, real etta, int face_num, int e, int xyz, int opposite) {

         for (int ip = 0; ip < 2; ip++)
            Jgrad[ip] = grad[ip]= 0.;
         
         for (int ip = 0; ip < 2; ip++)                                      
            for (int jp = 0; jp < 2; jp++)
               switch (xyz)
               {
               case 0: J2D[ip][jp] = prime_by_var(ip, jp, e, ksi, etta, (real)opposite); break;
               case 1:
                  J2D[0][0] = prime_by_var(0, 0, e, ksi, (real)opposite, etta);
                  J2D[0][1] = prime_by_var(0, 2, e, ksi, (real)opposite, etta);
                  J2D[1][0] = prime_by_var(2, 0, e, ksi, (real)opposite, etta);
                  J2D[1][1] = prime_by_var(2, 2, e, ksi, (real)opposite, etta); break;
               case 2: J2D[ip][jp] = prime_by_var(ip + 1, jp + 1, e, (real)opposite, ksi, etta); break;

               }
         reversedJ2D[0][0] = J2D[1][1];
         reversedJ2D[1][1] = J2D[0][0];
         reversedJ2D[1][0] = -J2D[1][0];
         reversedJ2D[0][1] = -J2D[0][1];
         
         for (int i = 0; i < 8; i++)
            calc_grad(i, mesh->hexas[e]->knots_num[i], ksi, etta, xyz, (real)opposite);
         
         real n[2];
         switch (xyz)
         {
         case 2: n[0] = mesh->faces[face_num]->normal.x;
                 n[1] = mesh->faces[face_num]->normal.y; break;
         case 1: n[0] = mesh->faces[face_num]->normal.x;
                 n[1] = mesh->faces[face_num]->normal.z; break;
         case 0: n[0] = mesh->faces[face_num]->normal.y;
                 n[1] = mesh->faces[face_num]->normal.z; break;

         }
         Jgrad[0] = reversedJ2D[0][0] * grad[0] + reversedJ2D[0][1] * grad[1];
         Jgrad[1] = reversedJ2D[1][0] * grad[0] + reversedJ2D[1][1] * grad[1];

         real res = Jgrad[0] * n[0] + Jgrad[1] * n[1];
         return res * abs((J2D[0][0]* J2D[1][1]) - (J2D[0][1] * J2D[1][0])); //J^-1 = 1/|J| * ...
      };
   }

   void Streams::FindStreams()
   {
      FindAverageStreams();
      BalanceStreams();
   }

   void Streams::FindAverageStreams()
   {
      real w = 1./2.;
      std::vector<bool> visited;
      visited.resize(mesh->faces.size(), false);

      for (int e = 0; e < mesh->hexas.size(); e++)
      {
         mesh_comps::hexahedron* hexa = mesh->hexas[e];
         for (int f = 0; f < 6; f++)
         {
            int xyz = f/2;
            int opposite = f%2;  
            int fn = hexa->faces_num[f];

            mesh_comps::face* face = mesh->faces[fn];

            real Qei = -Integrate2D(G2Dij, fn, e, xyz, opposite) * hexa->faces_sign[f] * mesh->hexas[e]->lam;
            if (visited[fn])
               Q_av[fn] = (1. - w) * Qei + w * Q_av[fn];
            else
            { 
               Q_av[fn] = Qei;
               visited[fn] = true;
            }
         }
      }
   }

   void Streams::BalanceStreams()
   {
      B = MakeSparseFormat(6, mesh->faces.size(), mesh->hexas.size(), false, mesh);
      int k = 0;
      AssembleRightPart();
      CheckStreams();
      do {

         AssembleMatrix();
         SolveSLAE(B, dQ, d);
         //adjust parameters
         CheckStreams();
         k++;
      } while (CheckAccuracy(1e-10));
   }
   
   void Streams::AssembleMatrix()
   {
      for (int i = 0; i < B->ig[B->dim]; i++)
      {
         B->l[i] = 0.;
         B->u[i] = 0.;
      }

      for (int i = 0; i < B->dim; i++)
         B->di[i] = 0.;

      for (int e = 0; e < mesh->hexas.size(); e++)
      {
         auto& el = mesh->hexas[e];
         for (int i = 0; i < 6; i++)
         {
            int gfi = el->faces_num[i];
            for (int j = 0; j < 6; j++)
            {
               int gfj = el->faces_num[j];
               if (i == j) continue;
               AddElement(B, i, j, beta[e] * el->faces_sign[i] * el->faces_sign[j]);
            }
         }
      }
      for (int i = 0; i < B->dim; i++)
         B->di[i] = 2. + alpha[i];

      for (auto cond : mesh->bounds2)
      {
         B->di[cond->face_num] = 1.;
         for (int j = B->ig[cond->face_num]; j < B->ig[cond->face_num + 1]; j++)
            B->l[j] = 0.;
         for (int j = 0; j < B->ig[B->dim]; j++)
            if (B->jg[j] == cond->face_num)
              B->u[j] = 0.;
         d[cond->face_num] = 0.;
      }
   }

   void Streams::AssembleRightPart()
   {
      for (int f = 0; f < mesh->faces.size(); f++)
      {
         auto& face = mesh->faces[f];
         bool noneigh = face->hexa_nums.size() != 1;
         auto& el1 = mesh->hexas[face->hexa_nums[0]];
         hexahedron* el2 = nullptr;
         if (noneigh) el2 = mesh->hexas[face->hexa_nums[1]];

         real sum1 = 0.;
         real sum2 = 0.;
         for (int i = 0; i < 6; i++)
         {
            int gfi1 = el1->faces_num[i];
            int gfi2 = noneigh ? el2->faces_num[i] : -1;
            sum1 += Q_av[gfi1];
            if (noneigh) sum2 += Q_av[gfi2];
         }
         sum1 *= el1->faces_sign[el1->containsFace(f)];
         if (noneigh) sum2 *= el2->faces_sign[el2->containsFace(f)];
         d[f] = -(sum1 + sum2);
      }
   }

   bool Streams::CheckAccuracy(real eps)
   {
      real sum = 0.;
      for (int e = 0; e < mesh->hexas.size(); e++)
      {
         auto& el = mesh->hexas[e];
         real sumQ = 0.;
         for (int i = 0; i < 6; i++)
            sumQ += el->faces_sign[i] * (abs(Q_av[el->faces_num[i]]) + dQ[el->faces_num[i]]);
         sum += abs(sumQ) * beta[e];
      }
      std::cout << "\nDisbalance = " << sum << '\n';
      return sum >= eps;
   }

   void Streams::CheckStreams()
   {
      std::cout << '\n' << "Summary streams:\n";
      for (int e = 0; e < mesh->hexas.size(); e++)
      {
         auto& el = mesh->hexas[e];
         real sum = 0.;
         for (int f = 0; f < 6; f++)
            sum += el->faces_sign[f] * (abs(Q[el->faces_num[f]]) + dQ[el->faces_num[f]]);
         std::cout <<  e << "\t:\t" << sum << '\n';
      }
      std::cout << '\n';
   }

   real Streams::Integrate2D(const std::function<real(real, real, int, int, int, int)> f, int face_num, int e, int xyz, int opposite)
   {
      const int nKnot = 3;//5; // Knots num

      const real xj[nKnot]
         = { .7745966692414833, 0., -.7745966692414833 }; 
      const real qj[nKnot]
         = { .55555555555555555, .8888888888888888, .55555555555555555 };

      real result = 0.;
      for (int ix = 0; ix < nKnot; ix++)
         for (int iy = 0; iy < nKnot; iy++)
            result += qj[ix] * qj[iy] * (f(.5 + xj[ix] * .5, .5 + xj[iy] * .5, face_num, e, xyz, opposite));
      return result / 4.;
   }

   int Streams::mu(int index)
   {
      return index % 2;
   }

   int Streams::v(int index)
   {
      return (index / 2) % 2;
   }

   int Streams::nu(int index)
   {
      return index / 4;
   }

   real Streams::W(int index, real alpha)
   {
      if (!index) return 1. - alpha;
      return alpha;
   }

   real Streams::d_phi(int index, int var, real ksi, real etta, real tetha)
   {
      real d_phi = 0.;
      switch (var)
      {
      case 0:    // ksi
      {
         d_phi = W(v(index), etta) * W(nu(index), tetha);
         if (!mu(index)) d_phi *= -1;
         break;
      }
      case 1:     // etha
      {
         d_phi = W(mu(index), ksi) * W(nu(index), tetha);
         if (!v(index)) d_phi *= -1;
         break;
      }
      case 2:     // theta
      {
         d_phi = W(mu(index), ksi) * W(v(index), etta);
         if (!nu(index)) d_phi *= -1;
         break;
      }
      }

      return d_phi;
   }

   real Streams::prime_by_var(int varOnCube, int varOnFE, int nFE, real ksi, real etta, real theta)
   {
      real var = 0.;
      for (int i = 0; i < 8; i++)
      {
         switch (varOnFE)
         {
         case 0: var += mesh->knots[mesh->hexas[nFE]->knots_num[i]]->x * d_phi(i, varOnCube, ksi, etta, theta); break;
         case 1: var += mesh->knots[mesh->hexas[nFE]->knots_num[i]]->y * d_phi(i, varOnCube, ksi, etta, theta); break;
         case 2: var += mesh->knots[mesh->hexas[nFE]->knots_num[i]]->z * d_phi(i, varOnCube, ksi, etta, theta); break;
         }
      }
      return var;
   }

   void Streams::calc_grad( int index, int Pindex, real ksi, real etta, int xyz, real opposite)
   {
      switch (xyz)
      {
      case 0: 
         grad[0] += P[Pindex] * d_phi(index, 0, ksi, etta, opposite);
         grad[1] += P[Pindex] * d_phi(index, 1, ksi, etta, opposite); break;
      case 1:
         grad[0] += P[Pindex] * d_phi(index, 0, ksi, opposite, etta);
         grad[1] += P[Pindex] * d_phi(index, 2, ksi, opposite, etta); break;
      case 2:
         grad[0] += P[Pindex] * d_phi(index, 1, opposite, ksi, etta);
         grad[1] += P[Pindex] * d_phi(index, 2, opposite, ksi, etta); break;
      }
         
   }
}