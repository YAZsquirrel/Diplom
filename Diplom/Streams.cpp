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
      alpha.resize(mesh->faces.size());
      beta.resize(mesh->hexas.size(), 1e+10);
      //disbalance.resize(mesh->hexas.size(), 0.0);
      //B.reserve(mesh->faces.size());
      //for (int i = 0; i < mesh->faces.size(); i++)
      //   B[i].reserve(mesh->faces.size());
      P = _P;

      G2Dij = [this](real u, real v, int face_num, int e, int xyz, int opposite) {

         for (int ip = 0; ip < 3; ip++)
            Jgrad[ip] = grad[ip] = 0.;

         for (int ip = 0; ip < 3; ip++)
            for (int jp = 0; jp < 3; jp++)
               switch (xyz)
               {
                  case 1: // x
                     J3D[ip][jp] = prime_by_var(ip, jp, e, opposite, u, v); break;
                  case 2: // y
                     J3D[ip][jp] = prime_by_var(ip, jp, e, u, opposite, v); break;
                  case 0: // z
                     J3D[ip][jp] = prime_by_var(ip, jp, e, u, v, opposite); break;
               }
         real x[4];
         real y[4];
         for (int i = 0; i < 4; i++)
         {
            switch (xyz)
            {
               case 1:
                  x[i] = mesh->knots[mesh->faces[face_num]->knots_num[i]]->y;
                  y[i] = mesh->knots[mesh->faces[face_num]->knots_num[i]]->z;
                  break;
               case 2:
                  x[i] = mesh->knots[mesh->faces[face_num]->knots_num[i]]->x;
                  y[i] = mesh->knots[mesh->faces[face_num]->knots_num[i]]->z; break;
               case 0:
                  x[i] = mesh->knots[mesh->faces[face_num]->knots_num[i]]->x;
                  y[i] = mesh->knots[mesh->faces[face_num]->knots_num[i]]->y; break;
            }
         }

         for (int ip = 0; ip < 2; ip++)
            for (int jp = 0; jp < 2; jp++)
               J2D[ip][jp] = prime_by_var2D(ip, jp, face_num, u, v, x, y);                  

         for (int ip = 0; ip < 3; ip++)
            for (int jp = 0; jp < 3; jp++)
            {
               real min[4]{};
               int k = 0;
               for (int im = 0; im < 3; im++)
                  for (int jm = 0; jm < 3; jm++)
                  {
                     if (im != ip && jm != jp)
                        min[k++] = J3D[im][jm];
                  }
               reversedJ3D[jp][ip] = ((ip + jp + 2) % 2 ? -1 : 1) * (min[0] * min[3] - min[1] * min[2]);
            }

         for (int i = 0; i < 8; i++)
            switch(xyz)
            {
               case 1: // x
                  calc_grad(i, mesh->hexas[e]->knots_num[i], opposite, u, v); break;
               case 2: // y
                  calc_grad(i, mesh->hexas[e]->knots_num[i], u, opposite, v); break;
               case 0: // z
                  calc_grad(i, mesh->hexas[e]->knots_num[i], u, v, opposite); break;
            }

         real n[3] = { mesh->faces[face_num]->normal.x, mesh->faces[face_num]->normal.y, mesh->faces[face_num]->normal.z };

         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
               Jgrad[i] += reversedJ3D[i][j] * grad[j];

         real res = 0.;
         for (int i = 0; i < 3; i++)
            res += Jgrad[i] * n[i];
         // ^3D


         return res * abs((J2D[0][0] * J2D[1][1]) - (J2D[0][1] * J2D[1][0])) / det_J3D(); //J^-1 = 1/|J| * ...
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
            int xyz;
            int opposite = f%2;  

            int fn = hexa->faces_num[f];
            mesh_comps::face* face = mesh->faces[fn];

            real xn = face->normal.x,
                 yn = face->normal.y,
                 zn = face->normal.z;
            
            if (abs(xn) > abs(zn))
               xyz = (abs(xn) > abs(yn)) ? 1 : 2;
            else
               xyz = (abs(zn) > abs(yn)) ? 0 : 2;
            
            //xyz = f/2;
            std::vector<real> norm = { xn, yn, zn };

            real Qei = -Integrate2D(G2Dij, fn, e, xyz, opposite) * hexa->faces_sign[f] * mesh->hexas[e]->lam;
            if (visited[fn])
               Q_av[fn] = (1. - w) * Qei + w * Q_av[fn] * hexa->faces_sign[f];
            else
            { 
               Q_av[fn] = Qei * hexa->faces_sign[f];
               visited[fn] = true;
            }
         }
      }
#ifdef DEBUG2
      for (int e = 0; e < mesh->hexas.size(); e++)
      {
         mesh_comps::hexahedron* hexa = mesh->hexas[e];
         std::cout << "\n-------------------\n";
         std::cout << e << " : \n";
         for (int f = 0; f < 6; f++)
         {
            int fn = hexa->faces_num[f];
            std::cout << "\t" << "(" << f << ") " << fn << " : " << Q_av[fn] << '\n';
         }

      }
#endif // DEBUG


   }

   void Streams::BalanceStreams()
   {
      B = MakeSparseFormat(6, mesh->faces.size(), mesh->hexas.size(), false, mesh);
      int k = 0;
#ifdef DEBUG0
      WriteMatrix(B);
      for (auto b : mesh->bounds2)
         std::cout << Q_av[b->face_num] << "\n";
#endif
      //SetKnownFlows();
      SetAlpha();
#ifdef DEBUG0
      for (auto b : mesh->bounds2)
         std::cout << Q_av[b->face_num] << "\n";
#endif     
      while (AdjustBeta(1e-9))
      {
         AssembleMatrix();
         AssembleRightPart();
         SolveSLAE(B, dQ, d);
         k++;
      } 

      for (int f = 0; f < mesh->faces.size(); f++)
         Q[f] = Q_av[f] + dQ[f];  
   }
   

   void Streams::SetKnownFlows()
   {
      for (auto w : mesh->w_info.wells)
      {
         
         real absQsum = 0.;
         real Qsum = 0.;
         real messum = 0.;
         for (int i = 0; i < w.bound_num.size(); i++)
         {
            if (abs(mesh->bounds2[w.bound_num[i]]->value) < 1e-12) continue;
            int f = mesh->bounds2[w.bound_num[i]]->face_num;
            real mes = 0.;
            absQsum += abs(Q_av[f]);
            Qsum += Q_av[f];
            knot* knots[4];
            for (int i = 0; i < 4; i++)
               knots[i] = mesh->knots[mesh->faces[f]->knots_num[i]];

            messum += mesh->faces[f]->GetArea(knots);
         }
         messum *= w.intake / 86400.;
         for (int i = 0; i < w.bound_num.size(); i++)
            Q_av[mesh->bounds2[w.bound_num[i]]->face_num] += 
               (messum - Qsum) * abs(Q_av[mesh->bounds2[w.bound_num[i]]->face_num]) / absQsum;
      }
   }

   void Streams::AssembleMatrix()
   {
      // A[i][j] = b_e * Sg[e][i] * Sg[e][j], (i,j) <- hexa.faces_num
      for (int i = 0; i < B->dim; i++)
         B->di[i] = 0.;
      for (int i = 0; i < B->ig[B->dim]; i++)
         B->l[i] = B->u[i] = 0.;

      for (int e = 0; e < mesh->hexas.size(); e++)
      {
         
         auto& el = mesh->hexas[e];
         for (int i = 0; i < 6; i++)
         {
            int gfi = el->faces_num[i];
            for (int j = 0; j < 6; j++)
            {
               int gfj = el->faces_num[j];
               real val = beta[e] * el->faces_sign[i] * el->faces_sign[j];
               AddElement(B, gfi, gfj, val);
            }
         }
      }

      for (int i = 0; i < B->dim; i++)
         B->di[i] += alpha[i];

      // fixate knowns
      for (auto cond : mesh->bounds2)
      {
         if (abs(cond->value) < 1e-12) 
            continue;
         B->di[cond->face_num] = 1.;
         for (int j = B->ig[cond->face_num]; j < B->ig[cond->face_num + 1]; j++)
            B->l[j] = 0.;
         for (int ii = 0; ii < B->dim; ii++)
            for (int j = B->ig[ii]; j < B->ig[ii + 1]; j++)
               if (B->jg[j] == cond->face_num)
                  B->u[j] = 0.;
         MatSymmetrisation(B, d, cond->face_num);

      }
   }

   void Streams::AssembleRightPart()
   {
      for (int i = 0; i < B->dim; i++)
         d[i] = 0.;

      for (int e = 0; e < mesh->hexas.size(); e++)
      {
         auto& el = mesh->hexas[e];
         real sumQ = 0.;
         for (int i = 0; i < 6; i++)
            sumQ += el->faces_sign[i] * Q_av[el->faces_num[i]];
         for (int i = 0; i < 6; i++)
            d[el->faces_num[i]] -= beta[e] * el->faces_sign[i] * sumQ; //d_local[i];
      }

      for (auto cond : mesh->bounds2)
         if (abs(cond->value) > 1e-12) 
            d[cond->face_num] = 0.;
   }

   void Streams::SetAlpha()
   {
      real max_flow = 0.;
      for (int i = 0; i < B->dim; i++)
         if (abs(Q_av[i]) > max_flow)
            max_flow = abs(Q_av[i]);
      for (int i = 0; i < B->dim; i++)
         alpha[i] = 1. / std::max(1e-8 * max_flow, abs(Q_av[i]));
      //for (int e = 0; e < mesh->hexas.size(); e++)
      //{
      //   auto& el = mesh->hexas[e];
      //   real disbalance = 0.;
      //   for (int f = 0; f < 6; f++)
      //      disbalance += el->faces_sign[f] * Q_av[el->faces_num[f]];
      //   for (int i = 0; i < 6; i++)
      //      if (mesh->faces[el->faces_num[i]]->hexa_nums.size() > 1) 
      //         alpha[el->faces_num[i]] += 1. / std::max(abs(disbalance), 1e-10);
      //}
      //for (int f = 0; f < mesh->faces.size(); f++)
      //   alpha[f] /= 2.;
   }

   bool Streams::AdjustBeta(real eps)
   {
      std::cout << '\n' << "Summary streams:\n";
      std::cout << std::scientific;

      real max_flow = 0.;
      for (int i = 0; i < B->dim; i++)
         max_flow = abs(Q_av[i]) < max_flow ? max_flow : abs(Q_av[i]);
      bool isBalanced = true;
      real dissum = 0.;
      int count = 0;
      for (int e = 0; e < mesh->hexas.size(); e++)
      {
         auto& el = mesh->hexas[e];
         real sum = 0.;
         for (int f = 0; f < 6; f++)
            sum += el->faces_sign[f] * (Q_av[el->faces_num[f]] + dQ[el->faces_num[f]]);
         //real maxQe = 0.;
         //for (int i = 0; i < 6; i++)
         //   maxQe = std::max(maxQe, abs(Q_av[el->faces_num[i]] + dQ[el->faces_num[i]]));
         if (abs(sum) / max_flow > eps)
         {
            beta[e] += 1e12*abs(sum)/eps;
            isBalanced = false;
            count++;
         }
         dissum += abs(sum);
      }
      int size = mesh->hexas.size();
      std::cout << "\t----Total Disbalance: " << abs(dissum) / mesh->hexas.size() 
                << ", Relative disbalance: " << abs(dissum) / max_flow / mesh->hexas.size() 
                << "\n\tNum of disbalanced cells " << count << "/" << size << '\n';

      return !isBalanced;
   }

   void Streams::OutputStreams(std::ofstream& outfull, std::ofstream& outone, std::ofstream& outdisbalance, int face_num_to_out)
   {
      outfull << std::scientific;
      outone << std::scientific;
      outdisbalance << std::scientific;

      real max_flow = 0.;
      for (int i = 0; i < B->dim; i++)
         max_flow = abs(Q_av[i]) < max_flow ? max_flow : abs(Q_av[i]);
      for (int i = 0; i < Q.size(); i++)
      {
         outfull << i << " " << Q_av[i] << " " << Q[i] << "\n";
         if (i == face_num_to_out)
            outone << i << " " << Q_av[i] << " " << Q[i] << "\n";
      }
      real av_dis = 0.;
      real av_dis_bal = 0.;

      for (int e = 0; e < mesh->hexas.size(); e++)
      {
         auto& el = mesh->hexas[e];
         real sum = 0.;
         real dissum = 0.;
         for (int f = 0; f < 6; f++){
            sum += el->faces_sign[f] * Q[el->faces_num[f]];
            dissum += el->faces_sign[f] * Q_av[el->faces_num[f]];}
         av_dis_bal += abs(sum) / max_flow;
         av_dis += abs(dissum) / max_flow;
         outdisbalance << e << " " << abs(dissum) << " " << abs(sum) << "\n";
      }
      outdisbalance << av_dis / mesh->hexas.size() << " " << av_dis_bal / mesh->hexas.size();
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

   real Streams::det_J3D()
   {
      return J3D[0][0] * J3D[1][1] * J3D[2][2] + J3D[2][0] * J3D[0][1] * J3D[1][2] + J3D[1][0] * J3D[2][1] * J3D[0][2]
         - (J3D[2][0] * J3D[1][1] * J3D[0][2] + J3D[0][0] * J3D[2][1] * J3D[1][2] + J3D[1][0] * J3D[0][1] * J3D[2][2]);
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

   real Streams::d_phi2D(int index, int var, real u, real v)
   {
      real d_phi = 0.;
      switch (var)
      {
         case 0:    // u
         {
            d_phi = W((index / 2), v);
            if (!(index % 2)) d_phi *= -1;
            break;
         }
         case 1:     // v
         {
            d_phi = W((index % 2), u);
            if (!(index / 2)) d_phi *= -1;
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

   real Streams::prime_by_var2D(int what, int varOnFace, int nFace, real u, real v, real x[4], real y[4])
   {
      real var = 0.;
      for (int i = 0; i < 4; i++)
      {
         switch (varOnFace)
         {
            case 0: var += x[i] * d_phi2D(i, what, u, v); break;
            case 1: var += y[i] * d_phi2D(i, what, u, v); break;
         }
      }
      return var;
   }

   void Streams::calc_grad(int index, int Pindex, real ksi, real etta, real tetha)
   {
      for (int i = 0; i < 3; i++)
         grad[i] += P[Pindex] * d_phi(index, i, ksi, etta, tetha);
   }
}