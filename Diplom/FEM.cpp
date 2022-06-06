#include "FEM.h"
#include <fstream>
//using namespace mesh_comps;

namespace FEMns
{
FEM::FEM(Mesh* _mesh)
{
#pragma region input
#if DEBUG
   mesh = new Mesh();
#else
   mesh = _mesh;
#endif // DEBUG


   std::ifstream fknots("Knots.txt");
   std::ifstream fhexas("Hexahedrons.txt");
   std::ifstream fbounds1("FirstBounds.txt");
   std::ifstream fbounds2("SecondBounds.txt");
   std::ifstream fparams("EnvParams.txt");
   fparams >> un;
   fparams.close();

   mesh->hexas.clear();
   mesh->knots.clear();

   fknots >> num_of_knots;
   mesh->knots.reserve(num_of_knots);
   
   for (int i = 0; i < num_of_knots; i++)
   {
      knot* k = new knot();
      fknots >> k->x >> k->y >> k->z;
      mesh->knots.push_back(k);
   }
   fknots.close();

   hexahedron* hexa;
   fhexas >> num_of_FE;
   mesh->hexas.reserve(num_of_FE);
   for (int i = 0; i < num_of_FE; i++)
   {
      mesh->hexas.push_back(hexa = new hexahedron());
      hexa->faces_num.reserve(6);
      for (int k = 0; k < 8; k++)
         fhexas >> hexa->knots_num[k];
   }
   fhexas.close();
   mesh->FindNeighborsAndFaces();
   mesh->FindFaceNormals();
   mesh->SetSignsForHexaNormals();

   int numOfBounds;
   fbounds1 >> numOfBounds;
   for (int i = 0; i < numOfBounds; i++)
   {
      bound* cond = new bound;
      for (int j = 0; j < 4; j++)
         fbounds1 >> cond->knots_num[j];

      fbounds1 >> cond->value;
      for (int f = 0; f < mesh->faces.size(); f++)
      {
         bool found[4]{};
         for (int fi = 0; fi < 4; fi++)
            for (int fj = 0; fj < 4 && !found[fi]; fj++)
               found[fi] = mesh->faces[f]->knots_num[fi] == cond->knots_num[fj];
         if ((found[0] && found[1] && found[2] && found[3]))
         {
            cond->face_num = f;
            break;
         }
      }
         
      mesh->bounds1.push_back(cond);
   }
   fbounds1.close();

   fbounds2 >> numOfBounds;
   for (int i = 0; i < numOfBounds; i++)
   {
      bound* cond = new bound;

      for (int j = 0; j < 4; j++)
         fbounds2 >> cond->knots_num[j];
      fbounds2 >> cond->value;
      for (int f = 0; f < mesh->faces.size(); f++)
      {
         bool found[4]{};
         for (int fi = 0; fi < 4; fi++)
            for (int fj = 0; fj < 4 && !found[fi]; fj++)
               found[fi] = mesh->faces[f]->knots_num[fi] == cond->knots_num[fj];
         if ((found[0] && found[1] && found[2] && found[3]))
         {
            cond->face_num = f;
            break;
         }
      }
      mesh->bounds2.push_back(cond);
   }
   fbounds2.close();

#pragma endregion

   A = MakeSparseFormat(8, num_of_knots, num_of_FE, true, mesh);
   q.resize(num_of_knots, 0.);
   b.resize(num_of_knots, 0.);

   Mij = [this](real ksi, real etta, real theta, int i, int j, int knot_num[8])
   {
      for (int ip = 0; ip < 3; ip++)                                  
         for (int jp = 0; jp < 3; jp++)                               
            J[ip][jp] = prime_by_var(ip, jp, knot_num, ksi, etta, theta);

      return phi(i, ksi, etta, theta) * phi(j, ksi, etta, theta) * abs(det_J());
   };

   Gij = [this](real ksi, real etta, real theta, int i, int j, int knot_num[8])
   {
      for (int ip = 0; ip < 3; ip++)
         Jgrad_i[ip] = Jgrad_j[ip] = 0.;
      
      for (int ip = 0; ip < 3; ip++)                                       // | dx/d(ksi)   | dy/d(ksi)   | dz/d(ksi)   |
         for (int jp = 0; jp < 3; jp++)                                    // | dx/d(etta)  | dy/d(etta)  | dz/d(etta)  |
            J[ip][jp] = prime_by_var(ip, jp, knot_num, ksi, etta, theta);  // | dx/d(theta) | dy/d(theta) | dz/d(theta) |
      
      // J^-1
      for (int ip = 0; ip < 3; ip++)
         for (int jp = 0; jp < 3; jp++)
         {
            real min[4]{};
            int k = 0;
            for (int im = 0; im < 3; im++)
               for (int jm = 0; jm < 3; jm++)
               {
                  if (im != ip && jm != jp)
                     min[k++] = J[im][jm];
               }
            reversed_J[jp][ip] = ((ip + jp + 2) % 2 ? -1 : 1) * (min[0] * min[3] - min[1] * min[2]);
         }
      
      // grad(phi(ksi, etta, theta))
      calc_grad(1, i, ksi, etta, theta);
      calc_grad(2, j, ksi, etta, theta); 

      // J^-1 * grad(phi)
      for (int ip = 0; ip < 3; ip++)
         for (int jp = 0; jp < 3; jp++)
         {
            Jgrad_i[ip] += reversed_J[ip][jp] * gradi[jp];
            Jgrad_j[ip] += reversed_J[ip][jp] * gradj[jp];
         }
 
      // Jgrad_i^T * Jgrad_j
      real res = 0;
      for (int ip = 0; ip < 3; ip++)
         res += Jgrad_i[ip] * Jgrad_j[ip];
      return res / abs(det_J());
   };

}

void FEM::SolveElliptic()
{
   CreateSLAE();
   SolveSLAE(A, q, b);
   std::ofstream out("Result.txt");

#ifdef DEBUG
   check_test();
#endif // DEBUG
   Output(out);
   out.close();
}

void FEM::GetSolutionOnPlane(real z)
{
   std::ofstream zout("ResultZ.txt");

   for (int i = 0; i < num_of_knots; i++)
   {
      //if (mesh->knots[i]->z <= z + 1e-13 && mesh->knots[i + 1]->z >= z + 1e-13)
      if (abs(mesh->knots[i]->z - z) <= 1e-12)
      {
         zout << mesh->knots[i]->x << " " 
               << mesh->knots[i]->y << " "
               //<< mesh->knots[i]->z << " "
               //<< q[i] * (1. - (z - mesh->knots[i]->z) / (mesh->knots[i + 1]->z - mesh->knots[i]->z)) +
               //   q[i + 1]     *       (z - mesh->knots[i]->z) / (mesh->knots[i + 1]->z - mesh->knots[i]->z)
               << q[i] << '\n';}
   }
   zout.close();
}

void FEM::Output(std::ofstream& out)
{
   //out.scientific;
   //out.precision(15);
   //std::cout.scientific;
   //std::cout.precision(15);
   out.setf(std::ios::right);
   out.width(15);
   out << "\n| x" << std::fixed;
   out.width(15);
   out << "| y";
   out.width(15);
   out << "| z";
   out.width(15);
   out << "| q";
   out.width(15);
   out << "|\n";
   //std::cout << title;

   for (int i = 0; i < num_of_knots; i++)
   {
      out << "|" << mesh->knots[i]->x;
      out.width(15);
      out << "|" << mesh->knots[i]->y;
      out.width(15);
      out << "|" << mesh->knots[i]->z;
      out.width(15);
      out << "|" << q[i];
      out.width(15);
      out << "|\n";
      //out  << "| " << "\t| " << mesh->knots[i]->x << "\t| " << mesh->knots[i]->y << mesh->knots[i]->z << "\t| "
      //   << std::scientific << ug(&mesh->knots[i]) << "\t| " << q[i] << "\t| "
      //   << abs(q[i] - ug(&mesh->knots[i])) << "\t|\n";
      // std::cout << std::scientific << "| " << ug(&knots[i]) << "\t| " << q[i] << "\t| "
      //    << abs(q[i] - ug(&knots[i])) << "\t|\n";
   }

}

void FEM::AddFirstBounds()
{
   for (auto cond : mesh->bounds1)
   {
      for (int i = 0; i < 4; i++)
      {
         A->di[cond->knots_num[i]] = 1e10;
         //for (int j = A->ig[cond->knots_num[i]]; j < A->ig[cond->knots_num[i] + 1]; j++)
         //   A->l[j] = 0.;
         //for (int j = 0; j < A->ig[num_of_knots]; j++)
         //   if (A->jg[j] == cond->knots_num[i])
         //      A->u[j] = 0.;
         #ifdef DEBUG
         b[cond->knots_num[i]] = 1e10 * ug(mesh->knots[cond->knots_num[i]]);//cond->value;// ;
         #else
         b[cond->knots_num[i]] = 1e10 * cond->value;
         #endif
      }
   }
}

void FEM::AddSecondBounds()
{  
   for (auto bound : mesh->bounds2)
   {
      real a0, a1, a2;
      real x[4], y[4];
      real xn = mesh->faces[bound->face_num]->normal.x,
           yn = mesh->faces[bound->face_num]->normal.y,
           zn = mesh->faces[bound->face_num]->normal.z;
      std::vector<real> planeNormal;
      planeNormal.resize(3, 0.);
      if (abs(xn) > abs(zn))
      {
         if (abs(xn) > abs(yn)) // max = x
         {
            x[0] = mesh->knots[bound->knots_num[0]]->y;
            x[1] = mesh->knots[bound->knots_num[1]]->y;
            x[2] = mesh->knots[bound->knots_num[2]]->y;
            x[3] = mesh->knots[bound->knots_num[3]]->y;
            y[0] = mesh->knots[bound->knots_num[0]]->z;
            y[1] = mesh->knots[bound->knots_num[1]]->z;
            y[2] = mesh->knots[bound->knots_num[2]]->z;
            y[3] = mesh->knots[bound->knots_num[3]]->z;
            planeNormal[0] = 1.;
         }
         else // max = y
         {
            x[0] = mesh->knots[bound->knots_num[0]]->x;
            x[1] = mesh->knots[bound->knots_num[1]]->x;
            x[2] = mesh->knots[bound->knots_num[2]]->x;
            x[3] = mesh->knots[bound->knots_num[3]]->x;
            y[0] = mesh->knots[bound->knots_num[0]]->z;
            y[1] = mesh->knots[bound->knots_num[1]]->z;
            y[2] = mesh->knots[bound->knots_num[2]]->z;
            y[3] = mesh->knots[bound->knots_num[3]]->z;
            planeNormal[1] = 1.;
         }
      }
      else
      {
         if (abs(zn) > abs(yn)) // max = z
         {
            x[0] = mesh->knots[bound->knots_num[0]]->x;
            x[1] = mesh->knots[bound->knots_num[1]]->x;
            x[2] = mesh->knots[bound->knots_num[2]]->x;
            x[3] = mesh->knots[bound->knots_num[3]]->x;
            y[0] = mesh->knots[bound->knots_num[0]]->y;
            y[1] = mesh->knots[bound->knots_num[1]]->y;
            y[2] = mesh->knots[bound->knots_num[2]]->y;
            y[3] = mesh->knots[bound->knots_num[3]]->y;
            planeNormal[2] = 1.;
         }
         else // max = y
         {
            x[0] = mesh->knots[bound->knots_num[0]]->x;
            x[1] = mesh->knots[bound->knots_num[1]]->x;
            x[2] = mesh->knots[bound->knots_num[2]]->x;
            x[3] = mesh->knots[bound->knots_num[3]]->x;
            y[0] = mesh->knots[bound->knots_num[0]]->z;
            y[1] = mesh->knots[bound->knots_num[1]]->z;
            y[2] = mesh->knots[bound->knots_num[2]]->z;
            y[3] = mesh->knots[bound->knots_num[3]]->z;
            planeNormal[1] = 1.;
         }
      }
      std::vector<real> bnormal = {xn, yn, zn};
      
      real Sfactor = abs(scalar(bnormal, planeNormal));
      a0 = (x[1] - x[0])*(y[2] - y[0]) - (y[1] - y[0])*(x[2] - x[0]);
      a1 = (x[1] - x[0])*(y[3] - y[2]) - (y[1] - y[0])*(x[3] - x[2]);
      a2 = (x[3] - x[1])*(y[2] - y[0]) - (y[3] - y[1])*(x[2] - x[0]);
      real s = 1;
      if (a0 < 0) s = -1;
      else if (abs(a0) < 1e-12) s = 0;
      for (int i = 0; i < 6; i++)
         if (mesh->hexas[mesh->faces[bound->face_num]->hexa_nums[0]]->faces_num[i] == bound->face_num) {
            s *= mesh->hexas[mesh->faces[bound->face_num]->hexa_nums[0]]->faces_sign[i]; break;}

      localM2d[0][0] = a0 / 9.  + a1 / 36. + a2 / 36.;  localM2d[0][1] = a0 / 18. + a1 / 36. + a2 / 72.;  localM2d[0][2] = a0 / 18. + a1 / 72. + a2 / 36.;  localM2d[0][3] = a0 / 36. + a1 / 72. + a2 / 72.;
      localM2d[1][0] = a0 / 18. + a1 / 36. + a2 / 72.;  localM2d[1][1] = a0 / 9.  + a1 / 12. + a2 / 36.;  localM2d[1][2] = a0 / 36. + a1 / 72. + a2 / 72.;  localM2d[1][3] = a0 / 18. + a1 / 24. + a2 / 36.;
      localM2d[2][0] = a0 / 18. + a1 / 72. + a2 / 36.;  localM2d[2][1] = a0 / 36. + a1 / 72. + a2 / 72.;  localM2d[2][2] = a0 / 9.  + a1 / 36. + a2 / 12.;  localM2d[2][3] = a0 / 18. + a1 / 36. + a2 / 24.;
      localM2d[3][0] = a0 / 36. + a1 / 72. + a2 / 72.;  localM2d[3][1] = a0 / 18. + a1 / 24. + a2 / 36.;  localM2d[3][2] = a0 / 18. + a1 / 36. + a2 / 24.;  localM2d[3][3] = a0 / 9.  + a1 / 12. + a2 / 12.;
      for (int i = 0; i < 4; i++)
         b[bound->knots_num[i]] += s * bound->value * (localM2d[i][0] + localM2d[i][1] + localM2d[i][2] + localM2d[i][3]) / Sfactor;
   }
}

void FEM::CreateSLAE()
{
   hexahedron* hexa;
   for (int i = 0; i < num_of_FE; i++)
   {
      hexa = mesh->hexas[i];
      CreateG(hexa);
      CreateM(hexa);
      AddToA(hexa);
      Createb(hexa);
   }
#ifdef DEBUG
   WriteMatrix(A);
#endif // DEBUG
   AddSecondBounds();
   AddFirstBounds();
}

void FEM::AddToA(hexahedron* hexa)
{
   for (int i = 0; i < 8; i++ )
      for (int j = 0; j < 8; j++ )
         AddElement(A, hexa->knots_num[i], hexa->knots_num[j], localG[i][j]);
         //localA[i][j] = localG[i][j] + localM[i][j];
   //AddLocal(A, hexa->knots_num, localA, 1);
}

void FEM::CreateM(hexahedron* hexa)
{
    for (int i = 0; i < 8; i++)
      for (int j = 0; j < 8; j++)
         localM[i][j] = Integrate(Mij, i, j, hexa->knots_num);
}

void FEM::CreateG(hexahedron* hexa)
{
#ifdef DEBUG
   for (int i = 0; i < 8; i++)
   {
      for (int j = 0; j < 8; j++)
      {
         localG[i][j] = Integrate(Gij, i, j, hexa->knots_num);
         std::cout << localG[i][j] << " ";
      }
      std::cout << '\n';
   }
   std::cout << "\n";
#else
   for (int i = 0; i < 8; i++)
      for (int j = 0; j < 8; j++)
         localG[i][j] = hexa->lam * Integrate(Gij, i, j, hexa->knots_num);
#endif // DEBUG



}

void FEM::Createb(hexahedron* hexa) 
{
   real localb[8]{};
   real f_[8]{};
#ifdef DEBUG
   for (int i = 0; i < 8; i++)
      f_[i] = f(mesh->knots[hexa->knots_num[i]]);
#else
   for (int i = 0; i < 8; i++)
      f_[i] = 0.0;//(mesh->knots[hexa->knots_num[i]]);
#endif // DEBUG
   

   for (int i = 0; i < 8; i++)
      for (int j = 0; j < 8; j++)
         localb[i] += localM[i][j] * f_[j];
   
   for (int i = 0; i < 8; i++)
      b[hexa->knots_num[i]] += localb[i];
}

int FEM::mu(int index)
{
   return index % 2;
}

int FEM::v(int index)
{
   return (index / 2) % 2;
}

int FEM::nu(int index)
{
   return index/ 4;
}

real FEM::W(int index, real alpha)
{
    if (!index) return 1. - alpha;
    return alpha;
}

real FEM::d_phi(int index, int var, real ksi, real etta, real tetha)
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

real FEM::prime_by_var(int varOnCube, int varOnFE, int knot_num[8], real ksi, real etta, real tetha)
{
   real var = 0.;
   for (int i = 0; i < 8; i++)
   {
      switch (varOnFE)
      {
         case 0: var += mesh->knots[knot_num[i]]->x * d_phi(i, varOnCube, ksi, etta, tetha); break;
         case 1: var += mesh->knots[knot_num[i]]->y * d_phi(i, varOnCube, ksi, etta, tetha); break;
         case 2: var += mesh->knots[knot_num[i]]->z * d_phi(i, varOnCube, ksi, etta, tetha); break;
      }
   }
   return var;
}

real FEM::phi(int index, real ksi, real etta, real tetha)
{
   return  W(mu(index), ksi) * W(v(index), etta) * W(nu(index), tetha);
}

real FEM::det_J()
{
   return J[0][0] * J[1][1] * J[2][2] + J[2][0] * J[0][1] * J[1][2] + J[1][0] * J[2][1] * J[0][2] 
       - (J[2][0] * J[1][1] * J[0][2] + J[0][0] * J[2][1] * J[1][2] + J[1][0] * J[0][1] * J[2][2]);
}

real FEM::Integrate(const std::function<real(real, real, real, int, int, int[8])> f, int i, int j, int knot_nums[8])
{
   const int nKnot = 3;//5; // Knots num

   const real xj[nKnot] 
   = { .7745966692414833, 0., -.7745966692414833 }; // sqrt(0.6)
   //= { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,	// Scales
   //           0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };

   const real qj[nKnot] 
   = { .55555555555555555, .8888888888888888, .55555555555555555 };
   //= { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,	// Weights
   //               (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };

   real result = 0.;
   for (int ix = 0; ix < nKnot; ix++)
      for (int iy = 0; iy < nKnot; iy++)
         for (int iz = 0; iz < nKnot; iz++)
            result += qj[ix] * qj[iy] * qj[iz] * (f(.5 + xj[ix] * .5, .5 + xj[iy] * .5, .5 + xj[iz] * .5, i, j, knot_nums));
   return result / 8.; 
}

void FEM::calc_grad(int ij, int index, real ksi, real etta, real tetha)
{
   switch (ij)
   {
      case 1: 
         for (int i = 0; i < 3; i++)
            gradi[i] = d_phi(index, i, ksi, etta, tetha);
         break;
      case 2:
         for (int i = 0; i < 3; i++)
            gradj[i] = d_phi(index, i, ksi, etta, tetha);
         break;
   }
}

void FEM::check_test()
{
   real qqt2 = 0;
   real ug2 = 0;
   for (int i = 0; i < num_of_knots; i++)
   {
      qqt2 += pow(q[i] - ug(mesh->knots[i]), 2);
      ug2 += pow(ug(mesh->knots[i]), 2);
   }
   std::cout << "\n |q - u| = " << sqrt(qqt2) / sqrt(ug2) << '\n';
}
}