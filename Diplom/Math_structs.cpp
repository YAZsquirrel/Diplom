#include "Math_structs.h"
#include <set>
#include <array>
#include <cmath>

namespace maths
{
   std::shared_ptr<Matrix> MakeSparseFormat(int localsize, size_t size, std::shared_ptr<Mesh> mesh)
   {
      const int N = localsize;
      // set connection table
      std::vector<std::set<int>> map;
      map.resize(size);
      for (auto& Elem : mesh->elems)
         for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
               if (Elem.knots_num[i] > Elem.knots_num[j])
                  map[Elem.knots_num[i]].insert(Elem.knots_num[j]);

      std::shared_ptr<Matrix> M = std::make_shared<Matrix>(SparseRowColumn);
      M->dim = size;
      M->ig.resize(size + 1, 0);

      for (size_t i = 0; i < size; i++)
         M->ig[i + 1] = M->ig[i] + (int)map[i].size();
      M->jg.resize(M->ig[size], 0);
      for (size_t i = 0; i < map.size(); i++)
      {
         std::vector<int> jind;
         jind.reserve(map[i].size());
         std::copy(map[i].begin(), map[i].end(), std::back_inserter(jind));
         for (int j = 0; j < jind.size(); j++)
            M->jg[M->ig[i] + j] = jind[j];
      }

      M->l.resize(M->ig[size], 0.);
      M->u.resize(M->ig[size], 0.);
      M->di.resize(size, 0.);

      return M;

   }

   void AddElement(std::shared_ptr<Matrix> M, int i, int j, real elem)
   {
      switch (M->format)
      {
      case Dense:
         M->dense[i][j] += elem;
         break;
      case SparseProfile:

         break;
      case SparseRowColumn:
         bool found = false;
         if (i == j)
            M->di[i] += elem;
         else if (i < j)
         {
            int m;
            for (m = M->ig[j]; m < M->ig[j + 1]; m++)
               if (M->jg[m] == i) { found = true; break; }
            if (found)
               M->u[m] += elem; // i-1?
         }
         else
         {
            int n;
            for (n = M->ig[i]; n < M->ig[i + 1]; n++)
               if (M->jg[n] == j) { found = true; break; }
            if (found)
               M->l[n] += elem; // i-1??
         }
         break;

      }


   }

   void MatxVec(std::vector<real>& v, std::shared_ptr<Matrix> M, std::vector<real>& b) // v = M_rcf*b
   {
      switch (M->format)
      {
      case Dense:
         for (size_t i = 0; i < M->dim; i++)
         {
            real sum = 0;
            for (int j = 0; j < M->dim; j++)
               sum += M->dense[i][j] * b[j];
            v[i] = sum;
         }
         break;
      case SparseProfile:

         break;
      case SparseRowColumn:
         for (int i = 0; i < M->dim; i++)
            v[i] = M->di[i] * b[i];

         for (int i = 0; i < M->dim; i++)
            for (int ii = M->ig[i]; ii < M->ig[i + 1]; ii++) // -1?
            {
               v[i] += M->l[ii] * b[M->jg[ii]];
               v[M->jg[ii]] += M->u[ii] * b[i];
            }
         break;
      case SparseRow:
         for (int i = 0; i < M->dim; i++)
         {
            real sum = 0.;
            for (int k = M->ig[i]; k < M->ig[i + 1]; k++)
               sum += M->gg[k] * b[M->jg[k]];

            v[i] = sum;
         }

         break;
      }
   }

   void copy(std::vector<real>& to, std::vector<real>& from)
   {
      for (int i = 0; i < to.size(); i++)
         to[i] = from[i];
   }
   void copy(std::vector<int>& to, std::vector<int>& from)
   {
      for (int i = 0; i < to.size(); i++)
         to[i] = from[i];
   }

   real scalar(std::vector<real>& v, std::vector<real>& u)
   {
      real sum = 0.;

      for (int i = 0; i < v.size(); i++)
         sum += v[i] * u[i];
      return sum;
   }

   real SLAEResidualOutput(std::vector<real>& q, std::shared_ptr<Matrix> M, std::vector<real>& b)
   {
      std::vector<real> y, & x = q;
      y.resize(M->dim, 0);
      MatxVec(y, M, x);
      for (size_t i = 0; i < M->dim; i++)
         y[i] -= b[i];
      return sqrt(scalar(y, y)) / sqrt(scalar(b, b));

      //std::cout << "res: " << res << '\n';
   }

   void SolveSLAE_LOS(std::shared_ptr<Matrix> M, std::vector<real>& q, std::vector<real>& b)
   {
      std::vector<real> z, r, p, Ar, & x = q;
      z.resize(M->dim);
      r.resize(M->dim);
      p.resize(M->dim);
      Ar.resize(M->dim);


      real res, alpha, beta, skp, eps = 1e-14;
      int i, k;
      //x = q;

      MatxVec(Ar, M, x);
      for (i = 0; i < M->dim; i++)
         z[i] = r[i] = b[i] - Ar[i];
      MatxVec(p, M, z);
      real b_norm = sqrt(scalar(b, b));
      res = sqrt(scalar(r, r)) / b_norm;

      for (k = 1; k < 100000 && res > eps; k++)
      {
         skp = scalar(p, p);
         alpha = scalar(p, r) / skp;
         for (i = 0; i < M->dim; i++)
         {
            x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
         }
         MatxVec(Ar, M, r);
         beta = -scalar(p, Ar) / skp;
         for (i = 0; i < M->dim; i++)
         {
            z[i] = r[i] + beta * z[i];
            p[i] = Ar[i] + beta * p[i];
         }
         res = sqrt(scalar(r, r)) / b_norm;
      }
      std::cout << "iter: " << k << " Residual: " << res << std::endl;
   }

   void SolveSLAE_predet_LOS(std::shared_ptr<Matrix> A, std::vector<real>& q, std::vector<real>& b)
   {
      std::vector<real> z, r, p, AQr, t, SAQr, & x = q;
      z.resize(A->dim);
      r.resize(A->dim);
      p.resize(A->dim);
      AQr.resize(A->dim);
      SAQr.resize(A->dim);
      t.resize(A->dim);
      real res, alpha, beta, pp, eps = 1e-14;
      int i, k;
      //x = q;

      std::shared_ptr<Matrix> SQ = MakeHolessky(A);

      MatxVec(t, A, x);  // A*x0                          
      real b_norm;

      for (int i = 0; i < A->dim; i++)    //b - A*x0     
         t[i] = b[i] - t[i];

      SolveForL(r, t, SQ); // r0 = 1/S * (b - Ax0) <-> S*r0 = b - A*x0 
      SolveForU(z, r, SQ);  // z0 = 1/Q * r0 <-> Q*z0 = r0

      MatxVec(t, A, z);    // A*z0                         
      SolveForL(p, t, SQ); // p0 = 1/S * A*z0 <-> S*p0 = A*z0

      b_norm = sqrt(scalar(b, b));
      real res0 = res = scalar(r, r); //sqrt(scalar(r, r)) / b_norm;

      for (k = 1; k < 100000 && res / res0 > eps * eps; k++)
      {
         for (int i = 0; i < A->dim; i++)
            SAQr[i] = AQr[i] = t[i] = 0.;

         pp = scalar(p, p);
         alpha = scalar(p, r) / pp;          // a = (pk0, rk0) / (pk0, pk0)

         for (int i = 0; i < A->dim; i++)     // xk1 = xk0 + a*zk0
         {                                    // rk1 = rk0 - a*pk0 
            x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
         }

         SolveForU(t, r, SQ);                 // t = 1/Q * rk <-> Q*t = rk
         MatxVec(AQr, A, t);                  // A * t1 = AQr -> AQr = A * 1/Q * rk
         SolveForL(SAQr, AQr, SQ);             // SAQr = 1/S AQr <-> S*SAQr = AQr -> SAQr = 1/S * A * 1/Q * rk

         beta = -scalar(SAQr, p) / pp;          // beta = (pk0, t) / (pk0, pk0)          

         //SolveForU(t2, r, SQ);                 // Q*zk1 = rk1
         for (int i = 0; i < A->dim; i++)
         {
            z[i] = t[i] + beta * z[i];        // zk1 += beta
            p[i] = SAQr[i] + beta * p[i];
         }

         res = scalar(r, r);//sqrt(scalar(r, r)) / b_norm;
         if (res != res)
            std::cerr << "Error: NaN detected!" << std::endl;

      }
      std::cout << "iter: " << k << " Residual: " << res << std::endl;
   }

   void SolveSLAE_Relax(std::shared_ptr<Matrix> M, std::vector<real>& q, std::vector<real>& b, real w)
   {
      w = w > 2. ? 2. : w;
      w = w <= 0. ? 0.5 : w;
      std::vector<real> Ar, q1, q2;
      Ar.resize(M->dim);
      q1.resize(M->dim);
      q2.resize(M->dim);

      MatxVec(Ar, M, q);
      copy(q1, q);

      real res = 0., prevres = 1e7;
      for (size_t i = 0; i < M->dim; i++)
         res += b[i] - Ar[i];
      int k = 0;

      for (; abs(res - prevres) < 1e-14 && k < 100000; k++)
      {
         copy(q2, q1);
         for (size_t i = 0; i < M->dim; i++)
         {
            Ar[i] = M->di[i] * q1[i];
            for (int j = M->ig[i]; j < M->ig[i + 1]; j++) // -1?
            {
               Ar[i] += M->l[j] * q1[M->jg[j]];
               Ar[M->jg[j]] += M->u[j] * q1[i];
            }

            q1[i] += w / M->di[i] * (b[i] - Ar[i]);
         }
         prevres = res;
         res = 0.;
         MatxVec(Ar, M, q1);

         for (size_t i = 0; i < M->dim; i++)
            res += b[i] - Ar[i];
      }
      copy(q, q2);

      SLAEResidualOutput(q, M, b);
   }

   void WriteMatrix(std::shared_ptr<Matrix> M, std::vector<real>& b)
   {
      double** mat = new double* [M->dim] {};
      for (int i = 0; i < M->dim; i++)
      {
         mat[i] = new double[M->dim] {};
      }

      switch (M->format)
      {
      case Dense:
         for (int i = 0; i < M->dim; i++)
         {
            for (int j = 0; j < M->dim; j++)
            {
               mat[i][j] = M->dense[i][j];
            }
         }
         break;
      case SparseRowColumn:

         for (int i = 0; i < M->dim; i++)
         {
            mat[i][i] = M->di[i];
            for (int j = M->ig[i]; j < M->ig[i + 1]; j++)
            {
               mat[i][M->jg[j]] = M->l[j];
               mat[M->jg[j]][i] = M->u[j];
            }
         }

         break;
      }

      std::ofstream out("matrix.txt");

      for (int i = 0; i < M->dim; i++)
      {
         for (int j = 0; j < M->dim; j++)
         {
            out.setf(std::ios::left);
            out.width(15);
            out << mat[i][j];
         }
         out << "\n";
      }
      out << "\n\n";
      for (int i = 0; i < M->dim; i++)
      {
         out.setf(std::ios::left);
         out.width(15);
         out << b[i];
      }

      out.close();
   }

   void MatSymmetrisation(std::shared_ptr<Matrix> M, std::vector<real>& b, int i)
   {
      switch (M->format)
      {
      case Dense:
         for (int j = 0; j < M->dim; j++)
         {
            b[j] -= b[i] * M->dense[i][j];
            M->dense[i][j] = 0;
         }
         M->dense[i][i] = 1.;
         break;
      case SparseProfile:

         break;
      case SparseRowColumn:
         for (int j = M->ig[i]; j < M->ig[i + 1]; j++)
         {
            b[M->jg[j]] -= b[i] * M->u[j];
            M->u[j] = 0;
         }

         for (int n = 0; n < M->dim; n++)
         {
            for (int j = M->ig[n]; j < M->ig[n + 1]; j++)
               if (M->jg[j] == i)
               {
                  b[n] -= b[i] * M->l[j];
                  M->l[j] = 0.0;
               }
         }
         break;

      }

   }

   std::shared_ptr<Matrix> MakeHolessky(std::shared_ptr<Matrix> A)
   {

      std::shared_ptr<Matrix> M = std::make_shared<Matrix>(SparseRowColumn);
      M->dim = A->dim;
      M->ig.resize(A->ig.size());
      M->jg.resize(A->jg.size());
      M->l.resize(A->l.size());
      M->u.resize(A->u.size());
      M->di.resize(A->di.size());

      copy(M->ig, A->ig);
      copy(M->jg, A->jg);

      for (int i = 0; i < M->dim; i++)
      {
         // Sij, Qji
         for (int ii = M->ig[i]; ii < M->ig[i + 1]; ii++) // ii -> index on row i
         {
            int j = M->jg[ii];

            real sumS = 0.;
            real sumQ = 0.;
            std::vector<int*> ks;
            int is = M->ig[i], ie = M->ig[i + 1];
            int js = M->ig[j], je = M->ig[j + 1];

            for (int ki = is; ki < ie; ki++)
               for (int kj = js; kj < je; kj++)
                  if (M->jg[ki] == M->jg[kj])
                     ks.push_back(new int[2] { ki, kj }); // find intersection of nonzero elem's ii

            for (auto k : ks)
            {
               int ki = k[0], kj = k[1];
               sumS += M->l[ki] * M->u[kj];
               sumQ += M->l[kj] * M->u[ki];
               delete[] k;
            }

            M->l[ii] = (A->l[ii] - sumS) / M->di[j];
            M->u[ii] = (A->u[ii] - sumQ) / M->di[j];

         }
         //Sii
         {
            real sumd = 0.;
            for (int ii = M->ig[i]; ii < M->ig[i + 1]; ii++)
               sumd += M->l[ii] * M->u[ii];

            M->di[i] = sqrt(A->di[i] - sumd);
         }
      }

      return M;
   }

   void SolveForL(std::vector<real>& q, std::vector<real>& b, std::shared_ptr<Matrix> M) // y = 1/L * b // q = 1/L * b // forward
   {
      switch (M->format)
      {
      case Dense:
         for (int i = 0; i < M->dim; i++)
         {
            real sum = 0;
            for (int j = 0; j < i; j++)
               sum += q[j] * M->dense[i][j];
            q[i] = b[i] - sum;
         }
         break;
      case SparseProfile:
         for (int i = 0; i < M->dim; i++)
         {
            int k = i - (M->ig[i + 1] - M->ig[i]);
            real sum = 0;
            for (int j = M->ig[i]; j < M->ig[i + 1]; j++, k++)
               sum += q[k] * M->l[j];
            //
            q[i] = (b[i] - sum) / M->di[i];
         }
         break;
      case SparseRowColumn:
         for (size_t k = 0; k < M->dim; k++)
         {
            real sum = 0.;
            for (size_t i = M->ig[k]; i < M->ig[k + 1]; i++)
               sum += M->l[i] * q[M->jg[i]];

            q[k] = (b[k] - sum) / M->di[k];
         }

         break;

      }
   }
   void SolveForU(std::vector<real>& q, std::vector<real>& b, std::shared_ptr<Matrix> M) // x = 1/U * y // q = 1/U * b // back
   {
      switch (M->format)
      {
      case Dense:
         for (int i = M->dim - 1; i > 0; i--)
         {
            real sum = b[i];
            for (int j = i + 1; j < M->dim; j++)
               sum -= q[j] * M->dense[i][j];
            q[i] = sum / M->dense[i][i];
         }
         break;

      case SparseProfile:
         for (int i = M->dim - 1; i > 0; i--)
         {
            //real sum = b[i];
            int j = i - (M->ig[i + 1] - M->ig[i]);
            q[i] += b[i];
            for (int k = M->ig[i]; k < M->ig[i + 1]; k++, j++)
               q[j] -= M->u[k] * q[i];
         }

         break;

      case SparseRowColumn:
         copy(q, b);
         for (int k = M->dim - 1; k >= 0; k--)
         {
            int ii = k;
            real qii = q[k] / M->di[k];
            q[k] = qii;

            for (int i = M->ig[k]; i < M->ig[k + 1]; i++)
               q[M->jg[i]] -= M->u[i] * qii;

         }
         break;

      }
   }
}