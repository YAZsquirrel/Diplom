#include "Mat_structs.h"

namespace mats
{
   Matrix* MakeSparseFormat(int localsize, int elemsize, int FEsize, bool isknots, mesh_comps::Mesh* mesh)
   {
      const int N = localsize;
      int* list1, * list2;
      int* listbeg = new int[elemsize];

      for (int i = 0; i < elemsize; i++)
         listbeg[i] = -1;

      list1 = new int[elemsize * elemsize]{};
      list2 = new int[elemsize * elemsize]{};
      int listsize = -1, iaddr, ind1, ind2, k;

      for (int iel = 0; iel < FEsize; iel++) // 
      {
         for (int i = 0; i < N; i++) // 
         {
            k = isknots ? mesh->hexas[iel]->knots_num[i] : mesh->hexas[iel]->faces_num[i]; //
            for (int j = i + 1; j < N; j++) // need to set N = ?
            {
               ind1 = k;
               ind2 = isknots ? mesh->hexas[iel]->knots_num[j] : mesh->hexas[iel]->faces_num[j];  //
               if (ind2 < ind1) //
               {
                  ind1 = ind2;
                  ind2 = k;
               }
               iaddr = listbeg[ind2];
               if (iaddr == -1) // 
               {
                  listsize++;
                  listbeg[ind2] = listsize;
                  list1[listsize] = ind1;
                  list2[listsize] = -1;
               }
               else // 
               {
                  while (list1[iaddr] < ind1 && list2[iaddr] >= 0)
                     iaddr = list2[iaddr];
                  if (list1[iaddr] > ind1)  // 
                  {                         // 
                     listsize++;
                     list1[listsize] = list1[iaddr];
                     list2[listsize] = list2[iaddr];
                     list1[iaddr] = ind1;
                     list2[iaddr] = listsize;
                  }
                  else if (list1[iaddr] < ind1) // 
                  {
                     listsize++;
                     list2[iaddr] = listsize;
                     list1[listsize] = ind1;
                     list2[listsize] = -1;
                  }
               }
            }
         }
      }

      Matrix* M = new Matrix;
      M->dim = elemsize;
      M->ig = new int[elemsize + 1]{};
      M->jg = new int[listsize + 1]{};  // +1???

      for (int i = 0; i < elemsize; i++)
      {
         M->ig[i + 1] = M->ig[i];

         for (iaddr = listbeg[i]; iaddr != -1; )
         {
            M->jg[M->ig[i + 1]] = list1[iaddr];
            M->ig[i + 1]++;
            iaddr = list2[iaddr];
         }
      }

      //for (int i = 0; i < num_of_knots + 1; i++)
      //   std::cout << A->ig[i] << " ";
      //std::cout << '\n';
      //for (int i = 0; i < listsize + 1; i++)
      //   std::cout << A->jg[i] << " ";
      //std::cout << '\n';

      delete[] listbeg;
      delete[] list1;
      delete[] list2;

      M->l = new real[listsize + 1]{};
      M->u = new real[listsize + 1]{};
      M->di = new real[elemsize]{};
      return M;

   }

   void AddElement(Matrix* M, int i, int j, real elem)
   {
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
   }

   void MatxVec(std::vector<real> &v, Matrix* M, std::vector<real> &b) // v = M*b
   {
      for (int i = 0; i < M->dim; i++)
         v[i] = M->di[i] * b[i];

      for (int i = 0; i < M->dim; i++)
         for (int j = M->ig[i]; j < M->ig[i + 1]; j++) // -1?
         {
            v[i] += M->l[j] * b[M->jg[j]];
            v[M->jg[j]] += M->u[j] * b[i];
         }
   }

   void copy(std::vector<real>& v, std::vector<real>& u)
   {
      for (int i = 0; i < v.size(); i++)
         v[i] = u[i];
   }

   real scalar(std::vector<real>& v, std::vector<real>& u)
   {
      real sum = 0.;
      for (int i = 0; i < v.size(); i++)
         sum += v[i] * u[i];
      return sum;
   }

   void SolveSLAE(Matrix* M, std::vector<real> &q, std::vector<real> &b)
   {
      std::vector<real> z, r, p, ff, &x = q;
      z.resize(M->dim);
      r.resize(M->dim);
      p.resize(M->dim);
      ff.resize(M->dim);

      
      real res, alpha, beta, skp, eps = 1e-17;
      int i, k;
      //x = q;
      
      real lastres;
      MatxVec(ff, M, x);
      for (i = 0; i < M->dim; i++)
         z[i] = r[i] = b[i] - ff[i];
      MatxVec(p, M, z);
      res = sqrt(scalar(r, r)) / sqrt(scalar(b, b));

      for (k = 1; k < 100000 && res > eps; k++)
      {
         lastres = res;
         skp = scalar(p, p);
         alpha = scalar(p, r) / skp;
         for (i = 0; i < M->dim; i++)
         {
            x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
         }
         MatxVec(ff, M, r);
         beta = -scalar(p, ff) / skp;
         for (i = 0; i < M->dim; i++)
         {
            z[i] = r[i] + beta * z[i];
            p[i] = ff[i] + beta * p[i];
         }
         res = sqrt(scalar(r, r)) / sqrt(scalar(b, b));
      }
      //std::cout << "iter: " << k << " Residual: " << res << std::endl;
   }

   void WriteMatrix(Matrix* M)
   {
      double** mat = new double* [M->dim] {};
      for (int i = 0; i < M->dim; i++)
      {
         mat[i] = new double[M->dim] {};
      }

      for (int i = 0; i < M->dim; i++)
      {
         mat[i][i] = M->di[i];
         for (int j = M->ig[i]; j < M->ig[i + 1]; j++)
         {
            mat[i][M->jg[j]] = M->l[j];
            mat[M->jg[j]][i] = M->u[j];
         }
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
   }

   void MatSymmetrisation(Matrix* M, std::vector<real>& b, int i)
   {
      for (int j = M->ig[i]; j < M->ig[i + 1]; j++)
      {
         b[M->jg[j]] -= b[i] * M->u[j];
         M->u[j] = 0;
      }

      for (int iIG = 0; iIG < M->dim; iIG++)
      {
         for (int j = M->ig[iIG]; j < M->ig[iIG + 1]; j++)
            if (M->jg[j] == i)
            {
               b[iIG] -= b[i] * M->l[j];
               M->l[j] = 0;
            }
      }
   }
}