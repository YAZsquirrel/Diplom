#include "Filtration.h"
#include <fstream>
#include "super_output.h"

using namespace mesh_comps;
//using namespace filtration;

namespace filtr {

   void Filtration::Start()
   {
      mesh = new Mesh();
      std::ifstream f_filtr("FiltrParams.txt");
      real trash;
      int comps_num;
      f_filtr >> trash;

      int por_num;
      f_filtr >> por_num;
      for (size_t i = 0; i < por_num; i++)
      {
         solid sol;
         f_filtr >> sol.K >> sol.Fi;
         f_filtr >> sol.corner0.x >> sol.corner0.y >> sol.corner0.z;
         f_filtr >> sol.corner1.x >> sol.corner1.y >> sol.corner1.z;
         pors.push_back(sol);
      }

      f_filtr >> comps_num;

      int phases_num;
      f_filtr >> phases_num;
      phases.reserve(phases_num);

      mesh->GenerateMesh();

      comps_in_phases.resize(phases_num);
      for (int i = 0; i < phases.size(); i++)
         comps_in_phases[i].reserve(comps_num);
      
      for (int i = 0; i < phases_num; i++)
      {
         phase ph;
         f_filtr >> ph.h >> ph.n >> ph.k;
         for (int j = 0; j < comps_num; j++)
         {
            real prop;
            f_filtr >> prop;
            comps_in_phases[i].push_back(prop);
         }
         phases.push_back(ph);

      }

      f_filtr.close();

      fem = new FEMns::FEM(mesh);
      SetDiffKoef();
      fem->SolveElliptic();
      Output2D(0, fem->GetKnots().data(), mesh);
      fem->GetSolutionOnPlane(60);
      fem->GetSolutionOnPlane(40);
      fem->GetSolutionOnPlane(20);
      fem->GetSolutionOnPlane(0);
      //return;

      str = new streams::Streams(mesh, fem->GetKnots());
      str->FindStreams();
      std::ofstream outall("StreamsResult.txt");
      std::ofstream outone("Stream1Result.txt");
      std::ofstream outdisb("HexaDisbalances.txt");
      str->OutputStreams(outall, outone, outdisb, 0);
      outall.close();
      outone.close();
      outdisb.close();
      delete mesh;
      delete fem;

   }

   void Filtration::SetDiffKoef()
   {
      std::sort(phases.begin(), phases.end(), [](phase& ph1, phase& ph2){ return ph1.h > ph2.h;});

      for (int i = mesh->hexas.size() - 1; i >= 0; i--)
      {
         hexahedron *hexa = mesh->hexas[i];
         hexa->Sm.resize(phases.size());
         real hh1 = mesh->knots[hexa->knots_num[4]]->z; // Верхняя грань
         real hh2 = mesh->knots[hexa->knots_num[0]]->z; // Нижняя грань
         real hh = hh1 - hh2;

         if (phases.size() == 1)
         {
            hexa->phases_num.push_back(0);
            hexa->Sm[0] = 1.;
         }
         else
            for (int j = 0; j < phases.size(); j++)
            {	
               if (abs(hh1 - phases[j].h) < 1e-12)
               {
                  if (j != phases.size() - 1) // 1-2
                  {
                     real phh = phases[j].h - phases[j + 1].h;

                     if (abs(hh2 - phases[j + 1].h) < 1e-12 || hh2 > phases[j + 1].h) // 2
                        hexa->Sm[j] = 1.;
                     else if (hh2 < phases[j + 1].h) // 1
                     {
                        hexa->Sm[j] = hh / phh;
                        hexa->Sm[j + 1] = 1. - hh / phh;
                        hexa->phases_num.push_back(j + 1);
                     }
                  }
                  else
                     hexa->Sm[j] = 1.;		// 3 (последний слой)
                  hexa->phases_num.push_back(j);
               }
               else if (hh1 < phases[j].h)
               {
                  if (j != phases.size() - 1)	// 2 (посл. слой)
                  {
                     hexa->Sm[j] = 1.;			
                     hexa->phases_num.push_back(j);
                  }
                  else
                  {
                     real phh1 = phases[j].h - phases[j + 1].h;

                     if (abs(hh2 - phases[j + 1].h) < 1e-12 || hh2 > phases[j + 1].h) // 2-6
                     {
                        hexa->Sm[j] = 1.;
                        hexa->phases_num.push_back(j);
                     }
                     else if (hh2 < phases[j + 1].h) // 4-5
                     {
                        if (j == phases.size() - 2) // 4
                        {
                           hexa->Sm[j] = hh / (hh1 - phases[j + 1].h);
                           hexa->Sm[j + 1] = 1. - hexa->Sm[j];
                           hexa->phases_num.push_back(j);
                           hexa->phases_num.push_back(j + 1);
                        }
                        else // 4-5
                        {
                           real phh2 = phases[j + 1].h - phases[j + 2].h;
                           if (hh2 > phases[j + 2].h) // 4
                           {
                              hexa->Sm[j] = hh / (hh1 - phases[j + 1].h);
                              hexa->Sm[j + 1] = 1. - hexa->Sm[j];
                              hexa->phases_num.push_back(j + 1);
                              hexa->phases_num.push_back(j);
                           }
                           else if (abs(hh2 - phases[j + 2].h) < 1e-12)
                           {
                              hexa->Sm[j] = 1. - hh / phh2;
                              hexa->Sm[j + 1] = hh / phh2;
                              hexa->phases_num.push_back(j + 1);
                              hexa->phases_num.push_back(j);
                           }
                           else if (hh2 < phases[j + 2].h) // 5
                           {
                              if (j == phases.size() - 3)
                              {
                                 hexa->Sm[j] = hh / (hh1 - phases[j + 1].h);
                                 hexa->Sm[j + 1] = hh / phh2;
                                 hexa->Sm[j + 2] = 1. - hexa->Sm[j] - hexa->Sm[j + 1];
                                 hexa->phases_num.push_back(j + 1);
                                 hexa->phases_num.push_back(j);
                                 hexa->phases_num.push_back(j + 2);
                              }
                              else
                              {
                                 if (abs(hh2 - phases[j + 3].h) < 1e-12)
                                 {
                                    real phh3 = phases[j + 2].h - phases[j + 3].h;
                                    hexa->Sm[j + 1] = hh / phh2;
                                    hexa->Sm[j + 2] = hh / phh3;
                                    hexa->Sm[j] = 1. - hexa->Sm[j + 1] - hexa->Sm[j + 2];
                                    hexa->phases_num.push_back(j + 1);
                                    hexa->phases_num.push_back(j);
                                    hexa->phases_num.push_back(j + 2);
                                 }
                                 else if(hh2 < phases[j + 3].h)
                                 {
                                    hexa->Sm[j] = hh / (hh1 - phases[j + 1].h);
                                    hexa->Sm[j + 1] = hh / phh2;
                                    hexa->Sm[j + 2] = 1. - hexa->Sm[j] - hexa->Sm[j + 1];
                                    hexa->phases_num.push_back(j + 1);
                                    hexa->phases_num.push_back(j);
                                    hexa->phases_num.push_back(j + 2);
                                 }
                                 //else if ...
                              }
                           }
                        }
                     }
                  }					
               }
            }
      }

      for (auto h : mesh->hexas)
      {
         for (auto ph : h->phases_num)
            h->lam += phases[ph].k / phases[ph].n;
         real K = 0.;
         for (int p = 0; p < pors.size(); p++)
         {
            bool isIn = true;
            for (int i = 0; i < 8 && isIn; i++)
               isIn = ( ( pors[p].corner0.x - 1e-12 < mesh->knots[h->knots_num[i]]->x && mesh->knots[h->knots_num[i]]->x + 1e-12 < pors[p].corner1.x) &&
                        ( pors[p].corner0.y - 1e-12 < mesh->knots[h->knots_num[i]]->y && mesh->knots[h->knots_num[i]]->y + 1e-12 < pors[p].corner1.y) &&
                        ( pors[p].corner0.z - 1e-12 < mesh->knots[h->knots_num[i]]->z && mesh->knots[h->knots_num[i]]->z + 1e-12 < pors[p].corner1.z) );
            if (isIn)
               h->lam *= pors[p].K;
         }
      }

   }

   
}