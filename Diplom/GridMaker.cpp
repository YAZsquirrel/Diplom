#include "GridMaker.h"
#include <set>
#include <algorithm>
#include <fstream>
#include <cmath>
namespace mesh_comps
{
   Mesh::Mesh()
   {
      std::ifstream fparam("EnvParams.txt");
      // get size of an env and steps
      {
         real x, y, z;
         fparam >> x;

         fparam >> x >> y >> z;
         env_corner1 = knot(x, y, z);
         fparam >> x >> y >> z;
         env_corner2 = knot(x, y, z);
         fparam >> x >> y >> z;
         step = knot(x, y, z);
      }

      // get info on wells
      {
         int well_num;
         real x, y, v, well_rad, h1, h2; // v = tetha,  
         fparam >> well_num;
         w_info.wells.reserve(well_num);

         for (size_t i = 0; i < well_num; i++)
         {
            fparam >> x >> y >> h1 >> h2 >> v >> well_rad;
            w_info.wells.push_back( well( x, y, well_rad, h1, h2, v ) );
         }
         fparam >> w_info.conc_rad >> w_info.rad_knots >> w_info.conc;
         //w_info.conc_rad = std::max(std::ceil(w_info.conc_rad / (0.9 * step.x)) * step.x , std::ceil(w_info.conc_rad / (0.9 * step.y)) * step.y);
      }

      fparam.close();
   }

   void Mesh::FindNeighborsAndFaces()
   {
      int lfn[6][4]{{0,1,2,3}, // z- // local face num
                    {4,5,6,7}, // z+
                    {2,0,6,4}, // x-
                    {3,1,7,5}, // x+
                    {0,1,4,5}, // y-
                    {2,3,6,7}};// y+  // низ, верх, перед, зад, права, лева  // если смотреть в (+x,0,0)
      std::vector<std::list<int>> arr;
      arr.resize(knots.size());
      for (int e = 0; e < hexas.size(); e++)
         for(int i = 0; i < 8; i++)
            arr[hexas[e]->knots_num[i]].push_back(e);

      for (int e = 0; e < hexas.size(); e++)
      {
         for (int f = 0; f < 6; f++)
         {
            short match = 0;
            int gf = 0;
            for (gf = 0; gf < faces.size() && match < 3; gf++)
            {
               match = 0;
               for (int i = 0; i < 4; i++)
                  for (int j = 0; j < 4; j++)
                     if (faces[gf]->knots_num[i] == hexas[e]->knots_num[lfn[f][j]]) match++;
            }
            gf--;
            if (match < 3)
            {
               face* _face;
               faces.push_back(_face = new face());
               for (int j = 0; j < 4; j++)
                  _face->knots_num[j] = hexas[e]->knots_num[lfn[f][j]];  
               _face->hexa_nums.push_back(e);
               hexas[e]->faces_num.push_back(faces.size() - 1);
            }
            else
            {
               // уже в списке
               faces[gf]->hexa_nums.push_back(e);
               hexas[e]->faces_num.push_back(gf);
            }
         }
      }

      neighbors.resize(hexas.size());
      for (int e = 0; e < hexas.size(); e++)
      {
         for (int f = 0; f < 6; f++)
         {
            std::set<int> ems1;
            std::set<int> ems2;
            std::set<int> ems3;
            std::set<int> ems12;
            std::set<int> finalems;
      
            for (auto es : arr[hexas[e]->knots_num[lfn[f][0]]])
               ems1.insert(es);
            for (auto es : arr[hexas[e]->knots_num[lfn[f][1]]])
               ems2.insert(es);
            for (auto es : arr[hexas[e]->knots_num[lfn[f][2]]])
               ems3.insert(es);
      
            std::set_intersection(ems1.begin(), ems1.end(),
                                  ems2.begin(), ems2.end(),
                                  std::inserter(ems12, ems12.begin()));
            std::set_intersection(ems12.begin(), ems12.end(),
                                  ems3.begin(), ems3.end(),
                                  std::inserter(finalems, finalems.begin()));
            for (auto el : finalems)
            {
               if (e != el)
               {
                  neighbors[e].insert(el);
                  neighbors[el].insert(e);
               }
            }
         }
      }
   }

   void Mesh::GenerateMesh()
   {
      set_layers();
      // env_corner1, env_corner2, step
      std::ofstream fknots("Knots.txt");
      fknots.precision(15);

      xn = (int)((env_corner2.x - env_corner1.x) / step.x) + 1;
      yn = (int)((env_corner2.y - env_corner1.y) / step.y) + 1;
      zn = (int)((env_corner2.z - env_corner1.z) / step.z) + 1;

      int plane_size = xn * yn;
      int total_knots = plane_size;// * zn; // 
      int* well_start_indecies = new int[w_info.wells.size()];
      std::vector<int>** inwell_indecies = new std::vector<int>*[w_info.wells.size()];
      knots.reserve(total_knots);
      real hx = step.x, hy = step.y, hz = step.z;
      bool is_near_well = false;
      real *zs = new real[zn]{};

      real z = env_corner1.z;
      for (int k = 0; k < zn; k++)
      {
         real z = k * hz;
         for (int l = 0; l < layers.size(); l++)            
            if (z + hz > layers[l] && z < layers[l])
            {
               z = layers[l];
               break;
            }
         zs[k] = z;



         if (k == 0)
         {
            std::vector<int> wells_area_x;
            std::vector<bool> wells_area_x_found;
            wells_area_x.resize(w_info.wells.size(), 0);
            wells_area_x_found.resize(w_info.wells.size(), false);
            bool is_prev_near_well = false;
            int last_well_num = -1;

            for (int j = 0; j < yn; j++)
            {
               real y = j * hy;
               for (int i = 0; i < xn; i++)
               {
                  is_near_well = false;
                  real x = i * hx;
                  int m = 0;
                  for (; m < w_info.wells.size() && !is_near_well; m++)
                     is_near_well = abs(w_info.wells[m].x - x) < w_info.conc_rad - 1e-12 && abs(w_info.wells[m].y - y) < w_info.conc_rad - 1e-12;
                  if (!is_near_well)
                  {
                     if (is_prev_near_well && !wells_area_x_found[last_well_num])
                        wells_area_x_found[last_well_num] = true;
                     knots.push_back(new knot(x, y, z));
                  }
                  else 
                  {
                     last_well_num = m - 1;
                     if (!wells_area_x_found[last_well_num]) wells_area_x[last_well_num]++;
                     plane_size--;
                  }
                  is_prev_near_well = is_near_well;
               }
            }

            for (int wk = 0; wk < w_info.wells.size(); wk++)
            {
               well *w = &w_info.wells[wk];
               real lxw = w->x - w_info.conc_rad, 
                    rxw = w->x + w_info.conc_rad, 
                    yuw = w->y + w_info.conc_rad, 
                    ydw = w->y - w_info.conc_rad; 
               real lx = lxw - (step.x - 1e-9), 
                    rx = rxw + (step.x - 1e-9), 
                    yu = yuw + (step.y - 1e-9), 
                    yd = ydw - (step.y - 1e-9);

               int well_area_x = wells_area_x[wk],
                   well_area_y = floor(2. * w_info.conc_rad / step.y);
               bool x1 = true, y1 = true;
               

               std::vector<int> inwell_inds;
               inwell_inds.reserve(2 * (well_area_y + well_area_x) + 4);

               for (int kn = 0; kn < knots.size() && (x1 || y1) && 
                        ( w->x + w_info.conc_rad + step.x > knots[kn]->x || w->y + w_info.conc_rad + step.y > knots[kn]->y); kn++)
               {
                  if (knots[kn]->x > rx && knots[kn]->y > yu)
                     break;
                  if (knots[kn]->x < rx && knots[kn]->x > lx && knots[kn]->y < yu && knots[kn]->y > yd)
                     inwell_inds.push_back(kn);
               }
               well_area_y = inwell_inds.size() / 2 - 2 - well_area_x;

               int rays = inwell_inds.size();
               // sort clockwise

               std::vector<int> *sorted_inds = new std::vector<int>;

               sorted_inds->reserve(inwell_inds.size());
               for (int iw = well_area_x, iw1 = 1; iw1 < well_area_y + 1; iw1++) // ^ +
                  sorted_inds->push_back(inwell_inds[iw + 2 * iw1]);

               for (int iw = well_area_x + 1; iw >= 0; iw--) // ->
                  sorted_inds->push_back(inwell_inds[inwell_inds.size() - 1 - iw]);

               for (int iw = well_area_x+1, iw1 = well_area_y; iw1 >= 1;  iw1--) // v +
                  sorted_inds->push_back(inwell_inds[iw + 2 * iw1]);

               for (int iw = well_area_x + 1; iw >= 0; iw--) // <-
                  sorted_inds->push_back(inwell_inds[iw]);

               inwell_indecies[wk] = sorted_inds;

               rays = sorted_inds->size();
               if (wk == 0) well_start_indecies[0] = knots.size();
               else well_start_indecies[wk] = well_start_indecies[wk-1] + w_info.rad_knots * rays;

               plane_size += w_info.rad_knots * rays;
               //knots.resize(plane_size);

               
               for (int iw = 0; iw < rays; iw++)
               {
                  real xiw = knots[sorted_inds->at(iw)]->x - w->x, yiw = knots[sorted_inds->at(iw)]->y - w->y;
                  real len = sqrt(xiw * xiw + yiw * yiw);
                  xiw /= len; yiw /= len;
                  //xiw += w->x;  yiw += w->y;

                  knots.push_back(new knot(xiw * w->radius + w->x, yiw * w->radius + w->y, z));

                  real length = sqrt(pow(xiw * (w->radius/* - step.x*/) + w->x - knots[sorted_inds->at(iw)]->x, 2)
                                   + pow(yiw * (w->radius/* - step.y*/) + w->y - knots[sorted_inds->at(iw)]->y, 2));
                  real stepr = length / w_info.rad_knots;

                  //for (float r = stepr, j = 1; j < num_of_knots_on_rad; r += stepr, j++)
                  // b1 = length * (1. - w_info.conc) / (1. - pow(w_info.conc, w_info.rad_knots)) 
                  for (real r = abs(w_info.conc - 1.) > 1e-12 ? length * (1. - w_info.conc) / (1. - pow(w_info.conc, w_info.rad_knots)) : stepr,
                  //for (real r = w_info.conc != 1. ? length * (1. - w_info.conc) * w_info.conc : stepr,
                     rr = 1, b = r; rr < w_info.rad_knots; rr++)
                  {
                     knots.push_back(new knot(xiw * (w->radius + r) + w->x, yiw * (w->radius + r) + w->y, z));
                     b *= w_info.conc; r += b;
                  }
               }
            }
         }
         else
            for (int i = 0; i < plane_size; i++)
               knots.push_back(new knot(knots[i]->x, knots[i]->y, z));
      }

      fknots << knots.size() << '\n';
      for (int i = 0; i < knots.size(); i++)
         fknots << knots[i]->x << " " << knots[i]->y << " " << knots[i]->z << '\n';
      fknots.close();
      if (well_start_indecies[0] < 0)
         well_start_indecies = new int[1]{ plane_size };
      FindAllHexasAndBounds(plane_size, well_start_indecies, inwell_indecies, zs);

      delete[] well_start_indecies;
   }

   void Mesh::SetSignsForHexaNormals()
   {
      int lfn[6][4]{ {0,1,2,3}, // local face num
                    {4,5,6,7},
                    {2,0,6,4},
                    {1,3,5,7},
                    {0,1,4,5},
                    {3,2,7,6} };

      for (auto el : hexas)
      {
         
         for (int j = 0; j < 6; j++)
         {

            knot *k1 = knots[el->knots_num[lfn[j][0]]],
                 *k2 = knots[el->knots_num[lfn[j][1]]],
                 *k3 = knots[el->knots_num[lfn[j][2]]];
            knot normal = knot();
            real x1 = k2->x - k1->x;
            real y1 = k2->y - k1->y;
            real z1 = k2->z - k1->z;
            real x2 = k3->x - k1->x;
            real y2 = k3->y - k1->y;
            real z2 = k3->z - k1->z;
            normal.x = y1 * z2 - z1 * y2;
            normal.y = z1 * x2 - x1 * z2;
            normal.z = x1 * y2 - y1 * x2;
            real len = sqrt((normal.x) * (normal.x) +
               (normal.y) * (normal.y) +
               (normal.z) * (normal.z));
            normal.x /= len;
            normal.y /= len;
            normal.z /= len;

            auto f = faces[el->faces_num[j]];
            el->faces_sign[j] = sign((normal.x) * (f->normal.x) + (normal.y) * (f->normal.y) + (normal.z) * (f->normal.z));
         }
      }
   }

   void Mesh::FindFaceNormals()
   {
      for (auto f : faces)
      {
         knot *k1 = knots[f->knots_num[0]], 
              *k2 = knots[f->knots_num[1]], 
              *k3 = knots[f->knots_num[2]];
         real x1 = k2->x - k1->x;
         real y1 = k2->y - k1->y;
         real z1 = k2->z - k1->z;
         real x2 = k3->x - k1->x;
         real y2 = k3->y - k1->y;
         real z2 = k3->z - k1->z;
         f->normal.x = y1*z2-z1*y2;
         f->normal.y = z1*x2-x1*z2;
         f->normal.z = x1*y2-y1*x2;
         real len = sqrt((f->normal.x)*(f->normal.x)+
                         (f->normal.y)*(f->normal.y)+
                         (f->normal.z)*(f->normal.z));
         f->normal.x /= len;
         f->normal.y /= len;
         f->normal.z /= len;
      }

   }

   void Mesh::FindAllHexasAndBounds(int plain_size, int* well_inds, std::vector<int>** inwell_inds, real* zs)
   {
      std::vector<int*> lower_faces;
      std::vector<int*> upper_faces;
      std::vector<std::vector<int*>> lower_bound2_edges;

      std::ofstream fhexas("Hexahedrons.txt");

      std::ofstream fbounds1("FirstBounds.txt");
      std::ofstream fbounds2("SecondBounds.txt");
      std::ifstream fFP("FiltrParams.txt");
      lower_bound2_edges.reserve(w_info.wells.size());
      for (int i = 0; i < w_info.wells.size(); i++)
         lower_bound2_edges.reserve(inwell_inds[i]->size());

      lower_faces.reserve(2 * plain_size);
      upper_faces.reserve(2 * plain_size);

      bool inwell_prevx = false;
      // add knots without well
      for (size_t i = 0, y = 0, x = 0; i < well_inds[0] && y < yn; i++, x++)
      {
         if (x > xn - 1)
         {
            x = 0;
            y++;
         }

         int *face = new int[4]{};
         bool found = false;

         face[0] = i;
         if (abs(knots[i]->x - knots[i + 1]->x) <= (step.x + 1e-10) && abs(knots[i]->y - knots[i + 1]->y) < 1e-10)  // o -> o
         {
            face[1] = i + 1;
         for (int j = 0; j < well_inds[0]; j++)
            if (abs(knots[i]->y - knots[j]->y) <= (step.y + 1e-10) && abs(knots[i]->x - knots[j]->x) < 1e-10 && i < j)
            {  
               face[2] = j;
               if (abs(knots[j]->x - knots[j + 1]->x) <= (step.x + 1e-10) && abs(knots[j]->y - knots[j + 1]->y) < 1e-10 && i < j)
               {
                  face[3] = j + 1;
                  found = true;
                  break;
               }
            }
         }
         if (found)
            lower_faces.push_back(face);
      }

      //add well knots
      for (int w = 0; w < w_info.wells.size(); w++)
      {
         std::vector<int*> bound_edges_well;
         for (int i = well_inds[w], iw = 0; iw < inwell_inds[w]->size() - 1; i += w_info.rad_knots, iw++) // r = inwell_inds[w]->size() - колво лучей, r * nr = колво в 
            for (int r = 0; r < w_info.rad_knots; r++)
            {
               int* face = new int[4]{};
               face[0] = i + r; face[1] = i + w_info.rad_knots + r;
               if (r < w_info.rad_knots - 1)
                  { face[2] = i + r + 1; face[3] = i + r + 1 + w_info.rad_knots; }
               else
                  { face[2] = inwell_inds[w]->at(iw); face[3] = inwell_inds[w]->at(iw + 1); }
               lower_faces.push_back(face);

               if (r == 0) bound_edges_well.push_back(new int[2]{face[0], face[1]});
            }
         for (int r = 0; r < w_info.rad_knots; r++)
         {
            int* face = new int[4]{};
            face[0] = well_inds[w] + w_info.rad_knots * (inwell_inds[w]->size() - 1) + r; 
            face[1] = well_inds[w] + r;
            if (r < w_info.rad_knots - 1)
            {
               face[2] = well_inds[w] + w_info.rad_knots * (inwell_inds[w]->size() - 1) + r + 1; 
               face[3] = well_inds[w] + r + 1;
            }
            else
            {
               face[2] = inwell_inds[w]->at(inwell_inds[w]->size() - 1); 
               face[3] = inwell_inds[w]->at(0);
            }
            lower_faces.push_back(face);
            if (r == 0) bound_edges_well.push_back(new int[2]{ face[0], face[1] });
         }

         lower_bound2_edges.push_back(bound_edges_well);
      }

      std::vector<int*> side_bound1;
      side_bound1.reserve(2 * (xn + yn - 2) * (zn - 1));

      // add sides (I)
      for (int i = 0; i < well_inds[0]; i++)
         for (int j = i + 1; j < well_inds[0]; j++)
         {
            if ((abs(knots[i]->y - env_corner1.y) < 1e-10 || abs(knots[i]->y - env_corner2.y) < 1e-10) 
                  && abs(knots[i]->x - knots[j]->x) < step.x + 1e-10 && abs(knots[i]->y - knots[j]->y) < 1e-10)
               {side_bound1.push_back(new int[2]{i, j}); continue;}
            if ((abs(knots[i]->x - env_corner1.x) < 1e-10 || abs(knots[i]->x - env_corner2.x) < 1e-10)
               && abs(knots[i]->y - knots[j]->y) < step.y + 1e-10 && abs(knots[i]->x - knots[j]->x) <  1e-10)
               {side_bound1.push_back(new int[2]{ i, j }); break;}
         }

      int plain_knot_size = plain_size;
      fbounds1 << 2 * (xn + yn - 2) * (zn - 1) << '\n';

      plain_size = lower_faces.size();
      fhexas << plain_size * (zn - 1) << '\n';
      int size = 0;
      for (int i = 0; i < w_info.wells.size(); i++)
         size += lower_bound2_edges[i].size();

      fbounds2 << size * (zn - 1) << '\n';// + plain_size * 2 << '\n';
      real P_plast;
      fFP >> P_plast;
      fFP.close();

      // add bottom layer (II)
      //for (int i = 0; i < plain_size; i++)
      //{
      //   fbounds2 << lower_faces[i][0] << " " 
      //            << lower_faces[i][1] << " " 
      //            << lower_faces[i][2] << " " 
      //            << lower_faces[i][3] << " " 
      //            << 0.0 << '\n';
      //}


      // add knots on upper plains 
      for (int k = 1; k < zn; k++)
      {

         for (int i = 0; i < lower_faces.size(); i++)
         {
            for (int j = 0; j < 4; j++)
               fhexas << lower_faces[i][j] + plain_knot_size * (k - 1) << " ";
            fhexas << " ";
            for (int j = 0; j < 4; j++)
               fhexas << lower_faces[i][j] + plain_knot_size * k << " ";
            fhexas << '\n';

         }
         //if (k == zn - 1)
         //{
         //   for (int i = 0; i < plain_size; i++)
         //   {
         //      fbounds2 << lower_faces[i][0] + plain_knot_size * k << " "
         //         << lower_faces[i][1] + plain_knot_size * k << " "
         //         << lower_faces[i][2] + plain_knot_size * k << " "
         //         << lower_faces[i][3] + plain_knot_size * k << " "
         //         << 0.0 << '\n';
         //   }
         //}

         for (int j = 0; j < side_bound1.size(); j++)
         {
            fbounds1 << side_bound1[j][0] + plain_knot_size * (k - 1) << " "
               << side_bound1[j][1] + plain_knot_size * (k - 1) << " "
               << side_bound1[j][0] + plain_knot_size * k << " "
               << side_bound1[j][1] + plain_knot_size * k << " "
               << P_plast << '\n';
         }

         for (int i = 0; i < w_info.wells.size(); i++)
         {
            real v = 0.0;
            real zmid = (zs[k - 1] + zs[k]) / 2.;
            if (w_info.wells[i].h2 < zmid && zmid < w_info.wells[i].h1)
               v = w_info.wells[i].intake;
            //else continue;

            for (int j = 0; j < lower_bound2_edges[i].size(); j++)
            {
               int *edge = lower_bound2_edges[i][j];
               int *knots = new int[4]{ edge[0] + plain_knot_size * (k - 1), edge[1] + plain_knot_size * (k - 1),
                           edge[0] + plain_knot_size * k, edge[1] + plain_knot_size * k };
               w_info.wells[i].faces_knots_num.push_back(knots);
               fbounds2 << knots[0] << " "
                        << knots[1] << " "
                        << knots[2] << " "
                        << knots[3] << " "
                        << v << '\n';

            }
         }
      }
      fhexas.close();
      fbounds1.close();
      fbounds2.close();

   }

}