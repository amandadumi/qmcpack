#include "SNADesc.h"
#include <iostream>    


namespace qmcplusplus
{
    SNADesc::SNADesc(){}

    SNADesc::~SNADesc(){}

      void SNADesc::init(int Npart, int input_twojmax,double input_rcut){
      twojmax = input_twojmax;
      nelements = 1;
      nmax = Npart -1;
      rmin0 = 0.0;
      compute_ncoeff();
      build_indexlist();
      
      create_arrays();
      init_clebsch_gordan();
      init_rootpqarray();

     rij = std::vector<std::vector<double>>(nmax,std::vector<double>(3,0.0));
     inside = std::vector<int>(Npart,Npart);
     wj = std::vector<double>(Npart,1.0);
     rcutij = std::vector<double>(Npart,1.0);
     element = std::vector<int>(Npart,1);
        
     bzero_flag = 0;
    };

  void SNADesc::create_arrays(){

  int jdimpq = twojmax + 2;
  rootpqarray = std::vector<std::vector<double>>(jdimpq,std::vector<double>(jdimpq,0));
  cglist = std::vector<double>(idxcg_max,0.0);
  ulisttot_r = std::vector<double>(idxu_max*nelements, 0.0);
  ulisttot_i = std::vector<double>(idxu_max*nelements, 0.0);
  dulist_r = std::vector<std::vector<double>>(idxu_max,std::vector<double>(3,0.0));
  dulist_i = std::vector<std::vector<double>>(idxu_max,std::vector<double>(3,0.0));
  zlist_r = std::vector<double>(idxz_max*ndoubles, 0.0);
  zlist_i = std::vector<double>(idxz_max*ndoubles, 0.0);
  blist = std::vector<double>(idxb_max*ntriples, 0.0);
  dblist = std::vector<std::vector<double>>(idxb_max*ntriples,std::vector<double>(3,0));
  ylist_r = std::vector<double>(idxu_max*nelements, 0.0);
  ylist_i = std::vector<double>(idxu_max*nelements, 0.0);
  ulist_r_ij = std::vector<std::vector<double>>(nmax, std::vector<double>(idxu_max,0.0));
  ulist_i_ij = std::vector<std::vector<double>>(nmax, std::vector<double>(idxu_max,0.0));
}

  void SNADesc::compute_ncoeff()
  {
    int ncount;

    ncount = 0;

    for (int j1 = 0; j1 <= twojmax; j1++)
      for (int j2 = 0; j2 <= j1; j2++)
        for (int j = j1 - j2;
             j <= std::min(twojmax, j1 + j2); j += 2)
          if (j >= j1) ncount++;

    ndoubles = nelements*nelements;
    ntriples = nelements*nelements*nelements;
    if (chem_flag)
      ncoeff = ncount*ntriples;
    else
      ncoeff = ncount;
  }

  double SNADesc::compute_sfac(double r, double rcut)
  {
    if (switch_flag == 0) return 1.0;
    if (switch_flag == 1) {
      if (r <= rmin0) return 1.0;
      else if (r > rcut) return 0.0;
      else {
        double rcutfac = mypi / (rcut - rmin0);
        return 0.5 * (std::cos((r - rmin0) * rcutfac) + 1.0);
      }
    }
    return 0.0;
  }

  double SNADesc::compute_dsfac(double r, double rcut)
  {
    if (switch_flag == 0) return 0.0;
    if (switch_flag == 1) {
      if (r <= rmin0) return 0.0;
      else if (r > rcut) return 0.0;
      else {
      double rcutfac = mypi / (rcut - rmin0);
      return -0.5 * std::sin((r - rmin0) * rcutfac) * rcutfac;
      }
    }
    return 0.0;
  }

  void SNADesc::compute_ui(int jnum, int ielem){
  double rsq, r, x, y, z, z0, theta0;

  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  zero_uarraytot(ielem);
  for (int j = 0; j < jnum; j++) {
    x = rij[j][0];
    y = rij[j][1];
    z = rij[j][2];
    rsq = x * x + y * y + z * z;
    r = std::sqrt(rsq);

    theta0 = (r - rmin0) * rfac0 * mypi / (rcutij[j] - rmin0);
    //    theta0 = (r - rmin0) * rscale0;
    z0 = r / tan(theta0);
    //std::cout <<"x "<< x<<" y "<<y<<" z " <<z<<" r " <<r<<" theta0 "<<theta0<<" z0 "<<z0<<std::endl;

    compute_uarray(x, y, z, z0, r, j);
//  std::cout<< "SNADesc::compute_ui after compute_uarray call" <<std::endl;
    if (chem_flag)
      add_uarraytot(r, wj[j], rcutij[j], j, element[j]);
    else
      add_uarraytot(r, wj[j], rcutij[j], j, 0);
  }
  // std::cout << "ulistot_i is " <<std::endl;
   //for (int i = 0; i < idxu_max*nelements; ++i){
   //  std::cout << ulisttot_i[i] << " " ;
  // }
  // std::cout << "ulistot_r is " <<std::endl;
  // for (int i = 0; i < idxu_max*nelements; ++i){
  //   std::cout << ulisttot_r[i] << " " ;
  // }

}

  void SNADesc::build_indexlist()
  {

  // index list for cglist

  int jdim = twojmax + 1;
  idxcg_block = std::vector<std::vector<std::vector<int>>>(jdim,std::vector<std::vector<int>>(jdim,std::vector<int>(jdim,0)));
  int idxcg_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
        // std::cout << j1 << j2 <<j << " "<<idxcg_count<<std::endl;
        idxcg_block[j1][j2][j] = idxcg_count;
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }



  idxcg_max = idxcg_count;
  // std::cout << "idxcg_max " << idxcg_max <<std::endl;
  // index list for uarray
  // need to include both halves

  idxu_block = std::vector<int>(jdim, 0);
  int idxu_count = 0;

  for (int j = 0; j <= twojmax; j++) {
    idxu_block[j] = idxu_count;
    // std::cout << "idxublock " << idxu_count <<"\n";
    for (int mb = 0; mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++)
        idxu_count++;


  }
  idxu_max = idxu_count;
  // index list for beta and B

  int idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  idxb_max = idxb_count;
  // std:: cout << "idxbmax is " << idxb_max << "n";
  idxb = new SNA_BINDICES[idxb_max];

  idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          idxb[idxb_count].j1 = j1;
          idxb[idxb_count].j2 = j2;
          idxb[idxb_count].j = j;
          idxb_count++;
        }

  // reverse index list for beta and b

  idxb_block = std::vector<std::vector<std::vector<int>>>(jdim,std::vector<std::vector<int>>(jdim,std::vector<int>(jdim,0)));
  idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          idxb_block[j1][j2][j] = idxb_count;
          idxb_count++;
          // std::cout <<"j1"<< j1<<"j2"<<j2 << "j"<< j << " idxb_block: " << idxb_block[j1][j2][j] << std::endl;
        }
      }

  // index list for zlist

  int idxz_count = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  idxz_max = idxz_count;
  // std::cout << "idxzmax is " << idxz_max << "\n";
  idxz = new SNA_ZINDICES[idxz_max];

  idxz_block = std::vector<std::vector<std::vector<int>>>(jdim,std::vector<std::vector<int>>(jdim,std::vector<int>(jdim,0)));

  idxz_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++){

      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
        idxz_block[j1][j2][j] = idxz_count;

        // find right beta[jjb] entry
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3

        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
          //  std::cout << "SNADesc::build_index in loop. ma "<< ma << " mb "<< mb<< std::endl;
            idxz[idxz_count].j1 = j1;
          //  std::cout << "SNADesc::build_index in loop. j1 "<< j1 <<std::endl;
          //  std::cout << "SNADesc::build_index in loop. idxz j1 val "<< idxz[idxz_count].j1 <<std::endl;
            idxz[idxz_count].j2 = j2;
            idxz[idxz_count].j = j;
            idxz[idxz_count].ma1min = std::max(0, (2 * ma - j - j2 + j1) / 2);
            idxz[idxz_count].ma2max = (2 * ma - j - (2 * idxz[idxz_count].ma1min - j1) + j2) / 2;
            idxz[idxz_count].na = std::min(j1, (2 * ma - j + j2 + j1) / 2) - idxz[idxz_count].ma1min + 1;
            idxz[idxz_count].mb1min = std::max(0, (2 * mb - j - j2 + j1) / 2);
            idxz[idxz_count].mb2max = (2 * mb - j - (2 * idxz[idxz_count].mb1min - j1) + j2) / 2;
            idxz[idxz_count].nb = std::min(j1, (2 * mb - j + j2 + j1) / 2) - idxz[idxz_count].mb1min + 1;
            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
            const int jju = idxu_block[j] + (j+1)*mb + ma;
            idxz[idxz_count].jju = jju;

            idxz_count++;
          }
      }
    }
}
    void SNADesc::compute_zi(){
     std::cout<< "SNADesc::compute_zi entered function" <<std::endl;

  int idouble = 0;
  //double * zptr_r;
  //double * zptr_i;
  for (int elem1 = 0; elem1 < nelements; elem1++)
    for (int elem2 = 0; elem2 < nelements; elem2++) {
   //  std::cout<< "SNADesc::compute_zi elem1 "<< elem1 << " elem2 "<< elem2 <<std::endl;

      int zidx = idouble*idxz_max;
      //zptr_i = &zlist_i[idouble*idxz_max];
   //  std::cout<< "SNADesc::compute_zi zidx is " << zidx <<std::endl;
   ////  std::cout<< "SNADesc::compute_zi idxz_max is " << idxz_max <<std::endl;

      for (int jjz = 0; jjz < idxz_max; jjz++) {
      ////  std::cout<< "SNADesc::compute_zi jjz" << jjz <<std::endl;
        const int j1 = idxz[jjz].j1;
      //  std::cout<< "SNADesc::compute_zi j1" << j1 <<std::endl;
        const int j2 = idxz[jjz].j2;
        const int j = idxz[jjz].j;
        const int ma1min = idxz[jjz].ma1min;
        const int ma2max = idxz[jjz].ma2max;
        const int na = idxz[jjz].na;
        const int mb1min = idxz[jjz].mb1min;
        const int mb2max = idxz[jjz].mb2max;
        const int nb = idxz[jjz].nb;
        //std::cout<< "SNADesc::compute_zi accessed the snaindex z struct okay"<<std::endl;
       //  std::cout << j1<< " "<< j2<< " "<< j<<" " <<ma1min<<" " <<ma2max<<" " <<na<<" " <<mb1min<<" " <<mb2max<<std::endl;

        int cgblock = idxcg_block[j1][j2][j];

        zlist_r[zidx + jjz] = 0.0;
        zlist_i[zidx + jjz] = 0.0;

        int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
        int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
        int icgb = mb1min * (j2 + 1) + mb2max;
        // std::cout << jju1<< " "<<jju2<< " "<<icgb<<std::endl;

        for (int ib = 0; ib < nb; ib++) {
          double suma1_r = 0.0;
          double suma1_i = 0.0;

          //const double *u1_r = &ulisttot_r[elem1*idxu_max+jju1];
          //const double *u1_i = &ulisttot_i[elem1*idxu_max+jju1];
          //const double *u2_r = &ulisttot_r[elem2*idxu_max+jju2];
          //const double *u2_i = &ulisttot_i[elem2*idxu_max+jju2];
          int u1 = elem1*idxu_max+jju1;
          int u2 = elem2*idxu_max+jju2;

          int ma1 = ma1min;
          int ma2 = ma2max;
          int icga = ma1min * (j2 + 1) + ma2max;
           //std::cout << icga <<std::endl;
          for (int ia = 0; ia < na; ia++) {
             //std::cout<< cglist[cgblock+icga]<<std::endl;
             //std::cout<< ulisttot_r[u1+ma1] << " "<< ulisttot_r[u2+ma2] <<" " << ulisttot_i[u1+ma1] <<" "<< ulisttot_i[u2+ma2]<<std::endl;
            suma1_r += cglist[cgblock+ icga] * (ulisttot_r[u1+ma1] * ulisttot_r[u2+ma2] - ulisttot_i[u1+ma1] * ulisttot_i[u2+ma2]);
            suma1_i += cglist[cgblock+ icga] * (ulisttot_r[u1+ma1] * ulisttot_i[u2+ma2] + ulisttot_i[u1+ma1] * ulisttot_r[u2+ma2]);
            ma1++;
            ma2--;
            icga += j2;
          } // end loop over ia
          //std::cout<< "sum_a1_i is " << suma1_i << std::endl;
          zlist_r[zidx+jjz] += cglist[cgblock + icgb] * suma1_r;
          zlist_i[zidx+jjz] += cglist[cgblock + icgb] * suma1_i;

          jju1 += j1 + 1;
          jju2 -= j2 + 1;
          icgb += j2;
        } // end loop over ib
        if (bnorm_flag) {
          zlist_r[zidx+jjz] /= (j+1);
          zlist_i[zidx+jjz] /= (j+1);
        }
      } // end loop over jjz
      idouble++;
    }
// std::cout << "zlist r is "<<std::endl;
// for (int i = 0 ; i <idxz_max*ndoubles; i ++){
//   std::cout << zlist_r[i]  << " " ;
//   }
//   std::cout <<std::endl;
  // std::cout << "zlist i is "<<std::endl;
  // for (int i = 0 ; i <idxz_max*ndoubles; i ++){
  //   std::cout << zlist_i[i] << " ";
  //   }
}
    void SNADesc::compute_yi(std::vector<double> beta){
  // std::cout << "computing yi"<< std::endl;
  int jju;
  double betaj;
  int itriple;

  for (int ielem1 = 0; ielem1 < nelements; ielem1++)
    for (int j = 0; j <= twojmax; j++) {
      jju = idxu_block[j];
      for (int mb = 0; 2*mb <= j; mb++)
        for (int ma = 0; ma <= j; ma++) {
          ylist_r[ielem1*idxu_max+jju] = 0.0;
          ylist_i[ielem1*idxu_max+jju] = 0.0;
          jju++;
        } // end loop over ma, mb
    } // end loop over j

  for (int elem1 = 0; elem1 < nelements; elem1++)
    for (int elem2 = 0; elem2 < nelements; elem2++) {
        for (int jjz = 0; jjz < idxz_max; jjz++) {
          // std::cout << "elem1 " <<elem1<<"elem2 "<<elem2<<"jjz "<<jjz<<std::endl;
          const int j1 = idxz[jjz].j1;
          // std::cout<<" j1 is "<<j1<<std::endl;
          const int j2 = idxz[jjz].j2;
          const int j = idxz[jjz].j;
          const int ma1min = idxz[jjz].ma1min;
          const int ma2max = idxz[jjz].ma2max;
          const int na = idxz[jjz].na;
          const int mb1min = idxz[jjz].mb1min;
          const int mb2max = idxz[jjz].mb2max;
          const int nb = idxz[jjz].nb;
          // std::cout<<" nb is "<<nb<<std::endl;

          int cgblock = idxcg_block[j1][j2][j];
          // std::cout << "idxcgblock is " << idxcg_block[j1][j1][j] << std::endl;
          // std::cout <<"cglist"<< *cglist <<std::endl;
          double ztmp_r = 0.0;
          double ztmp_i = 0.0;

          int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
          // std::cout <<"idxublock value"<< idxu_block[j1]  << "+"<< (j1+1) << "*" << mb1min<< " "<<jju1 <<std::endl;

          int jju2 = idxu_block[j2] + (j2 + 1) + mb2max;
          // std::cout <<"idxublock value"<< idxu_block[j2]  << "+"<< (j2+1) << "+" << mb2max <<" " << jju2<<std::endl;

          int icgb = mb1min * (j2 + 1) + mb2max;
          // std::cout<< "starting the sum! " << std::endl;
          for (int ib = 0; ib < nb; ib++) {
            // std::cout << "e1 " << elem1 << " e2 " << elem2 << " jjz " << jjz << " ib " << ib << std::endl;
            double suma1_r = 0.0;
            double suma1_i = 0.0;

            //const double *u1_r = &ulisttot_r[elem1*idxu_max+jju1];
            //const double *u1_i = &ulisttot_i[elem1*idxu_max+jju1];
            //const double *u2_r = &ulisttot_r[elem2*idxu_max+jju2];
            //const double *u2_i = &ulisttot_i[elem2*idxu_max+jju2];
            int u1 = elem1*idxu_max+jju1;
            int u2 = elem2*idxu_max+jju2;

            int ma1 = ma1min;
            int ma2 = ma2max;
            int icga = ma1min * (j2 + 1) + ma2max;

            for (int ia = 0; ia < na; ia++) {
              //suma1_r += cgblock[icga] * (u1_r[ma1] * u2_r[ma2] - u1_i[ma1] * u2_i[ma2]);
              //suma1_i += cgblock[icga] * (u1_r[ma1] * u2_i[ma2] + u1_i[ma1] * u2_r[ma2]);
              suma1_r += cglist[cgblock+ icga] * (ulisttot_r[u1+ma1] * ulisttot_r[u2+ma2] - ulisttot_i[u1+ma1] * ulisttot_i[u2+ma2]);
              suma1_i += cglist[cgblock+ icga] * (ulisttot_r[u1 + ma1] * ulisttot_i[u2+ma2] + ulisttot_i[u1+ ma1] * ulisttot_r[u2+ma2]);
              ma1++;
              ma2--;
              icga += j2;

            } // end loop over ia

            ztmp_r += cglist[cgblock+ icgb] * suma1_r;
            ztmp_i += cglist[cgblock+ icgb] * suma1_i;

            jju1 += j1 + 1;
            jju2 -= j2 + 1;
            icgb += j2;
          } // end loop over ib

          // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
          // find right y_list[jju] and beta[jjb] entries
          // multiply and divide by j+1 factors
          // account for multiplicity of 1, 2, or 3

        if (bnorm_flag) {
          ztmp_i /= j+1;
          ztmp_r /= j+1;
        }

        jju = idxz[jjz].jju;
        for (int elem3 = 0; elem3 < nelements; elem3++) {
        // pick out right beta value
          if (j >= j1) {
            const int jjb = idxb_block[j1][j2][j];
            itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb;
            // std::cout << "elem2 " <<elem2 << " nelements " << nelements <<" elem1 " <<elem1 <<" idxb_max " << idxb_max << " jjb " <<jjb <<std::endl;

            // std::cout<< "itriple is " << itriple << std::endl;
            if (j1 == j) {
              if (j2 == j) betaj = 3*beta[itriple];
              else betaj = 2*beta[itriple];
            } else betaj = beta[itriple];
          } else if (j >= j2) {
            const int jjb = idxb_block[j][j2][j1];
            itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb;
            if (j2 == j) betaj = 2*beta[itriple];
            else betaj = beta[itriple];
          } else {
            const int jjb = idxb_block[j2][j][j1];
            // std::cout <<j2 <<" "<< j << " " << j1 << " " << jjb <<std::endl;
            itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb;
            // std::cout <<"itriple is " <<itriple <<std::endl;
            // std::cout <<"itriple is " <<elem2 <<" "<< nelements << " " << elem3 << " " << elem1 <<" " << idxb_max << " " <<jjb <<std::endl;
            betaj = beta[itriple];
          }

          if (!bnorm_flag && j1 > j)
            betaj *= (j1 + 1) / (j + 1.0);

          ylist_r[elem3 * idxu_max + jju] += betaj * ztmp_r;
          ylist_i[elem3 * idxu_max + jju] += betaj * ztmp_i;
        }
      } // end loop over jjz
    }

}
    void SNADesc::compute_bi(int ielem){
   //  std::cout<< "SNADesc::compute_bi entered function" <<std::endl;
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

  int itriple = 0;
  int idouble = 0;
  for (int elem1 = 0; elem1 < nelements; elem1++)
    for (int elem2 = 0; elem2 < nelements; elem2++) {

      //double *zptr_r = &zlist_r[idouble*idxz_max];
      //double *zptr_i = &zlist_i[idouble*idxz_max];
      int zidx = idouble*idxz_max;

      for (int elem3 = 0; elem3 < nelements; elem3++) {
        for (int jjb = 0; jjb < idxb_max; jjb++) {
          const int j1 = idxb[jjb].j1;
          const int j2 = idxb[jjb].j2;
          const int j = idxb[jjb].j;

          int jjz = idxz_block[j1][j2][j];
          int jju = idxu_block[j];
          double sumzu = 0.0;
          for (int mb = 0; 2 * mb < j; mb++)
            for (int ma = 0; ma <= j; ma++) {
              sumzu += ulisttot_r[elem3*idxu_max+jju] * zlist_r[zidx+jjz] +
                       ulisttot_i[elem3*idxu_max+jju] * zlist_i[zidx+jjz];
              jjz++;
              jju++;
            } // end loop over ma, mb

          // For j even, handle middle column

          if (j % 2 == 0) {
            int mb = j / 2;
            for (int ma = 0; ma < mb; ma++) {
              sumzu += ulisttot_r[elem3*idxu_max+jju] * zlist_r[zidx+jjz] +
                       ulisttot_i[elem3*idxu_max+jju] * zlist_i[zidx+jjz];
              jjz++;
              jju++;
            }

            sumzu += 0.5 * (ulisttot_r[elem3*idxu_max+jju] * zlist_r[zidx+jjz] +
                            ulisttot_i[elem3*idxu_max+jju] * zlist_i[zidx+jjz]);
          } // end if jeven

          blist[itriple*idxb_max+jjb] = 2.0 * sumzu;

        }
        itriple++;
      }
      idouble++;
    }

  // apply bzero shift

  if (bzero_flag) {
    if (!wselfall_flag) {
      itriple = (ielem*nelements+ielem)*nelements+ielem;
      for (int jjb = 0; jjb < idxb_max; jjb++) {
        const int j = idxb[jjb].j;
        blist[itriple*idxb_max+jjb] -= bzero[j];
      } // end loop over JJ
    } else {
      int itriple = 0;
      for (int elem1 = 0; elem1 < nelements; elem1++)
        for (int elem2 = 0; elem2 < nelements; elem2++) {
          for (int elem3 = 0; elem3 < nelements; elem3++) {
            for (int jjb = 0; jjb < idxb_max; jjb++) {
              const int j = idxb[jjb].j;
              blist[itriple*idxb_max+jjb] -= bzero[j];
            } // end loop over JJ
            itriple++;
          } // end loop over elem3
        } // end loop over elem1,elem2
    }
  }
}
  
    void SNADesc::compute_dbidrj(){
   //  std::cout<< "SNADesc::compute_dbidrj entered function" <<std::endl;
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        zdb = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            zdb +=
  //              Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)
  //        dbdr(j1,j2,j) += 2*zdb
  //        zdb = 0
  //        for mb1 = 0,...,j1mid
  //          for ma1 = 0,...,j1
  //            zdb +=
  //              Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
  //        dbdr(j1,j2,j) += 2*zdb*(j+1)/(j1+1)
  //        zdb = 0
  //        for mb2 = 0,...,j2mid
  //          for ma2 = 0,...,j2
  //            zdb +=
  //              Conj(dudr(j2,ma2,mb2))*z(j1,j,j2,ma2,mb2)
  //        dbdr(j1,j2,j) += 2*zdb*(j+1)/(j2+1)

  int db_idx;
  std::vector<double> sumzdu_r =  std::vector<double>(3,0.0);
  int jjz, jju;

  int idouble;
  int itriple;

  // set all the derivatives to zero once

  for (int jjb = 0; jjb < idxb_max; jjb++) {

    for (int elem1 = 0; elem1 < nelements; elem1++)
      for (int elem2 = 0; elem2 < nelements; elem2++)
        for (int elem3 = 0; elem3 < nelements; elem3++) {

          itriple = (elem1 * nelements + elem2) * nelements + elem3;
          db_idx = itriple*idxb_max+jjb;

          //dbdr = dblist[itriple*idxb_max+jjb];
          dblist[db_idx][0] = 0.0;
          dblist[db_idx][1] = 0.0;
          dblist[db_idx][2] = 0.0;
        }

  }
//  std::cout << "dbdr accessed" << std::endl;

  int elem3 = elem_duarray;
//  std::cout << "elem_duarray"<< elem_duarray << std::endl;

  for (int jjb = 0; jjb < idxb_max; jjb++) {
    const int j1 = idxb[jjb].j1;
    const int j2 = idxb[jjb].j2;
    const int j = idxb[jjb].j;
  //  std::cout << "j1,j2,j3"<< j1 << " " << j2 <<" "<< j<< std::endl;


    // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)

    for (int elem1 = 0; elem1 < nelements; elem1++)
      for (int elem2 = 0; elem2 < nelements; elem2++) {

        jjz = idxz_block[j1][j2][j];
        jju = idxu_block[j];
        idouble = elem1*nelements+elem2;
        itriple = (elem1*nelements+elem2)*nelements+elem3;
        //dbdr = dblist[itriple*idxb_max+jjb];
        db_idx = itriple*idxb_max+jjb;
        //zptr_r = &zlist_r[idouble*idxz_max];
        //zptr_i = &zlist_i[idouble*idxz_max];
        int z_offset = idouble * idxz_max;
        for (int k = 0; k < 3; k++){
          sumzdu_r[k] = 0.0;
        //  std::cout << "sumzdu_r"<< sumzdu_r[k]<< std::endl;
        }
        for (int mb = 0; 2 * mb < j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            //dudr_r = dulist_r[jju];
            //dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++){
              sumzdu_r[k] +=
                  dulist_r[jju][k] * zlist_r[z_offset + jjz] +
                  dulist_i[jju][k] * zlist_i[z_offset + jjz];
            //  std::cout << "sumzdu_r"<< sumzdu_r[k]<< std::endl;
            }
            jjz++;
            jju++;
          } //end loop over ma mb
      ////  std::cout << "end ma/mb loop" << std::endl;
        // For j even, handle middle column

        if (j % 2 == 0) {
          int mb = j / 2;
          for (int ma = 0; ma < mb; ma++) {
            //dudr_r = dulist_r[jju];
            //dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++){
              sumzdu_r[k] +=
                  dulist_r[jju][k] * zlist_r[z_offset + jjz] +
                  dulist_i[jju][k] * zlist_i[z_offset + jjz];
            ////  std::cout << "sumzdu_r"<< sumzdu_r[k]<< std::endl;
            }
            jjz++;
            jju++;
          }
          //dudr_r = dulist_r[jju];
          //dudr_i = dulist_i[jju];
          for (int k = 0; k < 3; k++){
          //  std::cout << "dulist"<< dulist_r[jju][k]<< std::endl;
          //  std::cout << "dulist_i"<< dulist_i[jju][k]<< std::endl;
          //  std::cout << "zlist_i"<< zlist_i[z_offset+jjz]<< std::endl;
            sumzdu_r[k] +=
                (dulist_r[jju][k] * zlist_r[z_offset + jjz] +
                 dulist_i[jju][k] * zlist_i[z_offset + jjz]) * 0.5;
          //  std::cout << "sumzdu_r"<< sumzdu_r[k]<< std::endl;
          }
                //(dudr_r[k] * zlist_r[z_offset + jjz] +
                // dudr_i[k] * zlist_i[z_offset + jjz]) * 0.5;
          // jjz++;
          // jju++;
        } // end if jeven
      //  std::cout << "end if j even" << std::endl;

        for (int k = 0; k < 3; k++){
        //  std::cout << sumzdu_r[0] << std::endl;
        //  std::cout << dblist[db_idx][0] << std::endl;
          //dbdr[k] += 2.0 * sumzdu_r[k];
          dblist[db_idx][k]+= 2.0 * sumzdu_r[k];
        }
        // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
        
      //  std::cout << "end db assignment" << std::endl;

        double j1fac = (j + 1) / (j1 + 1.0);

        idouble = elem1*nelements+elem2;
        itriple = (elem3*nelements+elem2)*nelements+elem1;
        //dbdr = dblist[itriple*idxb_max+jjb];
        db_idx = itriple*idxb_max+jjb;
        jjz = idxz_block[j][j2][j1];
        jju = idxu_block[j1];
        //zptr_r = &zlist_r[idouble*idxz_max];
        //zptr_i = &zlist_i[idouble*idxz_max];
        z_offset = idouble * idxz_max;
        for (int k = 0; k < 3; k++)
          sumzdu_r[k] = 0.0;

        for (int mb = 0; 2 * mb < j1; mb++)
          for (int ma = 0; ma <= j1; ma++) {
            //dudr_r = dulist_r[jju];
            //dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dulist_r[jju][k] * zlist_r[z_offset + jjz] +
                  dulist_i[jju][k] * zlist_i[z_offset + jjz];
            jjz++;
            jju++;
          } //end loop over ma mb

        // For j1 even, handle middle column

        if (j1 % 2 == 0) {
          int mb = j1 / 2;
          for (int ma = 0; ma < mb; ma++) {
            //dudr_r = dulist_r[jju];
            //dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dulist_r[jju][k] * zlist_r[z_offset+ jjz] +
                  dulist_i[jju][k] * zlist_i[z_offset + jjz];
            jjz++;
            jju++;
          }
          //dudr_r = dulist_r[jju];
          //dudr_i = dulist_i[jju];
          for (int k = 0; k < 3; k++)
            sumzdu_r[k] +=
                (dulist_r[jju][k] * zlist_r[z_offset + jjz] +
                 dulist_i[jju][k] * zlist_i[z_offset + jjz]) * 0.5;
          // jjz++;
          // jju++;
        } // end if j1even
      //  std::cout << "end if j even middle" << std::endl;

        for (int k = 0; k < 3; k++)
          if (bnorm_flag)
           //dbdr[k] += 2.0 * sumzdu_r[k];
            dblist[db_idx][k] += 2.0 * sumzdu_r[k];
          else
            //dbdr[k] += 2.0 * sumzdu_r[k] * j1fac;
            dblist[db_idx][k] += 2.0 * sumzdu_r[k] * j1fac;

        // Sum over Conj(dudr(j2,ma2,mb2))*z(j,j1,j2,ma2,mb2)

        double j2fac = (j + 1) / (j2 + 1.0);

        idouble = elem2*nelements+elem1;
        itriple = (elem1*nelements+elem3)*nelements+elem2;
        //dbdr = dblist[itriple*idxb_max+jjb];
        db_idx = itriple*idxb_max+jjb;
        jjz = idxz_block[j][j1][j2];
        jju = idxu_block[j2];
        //zptr_r = &zlist_r[idouble*idxz_max];
        //zptr_i = &zlist_i[idouble*idxz_max];
        z_offset = idouble*idxz_max;
        for (int k = 0; k < 3; k++)
          sumzdu_r[k] = 0.0;

        for (int mb = 0; 2 * mb < j2; mb++)
          for (int ma = 0; ma <= j2; ma++) {
            //dudr_r = dulist_r[jju];
            //dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dulist_r[jju][k] * zlist_r[z_offset + jjz] +
                  dulist_i[jju][k] * zlist_i[z_offset + jjz];
            jjz++;
            jju++;
          } //end loop over ma mb

        // For j2 even, handle middle column

        if (j2 % 2 == 0) {
          int mb = j2 / 2;
          for (int ma = 0; ma < mb; ma++) {
            //dudr_r = dulist_r[jju];
            //dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dulist_r[jju][k] * zlist_r[z_offset + jjz] +
                  dulist_i[jju][k] * zlist_i[z_offset + jjz];
            jjz++;
            jju++;
          }
          //dudr_r = dulist_r[jju];
          //dudr_i = dulist_i[jju];
          for (int k = 0; k < 3; k++)
            sumzdu_r[k] +=
                (dulist_r[jju][k] * zlist_r[z_offset + jjz] +
                 dulist_i[jju][k] * zlist_i[z_offset + jjz]) * 0.5;
          // jjz++;
          // jju++;
        } // end if j2even
      //  std::cout << "end if j2 even " << std::endl;

        for (int k = 0; k < 3; k++)
          if (bnorm_flag)
            //dbdr[k] += 2.0 * sumzdu_r[k];
            dblist[db_idx][k] += 2.0 * sumzdu_r[k];
          else
            //dbdr[k] += 2.0 * sumzdu_r[k] * j2fac;
            dblist[db_idx][k] += 2.0 * sumzdu_r[k] * j2fac;
        
       // std::cout << " ddblist "<< dblist[db_idx][0] << std::endl;
       // std::cout << " ddblist "<< dblist[db_idx][1] << std::endl;
       // std::cout << " ddblist "<< dblist[db_idx][2] << std::endl;
      //  std::cout << "end assign to ddblist " << std::endl;

      }
  } //end loop over j1 j2 j

}

  void SNADesc::compute_duidrj(std::vector<double> rij, double wj, double rcut, int jj, int jelem){
    std::cout<< "SNADesc::compute_duidrj entered function" <<std::endl;
    double rsq, r, x, y, z, z0, theta0, cs, sn;
    double dz0dr;

    x = rij[0];
    y = rij[1];
    z = rij[2];
   std::cout<< "SNADesc::compute_duidrj xyz is "<<x<<" "<<y<<" "<<z<<std::endl;
    rsq = x * x + y * y + z * z;
    r = std::sqrt(rsq);
    double rscale0 = rfac0 * mypi / (rcut - rmin0);
    theta0 = (r - rmin0) * rscale0;
    cs = std::cos(theta0);
    sn = std::sin(theta0);
    z0 = r * cs / sn;
    dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;
  std::cout<< "SNADesc::compute_duidrj okay up to compute_duarray"<<std::endl;

  elem_duarray = jelem;
  compute_duarray(x, y, z, z0, r, dz0dr, wj, rcut, jj);
}

    void SNADesc::init_rootpqarray(){
      for (int p = 1; p <= twojmax; p++)
        for (int q = 1; q <= twojmax; q++){
          rootpqarray[p][q] = std::sqrt(static_cast<double>(p)/q);
        }
    }

    void SNADesc::zero_uarraytot(int ielem){
      for (int jelem = 0; jelem < nelements; jelem++)
      for (int j = 0; j <= twojmax; j++) {
        int jju = idxu_block[j];
        for (int mb = 0; mb <= j; mb++) {
          for (int ma = 0; ma <= j; ma++) {
            ulisttot_r[jelem*idxu_max+jju] = 0.0;
            ulisttot_i[jelem*idxu_max+jju] = 0.0;

            if (jelem == ielem || wselfall_flag){
              //std::cout << "we are in this wselfall_flag area"<< std::endl;
              //std::cout << "wself is "<< wself << std::endl;
              if (ma==mb)
                ulisttot_r[jelem*idxu_max+jju] = wself; ///// double check this
            }
            jju++;
          }
        }
      }
      //std::cout << std::endl;
    }

    void SNADesc::add_uarraytot(double r, double wj, double rcut, int jj, int jelem){

       std::cout << "SNA::add_uarraytot" << std::endl;
      double sfac;

      sfac = compute_sfac(r, rcut);

      sfac *= wj;
       std::cout << "sfac is " << sfac << std::endl;
      // std::cout << "ulist_r "<< std::endl;
      for (int j = 0; j <= twojmax; j++) {
        int jju = idxu_block[j];
        for (int mb = 0; mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            // std::cout << " " << ulist_r[jju] << std::endl;
            // std::cout << jelem*idxu_max+jju << std::endl;
            ulisttot_r[jelem*idxu_max+jju] +=
              sfac * ulist_r_ij[jj][jju];
            ulisttot_i[jelem*idxu_max+jju] +=
              sfac * ulist_i_ij[jj][jju];
            jju++;
            // std::cout << " " << ulist_r[jju] << std::endl;
            // std::cout << "jj " << jj <<" jju" << jju<<std::endl;

          }
      }
    }
    void SNADesc::compute_uarray(double x, double y, double z, double z0, double r, int jj){
    //std::cout<< "SNADesc::compute_uarray" <<std::endl;
    double r0inv;
    double a_r, b_r, a_i, b_i;
    double rootpq;

    // compute Cayley-Klein parameters for unit quaternion

    r0inv = 1.0 / std::sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;
    // std::cout<< "weird ar things:" <<//std::endl;
    std::cout << "a_r " << a_r <<" a_i " << a_i <<" b_r "<< b_r <<" b_i"<<b_i << std::endl;
    // VMK Section 4.8.2



    ulist_r_ij[jj][0] = 1.0;
    ulist_i_ij[jj][0] = 0.0;

    for (int j = 1; j <= twojmax; j++) {

      int jju = idxu_block[j];
      // std::cout<< "jju is " << jju <<std::endl;
      int jjup = idxu_block[j-1];
      // std::cout<< "jjup is " << jjup <<std::endl;

      // fill in left side of matrix layer from previous layer

      for (int mb = 0; 2*mb <= j; mb++) {
         //std::cout << "jju is " << jju << " !!"<<std::endl;
         //std::cout << "jjup is " << jjup << " !!"<<std::endl;
        ulist_r_ij[jj][jju] = 0.0;
        ulist_i_ij[jj][jju] = 0.0;

        for (int ma = 0; ma < j; ma++) {
          rootpq = rootpqarray[j - ma][j - mb];
          ulist_r_ij[jj][jju] +=
            rootpq *
            (a_r * ulist_r_ij[jj][jjup] +
            a_i * ulist_i_ij[jj][jjup]);
          // std::cout << "first term is " << a_r * ulist_r_ij[jj][jjup] << std::endl;
          // std::cout << "second term is " << a_i * ulist_i_ij[jj][jjup] << std::endl;
          // std::cout << "rootpq is " << rootpq << std::endl;
          // std::cout << "ulist_r of jjup " << ulist_r[jjup] <<std::endl;
          // std::cout << "ulist_i of jjup " << ulist_i[jjup] <<std::endl<<std::endl;
           //std::cout << "ulist_rij at jj " << jj << " and jju " <<jju << " is " << ulist_r_ij[jj][jju] <<std::endl<<std::endl;
          ulist_i_ij[jj][jju] +=
            rootpq *
            (a_r * ulist_i_ij[jj][jjup] -
            a_i * ulist_r_ij[jj][jjup]);
             //std::cout << "first term is " << a_r * ulist_i_ij[jj][jjup] << std::endl;
            // std::cout << "second term is " << a_i * ulist_r_ij[jj][jjup] << std::endl;
            // std::cout << "ulist_i_ij at jj " << jj << " and jju " <<jju << " is " << ulist_i_ij[jj][jju] <<std::endl;
          rootpq = rootpqarray[ma + 1][j - mb];
          ulist_r_ij[jj][jju+1] =
            -rootpq *
            (b_r * ulist_r_ij[jj][jjup] +
            b_i * ulist_i_ij[jj][jjup]);
            // std::cout << "rootpq is " << rootpq << std::endl;
            // std::cout << "first term is " << b_r * ulist_r_ij[jj][jjup] << std::endl;
            // std::cout << "second term is " << b_i * ulist_i_ij[jj][jjup] << std::endl;
            // std::cout << "ulist_r_ij at jj " << jj << " and jju (+1) " << jju << " is " << ulist_r_ij[jj][jju+1] <<std::endl;
          ulist_i_ij[jj][jju+1] =
            -rootpq *
            (b_r * ulist_i_ij[jj][jjup] -
            b_i * ulist_r_ij[jj][jjup]);
             //std::cout << "first term is " << b_r * ulist_i_ij[jj][jjup] << std::endl;
            // std::cout << "second term is " << b_i * ulist_r_ij[jj][jjup] << std::endl;
            //std::cout << "ulist_r_rij at jj " << jj << " and jju (+1) " << jju << " is " << ulist_i_ij[jj][jju+1] <<std::endl;
          jju++;
          jjup++;
        }
        jju++;
      }

      // copy left side to right side with inversion symmetry VMK 4.4(2)
      // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

      jju = idxu_block[j];
      jjup = jju+(j+1)*(j+1)-1;
      // std::cout << "jjup is " <<jjup << std::endl;
      int mbpar = 1;
      for (int mb = 0; 2*mb <= j; mb++) {
        int mapar = mbpar;
        for (int ma = 0; ma <= j; ma++) {
          // std::cout << "for ma " << ma << "  mapar is " << mapar <<std::endl;
          if (mapar == 1) {
            ulist_r_ij[jj][jjup] = ulist_r_ij[jj][jju];
            ulist_i_ij[jj][jjup] = -ulist_i_ij[jj][jju];
          } else {
            ulist_r_ij[jj][jjup] = -ulist_r_ij[jj][jju];
            ulist_i_ij[jj][jjup] = ulist_i_ij[jj][jju];
          }
          mapar = -mapar;
          jju++;
          jjup--;
        }
        mbpar = -mbpar;
      }
    }
    std::cout << "ulist_i_ij" <<std::endl;
     for (int i = 0; i <nmax; i++){
       for (int j = 0;j< idxu_max;j++){
         std::cout << ulist_i_ij[i][j] << " ";
       }
       std::cout << std::endl <<std::endl;
     }
    
    std::cout << "ulist_r_ij" <<std::endl;
    for (int i = 0; i <nmax; i++){
      for (int j = 0;j< idxu_max;j++){
        std::cout << ulist_r_ij[i][j] << " ";
      }
      std::cout << std::endl<< std::endl;
    }
  }
    void SNADesc::compute_duarray(double x, double y, double z,
                          double z0, double r, double dz0dr,
                          double wj, double rcut, int jj)
 {
  std::cout<< "SNADesc::compute_duarray" <<std::endl;
  std::cout<< "input_args: " << "z0 " << z0 <<   std::endl;
  std::cout<< "input_args: " << "r " << r <<   std::endl;
  std::cout<< "input_args: " << "dz0dr " << dz0dr   <<std::endl;
  std::cout<< "input_args: " << "wj " << wj<<  std::endl;
  std::cout<< "input_args: " << "rcut " <<rcut<<  std::endl;
  std::cout<< "input_args: " << "jj " <<jj<<   std::endl;
  double r0inv;
  double a_r, a_i, b_r, b_i;
  double da_r[3], da_i[3], db_r[3], db_i[3];
  double dz0[3], dr0inv[3], dr0invdr;
  double rootpq;

  double rinv = 1.0 / r;
  double ux = x * rinv;
  double uy = y * rinv;
  double uz = z * rinv;

  r0inv = 1.0 / std::sqrt(r * r + z0 * z0);
  a_r = z0 * r0inv;
  a_i = -z * r0inv;
  b_r = y * r0inv;
  b_i = -x * r0inv;

  dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);

  dr0inv[0] = dr0invdr * ux;
  dr0inv[1] = dr0invdr * uy;
  dr0inv[2] = dr0invdr * uz;

  dz0[0] = dz0dr * ux;
  dz0[1] = dz0dr * uy;
  dz0[2] = dz0dr * uz;

  for (int k = 0; k < 3; k++) {
    da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
    da_i[k] = -z * dr0inv[k];
  }

  da_i[2] += -r0inv;

  for (int k = 0; k < 3; k++) {
    db_r[k] = y * dr0inv[k];
    db_i[k] = -x * dr0inv[k];
  }

  db_i[0] += -r0inv;
  db_r[1] += r0inv;

  dulist_r[0][0] = 0.0;
  dulist_r[0][1] = 0.0;
  dulist_r[0][2] = 0.0;
  dulist_i[0][0] = 0.0;
  dulist_i[0][1] = 0.0;
  dulist_i[0][2] = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];
    std::cout<< " ulist_r is " << ulist_r_ij[jj][jjup]<<std::endl;
    std::cout<< " ulist_i is " << ulist_i_ij[jj][jjup]<<std::endl;
    for (int mb = 0; 2*mb <= j; mb++) {
      dulist_r[jju][0] = 0.0;
      dulist_r[jju][1] = 0.0;
      dulist_r[jju][2] = 0.0;
      dulist_i[jju][0] = 0.0;
      dulist_i[jju][1] = 0.0;
      dulist_i[jju][2] = 0.0;

      for (int ma = 0; ma < j; ma++) {
        rootpq = rootpqarray[j - ma][j - mb];
        for (int k = 0; k < 3; k++) {
          std::cout<< " da_k is " << da_r[k]<<std::endl;
          dulist_r[jju][k] +=
            rootpq * (da_r[k] * ulist_r_ij[jj][jjup] +
                      da_i[k] * ulist_i_ij[jj][jjup] +
                      a_r * dulist_r[jjup][k] +
                      a_i * dulist_i[jjup][k]);
          dulist_i[jju][k] +=
            rootpq * (da_r[k] * ulist_i_ij[jj][jjup] -
                      da_i[k] * ulist_r_ij[jj][jjup] +
                      a_r * dulist_i[jjup][k] -
                      a_i * dulist_r[jjup][k]);
        }

        rootpq = rootpqarray[ma + 1][j - mb];
        for (int k = 0; k < 3; k++) {
          dulist_r[jju+1][k] =
            -rootpq * (db_r[k] * ulist_r_ij[jj][jjup] +
                       db_i[k] * ulist_i_ij[jj][jjup] +
                       b_r * dulist_r[jjup][k] +
                       b_i * dulist_i[jjup][k]);
          dulist_i[jju+1][k] =
            -rootpq * (db_r[k] * ulist_i_ij[jj][jjup] -
                       db_i[k] * ulist_r_ij[jj][jjup] +
                       b_r * dulist_i[jjup][k] -
                       b_i * dulist_r[jjup][k]);
        }
        jju++;
        jjup++;
      }
      jju++;
    }

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    jju = idxu_block[j];
    jjup = jju+(j+1)*(j+1)-1;
    int mbpar = 1;
    for (int mb = 0; 2*mb <= j; mb++) {
      int mapar = mbpar;
      for (int ma = 0; ma <= j; ma++) {
        if (mapar == 1) {
          for (int k = 0; k < 3; k++) {
            dulist_r[jjup][k] = dulist_r[jju][k];
            dulist_i[jjup][k] = -dulist_i[jju][k];
          }
        } else {
          for (int k = 0; k < 3; k++) {
            dulist_r[jjup][k] = -dulist_r[jju][k];
            dulist_i[jjup][k] = dulist_i[jju][k];
          }
        }
        mapar = -mapar;
        jju++;
        jjup--;
      }
      mbpar = -mbpar;
    }
  }

  double sfac = compute_sfac(r, rcut);
  double dsfac = compute_dsfac(r, rcut);

  sfac *= wj;
  dsfac *= wj;
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];
    for (int mb = 0; 2*mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        dulist_r[jju][0] = dsfac * ulist_r_ij[jj][jju] * ux +
                                  sfac * dulist_r[jju][0];
        dulist_i[jju][0] = dsfac * ulist_i_ij[jj][jju] * ux +
                                  sfac * dulist_i[jju][0];
        dulist_r[jju][1] = dsfac * ulist_r_ij[jj][jju] * uy +
                                  sfac * dulist_r[jju][1];
        dulist_i[jju][1] = dsfac * ulist_i_ij[jj][jju] * uy +
                                  sfac * dulist_i[jju][1];
        dulist_r[jju][2] = dsfac * ulist_r_ij[jj][jju] * uz +
                                  sfac * dulist_r[jju][2];
        dulist_i[jju][2] = dsfac * ulist_i_ij[jj][jju] * uz +
                                  sfac * dulist_i[jju][2];
        jju++;
      }
  }
}  
    
    void SNADesc::init_clebsch_gordan(){
  std::cout<< "SNADesc:: in clebsch gordan" <<std::endl;
  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  int idxcg_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if (m < 0 || m > j) {
              cglist[idxcg_count] = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;

            for (int z = std::max(0, std::max(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= std::min((j1 + j2 - j) / 2,
                          std::min((j1 - aa2) / 2, (j2 + bb2) / 2));
                 z++) {
                   // std::cout<< "z is " <<z<<std::endl;
              ifac = z % 2 ? -1 : 1;
              // std::cout<< "ifac is " <<ifac<<std::endl;
              sum += ifac /
                (factorial(z) *
                 factorial((j1 + j2 - j) / 2 - z) *
                 factorial((j1 - aa2) / 2 - z) *
                 factorial((j2 + bb2) / 2 - z) *
                 factorial((j - j2 + aa2) / 2 + z) *
                 factorial((j - j1 - bb2) / 2 + z));
            }

            cc2 = 2 * m - j;
            dcg = deltacg(j1, j2, j);
            sfaccg = std::sqrt(factorial((j1 + aa2) / 2) *
                          factorial((j1 - aa2) / 2) *
                          factorial((j2 + bb2) / 2) *
                          factorial((j2 - bb2) / 2) *
                          factorial((j  + cc2) / 2) *
                          factorial((j  - cc2) / 2) *
                          (j + 1));
            // std::cout<< "cg initialize parts"<<std::endl;
            // std::cout << sum << " "<< dcg <<" "<<sfaccg <<std::endl;
            cglist[idxcg_count] = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
// for (int i = 0 ; i <= idxcg_max; i++){
//   std::cout <<cglist[i] << " " ;
//   }
// std::cout <<std::endl;

}

double SNADesc::factorial(int n){
  return std::tgamma(n+1);
}


double SNADesc::deltacg(int j1, int j2, int j){
  double sfaccg = factorial((j1 + j2 + j) / 2 + 1);
  return std::sqrt(factorial((j1 + j2 - j) / 2) *
              factorial((j1 - j2 + j) / 2) *
              factorial((-j1 + j2 + j) / 2) / sfaccg);
}

 
}
