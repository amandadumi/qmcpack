#ifndef QMCPLUSPLUS_SNADESC
#define QMCPLUSPLUS_SNADESC

#include <vector>
#include <cmath>
namespace qmcplusplus
{
struct SNA_ZINDICES {
  int j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, jju;
};

struct SNA_BINDICES {
  int j1, j2, j;
};

class SNADesc
{
    public:
     SNADesc();
     ~SNADesc();
      void build_indexlist();
      void init(int npart, int input_twojmax, double input_rcut);
      int ncoeff;
      std::vector<double> blist;
      std::vector<std::vector<double>> dblist;
      std::vector<int> inside;
      std::vector<double> wj;
      std::vector<double>rcutij;
      std::vector<std::vector<double>> rij;
      std::vector<int> element;    // index on [0,nelements)
      std::vector<int> nmax;
      int twojmax;
      std::vector<double> ylist_r, ylist_i;
      int idxcg_max, idxu_max, idxz_max, idxb_max;


      void compute_ncoeff();
      void create_arrays();
      double compute_sfac(double r, double rcut);
      double compute_dsfac(double r, double rcut);
      void compute_ui(int jnum,int ielem);
      void compute_zi();
      void compute_yi(std::vector<double> beta);
      void compute_yterm();
      void compute_bi(int ielem);
      void compute_duidrj(std::vector<double>rij,
                          double wj, double rcut, 
                          int jj, int jelem);
      void compute_dbidrj();
  
      SNA_ZINDICES *idxz;
      SNA_BINDICES *idxb;

      std::vector<std::vector<double>> rootpqarray;
      std::vector<double> cglist;
      std::vector<std::vector<std::vector<int>>> idxcg_block;

      std::vector<double> ulisttot_r, ulisttot_i;
      std::vector<std::vector<double>> dulist_r, dulist_i;
      std::vector<std::vector<double>> ulist_r_ij, ulist_i_ij;
      std::vector<int> idxu_block;

     std::vector<double> zlist_r, zlist_i;
      
      std::vector<std::vector<std::vector<int>>> idxz_block;

      std::vector<std::vector<std::vector<int>>> idxb_block;
      
      void create_twojmax_arrays();
      void destroy_twojmax_arrays();
      void init_clebsch_gordan();
      void print_clebsch_gordan();
      void init_rootpqarray();
      void zero_uarraytot(int);
      void add_uarraytot(double, double, double, int, int);
      void compute_uarray(double, double, double, double, double, int);
      void compute_duarray(double x, double y, double z,
                          double z0, double r, double dz0dr,
                          double wj, double rcut, int jj);
      double deltacg(int, int, int);
      double factorial(int n);

  
      double wself;

      int bzero_flag;       // 1 if bzero subtracted from barray
      int elem_duarray;       // 1 if bzero subtracted from barray
      std::vector<double> bzero;        // array of B values for isolated atoms
      int bnorm_flag=0;       // 1 if barray divided by j+1
      int chem_flag=0;        // 1 for multi-element bispectrum components
      int wselfall_flag=0;    // 1 for adding wself to all element labelings
      int nelements;        // number of elements
      int ndoubles;         // number of multi-element pairs
      int ntriples;         // number of multi-element triplets
      double mypi = 3.14159;
      double rfac0 = .99; // rfac0_in;
      double rmin0 = 0.0; //rmin0_in;
      int switch_flag = 1;//switch_flag_in;
};
}
#endif
