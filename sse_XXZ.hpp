#ifndef sse_XXZ_hpp
#define sse_XXZ_hpp

#include <iostream>
#include <iomanip>
#include <math.h>
#include <random>
#include <functional>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>

/// Stochastic Series Expansion for the Spin 1/2 Hamiltonian
/// H = \sum_{<i,j>}[ J_{i,j}/2 * ( S_i^+ S_j^- + S_i^- S_j^+) + V_{i,j} S_i^Z S_j^Z ]- \sum_{i} h_i S_i^Z
/// with arbitrary J_{i,j} and h_i

class sse_XXZ
{
  double Pi;
  enum class PATH {switch_reverse = 0, switch_continue = 1, straight = 2, bounce = 3};
  
  std::mt19937 MyGenerator;
  std::uniform_real_distribution<> uni_dist;
  std::function<double()> rnd;
  
  int Lx;
  int Ly;
  int Lz;
  int BCx;
  int BCy;
  int BCz;
  double T;
  double Beta;
  int Ntherm;
  int Nbins;
  int Nsweeps;
  double Jx;
  double Jy;
  double Delta;
  double V;
  double epsilon;
  double h;
  
  std::vector<int> spins;
  std::vector<int> spins_vertex;
  std::vector<int> spins_vertex_saved;
  std::vector<int> frstv;
  std::vector<int> lastv;
  std::vector<std::vector<int>> bonds;
  std::vector<double> t_bonds;
  std::vector<double> V_bonds;
  std::vector<double> c_bonds;
  std::vector<double> h_sites;
  std::vector<double> W_bonds;
  std::vector<std::vector<std::vector<double>>> path_prob;
  std::vector<std::vector<std::vector<std::vector<double>>>> path_prob_heatBath;
  std::vector<std::vector<std::vector<std::vector<double>>>> path_prob_nosymmetry;
  
  int Nsites;
  int Nbonds;
  int extra_bonds;
  
  int dim;
  double p0;
  int center_site;
  
  int Cutoff;
  std::vector<int> ops;
  std::vector<int> ops_saved;
  std::vector<int> vl;
  std::vector<int> vl_flipped;
  
  double constant_shift;
  
  std::vector<int> visited_vertex;
  int number_loops_per_sweep;
  
  // observables
  
  std::string file_extension;
  
  int NHam;
  double MC_NHam;
  std::vector<double> NHam_bins;
  
  std::vector<double> NHam_diagonal;
  std::vector<double> NHam_offdiagonal;
  double NHam_offdiag_x_plus;
  double NHam_offdiag_x_minus;
  double NHam_offdiag_y_plus;
  double NHam_offdiag_y_minus;
  
  double MC_energy;
  std::vector<double> energy_bins;
  
  double MC_sf_density;
  std::vector<double> sf_density_bins;
  
  double MC_wx;
  double MC_wy;
  double MC_wx2;
  double MC_wy2;
  std::vector<double> wx_bins;
  std::vector<double> wy_bins;
  std::vector<double> wx2_bins;
  std::vector<double> wy2_bins;
  
  double MC_mag;
  double MC_mag2;
  std::vector<double> mag_bins;
  std::vector<double> mag2_bins;
  
  std::vector<double> MC_energy_diag;
  std::vector<std::vector<double>> energy_diag_bins;
  
  std::vector<double> MC_energy_off;
  std::vector<std::vector<double>> energy_off_bins;
  
  std::vector<double> Sz_average;
  std::vector<double> Sz_average_beta;
  
public:
  
  sse_XXZ(std::map<std::string,double> params, std::string file_ext):
  Pi(std::acos(0)*2.),
  MyGenerator(10),
  uni_dist(0,1),
  rnd(bind(uni_dist,MyGenerator)),
  Lx(static_cast<int>(params["Lx"])),
  Ly(static_cast<int>(params["Ly"])),
  Lz(static_cast<int>(params["Lz"])),
  BCx(static_cast<int>(params["BCx"])),
  BCy(static_cast<int>(params["BCy"])),
  BCz(static_cast<int>(params["BCz"])),
  T(params["T"]),
  Beta(1./T),
  Ntherm(static_cast<int>(params["Ntherm"])),
  Nbins(static_cast<int>(params["Nbins"])),
  Nsweeps(static_cast<int>(params["Nsweeps"])),
  Jx(params["Jx"]),
  Jy(params["Jy"]),
  Delta(params["Delta"]),
  epsilon(params["epsilon"]),
  number_loops_per_sweep(10),
  file_extension(file_ext),
  NHam(0),
  MC_NHam(0),
  NHam_bins(std::vector<double>(Nbins,0)),
  MC_energy(0),
  energy_bins(std::vector<double>(Nbins,0)),
  MC_sf_density(0),
  sf_density_bins(std::vector<double>(Nbins,0)),
  MC_wx(0),
  MC_wy(0),
  MC_wx2(0),
  MC_wy2(0),
  wx_bins(std::vector<double>(Nbins,0)),
  wy_bins(std::vector<double>(Nbins,0)),
  wx2_bins(std::vector<double>(Nbins,0)),
  wy2_bins(std::vector<double>(Nbins,0)),
  MC_mag(0),
  MC_mag2(0),
  mag_bins(std::vector<double>(Nbins,0)),
  mag2_bins(std::vector<double>(Nbins,0))
  {
    if( BCy == 0 && BCz == 0)
      dim = 1;
    else if (BCz == 0)
      dim = 2;
    else
      dim = 3;
    
    std::vector<std::vector<double>> J_bonds;
    std::vector<std::vector<double>> Delta_bonds;
    
    if (dim ==1)
    {
      std::vector<double> J_loc = {Jx};
      J_bonds.resize(Lx*Ly*Lz,J_loc);
      Delta_bonds.resize(Lx*Ly*Lz,std::vector<double>(1,Delta));
    }
    else if(dim == 2)
    {
      std::vector<double> J_loc = {Jx,Jy};
      J_bonds.resize(Lx*Ly*Lz,J_loc);
      h_sites.resize(Lx*Ly*Lz,0);
      
      std::ifstream in_t("./qc_obs/t_hop_b");
      std::ifstream in_h("./qc_obs/h_sites");
      
      for(size_t i=0; i<J_bonds.size(); i++)
      {
        in_t >> J_loc[0] >> J_loc[1];
        J_bonds[i] = J_loc;
        in_h >> h_sites[i];
      }
      in_t.close();
      in_h.close();
      Delta_bonds.resize(Lx*Ly,std::vector<double>(2,Delta));
      for(size_t i=0; i<Delta_bonds.size(); i++)
        Delta_bonds[i][0] = Delta;
    }
    
    for(int z=0; z<Lz; z++)
    {
      for(int y=0; y<Ly; y++)
      {
        for(int x=0; x<Lx; x++)
        {
          spins.push_back(rnd()<0.5 ? -1 : 1);
          frstv.push_back(-1);
          lastv.push_back(-1);
          int s0 = x + y*Lx + z*Lx*Ly;
          if( (x<Lx-1) || (BCx==1) )
          {
            int s1 = (x+1)%Lx + y*Lx + z*Lx*Ly;
            bonds.push_back({s0,s1});
            t_bonds.push_back(J_bonds[s0][0]/2.);
            V_bonds.push_back(J_bonds[s0][0] * Delta_bonds[s0][0]);
          }
          if ( (y<Ly-1) || (BCy==1) )
          {
            int s1 = x + ((y+1)%Ly)*Lx + z*Lx*Ly;
            bonds.push_back({s0,s1});
            t_bonds.push_back(J_bonds[s0][1]/2.);
            V_bonds.push_back(J_bonds[s0][1] * Delta_bonds[s0][1]);
          }
          if ( (z<Lz-1) || (BCz==1) )
          {
            int s1 = x + y*Lx+((z+1)%Lz)*Lx*Ly;
            bonds.push_back({s0,s1});
          }
        }
      }
    }
    
    Nbonds = static_cast<int>(bonds.size());
    Nsites = static_cast<int>(spins.size());
    // get a exactly half filled initial state:
    int sum_m = std::accumulate(spins.begin(),spins.end(),0);
    double sum_h = std::accumulate(h_sites.begin(),h_sites.end(),0.);
    int number_of_flips = std::abs(sum_m/2);
    
    int counter = 0;
    while( std::accumulate(spins.begin(),spins.end(),0) != 0)
    {
      for(int i=0; i<Nsites; i++)
      {
        if( spins[i] == static_cast<int>(sum_m/std::abs(sum_m)) && rnd()<0.5)
        {
          spins[i] = -1 * spins[i];
          counter ++;
        }
        if(counter == number_of_flips)
          break;
      }
    }

    c_bonds.resize(Nbonds,0);
    path_prob.resize(Nbonds,std::vector<std::vector<double>>(3,std::vector<double>(4,0)));
    std::vector<std::vector<std::vector<double>>> path_prob_bond(6,std::vector<std::vector<double>>(4,std::vector<double>(4,0)));
    path_prob_heatBath.resize(Nbonds,path_prob_bond);
    path_prob_nosymmetry.resize(Nbonds,path_prob_bond);
    double number_n_n = 2*dim;
    for(int i=0; i<Nbonds; i++)
    {
      double W1 = t_bonds[i]; // off diagonal
      double W2 = V_bonds[i]/4. + h_sites[bonds[i][0]]/2./number_n_n - h_sites[bonds[i][1]]/2./number_n_n;  // diagonal unequal spins - up/down
      double W3 = V_bonds[i]/4. - h_sites[bonds[i][0]]/2./number_n_n + h_sites[bonds[i][1]]/2./number_n_n;  // diagonal unequal spins - down/up
      double W4 = -V_bonds[i]/4. + h_sites[bonds[i][0]]/2./number_n_n + h_sites[bonds[i][1]]/2./number_n_n;  // diagonal equal spins up
      double W5 = -V_bonds[i]/4. - h_sites[bonds[i][0]]/2./number_n_n - h_sites[bonds[i][1]]/2./number_n_n;  // diagonal equal spins down
      
      double weights[] = {W2,W3,W4,W5};
      c_bonds[i] = std::fabs(*std::min_element(weights,weights+4)) + epsilon;
      
      W_bonds.push_back( W1 ); W_bonds.push_back( W1 ); W_bonds.push_back( W2+c_bonds[i] ), W_bonds.push_back( W3+c_bonds[i] ); W_bonds.push_back( W4+c_bonds[i] ); W_bonds.push_back( W5+c_bonds[i] );
      
      std::vector<double>  W_b = {W1,W1,W2+c_bonds[i],W3+c_bonds[i],W4+c_bonds[i],W5+c_bonds[i]};
      set_path_prob_HeatBath(i,W_b);
      set_path_prob_dirloop_nosymmetry(i,W_b);
    }

    Cutoff = static_cast<int>(Nsites*Beta); // initial length of the operator string
    for(int i=0; i<Cutoff; i++)
    {
      spins_vertex.push_back(0); spins_vertex.push_back(0); spins_vertex.push_back(0); spins_vertex.push_back(0);
      spins_vertex_saved.push_back(0); spins_vertex_saved.push_back(0); spins_vertex_saved.push_back(0); spins_vertex_saved.push_back(0);
      ops.push_back(0);
      ops_saved.push_back(0);
      vl.push_back(-1); vl.push_back(-1); vl.push_back(-1); vl.push_back(-1);
    }
    vl_flipped.resize(vl.size(),0);
    constant_shift = 0;
    for(int i=0; i<Nbonds; i++)
      constant_shift +=  c_bonds[i];
    visited_vertex.resize(Cutoff);
    
    NHam_diagonal.resize(Nbonds,0);
    NHam_offdiagonal.resize(Nbonds,0);
    MC_energy_diag.resize(Nbonds,0);
    energy_diag_bins.resize(Nbins,MC_energy_diag);
    MC_energy_off.resize(Nbonds,0);
    energy_off_bins.resize(Nbins,MC_energy_diag);
    
    Sz_average.resize(Nsites,0);
    Sz_average_beta.resize(Nsites,0);
  }
  
  void adjust_cutoff(int step);
  
  void diagonal_update();
  
  void operator_loop_update_therm();
  
  void operator_loop_update();
  
  void build_vertex_list();
  
  double get_vertex_weight(int bond);
  
  int get_vertex_type(int vertex);
  
  int get_exit_leg_HeatBath(int vertex_type, int vertex, int entrance_leg, int bond);
  
  int get_exit_leg_dirloop_nosymmetry(int vertex_type, int vertex, int entrance_leg, int bond);
  
  void estimator();
  
  void get_Sz_average();
  
  void set_bin(int ibin);
  
  void set_path_prob_dirloop_nosymmetry(int bond, const std::vector<double>& W);
  
  void set_dirloop(int bond,const std::vector<double>& W_tilde, const std::vector<int>& index_in, const std::vector<int>& entrance, const std::vector<size_t>& idx);
  
  int get_exit(const std::vector<int>& trans);
  
  int get_path(const std::vector<size_t>& trans);
  
  void set_path_prob_HeatBath(int bond, const std::vector<double>& W);
  
  template<typename T>
  inline void print_vector(std::vector<T> v)
  {
    for(size_t i=0; i<v.size(); i++)
      std::cout << v[i] << " ";
    std::cout << std::endl;
  }
  
  template<typename T>
  inline void print_matrix(std::vector<std::vector<T>> v)
  {
    std::cout << std::endl;
    for(size_t i=0; i<v.size(); i++)
    {
      for(size_t j=0; j<v[0].size(); j++)
      {
        std::cout << v[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  
  template <typename T>
  inline std::vector<size_t> argsort(const std::vector<T>& v)
  {
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i)
      idx[i] = i;
    
    std::sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    
    return idx;
  }
  
};

#endif
