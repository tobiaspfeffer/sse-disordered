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
      
      std::vector<double> W_tilde = {W1,W2+c_bonds[i],W4+c_bonds[i]};
      set_path_prob_dirloop(i,W_tilde);
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
  
  void adjust_cutoff(int step)
  {
    int newCutoff = NHam + static_cast<int>(static_cast<double>(NHam)/3.);
    if(newCutoff > Cutoff)
    {
      while(static_cast<int>(ops.size())<newCutoff)
      {
        ops.push_back(0);
        ops_saved.push_back(0);
        vl.push_back(-1); vl.push_back(-1); vl.push_back(-1); vl.push_back(-1);
        vl_flipped.push_back(-1); vl_flipped.push_back(-1); vl_flipped.push_back(-1); vl_flipped.push_back(-1);
        spins_vertex.push_back(0); spins_vertex.push_back(0); spins_vertex.push_back(0); spins_vertex.push_back(0);
        spins_vertex_saved.push_back(0); spins_vertex_saved.push_back(0); spins_vertex_saved.push_back(0); spins_vertex_saved.push_back(0);
        visited_vertex.push_back(0);
      }
      Cutoff = static_cast<int>(ops.size());
      std::cout << "Adjusted cutoff to " << Cutoff << " after step " << step << std::endl;
    }
  }
  
  void diagonal_update()
  {
    std::fill(NHam_diagonal.begin(), NHam_diagonal.end(), 0);
    std::fill(NHam_offdiagonal.begin(), NHam_offdiagonal.end(), 0);
    NHam_offdiag_x_minus = 0; NHam_offdiag_x_plus = 0; NHam_offdiag_y_minus = 0; NHam_offdiag_y_plus = 0;
    for(size_t i=0; i<ops.size(); i++)
    {
      if(ops[i] == 0) // there is an identity
      {
        int b = static_cast<int>(rnd() * Nbonds);
        if( (rnd() * (Cutoff-NHam)) <= (get_vertex_weight(b) * Beta * Nbonds) ) // switching from the identity to a diagonal operator
        {
          NHam_diagonal[b] += 1;
          ops[i] = b+1;
          ops_saved[i] = b+1;
          NHam = NHam+1;
        }
      }
      else if(ops[i]>0) // there is a diagonal operator
      {
        NHam_diagonal[ops[i]-1] += 1;
        if( (rnd() * get_vertex_weight(ops[i]-1) * Beta * Nbonds) <= (Cutoff-NHam+1) )
        {
          NHam_diagonal[ops[i]-1] -= 1;
          ops[i] = 0;
          ops_saved[i] = 0;
          NHam = NHam-1;
        }
      }
      else // there is an offdiagonal operator
      {
        int b = -ops[i]-1; // map back to bond
        NHam_offdiagonal[b] += 1;
        if(b%2 == 0) // bond is in x-direction
        {
          if(spins[bonds[b][0]] == -1)
            NHam_offdiag_x_plus +=1;
          else
            NHam_offdiag_x_minus +=1;
        }
        else // bond is in y-direction
        {
          if(spins[bonds[b][0]] == -1)
            NHam_offdiag_y_plus +=1;
          else
            NHam_offdiag_y_minus +=1;
        }
        spins[bonds[b][0]] = -spins[bonds[b][0]];
        spins[bonds[b][1]] = -spins[bonds[b][1]];
      }
    }
  }
  
  void operator_loop_update_therm()
  {
    std::fill(vl_flipped.begin(), vl_flipped.end(), 0);
    std::fill(visited_vertex.begin(), visited_vertex.end(), 0);
    bool update;
    for(int l=0; l<number_loops_per_sweep ;l++)
    {
      int vertex=0;
      if(NHam != 0)
      {
        while(true)
        {
          vertex = static_cast<int> (rnd() * static_cast<double>(ops.size()) );
          if(ops[vertex] != 0 )
            break;
        }
      }
      int start_leg = 4*vertex + static_cast<int>(rnd()*4);
      int entrance_leg = start_leg;
      int loop_length = 0;
      int number_bounces = 0;
      while(true)
      {
        int vertex_type = get_vertex_type(vertex);
        int bond = 0;
        if( ops[vertex] > 0 ) // diagonal vertex
          bond = ops[vertex]-1;
        else // off_diagonal vertex
          bond = -ops[vertex]-1;

        int exit_leg = get_exit_leg_dirloop_nosymmetry(vertex_type,vertex,entrance_leg,bond);

        if(exit_leg == entrance_leg)
          number_bounces++;
        else
          visited_vertex[vertex] += 1;
        if(entrance_leg == exit_leg && exit_leg == start_leg)
        {
          update = true;
          break;
        }
        spins_vertex[entrance_leg] = -spins_vertex[entrance_leg];
        spins_vertex[exit_leg] = -spins_vertex[exit_leg];
        vertex_type = get_vertex_type(vertex); // after flipping the spins determine the vertex type
        if( vertex_type == 0 || vertex_type == 1 ) // after the flip the vertex is off diagonal
          ops[vertex] = -abs(ops[vertex]);
        else              // after the flip the vertex is diagonal
          ops[vertex] = abs(ops[vertex]);
        int next_leg = vl[exit_leg];

        vl_flipped[entrance_leg] -= 1;
        vl_flipped[exit_leg] -= 1;
        if(next_leg == start_leg)
        {
          update = true;
          break;
        }
        else
        {
          loop_length += 1;
          entrance_leg = next_leg;
          vertex = static_cast<int>(entrance_leg/4.);
        }
      }
    }
    int number_visited = 0;
    for(int i=0;i<NHam;i++)
      number_visited += visited_vertex[i];
    if (number_visited < NHam)
    {
      number_loops_per_sweep += number_loops_per_sweep/3;
      std::cout << "Adjusted number of operator loops in MC sweep to " << number_loops_per_sweep << std::endl;
    }
    if(update == true)
    {
      for(int i=0; i<Nsites; i++) // map back to base spin state
      {
        if(frstv[i] != -1)
        {
          if( abs(vl_flipped[frstv[i]])%2 != 0)
          {
            spins[i] = -spins[i];
            if( spins[i] != spins_vertex[frstv[i]])
            {
              std::cout << "NOT THE CORRECT SPIN" << std::endl;
              exit(1);
            }
          }
        }
        else // try to flip dangling spin
        {
          if(rnd()<0.5)
          {
            spins[i] = -spins[i];
          }
        }
      }
      ops_saved = ops;
      spins_vertex_saved = spins_vertex;
    }
    else
    {
      ops = ops_saved;
      spins_vertex = spins_vertex_saved;
    }
    
  }
  
  void operator_loop_update()
  {
    std::fill(vl_flipped.begin(), vl_flipped.end(), 0);
    bool update;
    for(int l=0; l<number_loops_per_sweep ;l++)
    {
      int vertex=0;
      if(NHam != 0)
        while(true)
        {
          vertex = static_cast<int> (rnd() * static_cast<double>(ops.size()) );
          if(ops[vertex] != 0 )
            break;
        }
      int start_leg = 4*vertex + static_cast<int>(rnd()*4);
      int entrance_leg = start_leg;
      int loop_length = 0;
      while(true)
      {
        int vertex_type = get_vertex_type(vertex);
        int bond = 0;
        if( ops[vertex] > 0 ) // diagonal vertex
          bond = ops[vertex]-1;
        else // off_diagonal vertex
          bond = -ops[vertex]-1;

        int exit_leg = get_exit_leg_dirloop_nosymmetry(vertex_type,vertex,entrance_leg,bond);
        
        if(entrance_leg == exit_leg && exit_leg == start_leg)
        {
          update = true;
          break;
        }
        spins_vertex[entrance_leg] = -spins_vertex[entrance_leg];
        spins_vertex[exit_leg] = -spins_vertex[exit_leg];
        vertex_type = get_vertex_type(vertex); // after flipping the spins determine the vertex type
        if( vertex_type == 0 || vertex_type == 1 ) // after the flip the vertex is off diagonal
          ops[vertex] = -abs(ops[vertex]);
        else              // after the flip the vertex is diagonal
          ops[vertex] = abs(ops[vertex]);
        int next_leg = vl[exit_leg];
        vl_flipped[entrance_leg] -= 1;
        vl_flipped[exit_leg] -= 1;
        if(next_leg == start_leg)
        {
          update = true;
          break;
        }
        else
        {
          loop_length += 1;
          entrance_leg = next_leg;
          vertex = static_cast<int>(entrance_leg/4.);
        }
      }
    }
    if(update == true)
    {
      for(int i=0; i<Nsites; i++) // map back to base spin state
      {
        if(frstv[i] != -1)
        {
          if( abs(vl_flipped[frstv[i]])%2 != 0)
          {
            spins[i] = -spins[i];
            if( spins[i] != spins_vertex[frstv[i]])
            {
              std::cout << "NOT THE CORRECT SPIN" << std::endl;
              exit(1);
            }
          }
        }
        else // try to flip dangling spin
        {
          if(rnd()<0.5)
          {
            spins[i] = -spins[i];
          }
        }
      }
      ops_saved = ops;
      spins_vertex_saved = spins_vertex;
    }
    else
    {
      ops = ops_saved;
      spins_vertex = spins_vertex_saved;
    }
    
  }
  
  void build_vertex_list()
  {
    for(int i=0; i<Nsites; i++)
    {
      frstv[i]=-1;
      lastv[i]=-1;
    }
    for(int vs=0; vs<static_cast<int>( vl.size()); vs+=4)
    {
      int op = ops[vs/4];
      if(op != 0)  // no identity
      {
        int b = abs(op)-1;
        int s0 = bonds[b][0];
        int s1 = bonds[b][1];
        int v0 = lastv[s0];
        int v1 = lastv[s1];
        if (v0 != -1)
        {
          vl[v0] = vs;
          vl[vs] = v0;
          spins_vertex[vs] = spins_vertex[v0];
          spins_vertex_saved[vs] = spins_vertex_saved[v0];
        }
        else
        {
          frstv[s0] = vs;
          spins_vertex[vs] = spins[s0];
          spins_vertex_saved[vs] = spins[s0];
        }
        if (v1 != -1)
        {
          vl[v1] = vs+1;
          vl[vs+1] = v1;
          spins_vertex[vs+1] = spins_vertex[v1];
          spins_vertex_saved[vs+1] = spins_vertex_saved[v1];

        }
        else
        {
          frstv[s1] = vs+1;
          spins_vertex[vs+1] = spins[s1];
          spins_vertex_saved[vs+1] = spins[s1];
        }
        lastv[s0] = vs+2;
        lastv[s1] = vs+3;
        if (op > 0) // there is a diagonal vertex
        {
          spins_vertex[vs+2] = spins_vertex[vs];
          spins_vertex[vs+3] = spins_vertex[vs+1];
          spins_vertex_saved[vs+2] = spins_vertex_saved[vs];
          spins_vertex_saved[vs+3] = spins_vertex_saved[vs+1];
        }
        else // there is a off diagonal vertex
        {
          spins_vertex[vs+2] = -spins_vertex[vs];
          spins_vertex[vs+3] = -spins_vertex[vs+1];
          spins_vertex_saved[vs+2] = -spins_vertex_saved[vs];
          spins_vertex_saved[vs+3] = -spins_vertex_saved[vs+1];
        }
      }
      else
      {
        for(int i=vs; i<vs+4; i++)
        {
          vl[i] = -1;
          spins_vertex[i] = 0;
          spins_vertex_saved[i] = 0;
        }
      }
    }
    for(int i=0; i<Nsites; i++)
    {
      int vfrst = frstv[i];
      if( vfrst != -1)
      {
        int vlast = lastv[i];
        vl[vlast] = vfrst;
        vl[vfrst] = vlast;
        if( spins_vertex[vfrst] != spins_vertex[vlast] )
        {
          std::cout << "There is sth. wrong in constructing the vertex list" << std::endl;
          exit(1);
        }
        
      }
    }
  }
  
  double get_vertex_weight(int bond)
  {
    if(spins[bonds[bond][0]] == spins[bonds[bond][1]]) // equal spins
    {
      if(spins[bonds[bond][0]] == 1) // up/up
        return W_bonds[6*bond+4];
      else
        return W_bonds[6*bond+5];  // down/down
    }
    else            // unequal spins
    {
      if(spins[bonds[bond][0]] == 1 ) // up/down
        return W_bonds[6*bond+2];
      else
        return W_bonds[6*bond+3];  // down/up
    }
  }
  
  int get_vertex_type(int vertex)
  {
    if( spins_vertex[4*vertex] != spins_vertex[4*vertex+2] )
    {
      if(spins_vertex[4*vertex] == 1)
        return 0; // off diag vertex - up/down
      else
        return 1; // off diag vertex - down/up
    }
    else if( spins_vertex[4*vertex] != spins_vertex[4*vertex+1] ) // diag vertex with unequal spins
    {
      if ( spins_vertex[4*vertex] == 1)  // up/down
        return 2;
      else  // down/up
        return 3;
    }
    else if( spins_vertex[4*vertex] == spins_vertex[4*vertex+1] ) // diag vertex with equal spins
    {
      if ( spins_vertex[4*vertex] == 1)  // up/up
        return 4;
      else  // down/down
        return 5;
    }
    else
    {
      std::cout << "Can't name vertex type " << vertex << std::endl;
      exit(1);
    }
    
  }
  
  int get_exit_leg_HeatBath(int vertex_type, int vertex, int entrance_leg, int bond)
  {
    int leg_number = entrance_leg > 3 ? entrance_leg%(4*vertex)  : entrance_leg;
    double r =rnd();
    double r0 = path_prob_heatBath[bond][vertex_type][leg_number][0];
    double r1 = path_prob_heatBath[bond][vertex_type][leg_number][1];
    double r2 = path_prob_heatBath[bond][vertex_type][leg_number][2];
    double r3 = path_prob_heatBath[bond][vertex_type][leg_number][3];
    if (std::fabs( r0+r1+r2+r3 - 1 ) > 1e-12)
    {
      std::cout << "The probability must sum up to 1 " << r0+r1+r2+r3 << " " << r0 << " " << r1 << " " << r2 << " " << r3  << std::endl;
      exit(1);
    }
    if (r < r0)
      return 4*vertex+0;
    else if (r < r0+r1)
      return 4*vertex+1;
    else if (r < r0+r1+r2)
      return 4*vertex+2;
    else
      return 4*vertex+3;
  }
  
  int get_exit_leg_dirloop(int vertex_type, int entrance_leg, int bond)
  {
    double r = rnd();
    double weight;
    if (vertex_type == 0 || vertex_type == 1)
    {
      weight = W_bonds[6*bond+vertex_type];
      vertex_type = 0;
    }
    if(vertex_type == 2 || vertex_type == 3)
    {
      weight = W_bonds[6*bond+vertex_type];
      vertex_type = 1;
    }
    else if(vertex_type == 4 || vertex_type == 5)
    {
      weight = W_bonds[6*bond+vertex_type];
      vertex_type = 2;
    }
    double r1 = path_prob[bond][vertex_type][static_cast<int>(PATH::switch_reverse)]/weight;
    double r2 = path_prob[bond][vertex_type][static_cast<int>(PATH::switch_continue)]/weight;
    double r3 = path_prob[bond][vertex_type][static_cast<int>(PATH::straight)]/weight;
    
    if(std::fabs(r1+r2+r3-1)>1e-12)
    {
      std::cout << "Sum must be 1 " << r1+r2+r3 << std::endl;
      std::cout << "vertex_type " << vertex_type << " " << weight << std::endl;
      exit(1);
    }
    
    if( r < r1 )
      return entrance_leg^1; // switch and reverse
    else if( r < r1+r2 )
    {
      return entrance_leg^3; // switch and continue
    }
    else if( r < r1+r2+r3 )
    {
      return entrance_leg^2; // straight
    }
    else
    {
      return entrance_leg;   //bounce
    }
  }
  
  int get_exit_leg_dirloop_nosymmetry(int vertex_type, int vertex, int entrance_leg, int bond)
  {
    double r = rnd();
    double weight = W_bonds[6*bond+vertex_type];
    int entrance = entrance_leg>3 ? entrance_leg%(4*vertex) : entrance_leg;
    double r1 = path_prob_nosymmetry[bond][vertex_type][entrance][static_cast<int>(PATH::switch_reverse)]/weight;
    double r2 = path_prob_nosymmetry[bond][vertex_type][entrance][static_cast<int>(PATH::switch_continue)]/weight;
    double r3 = path_prob_nosymmetry[bond][vertex_type][entrance][static_cast<int>(PATH::straight)]/weight;
    double r4 = path_prob_nosymmetry[bond][vertex_type][entrance][static_cast<int>(PATH::bounce)]/weight;
    
    if(std::fabs(r1+r2+r3+r4-1)>1e-12)
    {
      std::cout << " sum must be 1 " << r1 << " " << r2 << " " << r3 << " " << r4 << std::endl;
      std::cout << "vertex_type " << vertex_type << " " << weight << std::endl;
      exit(1);
    }
    
    if( r < r1 )
      return entrance_leg^1; // switch and reverse
    else if( r < r1+r2 )
    {
      return entrance_leg^3; // switch and continue
    }
    else if( r < r1+r2+r3 )
    {
      return entrance_leg^2; // straight
    }
    else
    {
      return entrance_leg;   //bounce
    }
  }
  
  void estimator()
  {
    MC_NHam += NHam;
    MC_energy += (-NHam/Beta+constant_shift)/Nsites; // energy density
    double W_x = (NHam_offdiag_x_plus - NHam_offdiag_x_minus);
    double W_y = (NHam_offdiag_y_plus - NHam_offdiag_y_minus);
    MC_wx += W_x/static_cast<double>(Lx);
    MC_wy += W_y/static_cast<double>(Ly);
    MC_wx2 += W_x/static_cast<double>(Lx)*W_x/static_cast<double>(Lx);
    MC_wy2 += W_y/static_cast<double>(Ly)*W_y/static_cast<double>(Ly);
    MC_sf_density += 1./static_cast<double>(Lx*Ly)*0.5*(W_x*W_x+W_y*W_y)/Beta; // this is in fact the sf-stiffness and 1./static_cast<double>(Lx*Ly) is due to the fact that W_x, W_y are not intensive in the above measurement
    double m = 0;
    for(int i=0; i<Nsites; i++)
      m += (0.5*spins[i]+0.5);
    MC_mag += m;
    MC_mag2 += m*m;
  }
  
  void get_Sz_average()
  {
    for(int i=0; i<Nsites; i++)
      Sz_average[i] = 0.5*spins[i];
      
  }
  
  void set_bin(int ibin)
  {
    NHam_bins[ibin] = MC_NHam/static_cast<double>(Nsweeps);
    MC_NHam = 0;
    energy_bins[ibin] = MC_energy/static_cast<double>(Nsweeps);
    MC_energy = 0;
    sf_density_bins[ibin] = MC_sf_density/static_cast<double>(Nsweeps);
    MC_sf_density = 0;
    wx_bins[ibin] = MC_wx/static_cast<double>(Nsweeps);
    MC_wx = 0;
    wy_bins[ibin] = MC_wy/static_cast<double>(Nsweeps);
    MC_wy = 0;
    wx2_bins[ibin] = MC_wx2/static_cast<double>(Nsweeps);
    MC_wx2 = 0;
    wy2_bins[ibin] = MC_wy2/static_cast<double>(Nsweeps);
    MC_wy2 = 0;
    mag_bins[ibin] = MC_mag/static_cast<double>(Nsweeps);
    MC_mag = 0;
    mag2_bins[ibin] = MC_mag2/static_cast<double>(Nsweeps);
    MC_mag2 = 0;
    for(int i=0; i<Nbonds; i++)
    {
      energy_diag_bins[ibin][i] = MC_energy_diag[i]/static_cast<double>(Nsweeps);
      energy_off_bins[ibin][i] = MC_energy_off[i]/static_cast<double>(Nsweeps);
      MC_energy_diag[i] = 0;
      MC_energy_off[i] = 0;
    }
    
    if(ibin%10 == 0 || ibin == Nbins-1)
    {
      std::ofstream out("./qc_obs/qc_energy_"+file_extension);
      for(int i=0; i<ibin+1; i++)
        out << energy_bins[i] << std::endl;
      out.close();
      out.open("./qc_obs/qc_sf_density_"+file_extension);
      for(int i=0; i<ibin+1; i++)
        out << sf_density_bins[i] << std::endl;
      out.close();
      out.open("./qc_obs/qc_wx_"+file_extension);
      for(int i=0; i<ibin+1; i++)
        out << wx_bins[i] << std::endl;
      out.close();
      out.open("./qc_obs/qc_wy_"+file_extension);
      for(int i=0; i<ibin+1; i++)
        out << wy_bins[i] << std::endl;
      out.close();
      out.open("./qc_obs/qc_wx2_"+file_extension);
      for(int i=0; i<ibin+1; i++)
        out << wx2_bins[i] << std::endl;
      out.close();
      out.open("./qc_obs/qc_wy2_"+file_extension);
      for(int i=0; i<ibin+1; i++)
        out << wy2_bins[i] << std::endl;
      out.close();
      out.open("./qc_obs/qc_Nparticles_"+file_extension);
      for(int i=0; i<ibin+1; i++)
        out << mag_bins[i] << std::endl;
      out.close();
      out.open("./qc_obs/qc_Nparticles2_"+file_extension);
      for(int i=0; i<ibin+1; i++)
        out << mag2_bins[i] << std::endl;
      out.close();
      
      std::vector<double> mean_expansion = {0,0};
      std::vector<double> mean_energy = {0,0};
      std::vector<double> mean_sfx = {0,0};
      std::vector<double> mean_sfy = {0,0};
      std::vector<double> mean_sf_density = {0,0};
      std::vector<double> mean_density = {0,0};
      std::vector<double> m = {0,0};
      std::vector<double> m2 = {0,0};
      std::vector<double> mean_compressibility = {0,0};
      std::vector<std::vector<double>> mean_diag(2,std::vector<double>(Nbonds,0));
      std::vector<std::vector<double>> mean_off(2,std::vector<double>(Nbonds,0));
      for(int i=0; i<=ibin; i++)
      {
        mean_expansion[0] += NHam_bins[i];
        mean_expansion[1] += NHam_bins[i]*NHam_bins[i];
        mean_energy[0] += energy_bins[i];
        mean_energy[1] += energy_bins[i]*energy_bins[i];
        mean_sfx[0] += wx2_bins[i]/Beta;
        mean_sfx[1] += wx2_bins[i]/Beta * wx2_bins[i]/Beta;
        mean_sfy[0] += wy2_bins[i]/Beta;
        mean_sfy[1] += wy2_bins[i]/Beta * wy2_bins[i]/Beta;
        mean_sf_density[0] += sf_density_bins[i];
        mean_sf_density[1] += sf_density_bins[i]*sf_density_bins[i];
        mean_density[0] += mag_bins[i]/Nsites;
        mean_density[1] += mag_bins[i]/Nsites*mag_bins[i]/Nsites;
        m[0] += mag_bins[i];
        m[1] += mag_bins[i]*mag_bins[i];
        m2[0] += mag2_bins[i];
        m2[1] += mag2_bins[i]*mag2_bins[i];
        mean_compressibility[1] += (mag2_bins[i] - mag_bins[i]*mag_bins[i])*Beta/Lx/Ly * (mag2_bins[i] - mag_bins[i]*mag_bins[i])*Beta/Lx/Ly;
        for(int b=0; b<Nbonds; b++)
        {
          mean_off[0][b] += energy_off_bins[i][b];
          mean_off[1][b] += energy_off_bins[i][b] * energy_off_bins[i][b];
          mean_diag[0][b] += energy_diag_bins[i][b];
          mean_diag[1][b] += energy_diag_bins[i][b] * energy_diag_bins[i][b];
        }
      }
      std::cout << "result after " << ibin << " bins" << std::endl;
      std::cout << "expansion_order = " << mean_expansion[0]/(ibin+1) << " +/- " << sqrt(mean_expansion[1]/(ibin+1)-pow(mean_expansion[0]/(ibin+1),2))/sqrt(ibin+1) << std::endl;
      std::cout << "energy = " << mean_energy[0]/(ibin+1) << " +/- " << sqrt(mean_energy[1]/(ibin+1)-pow(mean_energy[0]/(ibin+1),2))/sqrt(ibin+1) << std::endl;
      std::cout << "density = " << mean_density[0]/(ibin+1) << " +/- " << sqrt(mean_density[1]/(ibin+1)-pow(mean_density[0]/(ibin+1),2))/sqrt(ibin+1) << std::endl;
      std::cout << "sf_stifness_x = " << mean_sfx[0]/(ibin+1) << " +/- " << sqrt(mean_sfx[1]/(ibin+1)-pow(mean_sfx[0]/(ibin+1),2))/sqrt(ibin+1) << std::endl;
      std::cout << "sf_stifness_y = " << mean_sfy[0]/(ibin+1) << " +/- " << sqrt(mean_sfy[1]/(ibin+1)-pow(mean_sfy[0]/(ibin+1),2))/sqrt(ibin+1) << std::endl;
      std::cout << "sf_stifness = " << mean_sf_density[0]/(ibin+1) << " +/- " << sqrt(mean_sf_density[1]/(ibin+1)-pow(mean_sf_density[0]/(ibin+1),2))/sqrt(ibin+1) << std::endl;
      std::cout << "compressibility = " << (m2[0]/(ibin+1)-m[0]/(ibin+1)*m[0]/(ibin+1))*Beta/static_cast<double>(Lx*Ly) << std::endl;
      std::ofstream out_diag("./qc_obs/qc_energy_bonds_diag");
      for(int i=0; i<Nbonds; i++)
        out_diag << i << " " << mean_diag[0][i]/(ibin+1) << " +/- " << sqrt(mean_diag[1][i]/(ibin+1)-pow(mean_diag[0][i]/(ibin+1),2))/sqrt(ibin+1) << std::endl;
      out_diag.close();
      std::ofstream out_off("./qc_obs/qc_energy_bonds_off");
      for(int i=0; i<Nbonds; i++)
        out_off << i << " " << mean_off[0][i]/(ibin+1) << " +/- " << sqrt(mean_off[1][i]/(ibin+1)-pow(mean_off[0][i]/(ibin+1),2))/sqrt(ibin+1) << std::endl;
      out_off.close();
      std::cout << std::endl;
    }
  }
  
  template<typename T>
  void print_vector(std::vector<T> v)
  {
    for(size_t i=0; i<v.size(); i++)
      std::cout << v[i] << " ";
    std::cout << std::endl;
  }
  
  template<typename T>
  void print_matrix(std::vector<std::vector<T>> v)
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
  std::vector<size_t> argsort(const std::vector<T>& v)
  {
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i)
      idx[i] = i;
    
    std::sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
  }
    
  void set_path_prob_dirloop(int bond, const std::vector<double>& W_tilde)
  {
    std::vector<size_t> idx = argsort<double>(W_tilde);

    if(W_tilde[idx[2]]<=W_tilde[idx[1]]+W_tilde[idx[0]])
    {
      // for this bond there is a bounce free solution
      double a12 = -0.5*(W_tilde[idx[2]]-W_tilde[idx[1]]-W_tilde[idx[0]]);
      double a13 = 0.5*(W_tilde[idx[2]]-W_tilde[idx[1]]+W_tilde[idx[0]]);
      double a23 = 0.5*(W_tilde[idx[2]]+W_tilde[idx[1]]-W_tilde[idx[0]]);

      // determine the path type corresponding to a12, a13, a23
      std::vector<size_t> vertex_transition_a12 = {idx[0],idx[1]};
      std::vector<size_t> vertex_transition_a13 = {idx[0],idx[2]};
      std::vector<size_t> vertex_transition_a23 = {idx[1],idx[2]};

      path_prob[bond][idx[0]][get_path(vertex_transition_a12)] = a12;
      path_prob[bond][idx[0]][get_path(vertex_transition_a13)] = a13;
      path_prob[bond][idx[1]][get_path(vertex_transition_a23)] = a23;
      path_prob[bond][idx[1]][get_path(vertex_transition_a12)] = a12;
      path_prob[bond][idx[2]][get_path(vertex_transition_a13)] = a13;
      path_prob[bond][idx[2]][get_path(vertex_transition_a23)] = a23;
    }
    else
    {
      // bouncing for this bond can't be avoided, bounce on the vertex config. with the highest weight

      double b3 = W_tilde[idx[2]]-W_tilde[idx[1]]-W_tilde[idx[0]];
      double a13 = W_tilde[idx[0]];
      double a23 = W_tilde[idx[1]];
      
      // determine the path type corresponding to a12, a13, a23
      
      std::vector<size_t> vertex_transition_a13 = {idx[0],idx[2]};
      std::vector<size_t> vertex_transition_a23 = {idx[1],idx[2]};
      path_prob[bond][idx[2]][static_cast<int>(PATH::bounce)] = b3;
      path_prob[bond][idx[0]][get_path(vertex_transition_a13)] = a13;
      path_prob[bond][idx[1]][get_path(vertex_transition_a23)] = a23;
      path_prob[bond][idx[2]][get_path(vertex_transition_a13)] = a13;
      path_prob[bond][idx[2]][get_path(vertex_transition_a23)] = a23;
    }
  }
  
  void set_path_prob_dirloop_nosymmetry(int bond, const std::vector<double>& W)
  {
    // first field

    std::vector<double> W_tilde = {W[0],W[3],W[5]};
    std::vector<int> index_in = {0,3,5};
    std::vector<int> entrance = {0,1,3};
    std::vector<size_t> idx = argsort<double>(W_tilde);

    set_dirloop(bond,W_tilde,index_in,entrance,idx);
    
    // second field

    W_tilde = {W[0],W[3],W[4]};
    index_in = {0,3,4};
    entrance = {1,0,2};
    idx = argsort<double>(W_tilde);
    set_dirloop(bond,W_tilde,index_in,entrance,idx);
    
    // third field

    W_tilde = {W[0],W[2],W[4]};
    index_in = {0,2,4};
    entrance = {2,3,1};
    idx = argsort<double>(W_tilde);
    set_dirloop(bond,W_tilde,index_in,entrance,idx);
    
    // fourth field

    W_tilde = {W[0],W[2],W[5]};
    index_in = {0,2,5};
    entrance = {3,2,0};
    idx = argsort<double>(W_tilde);
    set_dirloop(bond,W_tilde,index_in,entrance,idx);
    
    // fiveth field

    W_tilde = {W[1],W[2],W[4]};
    index_in = {1,2,4};
    entrance = {0,1,3};
    idx = argsort<double>(W_tilde);
    set_dirloop(bond,W_tilde,index_in,entrance,idx);
    
    // sixth field

    W_tilde = {W[1],W[2],W[5]};
    index_in = {1,2,5};
    entrance = {1,0,2};
    idx = argsort<double>(W_tilde);
    set_dirloop(bond,W_tilde,index_in,entrance,idx);
    
    // seventh field

    W_tilde = {W[1],W[3],W[5]};
    index_in = {1,3,5};
    entrance = {2,3,1};
    idx = argsort<double>(W_tilde);
    set_dirloop(bond,W_tilde,index_in,entrance,idx);
    
    // eighth field

    W_tilde = {W[1],W[3],W[4]};
    index_in = {1,3,4};
    entrance = {3,2,0};
    idx = argsort<double>(W_tilde);
    set_dirloop(bond,W_tilde,index_in,entrance,idx);
    
  }
  
  void set_dirloop(int bond,const std::vector<double>& W_tilde, const std::vector<int>& index_in, const std::vector<int>& entrance, const std::vector<size_t>& idx)
  {
    if(W_tilde[idx[2]]<=W_tilde[idx[1]]+W_tilde[idx[0]])
    {

      // for this bond there is a bounce free solution
      double a12 = -0.5*(W_tilde[idx[2]]-W_tilde[idx[1]]-W_tilde[idx[0]]);
      double a13 = 0.5*(W_tilde[idx[2]]-W_tilde[idx[1]]+W_tilde[idx[0]]);
      double a23 = 0.5*(W_tilde[idx[2]]+W_tilde[idx[1]]-W_tilde[idx[0]]);

      // determine the path type corresponding to a12, a13, a23
      std::vector<int> vertex_transition_a12 = {entrance[idx[0]],index_in[idx[0]],index_in[idx[1]]};
      std::vector<int> vertex_transition_a13 = {entrance[idx[0]],index_in[idx[0]],index_in[idx[2]]};
      std::vector<int> vertex_transition_a23 = {entrance[idx[1]],index_in[idx[1]],index_in[idx[2]]};

      path_prob_nosymmetry[bond][index_in[idx[0]]][entrance[idx[0]]][get_exit(vertex_transition_a12)] = a12;
      path_prob_nosymmetry[bond][index_in[idx[0]]][entrance[idx[0]]][get_exit(vertex_transition_a13)] = a13;
      path_prob_nosymmetry[bond][index_in[idx[1]]][entrance[idx[1]]][get_exit(vertex_transition_a23)] = a23;
      path_prob_nosymmetry[bond][index_in[idx[1]]][entrance[idx[1]]][get_exit(vertex_transition_a12)] = a12;
      path_prob_nosymmetry[bond][index_in[idx[2]]][entrance[idx[2]]][get_exit(vertex_transition_a13)] = a13;
      path_prob_nosymmetry[bond][index_in[idx[2]]][entrance[idx[2]]][get_exit(vertex_transition_a23)] = a23;

    }
    else
    {
      // bouncing for this bond can't be avoided, bounce on the vertex config. with the highest weight

      double b3 = W_tilde[idx[2]]-W_tilde[idx[1]]-W_tilde[idx[0]];

      double a13 = W_tilde[idx[0]];
      double a23 = W_tilde[idx[1]];
      
      // determine the path type corresponding to a12, a13, a23
      std::vector<int> vertex_transition_a13 = {entrance[idx[0]],index_in[idx[0]],index_in[idx[2]]};
      std::vector<int> vertex_transition_a23 = {entrance[idx[1]],index_in[idx[1]],index_in[idx[2]]};

      path_prob_nosymmetry[bond][index_in[idx[2]]][entrance[idx[2]]][static_cast<int>(PATH::bounce)] = b3;
      path_prob_nosymmetry[bond][index_in[idx[0]]][entrance[idx[0]]][get_exit(vertex_transition_a13)] = a13;
      path_prob_nosymmetry[bond][index_in[idx[1]]][entrance[idx[1]]][get_exit(vertex_transition_a23)] = a23;
      path_prob_nosymmetry[bond][index_in[idx[2]]][entrance[idx[2]]][get_exit(vertex_transition_a13)] = a13;
      path_prob_nosymmetry[bond][index_in[idx[2]]][entrance[idx[2]]][get_exit(vertex_transition_a23)] = a23;
    }
  }
  
  int get_exit(const std::vector<int>& trans)
  {
    if(trans[1] == 0 ) // W0
    {
      if(trans[0] == 0)  // entrance_leg
      {
        if(trans[2] == 3) // to W3
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 5) // to W5
          return static_cast<int>(PATH::switch_continue);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 1)  // entrance_leg
      {
        if(trans[2] == 3) // to W3
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 4) // to W4
          return static_cast<int>(PATH::switch_continue);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 2)  // entrance_leg
      {
        if(trans[2] == 2) // to W2
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 4) // to W4
          return static_cast<int>(PATH::switch_continue);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 3)  // entrance_leg
      {
        if(trans[2] == 2) // to W2
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 5) // to W5
          return static_cast<int>(PATH::switch_continue);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else
      {
        std::cout << "can determine entrance leg " << trans[1] << std::endl;
        exit(1);
        return -1;
      }
    }
    
    else if(trans[1] == 1 ) // W1
    {
      if(trans[0] == 0)  // entrance_leg
      {
        if(trans[2] == 2) // to W2
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 4) // to W4
          return static_cast<int>(PATH::switch_continue);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 1)  // entrance_leg
      {
        if(trans[2] == 2) // to W2
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 5) // to W5
          return static_cast<int>(PATH::switch_continue);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 2)  // entrance_leg
      {
        if(trans[2] == 3) // to W3
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 5) // to W5
          return static_cast<int>(PATH::switch_continue);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 3)  // entrance_leg
      {
        if(trans[2] == 3) // to W3
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 4) // to W4
          return static_cast<int>(PATH::switch_continue);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else
      {
        std::cout << "can determine entrance leg " << trans[1] << std::endl;
        exit(1);
        return -1;
      }
    }
    
    else if(trans[1] == 2 ) // W2
    {
      if(trans[0] == 0)  // entrance_leg
      {
        if(trans[2] == 1) // to W1
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 5) // to W5
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 1)  // entrance_leg
      {
        if(trans[2] == 1) // to W1
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 4) // to W4
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 2)  // entrance_leg
      {
        if(trans[2] == 0) // to W0
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 5) // to W5
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 3)  // entrance_leg
      {
        if(trans[2] == 0) // to W0
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 4) // to W4
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else
      {
        std::cout << "can determine entrance leg " << trans[1] << std::endl;
        exit(1);
        return -1;
      }
    }
    
    else if(trans[1] == 3 ) // W3
    {
      if(trans[0] == 0)  // entrance_leg
      {
        if(trans[2] == 0) // to W0
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 4) // to W4
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 1)  // entrance_leg
      {
        if(trans[2] == 0) // to W0
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 5) // to W5
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 2)  // entrance_leg
      {
        if(trans[2] == 1) // to W1
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 4) // to W4
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 3)  // entrance_leg
      {
        if(trans[2] == 1) // to W1
          return static_cast<int>(PATH::switch_reverse);
        else if (trans[2] == 5) // to W5
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else
      {
        std::cout << "Can't identify entrance leg" << trans[1] << std::endl;
        exit(1);
        return -1;
      }
    }
    
    else if(trans[1] == 4 ) // W4
    {
      if(trans[0] == 0)  // entrance_leg
      {
        if(trans[2] == 1) // to W1
          return static_cast<int>(PATH::switch_continue);
        else if (trans[2] == 3) // to W3
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 1)  // entrance_leg
      {
        if(trans[2] == 0) // to W0
          return static_cast<int>(PATH::switch_continue);
        else if (trans[2] == 2) // to W2
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 2)  // entrance_leg
      {
        if(trans[2] == 0) // to W0
          return static_cast<int>(PATH::switch_continue);
        else if (trans[2] == 3) // to W3
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 3)  // entrance_leg
      {
        if(trans[2] == 1) // to W1
          return static_cast<int>(PATH::switch_continue);
        else if (trans[2] == 2) // to W2
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else
      {
        std::cout << "Can't identify entrance leg" << trans[1] << std::endl;
        exit(1);
        return -1;
      }
    }
    
    else if(trans[1] == 5 ) // W5
    {
      if(trans[0] == 0)  // entrance_leg
      {
        if(trans[2] == 0) // to W0
          return static_cast<int>(PATH::switch_continue);
        else if (trans[2] == 2) // to W2
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 1)  // entrance_leg
      {
        if(trans[2] == 1) // to W1
          return static_cast<int>(PATH::switch_continue);
        else if (trans[2] == 3) // to W3
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 2)  // entrance_leg
      {
        if(trans[2] == 1) // to W1
          return static_cast<int>(PATH::switch_continue);
        else if (trans[2] == 2) // to W2
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else if(trans[0] == 3)  // entrance_leg
      {
        if(trans[2] == 0) // to W0
          return static_cast<int>(PATH::switch_continue);
        else if (trans[2] == 3) // to W3
          return static_cast<int>(PATH::straight);
        else
        {
          std::cout << "no other path possible " << trans[1] << std::endl;
          exit(1);
          return -1;
        }
      }
      else
      {
        std::cout << "Can't identify entrance leg" << trans[1] << std::endl;
        exit(1);
        return -1;
      }
    }
    else
    {
      std::cout << "Can't identify the correct field" << trans[1] << std::endl;
      exit(1);
      return -1;
    }
  }
  
  int get_path(const std::vector<size_t>& trans)
  {
    if(trans[0]==0) // off diagonal
    {
      // going to
      if(trans[1]==1) // diagonal unequal
        return static_cast<int>(PATH::switch_reverse);
      else if(trans[1]==2) // diagonal equal
        return static_cast<int>(PATH::switch_continue);
      else
      {
        std::cout << "can't determine path" << std::endl;
        exit(1);
        return -1;
      }
    }
    else if(trans[0]==1) // diagonal unequal
    {
      // going to
      if(trans[1]==0) // off diagonal
        return static_cast<int>(PATH::switch_reverse);
      else if(trans[1]==2) // diagonal equal
        return static_cast<int>(PATH::straight);
      else
      {
        std::cout << "can't determine path" << std::endl;
        exit(1);
        return -1;
      }
    }
    else // diagonal equal spins
    {
      // going to
      if(trans[1]==0) // off-diagonal
        return static_cast<int>(PATH::switch_continue);
      else if(trans[1]==1) // diagonal unequal spins
        return static_cast<int>(PATH::straight);
      else
      {
        std::cout << "can't determine path" << std::endl;
        exit(1);
        return -1;
      }
    }
  }
  
  void set_path_prob_HeatBath(int bond, const std::vector<double>& W)
  {
    // off-diagonal vertex W0
    path_prob_heatBath[bond][0][0][0] = W[0]/(W[0]+W[3]+W[5]);
    path_prob_heatBath[bond][0][0][1] = W[3]/(W[0]+W[3]+W[5]);
    path_prob_heatBath[bond][0][0][2] = 0;
    path_prob_heatBath[bond][0][0][3] = W[5]/(W[0]+W[3]+W[5]);
    
    path_prob_heatBath[bond][0][1][0] = W[3]/(W[0]+W[3]+W[4]);
    path_prob_heatBath[bond][0][1][1] = W[0]/(W[0]+W[3]+W[4]);
    path_prob_heatBath[bond][0][1][2] = W[4]/(W[0]+W[3]+W[4]);
    path_prob_heatBath[bond][0][1][3] = 0;
    
    path_prob_heatBath[bond][0][2][0] = 0;
    path_prob_heatBath[bond][0][2][1] = W[4]/(W[0]+W[2]+W[4]);
    path_prob_heatBath[bond][0][2][2] = W[0]/(W[0]+W[2]+W[4]);
    path_prob_heatBath[bond][0][2][3] = W[2]/(W[0]+W[2]+W[4]);
    
    path_prob_heatBath[bond][0][3][0] = W[5]/(W[0]+W[2]+W[5]);
    path_prob_heatBath[bond][0][3][1] = 0;
    path_prob_heatBath[bond][0][3][2] = W[2]/(W[0]+W[2]+W[5]);
    path_prob_heatBath[bond][0][3][3] = W[0]/(W[0]+W[2]+W[5]);
    
    // off-diagonal vertex W1
    path_prob_heatBath[bond][1][0][0] = W[1]/(W[1]+W[2]+W[4]);
    path_prob_heatBath[bond][1][0][1] = W[2]/(W[1]+W[2]+W[4]);
    path_prob_heatBath[bond][1][0][2] = 0;
    path_prob_heatBath[bond][1][0][3] = W[4]/(W[1]+W[2]+W[4]);
    
    path_prob_heatBath[bond][1][1][0] = W[2]/(W[1]+W[2]+W[5]);
    path_prob_heatBath[bond][1][1][1] = W[1]/(W[1]+W[2]+W[5]);
    path_prob_heatBath[bond][1][1][2] = W[5]/(W[1]+W[2]+W[5]);
    path_prob_heatBath[bond][1][1][3] = 0;
    
    path_prob_heatBath[bond][1][2][0] = 0;
    path_prob_heatBath[bond][1][2][1] = W[5]/(W[1]+W[3]+W[5]);
    path_prob_heatBath[bond][1][2][2] = W[1]/(W[1]+W[3]+W[5]);
    path_prob_heatBath[bond][1][2][3] = W[3]/(W[1]+W[3]+W[5]);
    
    path_prob_heatBath[bond][1][3][0] = W[4]/(W[1]+W[3]+W[4]);
    path_prob_heatBath[bond][1][3][1] = 0;
    path_prob_heatBath[bond][1][3][2] = W[3]/(W[1]+W[3]+W[4]);
    path_prob_heatBath[bond][1][3][3] = W[1]/(W[1]+W[3]+W[4]);
    
    // diagonal - unequal spins vertex W2
    path_prob_heatBath[bond][2][0][0] = W[2]/(W[1]+W[2]+W[5]);
    path_prob_heatBath[bond][2][0][1] = W[1]/(W[1]+W[2]+W[5]);
    path_prob_heatBath[bond][2][0][2] = W[5]/(W[1]+W[2]+W[5]);
    path_prob_heatBath[bond][2][0][3] = 0;
    
    path_prob_heatBath[bond][2][1][0] = W[1]/(W[1]+W[2]+W[4]);
    path_prob_heatBath[bond][2][1][1] = W[2]/(W[1]+W[2]+W[4]);
    path_prob_heatBath[bond][2][1][2] = 0;
    path_prob_heatBath[bond][2][1][3] = W[4]/(W[1]+W[2]+W[4]);
    
    path_prob_heatBath[bond][2][2][0] = W[5]/(W[0]+W[2]+W[5]);
    path_prob_heatBath[bond][2][2][1] = 0;
    path_prob_heatBath[bond][2][2][2] = W[2]/(W[0]+W[2]+W[5]);
    path_prob_heatBath[bond][2][2][3] = W[0]/(W[0]+W[2]+W[5]);
    
    path_prob_heatBath[bond][2][3][0] = 0;
    path_prob_heatBath[bond][2][3][1] = W[4]/(W[0]+W[2]+W[4]);
    path_prob_heatBath[bond][2][3][2] = W[0]/(W[0]+W[2]+W[4]);
    path_prob_heatBath[bond][2][3][3] = W[2]/(W[0]+W[2]+W[4]);
    
    // diagonal - unequal spins vertex W3
    path_prob_heatBath[bond][3][0][0] = W[3]/(W[0]+W[3]+W[4]);
    path_prob_heatBath[bond][3][0][1] = W[0]/(W[0]+W[3]+W[4]);
    path_prob_heatBath[bond][3][0][2] = W[4]/(W[0]+W[3]+W[4]);
    path_prob_heatBath[bond][3][0][3] = 0;
    
    path_prob_heatBath[bond][3][1][0] = W[0]/(W[0]+W[3]+W[5]);
    path_prob_heatBath[bond][3][1][1] = W[3]/(W[0]+W[3]+W[5]);
    path_prob_heatBath[bond][3][1][2] = 0;
    path_prob_heatBath[bond][3][1][3] = W[5]/(W[0]+W[3]+W[5]);
    
    path_prob_heatBath[bond][3][2][0] = W[4]/(W[1]+W[3]+W[4]);
    path_prob_heatBath[bond][3][2][1] = 0;
    path_prob_heatBath[bond][3][2][2] = W[3]/(W[1]+W[3]+W[4]);
    path_prob_heatBath[bond][3][2][3] = W[1]/(W[1]+W[3]+W[4]);
    
    path_prob_heatBath[bond][3][3][0] = 0;
    path_prob_heatBath[bond][3][3][1] = W[5]/(W[1]+W[3]+W[5]);
    path_prob_heatBath[bond][3][3][2] = W[1]/(W[1]+W[3]+W[5]);
    path_prob_heatBath[bond][3][3][3] = W[3]/(W[1]+W[3]+W[5]);
    
    // diagonal - equal spins vertex W4
    path_prob_heatBath[bond][4][0][0] = W[4]/(W[1]+W[3]+W[4]);
    path_prob_heatBath[bond][4][0][1] = 0;
    path_prob_heatBath[bond][4][0][2] = W[3]/(W[1]+W[3]+W[4]);
    path_prob_heatBath[bond][4][0][3] = W[1]/(W[1]+W[3]+W[4]);
    
    path_prob_heatBath[bond][4][1][0] = 0;
    path_prob_heatBath[bond][4][1][1] = W[4]/(W[0]+W[2]+W[4]);
    path_prob_heatBath[bond][4][1][2] = W[0]/(W[0]+W[2]+W[4]);
    path_prob_heatBath[bond][4][1][3] = W[2]/(W[0]+W[2]+W[4]);
    
    path_prob_heatBath[bond][4][2][0] = W[3]/(W[0]+W[3]+W[4]);
    path_prob_heatBath[bond][4][2][1] = W[0]/(W[0]+W[3]+W[4]);
    path_prob_heatBath[bond][4][2][2] = W[4]/(W[0]+W[3]+W[4]);
    path_prob_heatBath[bond][4][2][3] = 0;
    
    path_prob_heatBath[bond][4][3][0] = W[1]/(W[1]+W[2]+W[4]);
    path_prob_heatBath[bond][4][3][1] = W[2]/(W[1]+W[2]+W[4]);
    path_prob_heatBath[bond][4][3][2] = 0;
    path_prob_heatBath[bond][4][3][3] = W[4]/(W[1]+W[2]+W[4]);
    
    // diagonal - equal spins vertex W5
    path_prob_heatBath[bond][5][0][0] = W[5]/(W[0]+W[2]+W[5]);
    path_prob_heatBath[bond][5][0][1] = 0;
    path_prob_heatBath[bond][5][0][2] = W[2]/(W[0]+W[2]+W[5]);
    path_prob_heatBath[bond][5][0][3] = W[0]/(W[0]+W[2]+W[5]);
    
    path_prob_heatBath[bond][5][1][0] = 0;
    path_prob_heatBath[bond][5][1][1] = W[5]/(W[1]+W[3]+W[5]);
    path_prob_heatBath[bond][5][1][2] = W[1]/(W[1]+W[3]+W[5]);
    path_prob_heatBath[bond][5][1][3] = W[3]/(W[1]+W[3]+W[5]);
    
    path_prob_heatBath[bond][5][2][0] = W[2]/(W[1]+W[2]+W[5]);
    path_prob_heatBath[bond][5][2][1] = W[1]/(W[1]+W[2]+W[5]);
    path_prob_heatBath[bond][5][2][2] = W[5]/(W[1]+W[2]+W[5]);
    path_prob_heatBath[bond][5][2][3] = 0;
    
    path_prob_heatBath[bond][5][3][0] = W[0]/(W[0]+W[3]+W[5]);
    path_prob_heatBath[bond][5][3][1] = W[3]/(W[0]+W[3]+W[5]);
    path_prob_heatBath[bond][5][3][2] = 0;
    path_prob_heatBath[bond][5][3][3] = W[5]/(W[0]+W[3]+W[5]);
  }
  
};

#endif
