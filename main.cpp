#include "sse_XXZ.hpp"


int main (int argc, char *argv[])
{
  
  std::map<std::string,double> params;
  std::ifstream in_config("params");
  std::string buffer;
  in_config >> params["Lx"]; getline(in_config,buffer); // extend in x direction
  in_config >> params["Ly"]; getline(in_config,buffer); // extend in y direction
  in_config >> params["Lz"]; getline(in_config,buffer); // extend in z direction
  in_config >> params["BCx"]; getline(in_config,buffer); // 1: PBC, 0: PBC in x direction
  in_config >> params["BCy"]; getline(in_config,buffer); // 1: PBC, 0: PBC in y direction
  in_config >> params["BCz"]; getline(in_config,buffer); // 1: PBC, 0: PBC in z direction
  in_config >> params["T"]; getline(in_config,buffer); 
  in_config >> params["Ntherm"]; getline(in_config,buffer); 
  in_config >> params["Nbins"]; getline(in_config,buffer); 
  in_config >> params["Nsweeps"]; getline(in_config,buffer); // number of sweeps per bin
  in_config >> params["Jx"]; getline(in_config,buffer); 
  in_config >> params["Jy"]; getline(in_config,buffer); 
  in_config >> params["Delta"]; getline(in_config,buffer); // anisotropy
  in_config >> params["V"]; getline(in_config,buffer);  // hopping disorder
  in_config >> params["epsilon"]; getline(in_config,buffer); // to ensure positive weights
  in_config >> params["h"]; getline(in_config,buffer);  // diagonal disorder
  
  std::string file_extension = "";
  sse_XXZ disordered_XXZ(params,file_extension);
  
  std::cout << std::endl << "# Starting simulation with parameters" << std::endl << std::endl;
  for(auto iter_params : params)
    std::cout << "# " << iter_params.first << " " << iter_params.second << std::endl;
  
  for(int step=0; step<params["Ntherm"]; step++)
  {
    disordered_XXZ.diagonal_update();
    disordered_XXZ.adjust_cutoff(step);
    disordered_XXZ.build_vertex_list();
    disordered_XXZ.operator_loop_update_therm();
  }
  std::cout << std::endl <<  "# Finished thermalization phase" << std::endl << std::endl;
  for(int ibin = 0; ibin<params["Nbins"]; ibin++)
  {
    for(int step = 0; step<params["Nsweeps"]; step++)
    {
      disordered_XXZ.diagonal_update();
      disordered_XXZ.build_vertex_list();
      disordered_XXZ.operator_loop_update();
      disordered_XXZ.estimator();
    }
    disordered_XXZ.set_bin(ibin);
  }
}

