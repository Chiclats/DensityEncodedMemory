#include<bits/stdc++.h>
#include"./StateProcessing_new.hpp"
#include"./ActiveNematics_new.hpp"

int main(){

  const BoxMeshing_2D::BoxShape box_shape = BoxMeshing_2D::create_box_shape(4.0, 3.0, 1.0);
  
  const ActiveNematics_2D::Active_Nematics_Parameters params={
    20,
    box_shape,
    0.1,
    0.5,
    0.01
  };

  std::vector<ActiveNematics_2D::Particle> particles=ActiveNematics_2D::initialize_particles(params);

  std::cout<<"==========\n";
  
  for(const auto& p : particles)
    p.output();

  std::cout<<"==========\n";

  std::vector<std::vector<unsigned int>> meshed_index_mat=ActiveNematics_2D::initialize_meshed_index_mat(particles, box_shape);

  for(unsigned int i_box = 0; i_box < std::get<5>(box_shape); i_box++) {
    std::cout << "Box " << i_box << ": ";
    for (unsigned int particle_index : meshed_index_mat[i_box]) {
      std::cout << particle_index << " ";
    }
    std::cout << "\n";
  }
  
}
