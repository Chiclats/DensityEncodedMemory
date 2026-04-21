/*--------------------
ver 260224
--------------------*/

#ifndef ActiveNematics_new_hpp
#define ActiveNematics_new_hpp

#include<bits/stdc++.h>
#include<omp.h>
#include"./BasicTools_new.hpp"
#include"./StateProcessing_new.hpp"

namespace ActiveNematics_2D{

  struct Particle{
    
    double x,y,vx,vy;
    // position and velocity of the particle
    // WARNING: velocity is a unit vector
    
    unsigned int index;
    // to store the index of the particle in the original PG

    unsigned int box_index;
    // to store the index of the box to which the particle belongs, for updating MeshedIndexMat

    void update_position_periodic(const double velocity,const BoxMeshing_2D::BoxShape& box_shape)
    {// to update the position of the particle based on its velocity and the time step, and apply periodic boundary condition
      x += vx * velocity;
      y += vy * velocity;
      
      auto [system_size_x, system_size_y, _, _, _, _] = box_shape;
      
      if (x < 0) x += system_size_x;
      else if (x >= system_size_x) x -= system_size_x;
      
      if (y < 0) y += system_size_y;
      else if (y >= system_size_y) y -= system_size_y;
    }
    
    void decide_box_index(const BoxMeshing_2D::BoxShape& box_shape)
    {// to decide the box index for the particle and update `box_index`
      box_index = BoxMeshing_2D::decide_box_index(x, y, box_shape);
    }

    void output(std::ostream& os = std::cout) const
    {// to output the particle information to an output stream
      os << x <<' ' << y <<' ' << vx <<' ' << vy <<' '<< index <<' ' << box_index << std::endl;
    }
  };

  struct Active_Nematics_Parameters{
    unsigned int num_particles;
    BoxMeshing_2D::BoxShape box_shape; // interaction range is always equal to char_length
    double velocity;       // self-propulsion velocity
    double sigma;          // noise strength
    double alpha;          // flipping rate
    double guiding_coeff;  // tetratic field strength
  };

  double distance_square_periodic(const Particle& p1,const Particle& p2,const BoxMeshing_2D::BoxShape& box_shape)
  {// to calculate the distance^2 between two particles under periodic boundary condition
    auto [system_size_x, system_size_y, _, _, _, _] = box_shape;
    
    double dx = std::abs(p1.x - p2.x);
    double dy = std::abs(p1.y - p2.y);
    
    if (dx > system_size_x / 2) dx = system_size_x - dx;
    if (dy > system_size_y / 2) dy = system_size_y - dy;
    
    return dx * dx + dy * dy;
  }

  std::vector<Particle> initialize_particles(const Active_Nematics_Parameters& params)
  {// to initialize the particles with random positions and velocities
    std::vector<Particle> particles(params.num_particles);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_x(0, std::get<0>(params.box_shape));
    std::uniform_real_distribution<> dis_y(0, std::get<1>(params.box_shape));
    std::uniform_real_distribution<> dis_angle(0, 2 * std::numbers::pi);
		
    for (unsigned int i = 0; i < params.num_particles; i++) {
      double angle = dis_angle(gen);
      particles[i] = {dis_x(gen), dis_y(gen), cos(angle), sin(angle), i, 0};
      particles[i].decide_box_index(params.box_shape);
    }
   
    return particles;
  }

  std::vector<std::vector<unsigned int>> initialize_meshed_index_mat(const std::vector<Particle>& particles,const BoxMeshing_2D::BoxShape& box_shape)
  {// to initialize the meshed index matrix based on the initial particle positions
    auto [_, _, _, _, _, total_box_num] = box_shape;
    std::vector<std::vector<unsigned int>> meshed_index_mat(total_box_num);
		
    for (const auto& p : particles) {
      meshed_index_mat[p.box_index].push_back(p.index);
    }
    
    return meshed_index_mat;
  }

  void Evolve_Turned1FoldMean(std::vector<Particle>& PG, std::vector<std::vector<unsigned int>>& meshed_index_mat, const Active_Nematics_Parameters& params, std::mt19937_64& rd_gen)
  {// Dir(t+1) = (\pm 1) * arg( mean( sgn( cos( NeighborDir - Dir ) ) * exp( i NeighborDir ) ) )
    // Each step the vector have a probability Alpha of reversing (+1 to -1)

    const auto [num_particles, box_shape, velocity, sigma, alpha, guiding_coeff] = params;

    std::vector<double> neighbor_dir_list_x(num_particles, 0), neighbor_dir_list_y(num_particles, 0);
    // to store sum(vx) and sum(vy)

    // calculate the mean velocity dir
    for( unsigned int i_box=0; i_box<std::get<5>(box_shape); i_box++)
      for(unsigned int i_p=0; i_p<meshed_index_mat[i_box].size(); i_p++){

	const unsigned int index_i = meshed_index_mat[i_box][i_p]; // index of the particle in PG

	// itself is taken into account
	neighbor_dir_list_x[index_i] += PG[index_i].vx; // cos Dir
	neighbor_dir_list_y[index_i] += PG[index_i].vy; // sin Dir

	// in i_box
	for(unsigned int j_p=i_p+1; j_p<meshed_index_mat[i_box].size(); j_p++){
	  
	  const unsigned int index_j = meshed_index_mat[i_box][j_p]; // index of the particle in PG

	  if(distance_square_periodic(PG[index_i], PG[index_j], box_shape) <= std::get<2>(box_shape) * std::get<2>(box_shape)){
	    // if distance <= char_length
	    if(PG[index_j].vx * PG[index_i].vx + PG[index_j].vy * PG[index_i].vy > 0){ // cos( Dir_i - Dir_j ) > 0
	      neighbor_dir_list_x[index_i] += PG[index_j].vx;
	      neighbor_dir_list_y[index_i] += PG[index_j].vy;
	      neighbor_dir_list_x[index_j] += PG[index_i].vx;
	      neighbor_dir_list_y[index_j] += PG[index_i].vy;
	    }
	    else{ // cos( Dir_i - Dir_j ) < 0
	      neighbor_dir_list_x[index_i] += -PG[index_j].vx;
	      neighbor_dir_list_y[index_i] += -PG[index_j].vy;
	      neighbor_dir_list_x[index_j] += -PG[index_i].vx;
	      neighbor_dir_list_y[index_j] += -PG[index_i].vy;
	    }
	  }
	}

	// inter-box
	for(unsigned int i_adj_box=0; i_adj_box<4; i_adj_box++){
	  // only upper and left boxes are counted to avoid double counting

	  // absolute box index
	  const unsigned int adj_box_index = BoxMeshing_2D::decide_box_neighbor_index_ordered(i_box, i_adj_box, box_shape);

	  for(const auto index_j : meshed_index_mat[adj_box_index]){

	    if(distance_square_periodic(PG[index_i],PG[index_j],box_shape) <= std::get<2>(box_shape) * std::get<2>(box_shape)){
	      // if distance <= char_length
	      if(PG[index_j].vx * PG[index_i].vx + PG[index_j].vy * PG[index_i].vy > 0){ // cos( Dir_i - Dir_j ) > 0
		neighbor_dir_list_x[index_i] += PG[index_j].vx;
		neighbor_dir_list_y[index_i] += PG[index_j].vy;
		neighbor_dir_list_x[index_j] += PG[index_i].vx;
		neighbor_dir_list_y[index_j] += PG[index_i].vy;
	      }
	      else{ // cos( Dir_i - Dir_j ) < 0
		neighbor_dir_list_x[index_i] += -PG[index_j].vx;
		neighbor_dir_list_y[index_i] += -PG[index_j].vy;
		neighbor_dir_list_x[index_j] += -PG[index_i].vx;
		neighbor_dir_list_y[index_j] += -PG[index_i].vy;
	      }
	    }
	  }
	}
      }

    //update particle directions and positions
    double new_dir;
    std::uniform_real_distribution<double> noise_dist(-std::numbers::pi * sigma, std::numbers::pi * sigma);
    std::uniform_real_distribution<double> flip_dist(0, 1);
    
    for(unsigned int i_p=0; i_p<num_particles; i_p++){
      // to avoid everything cancelled out
      if(neighbor_dir_list_y[i_p] == 0 && neighbor_dir_list_x[i_p] == 0)
	new_dir = atan2(PG[i_p].vy, PG[i_p].vx);
      else new_dir = atan2(neighbor_dir_list_y[i_p], neighbor_dir_list_x[i_p]);

      //noise
      new_dir += noise_dist(rd_gen);
      //tetratic field
      new_dir += guiding_coeff * (4 * PG[i_p].vx * PG[i_p].vx * PG[i_p].vx * PG[i_p].vy
				  - 4 * PG[i_p].vx * PG[i_p].vy * PG[i_p].vy * PG[i_p].vy); // GuidingCoeff*sin(4*Dir)

      //probable flipping
      if(flip_dist(rd_gen) <= alpha)
	new_dir += std::numbers::pi;

      //update velocity
      PG[i_p].vx = cos(new_dir);
      PG[i_p].vy = sin(new_dir);

      //update position
      PG[i_p].update_position_periodic(velocity, box_shape);
    }

    //update meshed index mat
    std::vector<std::vector<unsigned int>> new_meshed_index_mat(std::get<5>(box_shape));

    for(auto& p : PG) {
      p.decide_box_index(box_shape);
      new_meshed_index_mat[p.box_index].push_back(p.index);
    }

    meshed_index_mat = std::move(new_meshed_index_mat);
  }
  
}

#endif//ActiveNematics_new_hpp
