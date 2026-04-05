#include<bits/stdc++.h>
#include"./StateProcessing_new.hpp"

int main(){

  double
    system_size_x = 1000.5,
    system_size_y = 2000.3,
    char_length = 1.0;

  auto box_shape = BoxMeshing_2D::create_box_shape(system_size_x, system_size_y, char_length);

  std::cout<<"system_size_x = "<<std::get<0>(box_shape)<<std::endl
	   <<"system_size_y = "<<std::get<1>(box_shape)<<std::endl
	   <<"char_length = "<<std::get<2>(box_shape)<<std::endl
	   <<"box_num_x = "<<std::get<3>(box_shape)<<std::endl
	   <<"box_num_y = "<<std::get<4>(box_shape)<<std::endl
	   <<"total_box_num = "<<std::get<5>(box_shape)<<std::endl;

  auto adjacent_box_index_vec = BoxMeshing_2D::all_adjacent_box_index(box_shape);

  //----------

  const int ITERATIONS = 1e7;
  std::vector<int> indices(ITERATIONS);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, 2000000 - 1);
  for (int i = 0; i < ITERATIONS; ++i) {
    indices[i] = dis(gen);
  }

  double sum = 0;

  //---------- calculate

  auto start_time = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < ITERATIONS; ++i)
    sum += BoxMeshing_2D::decide_box_neighbor_index(indices[i], -1, -1, box_shape);

  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

  std::cout << "(1)Time taken: " << duration << " milliseconds"
	    << " Sum: " << sum << std::endl;
  
  //---------- read from table

  sum = 0;

  start_time = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < ITERATIONS; ++i)
    sum += adjacent_box_index_vec[indices[i]][0];

  end_time = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

  std::cout << "(2)Time taken: " << duration << " milliseconds"
	    << " Sum: " << sum << std::endl;
  
  
}
