/*--------------------
ver 260224
--------------------*/

#ifndef StateProcessing_new_hpp
#define StateProcessing_new_hpp

namespace BoxMeshing_2D{

  using BoxShape=std::tuple<double,double,double,unsigned int,unsigned int,unsigned int>;
  // to store the shape of the box, (system_size_x, system_size_y, char_length, box_num_x, box_num_y, total_box_num);

  BoxShape create_box_shape(const double system_size_x,const double system_size_y,const double char_length)
  {// to create a BoxShape based on the system size and characteristic length
    // last box may be larger than char_length
    unsigned int box_num_x = static_cast<unsigned int>(floor(system_size_x / char_length));
    unsigned int box_num_y = static_cast<unsigned int>(floor(system_size_y / char_length));
    unsigned int total_box_num = box_num_x * box_num_y;
    return std::make_tuple(system_size_x, system_size_y, char_length, box_num_x, box_num_y, total_box_num);
  }

  unsigned int decide_box_index(const double x,const double y,const BoxShape& box_shape)
  {// to decide the box index for a particle at (x,y) based on the box shape
    auto [system_size_x, system_size_y, char_length, box_num_x, box_num_y, _] = box_shape;
    if(x < 0 || x >= system_size_x || y < 0 || y >= system_size_y)
      throw std::runtime_error("In BoxMeshing_2D::decide_box_index : the coordinate of particle should not exceed the system size.");
		
    unsigned int index_x = static_cast<unsigned int>(floor(x / char_length));
    unsigned int index_y = static_cast<unsigned int>(floor(y / char_length));
    if(index_x >= box_num_x) index_x = box_num_x - 1;
    if(index_y >= box_num_y) index_y = box_num_y - 1;
    
    return index_x + box_num_x * index_y;
  }

  unsigned int decide_box_neighbor_index(const unsigned int box_index,const int offset_x,const int offset_y,const BoxShape& box_shape)
  {// to decide the neighbor box index based on the box index and the offset (offset_x, offset_y)
    // WARNING: in order to decide the neighbor box index, calculating directly is faster than reading from a table
    auto [system_size_x, system_size_y, char_length, box_num_x, box_num_y, _] = box_shape;
    
    unsigned int index_x = box_index % box_num_x;
    unsigned int index_y = box_index / box_num_x;
    
    int neighbor_index_x = static_cast<int>(index_x) + offset_x;
    int neighbor_index_y = static_cast<int>(index_y) + offset_y;
    
    if(neighbor_index_x < 0) neighbor_index_x += box_num_x;
    if(neighbor_index_x >= static_cast<int>(box_num_x)) neighbor_index_x -= box_num_x;
    if(neighbor_index_y < 0) neighbor_index_y += box_num_y;
    if(neighbor_index_y >= static_cast<int>(box_num_y)) neighbor_index_y -= box_num_y;
    
    return neighbor_index_x + box_num_x * neighbor_index_y;
  }

  unsigned int decide_box_neighbor_index_ordered(const unsigned int box_index,const int nbr_index, const BoxShape& box_shape)
  {// to decide the neighbor box index based on the box index and the neighbor index (0-7)
    // neighbor index order: 0-bottom-left, 1-bottom, 2-bottom-right, 3-left, 4-right, 5-top-left, 6-top, 7-top-right
    static const std::array<std::pair<int,int>,8> offset_array = {{
	{-1,-1}, {0,-1}, {1,-1},
	{-1,0},           {1,0},
	{-1,1},  {0,1},  {1,1}
      }};
		
    if(nbr_index < 0 || nbr_index >= 8)
      throw std::runtime_error("In BoxMeshing_2D::decide_box_neighbor_index_ordered : neighbor index should be between 0 and 7.");
		
    const auto [offset_x, offset_y] = offset_array[nbr_index];
    return decide_box_neighbor_index(box_index, offset_x, offset_y, box_shape);
  }

  std::vector<std::array<unsigned int, 8>> all_adjacent_box_index(const BoxShape& box_shape)
  {// to return the indices of all adjacent boxes for each box
    // WARNING: in order to decide the neighbor box index, calculating directly is faster than reading from a table
    auto [_, _, _, box_num_x, box_num_y, total_box_num] = box_shape;
    std::vector<std::array<unsigned int, 8>> ans_vec(total_box_num);
		
    for(unsigned int i_box = 0; i_box < total_box_num; i_box++) {
      ans_vec[i_box] = {
	decide_box_neighbor_index(i_box, -1, -1, box_shape), // bottom-left
	decide_box_neighbor_index(i_box, 0, -1, box_shape),  // bottom
	decide_box_neighbor_index(i_box, 1, -1, box_shape),  // bottom-right
	decide_box_neighbor_index(i_box, -1, 0, box_shape),  // left
	decide_box_neighbor_index(i_box, 1, 0, box_shape),   // right
	decide_box_neighbor_index(i_box, -1, 1, box_shape),  // top-left
	decide_box_neighbor_index(i_box, 0, 1, box_shape),   // top
	decide_box_neighbor_index(i_box, 1, 1, box_shape)    // top-right
      };
    }
		
    return ans_vec;
  }

}

#endif//StateProcessing_new_hpp
