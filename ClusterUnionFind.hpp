/*--------------------
ver 250417
--------------------*/

#ifndef ClusterUnionFind_hpp
#define ClusterUnionFind_hpp

#include<bits/stdc++.h>
#include"BasicTools.hpp"
#include"AgentRealization.hpp"
#include"StateProcessing.hpp"
#include"DataProcessing.hpp"

#ifndef LocalProgram
  #include<omp.h>
#endif//LocalProgram

namespace DataProcessing::Cluster{
  
  template<typename Type>
  vector<int> ClusterDivide(const vector<Type>& PG,const double& R,
			    const double& SystemSize_X,const double& SystemSize_Y,
			    const function<double(Type,Type)>& Distance)
  {// 

    // Step 1: grids meshing
    const auto MeshedIndexMat = BoxMeshing_2D::MeshedIndex(PG, SystemSize_X, SystemSize_Y, R);
    const BoxInfoT BoxInfo = BoxMeshing_2D::BoxInfo(SystemSize_X,SystemSize_Y,R);
    int Num_x = floor(SystemSize_X / R);
    int Num_y = floor(SystemSize_Y / R);

    // Step 2: initialize union find
    vector<int> parent(PG.size());
    vector<int> rank(PG.size(), 1);
    for (int i = 0; i < PG.size(); ++i) parent[i] = i;

    // path compression and find
    auto find = [&](int x) {
      while (parent[x] != x) {
	parent[x] = parent[parent[x]];
	x = parent[x];
      }
      return x;
    };

    // unite
    auto unite = [&](int x, int y) {
      int rx = find(x), ry = find(y);
      if (rx == ry) return;
      if (rank[rx] > rank[ry]) swap(rx, ry);
      parent[rx] = ry;
      if (rank[rx] == rank[ry]) rank[ry]++;
    };

    // Step 3: deal with possible pairs
    for (int box_idx = 0; box_idx < MeshedIndexMat.size(); ++box_idx) {
      const auto& curr = MeshedIndexMat[box_idx];

      // in the same box (i < j)
      for (int i = 0; i < curr.size(); ++i)
	for (int j = i+1; j < curr.size(); ++j)
	  if (Distance(PG[curr[i]], PG[curr[j]]) <= R)
	    unite(curr[i], curr[j]);

      const array<int,4> LDAdjBox={
	BoxInfo[box_idx].second[0],
	BoxInfo[box_idx].second[2],
	BoxInfo[BoxInfo[box_idx].second[2]].second[0],
	BoxInfo[BoxInfo[box_idx].second[2]].second[1]
      };

      // adjacent boxes
      for (auto adj_idx : LDAdjBox) {
	const auto& adj = MeshedIndexMat[adj_idx];

	for (int i : curr)
	  for (int j : adj)
	    if (Distance(PG[i], PG[j]) <= R)
	      unite(i, j);
      }
    }

    // Step 4: generate continuous cluster index
    vector<int> clusters(PG.size());
    unordered_map<int, int> root_map;
    int cluster_id = 1;
    
    for (int i = 0; i < PG.size(); ++i) {
      int root = find(i);
      if (!root_map.count(root)) {
	root_map[root] = cluster_id;
	cluster_id++;
      }
      clusters[i] = root_map[root];
    }

    return clusters;    
  }

  
  


};

#endif//ClusterUnionFind_hpp
