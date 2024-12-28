/*--------------------
  ver. 240802
  --------------------*/

#ifndef ExtendedCollectiveMotionRealization_hpp
#define ExtendedCollectiveMotionRealization_hpp

#include<bits/stdc++.h>
#include"BasicTools.hpp"
#include"StateProcessing.hpp"

using namespace std;

template<typename AgentT>
using EvolutionT=function<AgentT(AgentT,vector<AgentT>)>;

template<typename AgentT>
using DistanceExtendedT=function<double(AgentT,AgentT)>;

template<typename AgentT>
using BoundaryExtendedT=function<AgentT(AgentT)>;

template<typename AgentT>
vector<vector<AgentT>> EvolveMeshedShortRangeSystem(const vector<vector<AgentT>>& MeshedSystem,
						    const BoxInfoT& BoxInfo,
						    const EvolutionT<AgentT>& Evolution,
						    const BoundaryExtendedT<AgentT>& Boundary,
						    const DistanceExtendedT<AgentT>& Distance)
{/*Extended ver of EvolveMeshedList*/
  auto [SystemSize_X,SystemSize_Y,CharLength,BoxNum_X,BoxNum_Y]=BoxMeshing_2D::UnzipBoxInfo(BoxInfo);
  vector<vector<AgentT>> AnsSystem;

  for( int i_Box=0 ; i_Box<BoxInfo.size() ; i_Box++ ){
    const array<int,8> AdjacentBoxList={
      BoxInfo[i_Box].second[0],BoxInfo[i_Box].second[1],
      BoxInfo[i_Box].second[2],BoxInfo[i_Box].second[3],
      BoxInfo[BoxInfo[i_Box].second[0]].second[2],
      BoxInfo[BoxInfo[i_Box].second[0]].second[3],
      BoxInfo[BoxInfo[i_Box].second[1]].second[2],
      BoxInfo[BoxInfo[i_Box].second[1]].second[3],
    };
    //[Warning] Assuming the box matrix has more than 3 rows and 3 collumns.

    for( int ip=0 ; ip<MeshedSystem[i_Box].size() ; ip++ ){
      AgentT p=MeshedSystem[i_Box][ip];
      
      vector<AgentT> Neighbors;
      
      for( auto i_AdjBox : AdjacentBoxList )
	for( AgentT np : MeshedSystem[i_AdjBox] )
	  if(Distance(p,np)<CharLength)
	    Neighbors.push_back(np);

      for( AgentT np : MeshedSystem[i_Box] )
	if(Distance(p,np)<CharLength and p!=np)//[Warning] Class AgentT should have !=.
	  Neighbors.push_back(np);

      AgentT ap=Boundary(Evolution(p,Neighbors));

      AnsSystem[BoxMeshing_2D::DecideBoxIndex(ap,SystemSize_X,SystemSize_Y,CharLength)].push_back(ap);
      //[Warning] Class AgentT should have attribute complex<double> Pos.
    }
  }

  return AnsSystem;
}

#ifdef CheckingHeader
HeaderChecker(ExtendedCollectiveMotionRealization_hpp);
#endif //CheckingHeader

#endif//Extendedcollectivemotionrealization_hpp
