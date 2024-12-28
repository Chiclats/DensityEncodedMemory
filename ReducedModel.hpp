/*----------
  ver 240912
  ----------*/
#ifndef ReducedModel_hpp
#define ReducedModel_hpp

#include<bits/stdc++.h>
#include"BasicTools.hpp"

using namespace std;

using R1D2VDSys=vector<pair<int,int>>;

struct Reduced1D2ValuedDiscreteSystem{
  
  R1D2VDSys Profile;
  //first: rightward, second: leftward

  int size(){
    return Profile.size();
  }

  function<bool(pair<int,int>)> PolarityEvolutionFunc;
  // pair<bool,bool> first for R-particles, second for L-particles, T rightward and F leftward

  function<R1D2VDSys(R1D2VDSys)> PositionEvolutionFunc;

  void EvolveSystem(){

    R1D2VDSys TempProfile(size(),{0,0});

    //polarity evolution
    for( int i_Site=0 ; i_Site<size() ; i_Site++ ){

      //Rightward particles
      for( int i_Rp=1 ; i_Rp<=Profile[i_Site].first ; i_Rp++ ){
	auto NextPolarity=PolarityEvolutionFunc(Profile[i_Site]);
	if(NextPolarity) TempProfile[i_Site].first++;
	else TempProfile[i_Site].second++;
      }

      //Leftward particles
      for( int i_Lp=1 ; i_Lp<=Profile[i_Site].second ; i_Lp++ ){
	auto NextPolarity=PolarityEvolutionFunc(Profile[i_Site]);
	if(NextPolarity) TempProfile[i_Site].first++;
	else TempProfile[i_Site].second++;
      }
    }

    //position evolution
    Profile=PositionEvolutionFunc(TempProfile);

    return;
  }
  
};

namespace R1D2VDSys_Func{

  function<bool(pair<int,int>)> PolarityEvolution_ReducedVMVecFluc(double Sigma,int i_MP=0){
    return [Sigma,i_MP](pair<int,int> SiteParticles)->bool{
      int RN=SiteParticles.first,LN=SiteParticles.second;
      if((RN-LN)*1.0/(RN+LN)+Sigma*cos(RandGen[i_MP].RandomDouble()*2*pi)>0) return true;
      else return false;
    };
  }

  function<R1D2VDSys(R1D2VDSys)> PositionEvolution_Periodic(int SystemSize){
    return [SystemSize](R1D2VDSys TempProfile)->R1D2VDSys{
      R1D2VDSys AnsProfile(SystemSize,{0,0});

      for( int i=1 ; i<SystemSize-1 ; i++ ){
	AnsProfile[i].first=TempProfile[i-1].first;
	AnsProfile[i].second=TempProfile[i+1].second;
      }

      //for boundary
      AnsProfile[0].first=TempProfile[SystemSize-1].first;
      AnsProfile[0].second=TempProfile[1].second;
      AnsProfile[SystemSize-1].first=TempProfile[SystemSize-2].first;
      AnsProfile[SystemSize-1].second=TempProfile[0].second;

      return AnsProfile;
    };
  }
  
};

#ifdef CheckingHeader
HeaderChecker(ReducedModel_hpp);
#endif //CheckingHeader

#endif //ReducedModel
