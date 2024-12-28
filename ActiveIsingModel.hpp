/*----------
  ver 240719
  ----------*/
#ifndef ActiveIsingModel_hpp
#define ActiveIsingModel_hpp

#include<bits/stdc++.h>
#include"BasicTools.hpp"

using namespace std;

struct ActiveIsingModel2D{

private:
  using ProfileT=vector<vector<pair<int,int>>>;
  using CoorT=pair<int,int>;

public:
  ProfileT Profile;

  double Bias; //i.e. epsilon
  double Diffusion; //i.e. D
  double Beta;
  vector<vector<vector<CoorT>>> NeighborTable;//similar as BoxInfo, {LCoor,RCoor,DCoor,UCoor}
  
  double TimeStep;

  void EvolveTimeStep()
  {//Evolve approximately with discrete time
    ProfileT AnsProfile=Profile;
    int SystemSize_X=Profile.size(),SystemSize_Y=Profile[0].size();

    int NR,NL,NR_RMove,NR_LMove,NR_UMove,NR_DMove,NR_Flip,
      NL_RMove,NL_LMove,NL_UMove,NL_DMove,NL_Flip;
    double p_D,pDP,pDA,pFR,pFL;

    p_D=TimeStep*Diffusion;
    pDP=TimeStep*Diffusion*(1+Bias); //go along the polarity
    pDA=TimeStep*Diffusion*(1-Bias); //go against the polarity

    for( int iX=0 ; iX<SystemSize_X ; iX++ )
      for( int iY=0 ; iY<SystemSize_Y ; iY++ ){
	if(Profile[iX][iY].first==0 and Profile[iX][iY].second==0) continue;

	NR=Profile[iX][iY].first;
	NL=Profile[iX][iY].second;

	pFR=TimeStep*exp(-Beta*(NR-NL)/(NR+NL));  //flip R -> L
	pFL=TimeStep*exp(Beta*(NR-NL)/(NR+NL));   //flip L -> R
	NR_RMove=RandGen[0].RandomBinomialInt(NR,pDP);
	NR_LMove=RandGen[0].RandomBinomialInt(NR,pDA);
	NR_UMove=RandGen[0].RandomBinomialInt(NR,p_D);
	NR_DMove=RandGen[0].RandomBinomialInt(NR,p_D);
	NR_Flip=RandGen[0].RandomBinomialInt(NR,pFR);
	NL_RMove=RandGen[0].RandomBinomialInt(NL,pDA);
	NL_LMove=RandGen[0].RandomBinomialInt(NL,pDP);
	NL_UMove=RandGen[0].RandomBinomialInt(NL,p_D);
	NL_DMove=RandGen[0].RandomBinomialInt(NL,p_D);
	NL_Flip=RandGen[0].RandomBinomialInt(NL,pFL);

	// cout<<iX<<','<<iY<<' '<<NR<<' '<<NL<<endl;
	// cout<<p_D<<' '<<pDP<<' '<<pDA<<' '<<pFR<<' '<<pFL<<endl;//
	// cout<<NR_RMove<<' '<<NR_LMove<<' '<<NR_UMove<<' '<<NR_DMove<<' '<<NR_Flip<<' '
	//     <<NL_RMove<<' '<<NL_LMove<<' '<<NL_UMove<<' '<<NL_DMove<<' '<<NL_Flip<<endl;//

	// if(NR_RMove+NR_LMove+NL_UMove+NL_DMove+NL_RMove+NL_LMove+NL_UMove+NL_DMove>0) cout<<"changed!\n";//
	
	if(NR_RMove+NR_LMove+NR_UMove+NR_DMove+NR_Flip>NR or NL_RMove+NL_LMove+NL_UMove+NL_DMove+NL_Flip>NL)
	{//to avoid that the changed number exceed the total number 
	  iY--;
	  //cout<<"Warning "<<iX<<' '<<iY<<endl;//
	  continue;
	}

	AnsProfile[iX][iY].first-=(NR_RMove+NR_LMove+NR_UMove+NR_DMove+NR_Flip);
	AnsProfile[iX][iY].second-=(NL_RMove+NL_LMove+NL_UMove+NL_DMove+NL_Flip);
	AnsProfile[iX][iY].first+=NL_Flip;
	AnsProfile[iX][iY].second+=NR_Flip;
	AnsProfile[NeighborTable[iX][iY][0].first][NeighborTable[iX][iY][0].second].first+=NR_LMove;
	AnsProfile[NeighborTable[iX][iY][1].first][NeighborTable[iX][iY][1].second].first+=NR_RMove;
	AnsProfile[NeighborTable[iX][iY][2].first][NeighborTable[iX][iY][2].second].first+=NR_DMove;
	AnsProfile[NeighborTable[iX][iY][3].first][NeighborTable[iX][iY][3].second].first+=NR_UMove;
	AnsProfile[NeighborTable[iX][iY][0].first][NeighborTable[iX][iY][0].second].second+=NL_LMove;
	AnsProfile[NeighborTable[iX][iY][1].first][NeighborTable[iX][iY][1].second].second+=NL_RMove;
	AnsProfile[NeighborTable[iX][iY][2].first][NeighborTable[iX][iY][2].second].second+=NL_DMove;
	AnsProfile[NeighborTable[iX][iY][3].first][NeighborTable[iX][iY][3].second].second+=NL_UMove;
      }

    Profile=AnsProfile;
  }
  
  //generating profile-------------------
  
  static ProfileT GenerateHomogeneousRandomProfile(int SystemSize_X,int SystemSize_Y,int N)
  {//randomly distributed with random polarity
    ProfileT AnsMat(SystemSize_X,vector<pair<int,int>>(SystemSize_Y,{0,0}));
    for( int i=0 ; i<N ; i++ ){
      int
	ti1=RandGen[0].RandomInt(0,SystemSize_X-1),
	ti2=RandGen[0].RandomInt(0,SystemSize_Y-1),
	ti3=RandGen[0].RandomInt(0,1);
      if(ti3) AnsMat[ti1][ti2].first++;
      else AnsMat[ti1][ti2].second++;
    }

    return AnsMat;
  }

  //NeighborTable-------------------
  
  static vector<vector<vector<CoorT>>> NeighborTable_Periodic(int SystemSize_X,int SystemSize_Y)
  {
    vector<vector<vector<CoorT>>> AnsMat(SystemSize_X,vector<vector<CoorT>>(SystemSize_Y));

    for( int iX=0 ; iX<SystemSize_X ; iX++ ){
      for( int iY=0 ; iY<SystemSize_Y ; iY++ ){
	AnsMat[iX][iY]={{iX-1,iY},{iX+1,iY},{iX,iY-1},{iX,iY+1}};
	if(iX==0) AnsMat[iX][iY][0].first=SystemSize_X-1;
	if(iX==SystemSize_X-1) AnsMat[iX][iY][1].first=0;
	if(iY==0) AnsMat[iX][iY][2].second=SystemSize_Y-1;
	if(iY==SystemSize_Y-1) AnsMat[iX][iY][3].second=0;
      }
    }

    return AnsMat;
  }

};

struct ActiveIsingModel1D{

private:
  using ProfileT=vector<pair<int,int>>;

public:
  ProfileT Profile;

  double Bias;
  double Diffusion;
  double Beta;
  vector<vector<int>> NeighborTable; //{{iL,iR},...}

  double TimeStep;

  void EvolveTimeStep()
  {//Evolve approximately with discrete time
    ProfileT AnsProfile=Profile;
    int SystemSize_X=Profile.size();

    int NR,NL,NR_RMove,NR_LMove,NR_UMove,NR_DMove,NR_Flip,
      NL_RMove,NL_LMove,NL_UMove,NL_DMove,NL_Flip;
    double
      pDP=TimeStep*Diffusion*(1+Bias), //go along the polarity
      pDA=TimeStep*Diffusion*(1-Bias), //go against the polarity
      pFR,pFL;

    for( int iX=0 ; iX<SystemSize_X ; iX++ ){
      if(Profile[iX].first==0 and Profile[iX].second==0) continue;

      NR=Profile[iX].first;
      NL=Profile[iX].second;

      pFR=TimeStep*exp(-Beta*(NR-NL)/(NR+NL));  //flip R -> L
      pFL=TimeStep*exp(Beta*(NR-NL)/(NR+NL));   //flip L -> R
      NR_RMove=RandGen[0].RandomBinomialInt(NR,pDP);
      NR_LMove=RandGen[0].RandomBinomialInt(NR,pDA);
      NR_Flip=RandGen[0].RandomBinomialInt(NR,pFR);
      NL_RMove=RandGen[0].RandomBinomialInt(NL,pDA);
      NL_LMove=RandGen[0].RandomBinomialInt(NL,pDP);
      NL_Flip=RandGen[0].RandomBinomialInt(NL,pFL);
	
      if(NR_RMove+NR_LMove+NR_Flip>NR or NL_RMove+NL_LMove+NL_Flip>NL)
	{//to avoid that the changed number exceed the total number 
	  iX--;
	  continue;
	}

      AnsProfile[iX].first-=(NR_RMove+NR_LMove+NR_Flip);
      AnsProfile[iX].second-=(NL_RMove+NL_LMove+NL_Flip);
      AnsProfile[iX].first+=NL_Flip;
      AnsProfile[iX].second+=NR_Flip;
      AnsProfile[NeighborTable[iX][0]].first+=NR_LMove;
      AnsProfile[NeighborTable[iX][1]].first+=NR_RMove;
      AnsProfile[NeighborTable[iX][0]].second+=NL_LMove;
      AnsProfile[NeighborTable[iX][1]].second+=NL_RMove;
    }

    Profile=AnsProfile;
  }

  //generating profile-------------------

  static ProfileT GenerateHomogeneousRandomProfile(int SystemSize_X,int N)
  {//randomly distributed with random polarity
    ProfileT AnsMat(SystemSize_X,{0,0});
    for( int i=0 ; i<N ; i++ ){
      int
	ti1=RandGen[0].RandomInt(0,SystemSize_X-1),
	ti3=RandGen[0].RandomInt(0,1);
      if(ti3) AnsMat[ti1].first++;
      else AnsMat[ti1].second++;
    }

    return AnsMat;
  }

  //neighbortable-------------------

  static vector<vector<int>> NeighborTable_Periodic(int SystemSize_X)
  {
    vector<vector<int>> AnsList;
    for( int i=0 ; i<SystemSize_X ; i++ ){
      AnsList.push_back({i-1,i+1});
    }
    AnsList[0][0]=SystemSize_X-1;
    AnsList[SystemSize_X-1][1]=0;
    
    return AnsList;
  }
  
};

namespace FileIO{

  void OutAIM2D(const string& Filename,const ActiveIsingModel2D& Sys){
    OutMatPair(Filename,Sys.Profile);
  }
  
  void IntoAIM2D(string Filename,ActiveIsingModel2D& Sys,int rows,int cols){
    Sys.Profile=vector<vector<pair<int,int>>>(rows,vector<pair<int,int>>(cols));
    
    int Rt,Lt;
    
    ifstream fin;
    fin.open(Filename);
    
    for( int i=0 ; i<rows ; i++)
      for( int j=0 ; j<cols ; j++){
	fin>>Rt>>Lt;
	Sys.Profile[i][j]={Rt,Lt};
      }
    
    fin.close();
  }
  
};

namespace DirectlyOutput{

  void AIM2D(ActiveIsingModel2D Sys){
    auto Profile=Sys.Profile;
    for( auto v : Profile ){
      for( auto p : v )
	cout<<'('<<p.first<<','<<p.second<<")\t";
      cout<<endl;
    }
  }

};

#ifdef CheckingHeader
HeaderChecker(ActiveIsingModel_hpp);
#endif //CheckingHeader

#endif //ActiveIsingModel_hpp
