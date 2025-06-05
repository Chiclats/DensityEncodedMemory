/*--------------------
ver 250505
--------------------*/

#ifndef VicsekModel_VectorialNoise_hpp
#define VicsekModel_VectorialNoise_hpp

#include<bits/stdc++.h>
#include<omp.h>
#include"BasicTools.hpp"
#include"AgentRealization.hpp"
#include"StateProcessing.hpp"
#include"VicsekModel.hpp"

namespace VicsekModel_2D{

  void Evolve_VicsekModel_VecNoise(vector<Particle>& PG,
				   vector<vector<int>>& MeshedIndexMat,
				   const vector<array<int,4>>& LDAdjBox,
				   const vector<double>& Parameters,// Systemsize_x Systemsize_y Charlength Velocity Sigma Guidingcoeff
				   const function<double(Particle,Particle)>& Distance,
				   const int& i_MP=0)
  {
    double
      SystemSize_X=Parameters[0],
      SystemSize_Y=Parameters[1],
      CharLength=Parameters[2],
      Velocity=Parameters[3],
      Sigma=Parameters[4],
      GuidingCoeff=0;
    if(Parameters.size()==6) GuidingCoeff=Parameters[5];
    
    vector<double> NeighborDirList_x(PG.size(),0),NeighborDirList_y(PG.size(),0);
    vector<int> NeighborNumList(PG.size(),0);// to store numbers of neighbors of each particle

    const int N_x=SystemSize_X/CharLength, N_y=SystemSize_Y/CharLength;
    int i_i,i_j;
    double x_i,y_i,x_j,y_j;
    //calculate the mean velocity dir
    for( int i_Box=0 ; i_Box<LDAdjBox.size() ; i_Box++ )
      for( int i=0 ; i<MeshedIndexMat[i_Box].size() ; i++ ){
	i_i=MeshedIndexMat[i_Box][i];
	x_i=PG[i_i].x;
	y_i=PG[i_i].y;

	//itself is taken into account
	NeighborDirList_x[i_i]+=PG[i_i].vx;
	NeighborDirList_y[i_i]+=PG[i_i].vy;
	NeighborNumList[i_i]++;
	
	//in i_Box
	for( int j=i+1 ; j<MeshedIndexMat[i_Box].size() ; j++ ){
	  i_j=MeshedIndexMat[i_Box][j];
	  x_j=PG[i_j].x;
	  y_j=PG[i_j].y;
	  
	  if((x_i-x_j)*(x_i-x_j)+(y_i-y_j)*(y_i-y_j)<=CharLength*CharLength){
	    NeighborDirList_x[i_i]+=PG[i_j].vx;
	    NeighborDirList_y[i_i]+=PG[i_j].vy;
	    NeighborNumList[i_i]++;
	    NeighborDirList_x[i_j]+=PG[i_i].vx;
	    NeighborDirList_y[i_j]+=PG[i_i].vy;
	    NeighborNumList[i_j]++;
	  }
	}

	//inter-box
	for( int i_AdjBox : LDAdjBox[i_Box] )
	  for( int i_j : MeshedIndexMat[i_AdjBox] )
	    if(Distance(PG[i_i],PG[i_j])<=CharLength*CharLength){
	      NeighborDirList_x[i_i]+=PG[i_j].vx;
	      NeighborDirList_y[i_i]+=PG[i_j].vy;
	      NeighborNumList[i_i]++;
	      NeighborDirList_x[i_j]+=PG[i_i].vx;
	      NeighborDirList_y[i_j]+=PG[i_i].vy;
	      NeighborNumList[i_j]++;
	    }
      }

    //move particles
    double NewDir,RandDir;
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      RandDir=RandGen[i_MP].RandomDouble()*2*pi;//add noise
      NeighborDirList_x[i_Particle]+=Sigma*NeighborNumList[i_Particle]*cos(RandDir);
      NeighborDirList_y[i_Particle]+=Sigma*NeighborNumList[i_Particle]*sin(RandDir);
      
      NewDir=atan2(NeighborDirList_y[i_Particle],NeighborDirList_x[i_Particle]);
      NewDir+=GuidingCoeff*(2*PG[i_Particle].vx*PG[i_Particle].vy);// GuidingCoeff*sin(2*Dir)
      //WARNING: a little difference between the old version: Nematic force is on the old value of Dir
      PG[i_Particle].vx=cos(NewDir);
      PG[i_Particle].vy=sin(NewDir);
      PG[i_Particle].x+=Velocity*PG[i_Particle].vx;
      PG[i_Particle].y+=Velocity*PG[i_Particle].vy;

      //WARNING: for periodic only
      if(PG[i_Particle].x>=SystemSize_X) PG[i_Particle].x-=SystemSize_X;
      if(PG[i_Particle].x<=0) PG[i_Particle].x+=SystemSize_X;
      if(PG[i_Particle].y>=SystemSize_Y) PG[i_Particle].y-=SystemSize_Y;
      if(PG[i_Particle].y<=0) PG[i_Particle].y+=SystemSize_Y;
    }

    //update MeshedIndexMat
    vector<int> Index_MeshedIndexMat_List(PG.size()); // store corresponding box index for each particle
    vector<int> Num_MeshedIndexMat_List(MeshedIndexMat.size(),0); // store particle numbers for each box
    int Index_MeshedIndexMat;
    //calculate corresponding box index
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      int Index_x=PG[i_Particle].x/CharLength,Index_y=PG[i_Particle].y/CharLength;
      if(Index_x>=N_x) Index_x=N_x-1;
      if(Index_y>=N_y) Index_y=N_y-1;
      
      Index_MeshedIndexMat=Index_x+N_x*Index_y;
      Index_MeshedIndexMat_List[i_Particle]=Index_MeshedIndexMat;
      Num_MeshedIndexMat_List[Index_MeshedIndexMat]++;
    }

    vector<vector<int>> NewMeshedIndexMat(MeshedIndexMat.size());
    //resize the Mat
    for( int i_Box=0 ; i_Box<MeshedIndexMat.size() ; i_Box++ )
      NewMeshedIndexMat[i_Box].resize(Num_MeshedIndexMat_List[i_Box]);
    //put particle index into the Mat
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      Index_MeshedIndexMat=Index_MeshedIndexMat_List[i_Particle];
      Num_MeshedIndexMat_List[Index_MeshedIndexMat]--;
      NewMeshedIndexMat[Index_MeshedIndexMat][Num_MeshedIndexMat_List[Index_MeshedIndexMat]]=i_Particle;
    }
    
    MeshedIndexMat=NewMeshedIndexMat;
    
  }

#ifndef LocalProgram

  void Evolve_VicsekModel_VecNoise_MP(vector<Particle>& PG,
				      vector<vector<int>>& MeshedIndexMat,
				      const vector<array<int,8>>& AllAdjBox,
				      const vector<double>& Parameters,// Systemsize_x Systemsize_y Charlength Velocity Sigma Guidingcoeff
				      const function<double(Particle,Particle)>& Distance)
  {
    double
      SystemSize_X=Parameters[0],
      SystemSize_Y=Parameters[1],
      CharLength=Parameters[2],
      Velocity=Parameters[3],
      Sigma=Parameters[4],
      GuidingCoeff=0;
    if(Parameters.size()==6) GuidingCoeff=Parameters[5];
    
    vector<double> NeighborDirList_x(PG.size(),0),NeighborDirList_y(PG.size(),0);
    vector<int> NeighborNumList(PG.size(),0);

    //calculate the mean velocity dir
#pragma omp parallel for
    for( int i_Box=0 ; i_Box<AllAdjBox.size() ; i_Box++ ){

      int i_i,i_j;
      double x_i,y_i,x_j,y_j;
      
      for( int i=0 ; i<MeshedIndexMat[i_Box].size() ; i++ ){
	i_i=MeshedIndexMat[i_Box][i];
	x_i=PG[i_i].x;
	y_i=PG[i_i].y;

	//itself is taken into account
	NeighborDirList_x[i_i]+=PG[i_i].vx;
	NeighborDirList_y[i_i]+=PG[i_i].vy;
	NeighborNumList[i_i]++;
	
	//in i_Box
	for( int j=i+1 ; j<MeshedIndexMat[i_Box].size() ; j++ ){
	  i_j=MeshedIndexMat[i_Box][j];
	  x_j=PG[i_j].x;
	  y_j=PG[i_j].y;
	  
	  if((x_i-x_j)*(x_i-x_j)+(y_i-y_j)*(y_i-y_j)<=CharLength*CharLength){
	    NeighborDirList_x[i_i]+=PG[i_j].vx;
	    NeighborDirList_y[i_i]+=PG[i_j].vy;
	    NeighborNumList[i_i]++;
	    NeighborDirList_x[i_j]+=PG[i_i].vx;
	    NeighborDirList_y[i_j]+=PG[i_i].vy;
	    NeighborNumList[i_j]++;
	  }
	}

	//inter-box
	for( int i_AdjBox : AllAdjBox[i_Box] )
	  for( int i_j : MeshedIndexMat[i_AdjBox] )
	    if(Distance(PG[i_i],PG[i_j])<=CharLength*CharLength){
	      NeighborDirList_x[i_i]+=PG[i_j].vx;
	      NeighborDirList_y[i_i]+=PG[i_j].vy;
	      NeighborNumList[i_i]++;
	    }
      }
    }

    //move particles
#pragma omp parallel for
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      double RandDir=RandGen[omp_get_thread_num()].RandomDouble()*2*pi;
      NeighborDirList_x[i_Particle]+=Sigma*NeighborNumList[i_Particle]*cos(RandDir);
      NeighborDirList_y[i_Particle]+=Sigma*NeighborNumList[i_Particle]*sin(RandDir);// add vectorial noise
      
      double NewDir=atan2(NeighborDirList_y[i_Particle],NeighborDirList_x[i_Particle]);
      NewDir+=GuidingCoeff*(2*PG[i_Particle].vx*PG[i_Particle].vy);// GuidingCoeff*sin(2*Dir)
      //WARNING: a little difference between the old version: Nematic force is on the old value of Dir
      PG[i_Particle].vx=cos(NewDir);
      PG[i_Particle].vy=sin(NewDir);
      PG[i_Particle].x+=Velocity*PG[i_Particle].vx;
      PG[i_Particle].y+=Velocity*PG[i_Particle].vy;

      //WARNING: for periodic only
      if(PG[i_Particle].x>=SystemSize_X) PG[i_Particle].x-=SystemSize_X;
      if(PG[i_Particle].x<=0) PG[i_Particle].x+=SystemSize_X;
      if(PG[i_Particle].y>=SystemSize_Y) PG[i_Particle].y-=SystemSize_Y;
      if(PG[i_Particle].y<=0) PG[i_Particle].y+=SystemSize_Y;
    }

    //update MeshedIndexMat
    vector<int> Index_MeshedIndexMat_List(PG.size()); // store corresponding box index for each particle
    vector<int> Num_MeshedIndexMat_List(MeshedIndexMat.size(),0); // store particle numbers for each box
    int Index_MeshedIndexMat;
    const int N_x=SystemSize_X/CharLength, N_y=SystemSize_Y/CharLength;
    //calculate corresponding box index
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      int Index_x=PG[i_Particle].x/CharLength,Index_y=PG[i_Particle].y/CharLength;
      if(Index_x>=N_x) Index_x=N_x-1;
      if(Index_y>=N_y) Index_y=N_y-1;
      
      Index_MeshedIndexMat=Index_x+N_x*Index_y;
      Index_MeshedIndexMat_List[i_Particle]=Index_MeshedIndexMat;
      Num_MeshedIndexMat_List[Index_MeshedIndexMat]++;
    }

    vector<vector<int>> NewMeshedIndexMat(MeshedIndexMat.size());
    //resize the Mat
    for( int i_Box=0 ; i_Box<MeshedIndexMat.size() ; i_Box++ )
      NewMeshedIndexMat[i_Box].resize(Num_MeshedIndexMat_List[i_Box]);
    //put particle index into the Mat
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      Index_MeshedIndexMat=Index_MeshedIndexMat_List[i_Particle];
      Num_MeshedIndexMat_List[Index_MeshedIndexMat]--;
      NewMeshedIndexMat[Index_MeshedIndexMat][Num_MeshedIndexMat_List[Index_MeshedIndexMat]]=i_Particle;
    }
    
    MeshedIndexMat=NewMeshedIndexMat;
    
  }
  
#endif

  PolarParticle2D_GT Evolve_VicsekModel_VecNoise_Periodic(const PolarParticle2D_GT& PG,
							  const BoxInfoT& BoxInfo,
							  const vector<double>& Parameters,
							  const int& StepNumber,
							  const int& i_MP=0)
  {//an integrated version of Evolve_VicsekModel, evolving for StepNumber steps
    const double
      SystemSize_X=Parameters[0],
      SystemSize_Y=Parameters[1],
      CharLength=Parameters[2];
    
    vector<Particle> vP=PolarParticle2D_GT_To_vParticle(PG);
    vector<vector<int>> MeshedIndex=MeshingIndex(PG,SystemSize_X,SystemSize_Y,CharLength);
    const vector<array<int,4>> LDAdjBox=LowerRightAdjacentBoxIndex(BoxInfo);
    const function<double(Particle,Particle)> Distance_=Distance::Periodic(SystemSize_X,SystemSize_Y,CharLength);

    for( int i_Step=0 ; i_Step<StepNumber ; i_Step++ )
      Evolve_VicsekModel_VecNoise(vP,MeshedIndex,LDAdjBox,Parameters,Distance_,i_MP);

    auto ansPG=PG;
    UpdatePGFromVParticle(ansPG,vP);
    return ansPG;
  }

#ifndef LocalProgram
  PolarParticle2D_GT Evolve_VicsekModel_VecNoise_Periodic_MP(const PolarParticle2D_GT& PG,
							     const BoxInfoT& BoxInfo,
							     const vector<double>& Parameters,
							     const int& StepNumber)
  {//an integrated version of Evolve_VicsekModel, evolving for StepNumber steps
    const double
      SystemSize_X=Parameters[0],
      SystemSize_Y=Parameters[1],
      CharLength=Parameters[2];
    
    vector<Particle> vP=PolarParticle2D_GT_To_vParticle(PG);
    vector<vector<int>> MeshedIndex=MeshingIndex(PG,SystemSize_X,SystemSize_Y,CharLength);
    const vector<array<int,8>> AllAdjBox=AllAdjacentBoxIndex(BoxInfo);
    const function<double(Particle,Particle)> Distance_=Distance::Periodic(SystemSize_X,SystemSize_Y,CharLength);

    for( int i_Step=0 ; i_Step<StepNumber ; i_Step++ )
      Evolve_VicsekModel_VecNoise_MP(vP,MeshedIndex,AllAdjBox,Parameters,Distance_);

    auto ansPG=PG;
    UpdatePGFromVParticle(ansPG,vP);
    return ansPG;
  }
#endif

}
			
#endif//VicsekModel_VectorialNoise_hpp
