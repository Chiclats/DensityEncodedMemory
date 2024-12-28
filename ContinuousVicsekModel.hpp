/*--------------------
ver 241103
--------------------*/

#ifndef ContinuousVicsekModel_hpp
#define ContinuousVicsekModel_hpp

#include<bits/stdc++.h>
#include<omp.h>
#include"BasicTools.hpp"
#include"AgentRealization.hpp"
#include"StateProcessing.hpp"
#include"VicsekModel.hpp"

using namespace std;

namespace VicsekModel_2D
{//Continuous Pairwise Vicsek Model

  int TestContinuity;

  void Evolve_ContinuousVicsekModel(vector<Particle>& PG,
				    vector<vector<int>>& MeshedIndexMat,
				    const vector<array<int,4>>& LDAdjBox,
				    const vector<double>& Parameters,
				    // SystemSize_X SystemSize_Y Charlength Velocity AligningCoeff Sigma GuidingCoeff TimeInterval SqrtTimeinterval
				    const function<double(Particle,Particle)>& Distance,
				    const int& i_MP=0)
  {
    const double
      SystemSize_X=Parameters[0],
      SystemSize_Y=Parameters[1],
      CharLength=Parameters[2],
      Velocity=Parameters[3],
      AligningCoeff=Parameters[4],
      Sigma=Parameters[5],
      GuidingCoeff=Parameters[6],
      TimeInterval=Parameters[7],
      SqrtTimeInterval=Parameters[8];
    
    vector<double> NeighborDirList_x(PG.size(),0),NeighborDirList_y(PG.size(),0);
    vector<double> DirIncrementList(PG.size(),0);

    int i_i,i_j;
    double x_i,y_i,x_j,y_j;
    double Increment;
    //calculate the mean velocity dir
    for( int i_Box=0 ; i_Box<LDAdjBox.size() ; i_Box++ )
      for( int i=0 ; i<MeshedIndexMat[i_Box].size() ; i++ ){
	i_i=MeshedIndexMat[i_Box][i];
	x_i=PG[i_i].x;
	y_i=PG[i_i].y;

	//itself is taken into account
	NeighborDirList_x[i_i]+=PG[i_i].vx;
	NeighborDirList_y[i_i]+=PG[i_i].vy;
	
	//in i_Box
	for( int j=i+1 ; j<MeshedIndexMat[i_Box].size() ; j++ ){
	  i_j=MeshedIndexMat[i_Box][j];
	  x_j=PG[i_j].x;
	  y_j=PG[i_j].y;
	  
	  if((x_i-x_j)*(x_i-x_j)+(y_i-y_j)*(y_i-y_j)<=CharLength*CharLength){
	    NeighborDirList_x[i_i]+=PG[i_j].vx;
	    NeighborDirList_y[i_i]+=PG[i_j].vy;
	    NeighborDirList_x[i_j]+=PG[i_i].vx;
	    NeighborDirList_y[i_j]+=PG[i_i].vy;

	    Increment=PG[i_j].vy*PG[i_i].vx-PG[i_j].vx*PG[i_i].vy;
	    DirIncrementList[i_i]+=Increment;
	    DirIncrementList[i_j]+=(-Increment);
	  }
	}

	//inter-box
	for( int i_AdjBox : LDAdjBox[i_Box] )
	  for( int i_j : MeshedIndexMat[i_AdjBox] )
	    if(Distance(PG[i_i],PG[i_j])<=CharLength*CharLength){
	      NeighborDirList_x[i_i]+=PG[i_j].vx;
	      NeighborDirList_y[i_i]+=PG[i_j].vy;
	      NeighborDirList_x[i_j]+=PG[i_i].vx;
	      NeighborDirList_y[i_j]+=PG[i_i].vy;

	      Increment=PG[i_j].vy*PG[i_i].vx-PG[i_j].vx*PG[i_i].vy;
	      DirIncrementList[i_i]+=Increment;
	      DirIncrementList[i_j]+=(-Increment);
	    }
      }

    //move particles
    double OldDir,NewDir,NewDirLim,IncrementLim; //NewDirLim is to store the mean direction. NewDir will not exceed NewDirLim
    vector<vector<int>> NewMeshedIndexMat(MeshedIndexMat.size());
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      //ensure that NewDir doesn't exceed NewDirLim
      NewDirLim=atan2(NeighborDirList_y[i_Particle],NeighborDirList_x[i_Particle]);
      OldDir=atan2(PG[i_Particle].vy,PG[i_Particle].vx);
      Increment=AligningCoeff*DirIncrementList[i_Particle]*TimeInterval; //Increment reused
      IncrementLim=NewDirLim-OldDir;
      if(Increment>0 and IncrementLim<0) IncrementLim+=2*pi;
      if(0>Increment and 0<IncrementLim) IncrementLim-=2*pi;

      if(Increment>pi/2 or Increment<-pi/2)
	throw runtime_error("Too large AligningCoeff / TimeInterval / particle density causes increment out of bound!");
      if(abs(Increment)>abs(IncrementLim))
	//NewDir=NewDirLim;
	NewDir=NewDirLim, TestContinuity++;
      else NewDir=OldDir+Increment;
      
      NewDir+=Sigma*RandGen[i_MP].GaussianRandomDouble()*SqrtTimeInterval; 
      NewDir+=GuidingCoeff*(2*PG[i_Particle].vx*PG[i_Particle].vy)*TimeInterval;// GuidingCoeff*sin(2*Dir)
      //WARNING: a little difference between the old version: Nematic force is on the old value of Dir
      PG[i_Particle].vx=cos(NewDir);
      PG[i_Particle].vy=sin(NewDir);
      PG[i_Particle].x+=Velocity*PG[i_Particle].vx*TimeInterval;
      PG[i_Particle].y+=Velocity*PG[i_Particle].vy*TimeInterval;

      //WARNING: for periodic only
      if(PG[i_Particle].x>=SystemSize_X) PG[i_Particle].x-=SystemSize_X;
      if(PG[i_Particle].x<=0) PG[i_Particle].x+=SystemSize_X;
      if(PG[i_Particle].y>=SystemSize_Y) PG[i_Particle].y-=SystemSize_Y;
      if(PG[i_Particle].y<=0) PG[i_Particle].y+=SystemSize_Y;

      NewMeshedIndexMat[int(PG[i_Particle].x/CharLength)+int(SystemSize_X/CharLength)*int(PG[i_Particle].y/CharLength)].push_back(i_Particle);
    }
    MeshedIndexMat=NewMeshedIndexMat;
    
  }

  PolarParticle2D_GT Evolve_ContinuousVicsekModel_Periodic(const PolarParticle2D_GT& PG,
							   const BoxInfoT& BoxInfo,
							   const vector<double>& Parameters,
							   // SystemSize_X SystemSize_Y Charlength Velocity AligningCoeff Sigma GuidingCoeff TimeInterval SqrtTimeinterval
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
      Evolve_ContinuousVicsekModel(vP,MeshedIndex,LDAdjBox,Parameters,Distance_,i_MP);

    auto ansPG=PG;
    UpdatePGFromVParticle(ansPG,vP);
    return ansPG;
  }
  
}

namespace VicsekModel_2D
{//Continuous Normalized Vicsek Model

  void Evolve_CntNml_VicsekModel(vector<Particle>& PG,
				 vector<vector<int>>& MeshedIndexMat,
				 const vector<array<int,4>>& LDAdjBox,
				 const vector<double>& Parameters,
				 // Systemsize_x Systemsize_y Charlength Velocity AligningCoeff Sigma Guidingcoeff Timeinterval SqrttTimeInterval
				 const function<double(Particle,Particle)>& Distance,
				 const int& i_MP=0)
  {
    const double
      SystemSize_X=Parameters[0],
      SystemSize_Y=Parameters[1],
      CharLength=Parameters[2],
      Velocity=Parameters[3],
      AligningCoeff=Parameters[4],
      Sigma=Parameters[5],
      GuidingCoeff=Parameters[6],
      TimeInterval=Parameters[7],
      SqrtTimeInterval=Parameters[8];
    
    vector<double> NeighborDirList_x(PG.size(),0),NeighborDirList_y(PG.size(),0);

    int i_i,i_j;
    double x_i,y_i,x_j,y_j;
    //calculate the mean velocity dir
    for( int i_Box=0 ; i_Box<LDAdjBox.size() ; i_Box++ )
      for( int i=0 ; i<MeshedIndexMat[i_Box].size() ; i++ ){
	i_i=MeshedIndexMat[i_Box][i];
	x_i=PG[i_i].x;
	y_i=PG[i_i].y;

	//itself is NOT taken into account
	//NeighborDirList_x[i_i]+=PG[i_i].vx;
	//NeighborDirList_y[i_i]+=PG[i_i].vy;
	
	//in i_Box
	for( int j=i+1 ; j<MeshedIndexMat[i_Box].size() ; j++ ){
	  i_j=MeshedIndexMat[i_Box][j];
	  x_j=PG[i_j].x;
	  y_j=PG[i_j].y;
	  
	  if((x_i-x_j)*(x_i-x_j)+(y_i-y_j)*(y_i-y_j)<=CharLength*CharLength){
	    NeighborDirList_x[i_i]+=PG[i_j].vx;
	    NeighborDirList_y[i_i]+=PG[i_j].vy;
	    NeighborDirList_x[i_j]+=PG[i_i].vx;
	    NeighborDirList_y[i_j]+=PG[i_i].vy;
	  }
	}

	//inter-box
	for( int i_AdjBox : LDAdjBox[i_Box] )
	  for( int i_j : MeshedIndexMat[i_AdjBox] )
	    if(Distance(PG[i_i],PG[i_j])<=CharLength*CharLength){
	      NeighborDirList_x[i_i]+=PG[i_j].vx;
	      NeighborDirList_y[i_i]+=PG[i_j].vy;
	      NeighborDirList_x[i_j]+=PG[i_i].vx;
	      NeighborDirList_y[i_j]+=PG[i_i].vy;
	    }
      }

    //move particles
    double NewDir,ModuNeighborDir;
    vector<vector<int>> NewMeshedIndexMat(MeshedIndexMat.size());
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      ModuNeighborDir=sqrt(NeighborDirList_x[i_Particle]*NeighborDirList_x[i_Particle]+NeighborDirList_y[i_Particle]*NeighborDirList_y[i_Particle]);
      NewDir=atan2(PG[i_Particle].vy,PG[i_Particle].vx);

      if(ModuNeighborDir>0)// Alignment
	NewDir+=AligningCoeff*(NeighborDirList_y[i_Particle]*PG[i_Particle].vx/ModuNeighborDir
			       -NeighborDirList_x[i_Particle]*PG[i_Particle].vy/ModuNeighborDir)*TimeInterval;
      NewDir+=Sigma*RandGen[i_MP].GaussianRandomDouble()*SqrtTimeInterval;// Noise
      NewDir+=GuidingCoeff*(2*PG[i_Particle].vx*PG[i_Particle].vy)*TimeInterval;// GuidingCoeff*sin(2*Dir): Nematic field
      //WARNING: a little difference between the old version: Nematic force is on the old value of Dir
      
      PG[i_Particle].vx=cos(NewDir);
      PG[i_Particle].vy=sin(NewDir);
      PG[i_Particle].x+=Velocity*PG[i_Particle].vx*TimeInterval;
      PG[i_Particle].y+=Velocity*PG[i_Particle].vy*TimeInterval;

      //WARNING: for periodic only
      if(PG[i_Particle].x>=SystemSize_X) PG[i_Particle].x-=SystemSize_X;
      if(PG[i_Particle].x<=0) PG[i_Particle].x+=SystemSize_X;
      if(PG[i_Particle].y>=SystemSize_Y) PG[i_Particle].y-=SystemSize_Y;
      if(PG[i_Particle].y<=0) PG[i_Particle].y+=SystemSize_Y;

      NewMeshedIndexMat[int(PG[i_Particle].x/CharLength)+int(SystemSize_X/CharLength)*int(PG[i_Particle].y/CharLength)].push_back(i_Particle);
    }
    MeshedIndexMat=NewMeshedIndexMat;
    
  }

  PolarParticle2D_GT Evolve_CtnNml_VicsekModel_Periodic(const PolarParticle2D_GT& PG,
							const BoxInfoT& BoxInfo,
							const vector<double>& Parameters,
							// SystemSize_X SystemSize_Y Charlength Velocity AligningCoeff Sigma GuidingCoeff TimeInterval SqrtTimeinterval
							const int& StepNumber,
							const int& i_MP=0)
  {//an integrated version of Evolve_CtnNml_VicsekModel, evolving for StepNumber steps
    const double
      SystemSize_X=Parameters[0],
      SystemSize_Y=Parameters[1],
      CharLength=Parameters[2];
    
    vector<Particle> vP=PolarParticle2D_GT_To_vParticle(PG);
    vector<vector<int>> MeshedIndex=MeshingIndex(PG,SystemSize_X,SystemSize_Y,CharLength);
    const vector<array<int,4>> LDAdjBox=LowerRightAdjacentBoxIndex(BoxInfo);
    const function<double(Particle,Particle)> Distance_=Distance::Periodic(SystemSize_X,SystemSize_Y,CharLength);

    for( int i_Step=0 ; i_Step<StepNumber ; i_Step++ )
      Evolve_CntNml_VicsekModel(vP,MeshedIndex,LDAdjBox,Parameters,Distance_,i_MP);

    auto ansPG=PG;
    UpdatePGFromVParticle(ansPG,vP);
    return ansPG;
  }
  
}

#ifdef CheckingHeader
HeaderChecker(ContinuousVicsekModel_hpp);
#endif //CheckingHeader

#endif //ContinuousVicsekModel_hpp
