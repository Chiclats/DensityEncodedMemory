/*--------------------
ver 241229
--------------------*/

#ifndef NematicAlgnPP_hpp
#define NematicAlgnPP_hpp

#include<bits/stdc++.h>
#include<omp.h>
#include"BasicTools.hpp"
#include"AgentRealization.hpp"
#include"StateProcessing.hpp"

namespace NematicAlgnPolarParticles_2D{

  struct Particle{
    double x,y,vx,vy;
  };

  Particle PolarParticle2D_To_Particle(PolarParticle2D p){
    double Dir=p.Dir;
    Particle ans;
    ans.x=p.Pos.real();
    ans.y=p.Pos.imag();
    ans.vx=cos(Dir);
    ans.vy=sin(Dir);
    return ans;
  }

  vector<Particle> PolarParticle2D_GT_To_NAPP(PolarParticle2D_GT PG){
    vector<Particle> ansPG;
    for( auto p : PG ){
      ansPG.push_back(PolarParticle2D_To_Particle(p));
    }
    return ansPG;
  }

  void UpdatePGFromNAPP(PolarParticle2D_GT& PG,vector<Particle> NAPP)
  {//Update a PolarParticle2D_GT based on a vector<Particle>, only Pos and Dir changes
    for( int i=0 ; i<PG.size() ; i++ ){
      PG[i].Pos=NAPP[i].x+NAPP[i].y*ii;
      PG[i].Dir=atan2(NAPP[i].vy,NAPP[i].vx);
      if(PG[i].Dir<0) PG[i].Dir+=2*pi;
    }
  }

  vector<array<int,4>> LowerRightAdjacentBoxIndex(const BoxInfoT& BoxInfo)
  {//to return the indices of lower and right adjacent boxes
    vector<array<int,4>> AnsVec(BoxInfo.size());

    for( int i_Box=0 ; i_Box<BoxInfo.size() ; i_Box++ )
      AnsVec[i_Box]={
	BoxInfo[i_Box].second[1],
	BoxInfo[i_Box].second[2],
	BoxInfo[BoxInfo[i_Box].second[2]].second[0],
	BoxInfo[BoxInfo[i_Box].second[2]].second[1]
      };
    
    return AnsVec;
  }

  vector<vector<int>> MeshingIndex(const vector<PolarParticle2D>& PG,
				   double SystemSize_X,double SystemSize_Y,double CharLength)
  {
    vector<vector<int>> AnsMat((int) floor(SystemSize_X/CharLength)*floor(SystemSize_Y/CharLength));
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      AnsMat[BoxMeshing_2D::DecideBoxIndex(PG[i_Particle],SystemSize_X,SystemSize_Y,CharLength)].push_back(i_Particle);
    }
    return AnsMat;
  }

  void Evolve_NAPP(vector<Particle>& PG,
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
    //to store sum(exp(i 2Dir))

    int i_i,i_j;
    double x_i,y_i,x_j,y_j;
    //calculate the mean velocity dir
    for( int i_Box=0 ; i_Box<LDAdjBox.size() ; i_Box++ )
      for( int i=0 ; i<MeshedIndexMat[i_Box].size() ; i++ ){
	i_i=MeshedIndexMat[i_Box][i];
	x_i=PG[i_i].x;
	y_i=PG[i_i].y;

	//itself is taken into account
	NeighborDirList_x[i_i]+=2*PG[i_i].vx*PG[i_i].vx-1; //2 * cos^2 Dir - 1
	NeighborDirList_y[i_i]+=2*PG[i_i].vy*PG[i_i].vx;   //2 * cos Dir * sin Dir
	
	//in i_Box
	for( int j=i+1 ; j<MeshedIndexMat[i_Box].size() ; j++ ){
	  i_j=MeshedIndexMat[i_Box][j];
	  x_j=PG[i_j].x;
	  y_j=PG[i_j].y;

	  if((x_i-x_j)*(x_i-x_j)+(y_i-y_j)*(y_i-y_j)<=CharLength*CharLength){
	    NeighborDirList_x[i_i]+=2*PG[i_j].vx*PG[i_j].vx-1;
	    NeighborDirList_y[i_i]+=2*PG[i_j].vy*PG[i_j].vx;
	    NeighborDirList_x[i_j]+=2*PG[i_i].vx*PG[i_i].vx-1;
	    NeighborDirList_y[i_j]+=2*PG[i_i].vy*PG[i_i].vx;
	  }
	}

	//inter-box
	for( int i_AdjBox : LDAdjBox[i_Box] )
	  for( int i_j : MeshedIndexMat[i_AdjBox] )
	    if(Distance(PG[i_i],PG[i_j])<=CharLength*CharLength){
	      NeighborDirList_x[i_i]+=2*PG[i_j].vx*PG[i_j].vx-1;
	      NeighborDirList_y[i_i]+=2*PG[i_j].vy*PG[i_j].vx;
	      NeighborDirList_x[i_j]+=2*PG[i_i].vx*PG[i_i].vx-1;
	      NeighborDirList_y[i_j]+=2*PG[i_i].vy*PG[i_i].vx;
	    }
      }

    //move particles
    double NewDir,New2Dir;
    vector<vector<int>> NewMeshedIndexMat(MeshedIndexMat.size());
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){

      if(NeighborDirList_x[i_Particle]==0 and NeighborDirList_y[i_Particle]==0){//to avoid everything cancelled out
	NewDir=atan2(PG[i_Particle].vy,PG[i_Particle].vx);
	cerr<<"Attention: zero vector sum causing dir decision failure."<<endl;
      }
      else{
	New2Dir=atan2(NeighborDirList_y[i_Particle],NeighborDirList_x[i_Particle]); //(-pi,pi]
	if(New2Dir<0) New2Dir+=2*pi; //2-time NewDir
	NewDir=New2Dir/2;
	if(PG[i_Particle].vx*cos(NewDir)+PG[i_Particle].vy*sin(NewDir)<0)
	  NewDir+=pi; //choose from NewDir or NewDir+pi the one closer to original Dir
      }
            
      NewDir+=Sigma*(RandGen[i_MP].RandomDouble()*2*pi-pi); 
      NewDir+=GuidingCoeff*(4*PG[i_Particle].vx*PG[i_Particle].vx*PG[i_Particle].vx*PG[i_Particle].vy-
			    4*PG[i_Particle].vx*PG[i_Particle].vy*PG[i_Particle].vy*PG[i_Particle].vy);// GuidingCoeff*sin(4*Dir)
      //WARNING: Nematic force is on the old value of Dir
      PG[i_Particle].vx=cos(NewDir);
      PG[i_Particle].vy=sin(NewDir);
      PG[i_Particle].x+=Velocity*PG[i_Particle].vx;
      PG[i_Particle].y+=Velocity*PG[i_Particle].vy;

      //WARNING: for periodic only
      if(PG[i_Particle].x>=SystemSize_X) PG[i_Particle].x-=SystemSize_X;
      if(PG[i_Particle].x<=0) PG[i_Particle].x+=SystemSize_X;
      if(PG[i_Particle].y>=SystemSize_Y) PG[i_Particle].y-=SystemSize_Y;
      if(PG[i_Particle].y<=0) PG[i_Particle].y+=SystemSize_Y;

      NewMeshedIndexMat[int(PG[i_Particle].x/CharLength)+int(SystemSize_X/CharLength)*int(PG[i_Particle].y/CharLength)].push_back(i_Particle);
    }
    MeshedIndexMat=NewMeshedIndexMat;
    
  }

  namespace Distance
  {//Distance is the square value

    function<double(Particle,Particle)> Periodic(double SystemSize_X,double SystemSize_Y,double CharLength){
      return [SystemSize_X,SystemSize_Y,CharLength](Particle p1,Particle p2)->double{
	double x1=p1.x,x2=p2.x,y1=p1.y,y2=p2.y;
	double dx=abs(x1-x2),dy=abs(y1-y2);
	if(dx>=SystemSize_X-2*CharLength) dx=SystemSize_X-dx;
	if(dy>=SystemSize_Y-2*CharLength) dy=SystemSize_Y-dy;
	return dx*dx+dy*dy;
      };
    }
    
  }

  PolarParticle2D_GT Evolve_NAPP_Periodic(const PolarParticle2D_GT& PG,
					  const BoxInfoT& BoxInfo,
					  const vector<double>& Parameters,// Systemsize_x Systemsize_y Charlength Velocity Sigma Guidingcoeff
					  const int& StepNumber,
					  const int& i_MP=0)
  {//an integrated version of Evolve_VicsekModel, evolving for StepNumber steps
    const double
      SystemSize_X=Parameters[0],
      SystemSize_Y=Parameters[1],
      CharLength=Parameters[2];
    
    vector<Particle> NAPP=PolarParticle2D_GT_To_NAPP(PG);
    vector<vector<int>> MeshedIndex=MeshingIndex(PG,SystemSize_X,SystemSize_Y,CharLength);
    const vector<array<int,4>> LDAdjBox=LowerRightAdjacentBoxIndex(BoxInfo);
    const function<double(Particle,Particle)> Distance_=Distance::Periodic(SystemSize_X,SystemSize_Y,CharLength);

    for( int i_Step=0 ; i_Step<StepNumber ; i_Step++ )
      Evolve_NAPP(NAPP,MeshedIndex,LDAdjBox,Parameters,Distance_,i_MP);
    
    auto ansPG=PG;
    UpdatePGFromNAPP(ansPG,NAPP);
    return ansPG;
  }
  
}

#endif//Nematicalgnpp_hpp
