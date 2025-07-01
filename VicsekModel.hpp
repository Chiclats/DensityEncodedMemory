/*--------------------
ver 250701
--------------------*/

#ifndef VicsekModel_hpp
#define VicsekModel_hpp

#ifndef vs
#include<bits/stdc++.h>
#endif

#include<omp.h>
#include"BasicTools.hpp"
#include"AgentRealization.hpp"
#include"StateProcessing.hpp"

namespace VicsekModel_2D{

  struct Particle{
    double x,y,vx,vy;
  };

  Particle PolarParticle2D_To_Particle(const PolarParticle2D& p){
    double Dir=p.Dir;
    Particle ans;
    ans.x=p.Pos.real();
    ans.y=p.Pos.imag();
    ans.vx=cos(Dir);
    ans.vy=sin(Dir);
    return ans;
  }

  vector<Particle> PolarParticle2D_GT_To_vParticle(const PolarParticle2D_GT& PG){
    vector<Particle> ansPG;
    for( auto p : PG ){
      ansPG.push_back(PolarParticle2D_To_Particle(p));
    }
    return ansPG;
  }

  void UpdatePGFromVParticle(PolarParticle2D_GT& PG,const vector<Particle>& vP)
  {//Update a PolarParticle2D_GT based on a vector<Particle>, only Pos and Dir changes
    for( int i=0 ; i<PG.size() ; i++ ){
      PG[i].Pos=vP[i].x+vP[i].y*ii;
      PG[i].Dir=atan2(vP[i].vy,vP[i].vx);
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

  vector<array<int,8>> AllAdjacentBoxIndex(const BoxInfoT& BoxInfo)
  {//to return the indices of all adjacent boxes
    vector<array<int,8>> AnsVec(BoxInfo.size());

    for( int i_Box=0 ; i_Box<BoxInfo.size() ; i_Box++ )
      AnsVec[i_Box]={
	BoxInfo[i_Box].second[0],
	BoxInfo[i_Box].second[1],
	BoxInfo[i_Box].second[2],
	BoxInfo[i_Box].second[3],
	BoxInfo[BoxInfo[i_Box].second[2]].second[0],
	BoxInfo[BoxInfo[i_Box].second[2]].second[1],
	BoxInfo[BoxInfo[i_Box].second[3]].second[0],
	BoxInfo[BoxInfo[i_Box].second[3]].second[1],
      };
    
    return AnsVec;
  }

  vector<vector<int>> MeshingIndex(const vector<PolarParticle2D>& PG,
				   double SystemSize_X,double SystemSize_Y,double CharLength)
  {
    auto BoxInfo=BoxMeshing_2D::BoxInfo(SystemSize_X,SystemSize_Y,CharLength);
    vector<vector<int>> AnsMat(BoxInfo.size());
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      AnsMat[BoxMeshing_2D::DecideBoxIndex(PG[i_Particle],SystemSize_X,SystemSize_Y,CharLength)].push_back(i_Particle);
    }
    return AnsMat;
  }

  void Evolve_VicsekModel(vector<Particle>& PG,
			  vector<vector<int>>& MeshedIndexMat,
			  const vector<array<int,4>>& LDAdjBox,
			  const vector<double>& Parameters,// Systemsize_x Systemsize_y Charlength Velocity Sigma Guidingcoeff
			  const function<double(Particle,Particle)>& Distance,
			  const int& i_MP=0)
  {
    const double
      SystemSize_X=Parameters[0],
      SystemSize_Y=Parameters[1],
      CharLength=Parameters[2],
      Velocity=Parameters[3],
      Sigma=Parameters[4];
    double GuidingCoeff=0;
    if(Parameters.size()==6) GuidingCoeff=Parameters[5];
    
    vector<double> NeighborDirList_x(PG.size(),0),NeighborDirList_y(PG.size(),0);

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
    double NewDir;
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      if(NeighborDirList_y[i_Particle]==0 and NeighborDirList_x[i_Particle]==0){//to avoid everything cancelled out
	NewDir=atan2(PG[i_Particle].vy,PG[i_Particle].vx);
	cerr<<"Attention: zero vector sum causing dir decision failure."<<endl;
      }
      else NewDir=atan2(NeighborDirList_y[i_Particle],NeighborDirList_x[i_Particle]);
      
      NewDir+=Sigma*(RandGen[i_MP].RandomDouble()*2*pi-pi); 
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

  void Evolve_VicsekModel_MP(vector<Particle>& PG,
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
	for( int i_AdjBox : AllAdjBox[i_Box] )
	  for( int i_j : MeshedIndexMat[i_AdjBox] )
	    if(Distance(PG[i_i],PG[i_j])<=CharLength*CharLength){
	      NeighborDirList_x[i_i]+=PG[i_j].vx;
	      NeighborDirList_y[i_i]+=PG[i_j].vy;
	    }
      }
    }

    //move particles
    #pragma omp parallel for
    for( int i_Particle=0 ; i_Particle<PG.size() ; i_Particle++ ){
      double NewDir=atan2(NeighborDirList_y[i_Particle],NeighborDirList_x[i_Particle]);
      NewDir+=Sigma*(RandGen[omp_get_thread_num()].RandomDouble()*2*pi-pi); 
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
    
  };

  PolarParticle2D_GT Evolve_VicsekModel_Periodic(const PolarParticle2D_GT& PG,
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

    for( int i_Step=0 ; i_Step<StepNumber ; i_Step++ ){
      Evolve_VicsekModel(vP,MeshedIndex,LDAdjBox,Parameters,Distance_,i_MP);

      #ifdef EVOLVE_SHRINK_TO_FIT_TURN
      //To avoid the expansion of MeshedIndex, which sometimes leads to OOM
      //EVOLVE_SHRINK_TO_FIT_TURN, when defined, is an integer
      if(i_Step%EVOLVE_SHRINK_TO_FIT_TURN==0){
	for( int i_Box=0 ; i_Box<BoxInfo.size() ; i_Box++ )
	  MeshedIndex[i_Box].shink_to_fit();
	MeshedIndex.shink_to_fit();
      }	
      #endif
    }

    auto ansPG=PG;
    UpdatePGFromVParticle(ansPG,vP);
    return ansPG;
  }

#ifndef LocalProgram
  PolarParticle2D_GT Evolve_VicsekModel_Periodic_MP(const PolarParticle2D_GT& PG,
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
      Evolve_VicsekModel_MP(vP,MeshedIndex,AllAdjBox,Parameters,Distance_);

    auto ansPG=PG;
    UpdatePGFromVParticle(ansPG,vP);
    return ansPG;
  }
#endif

}

namespace DirectlyOutput{

  void Vec(vector<VicsekModel_2D::Particle> PG){

    for( auto p : PG ){
      cout<<'('<<p.x<<' '<<p.y<<") "<<atan2(p.vy,p.vx)<<" ("<<p.vx<<' '<<p.vy<<")\n";
    }

  }

};

namespace FileIO{

  void OutParticleVec(const string& Filename,const vector<VicsekModel_2D::Particle>& PG,
		      const double& PosPrecession=DoublePrecession,const double& DirPrecession=DoublePrecession){
    ofstream fout;
    fout.open(Filename);
    for( int i=0 ; i<PG.size() ; i++){
      fout<<setprecision(PosPrecession)<<PG[i].x<<' '<<PG[i].y<<' ';
      fout<<setprecision(DirPrecession)<<atan2(PG[i].vy,PG[i].vx)<<endl;
    }
    
    fout.close();
    return;
  }
  
}
			
#endif
