/*--------------------
ver 250204
--------------------*/

#ifndef StateProcessing_hpp
#define StateProcessing_hpp

#include<bits/stdc++.h>
#include"BasicTools.hpp"
#include"AgentRealization.hpp"

using namespace std;

namespace StateProcessing{

  vector<PolarParticle2D> GenerateUniformRandomParticleGroup(double SystemSize_X,double SystemSize_Y,int N,int i_MP=0)
  {//Uniformly distributed in space and orientation, MP ver: i_MP=omp_get_thread_num
    vector<PolarParticle2D> PG;
    for( int i=0 ; i<N ; i++ ){
      #ifdef LocalProgram
      PG.push_back(Make_PolarParticle2D(RandGen[i_MP].RandomDouble()*SystemSize_X,RandGen[i_MP].RandomDouble()*SystemSize_Y,RandGen[i_MP].RandomDouble()*2*pi));
      #else
      PG.push_back(Make_PolarParticle2D(RandGen[i_MP].RandomDouble()*SystemSize_X,
					RandGen[i_MP].RandomDouble()*SystemSize_Y,
					RandGen[i_MP].RandomDouble()*2*pi));
      #endif
    }
    return PG;
  }

  vector<PolarParticle2D> GenerateDensityProfileRandomParticleGroup(double SystemSize_X,double SystemSize_Y,int N,
								       function<pair<double,double>(pair<double,double>)> DensityProfileMap,
								       int i_MP=0)
  {//uniformly distributed in orientation, and based on DensityProfileMap to generate a density profile
   //MP ver: i_MP=omp_get_thread_num
    vector<PolarParticle2D> PG;
    for( int i=0 ; i<N ; i++ ){
      double
	x0=RandGen[i_MP].RandomDouble()*SystemSize_X,
	y0=RandGen[i_MP].RandomDouble()*SystemSize_Y;
      auto pair_x_y=DensityProfileMap({x0,y0});
      double
	x=pair_x_y.first,
	y=pair_x_y.second;
      
      PG.push_back(Make_PolarParticle2D(x,y,RandGen[i_MP].RandomDouble()*2*pi));
    }
    return PG;
  }

  namespace DensityProfileMap
  {//some default DensityProfileMap
    function<pair<double,double>(pair<double,double>)> LinearX(double SystemSize_X,double SystemSize_Y)
    {//Density propto x
      return
	[SystemSize_X,SystemSize_Y](pair<double,double> pair_x0_y0)->pair<double,double>
	{
	  double x0=pair_x0_y0.first,y0=pair_x0_y0.second;
	  double y=y0,x=sqrt(SystemSize_X*x0);
	  return {x,y};
	};
    }

    function<pair<double,double>(pair<double,double>)> LinearY(double SystemSize_X,double SystemSize_Y)
    {//Density propto y
      return
	[SystemSize_X,SystemSize_Y](pair<double,double> pair_x0_y0)->pair<double,double>
	{
	  double x0=pair_x0_y0.first,y0=pair_x0_y0.second;
	  double x=x0,y=sqrt(SystemSize_Y*y0);
	  return {x,y};
	};
    }
    
  };

  vector<PolarParticle2D> RandomlySetOrientation(vector<PolarParticle2D> PG,int i_MP=0)
  {//randomly set the orientation of particles in PG but keep their positions
    vector<PolarParticle2D> ansPG=PG;
    for( int i=0 ; i<PG.size() ; i++ )
      ansPG[i].Dir=RandGen[i_MP].RandomDouble()*2*pi;
    return ansPG;
  }

  Meshed_PolarParticle2D_GT RandomlySetOrientation(Meshed_PolarParticle2D_GT MeshedPG,int i_MP=0)
  {//randomly set the orientation of particles in PG but keep their positions
    for( auto& PG : MeshedPG )
      for( auto& p : PG)
	p.Dir=RandGen[i_MP].RandomDouble()*2*pi;

    return MeshedPG;
  }

  vector<PolarParticle2D> FlipParticleGroupLeftRight(const vector<PolarParticle2D>& PG,double SystemSize_X,double SystemSize_Y)
  {//Turn an polar state leftward right

    vector<PolarParticle2D> AnsPG;
    double NewX,NewY,NewDir;
    for( int i=0 ; i<PG.size() ; i++ ){
      NewX=SystemSize_X-PG[i].Pos.real();
      NewY=PG[i].Pos.imag();
      NewDir=pi-PG[i].Dir;
      NewDir=(NewDir<0)?(NewDir+2*pi):NewDir;
      AnsPG.push_back(PolarParticle2D(NewX+ii*NewY,NewDir));
    }

    return AnsPG;
  }

  vector<PolarParticle2D> FlipParticleGroupRightUp(const vector<PolarParticle2D>& PG,double SystemSize_X,double SystemSize_Y)
  {// (1,2) pi/6 => (2,1) pi/3, requiring SystemSize_X==SystemSize_Y

    if(SystemSize_X!=SystemSize_Y)
      throw runtime_error("In StateProcessing::FlipParticleGroupRightUp, the system should be square.");
    
    vector<PolarParticle2D> AnsPG;
    double NewX,NewY,NewDir;
    for( int i=0 ; i<PG.size() ; i++ ){
      NewX=PG[i].Pos.imag();
      NewY=PG[i].Pos.real();
      NewDir=pi/2-PG[i].Dir;
      while(NewDir<0) NewDir+=2*pi;
      while(NewDir>=2*pi) NewDir-=2*pi;
      AnsPG.push_back(PolarParticle2D(NewX+ii*NewY,NewDir));
    }

    return AnsPG;
  }

  vector<PolarParticle2D> MovingParticleGroup(vector<PolarParticle2D> PG,complex<double> displacement,
					      function<PolarParticle2D(PolarParticle2D)> BoundaryCondFunc=[](PolarParticle2D p)->PolarParticle2D{return p;})
  {//to move a PG collectively
    for( int i=0 ; i<PG.size() ; i++ ){
      PG[i].Pos+=displacement;
      PG[i]=BoundaryCondFunc(PG[i]);
    }
    return PG;
  }

};

using BoxInfoT=vector<pair<vector<double>,vector<int>>>;

namespace BoxMeshing_2D{

  template<typename Type>
  vector<Type> Gathering(vector<vector<Type>> MeshedPG)
  {//gather all particles in each small boxes
    vector<PolarParticle2D> AnsVec;

    for( auto p : MeshedPG  )
      AnsVec.insert(AnsVec.end(),p.begin(),p.end());
  
    return AnsVec;
  }

  BoxInfoT BoxInfo(double SystemSize_X,double SystemSize_Y,double CharLength)
  {//Divide the system into many small boxes, return the information of boxes; {{LEnd,REnd,DEnd,UEnd},{LIndex,RIndex,DIndex,UIndex}}

    vector<pair<vector<double>,vector<int>>> AnsList;
  
    int Num_x=floor(SystemSize_X/CharLength),Num_y=floor(SystemSize_Y/CharLength);
    for( int j=0 ; j<Num_y ; j++ )
      for( int i=0 ; i<Num_x ; i++ ){
	vector<double> BdrList={};
	vector<int> IndexList={};

	BdrList.push_back(i*CharLength);
	if(i==0) IndexList.push_back(Num_x*j+Num_x-1);
	else IndexList.push_back(Num_x*j+i-1);

	if(i==Num_x-1){
	  BdrList.push_back(SystemSize_X);
	  IndexList.push_back(Num_x*j);
	}
	else{
	  BdrList.push_back((i+1)*CharLength);
	  IndexList.push_back(Num_x*j+i+1);
	}

	BdrList.push_back(j*CharLength);
	if(j==0) IndexList.push_back((Num_y-1)*Num_x+i);
	else IndexList.push_back(Num_x*(j-1)+i);

	if(j==Num_y-1){
	  BdrList.push_back(SystemSize_Y);
	  IndexList.push_back(i);
	}
	else{
	  BdrList.push_back((j+1)*CharLength);
	  IndexList.push_back(Num_x*(j+1)+i);
	}

	AnsList.push_back({BdrList,IndexList});
      }

    return AnsList;
  }

  tuple<double,double,double,int,int> UnzipBoxInfo(BoxInfoT BoxInfo_)
  {//Derive parameters from BoxInfo, return (SystemSize_X,SystemSize_Y,CharLength,BoxNum_X,BoxNum_Y)
    double SystemSize_X=BoxInfo_[BoxInfo_.size()-1].first[1];
    double SystemSize_Y=BoxInfo_[BoxInfo_.size()-1].first[3];
    double CharLength=BoxInfo_[0].first[1]-BoxInfo_[0].first[0];
    int BoxNum_X=BoxInfo_[0].second[0]+1;
    int BoxNum_Y=BoxInfo_[0].second[2]/BoxNum_X+1;
    return make_tuple(SystemSize_X,SystemSize_Y,CharLength,BoxNum_X,BoxNum_Y);
  }

  template<typename Type>
  int DecideBoxIndex(Type p,double SystemSize_X,double SystemSize_Y,double CharLength)
  {//decide the index of the box to which the particle p will be sent
    if(p.Pos.real()<0 or p.Pos.real()>SystemSize_X or p.Pos.imag()<0 or p.Pos.imag()>SystemSize_Y)
      throw runtime_error("In BoxMeshing_2D::DecideBoxIndex : the coordinate of particle should not exceed the system size.");
    
    int
      Num_x=floor(SystemSize_X/CharLength),
      Num_y=floor(SystemSize_Y/CharLength),
      XIndex=floor(p.Pos.real()/CharLength),
      YIndex=floor(p.Pos.imag()/CharLength);
    if(XIndex>=Num_x)XIndex--;
    if(YIndex>=Num_y)YIndex--;
    int BoxIndex=XIndex+Num_x*YIndex;
    return BoxIndex;
  }

  template<typename Type>
  vector<vector<Type>> Meshing(vector<Type> list,double SystemSize_X,double SystemSize_Y,double CharLength)
  {//Divide the system into many small boxes
  
    int Num_x=floor(SystemSize_X/CharLength);
    int Num_y=floor(SystemSize_Y/CharLength);
    vector<vector<Type>> AnsList(Num_x*Num_y);
  
    for( int i=0 ; i<list.size() ; i++ )
      AnsList[DecideBoxIndex(list[i],SystemSize_X,SystemSize_Y,CharLength)].push_back(list[i]);
  
    return AnsList;
  }

  template<typename Type>
  pair<vector<vector<Type>>,vector<pair<vector<double>,vector<int>>>>
  Coarsening(const vector<vector<Type>>& MeshedPG,const vector<pair<vector<double>,vector<int>>>& BoxInfo_,int TimeCoarsening)
  {//to get a MeshedPG with originally TimeCoarsening*TimeCoarsening boxes into one box, return a pair {NewMeshedPG,NewBoxInfo}
    double SystemSize_X=BoxInfo_[BoxInfo_.size()-1].first[1];
    double SystemSize_Y=BoxInfo_[BoxInfo_.size()-1].first[3];
    double CharLength=BoxInfo_[0].first[1]-BoxInfo_[0].first[0];
    int N_X=BoxInfo_[0].second[0]+1;
    int N_Y=BoxInfo_[0].second[2]/N_X+1;

    //WARNING: Timecoarsening should be a factor of both N_X and N_Y

    vector<vector<Type>> NewMeshedPG(BoxInfo_.size()/TimeCoarsening/TimeCoarsening);

    for( int i=0 ; i<BoxInfo_.size() ; i++ ){
      int i_X=i%N_X, i_Y=i/N_X;
      int New_i_X=i_X/TimeCoarsening, New_i_Y=i_Y/TimeCoarsening;
      int New_i=New_i_X+New_i_Y*(N_X/TimeCoarsening);

      NewMeshedPG[New_i].insert(NewMeshedPG[New_i].end(),MeshedPG[i].begin(),MeshedPG[i].end());
    }
    
    auto NewBoxInfo=BoxInfo(SystemSize_X,SystemSize_Y,CharLength*TimeCoarsening);
    return {NewMeshedPG,NewBoxInfo};
  }
  
};


#ifdef CheckingHeader
HeaderChecker(StateProcessing_hpp);
#endif //CheckingHeader

#endif //StateProcessing_hpp
