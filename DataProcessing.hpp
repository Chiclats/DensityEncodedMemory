/*--------------------
ver 240805
--------------------*/

#ifndef DataProcessing_hpp
#define DataProcessing_hpp

#include<bits/stdc++.h>
#include"BasicTools.hpp"
#include"AgentRealization.hpp"
#include"StateProcessing.hpp"

#ifndef LocalProgram
#include<omp.h>
#endif//LocalProgram

using namespace std;

namespace DataProcessing
{//a collection of data processing functions

  vector<vector<PolarParticle2D>> ParticleGroupDivision(vector<PolarParticle2D> PG,int TagNum)
  {//Divide a PG with particles of different tags into several PG, Tag should always start from 1 and be integers
    if(TagNum==1)
      return {PG};
    
    vector<vector<PolarParticle2D>> ansVecPG(TagNum);
    for( auto p : PG )
      ansVecPG[p.Tag-1].push_back(p);
      
    return ansVecPG;
  }//checked 240502

  pair<double,double> NFoldOrderParameter(vector<PolarParticle2D> PG,double fold)
  {//calculate the n-fold order parameter of a PG
    double n=PG.size();
    complex<double> sum(0,0);

    for( auto p : PG )
      sum+=exp(ii*p.Dir*fold);

    sum=sum/n;
    double modu=abs(sum),dir=arg(sum);
    
    while(dir>2*pi) dir-=(2*pi);
    while(dir<0) dir+=(2*pi);
    //draw dir back to the range (0,2*pi)
    
    return {modu,dir};
  }//checked 240502

  tuple<double,double> NFoldOrderParameter(Meshed_PolarParticle2D_GT MeshedPG,double fold)
  {//calculate the n-fold order parameter of a PG
    double n=0;
    complex<double> sum(0,0);

    for( auto ps : MeshedPG ){
      for( auto p : ps)
	sum+=exp(ii*p.Dir*fold);
      n+=ps.size();
    }

    sum=sum/n;
    double modu=abs(sum),dir=arg(sum);
    
    while(dir>2*pi) dir-=(2*pi);
    while(dir<0) dir+=(2*pi);
    //draw dir back to the range (0,2*pi)
    
    return make_tuple(modu,dir);
  }

  vector<PolarParticle2D> PutIndexToTag(vector<PolarParticle2D> PG)
  {//Set the Tag slot to be the index of particles in PG
    for( int i=0 ; i<PG.size() ; i++ )
      PG[i].Tag=i;
    return PG;
  }

  vector<vector<PolarParticle2D>> PutIndexToTag(vector<vector<PolarParticle2D>> MeshedPG)
  {//Set the Tag slot to be the index of particles in MeshedPG
    int cnt=0;
    for( int i=0 ; i<MeshedPG.size() ; i++ )
      for( int j=0 ; j<MeshedPG[i].size() ; j++ ){
	MeshedPG[i][j].Tag=cnt;
	cnt++;
      }
    return MeshedPG;
  }
  
  vector<pair<complex<int>,double>> _CoorPairList_Initialization(){//initialize CoorPairList
    vector<pair<complex<int>,double>> CoorPairList;
    for( int i=-50 ; i<=50 ; i++ )
      for( int j=-50 ; j<=50 ; j++ ){
	double min_dist;
	if(i==0 and j==0) min_dist=0;
	else if(i==0) min_dist=(j>0)?(j-1):(-j-1);
	else if(j==0) min_dist=(i>0)?(i-1):(-i-1);
	else min_dist=sqrt(((j>0)?(j-1):(j+1))*((j>0)?(j-1):(j+1))+((i>0)?(i-1):(i+1))*((i>0)?(i-1):(i+1)));
      
	CoorPairList.push_back({complex<int>(i,j),min_dist});
      }
    sort(CoorPairList.begin(),CoorPairList.end(),
	 [](pair<complex<int>,double> p1,pair<complex<int>,double> p2){
	   return p1.second<p2.second;
	 });
    return CoorPairList;
  }
  vector<pair<complex<int>,double>> _CoorPairList=_CoorPairList_Initialization(); // {(x,y), min_dist}
  
  vector<vector<PolarParticle2D>> ClosestSeveralNeighborsList(Meshed_PolarParticle2D_GT MeshedPG,
							      BoxInfoT BoxInfo,
							      function<double(PolarParticle2D,PolarParticle2D)> DistanceOfTwoParticles,
							      int NeighborNum=1)
  {//List out the closest particle of each particle: {{PolarParticle2D,Closest PolarParticle2D,2nd Closest PolarParticle2D},...}

    vector<vector<PolarParticle2D>> AnsVec;

    double CharLength;
    int N_X,N_Y;
    tie(ignore,ignore,CharLength,N_X,N_Y)=BoxMeshing_2D::UnzipBoxInfo(BoxInfo);

    auto IndexShift=[N_X,N_Y](int i,complex<int> Coor)->int{
      int
	n_X=i%N_X,
	n_Y=i/N_X,
	new_n_X=(n_X+Coor.real())%N_X,
	new_n_Y=(n_Y+Coor.imag())%N_Y;
      if(new_n_X<0) new_n_X+=N_X;
      if(new_n_Y<0) new_n_Y+=N_Y;
      return new_n_X+new_n_Y*N_X;
    };

    auto InsertParticle=[NeighborNum,DistanceOfTwoParticles](const PolarParticle2D& p,const PolarParticle2D& np,
							     vector<PolarParticle2D>& Neighbors)->void{
      //insert a new neighbor particle into the list `Neighbors`, `p` itself will be also contained
      double nd=DistanceOfTwoParticles(np,p);
      for( auto it=Neighbors.begin() ; it<Neighbors.end() ; it++ )
	//particles with tag=-1 are placeholders
	if((*it).Tag==-1 or DistanceOfTwoParticles(p,*it)>nd)
	  if(it==Neighbors.begin() or (*(it-1)).Pos!=np.Pos or (*(it-1)).Dir!=np.Dir){
	    Neighbors.insert(it,np);
	    break;
	  }
      Neighbors.resize(NeighborNum+1);
    };

    for( int i_Box=0 ; i_Box<BoxInfo.size() ; i_Box++ )
      for( auto p : MeshedPG[i_Box] ){

	vector<PolarParticle2D> Neighbors(NeighborNum+1,PolarParticle2D(0,0,-1));
	//particles with tag=-1 are placeholders

	//find near enough particles
	int i_Coor;
	for( i_Coor=0 ; Neighbors.back().Tag==-1 ; i_Coor++ ){
	  //repeat until `Neighbors` is full
	  auto Coor=_CoorPairList[i_Coor].first;
	  
	  for( auto np : MeshedPG[IndexShift(i_Box,Coor)] )
	    InsertParticle(p,np,Neighbors);
	}

	double dist=DistanceOfTwoParticles(p,Neighbors.back());
	
	//find closest particles

	for( int j_Coor=i_Coor ; _CoorPairList[j_Coor].second*CharLength<dist ; j_Coor++ ){
	  auto Coor=_CoorPairList[j_Coor].first;
	  
	  for( auto np : MeshedPG[IndexShift(i_Box,Coor)] )
	    InsertParticle(p,np,Neighbors);

	  dist=DistanceOfTwoParticles(p,Neighbors.back());
	}

	AnsVec.push_back(Neighbors);
      }
    
    return AnsVec;
  }

  vector<complex<double>> RelativePosList(const vector<vector<PolarParticle2D>>& MeshedPG,
					  const vector<pair<vector<double>,vector<int>>>& BoxInfo,
					  int GridNumConsidered, /* an even number */
					  int SampleRatio /*a sample per SampleRatio particles*/) 
  {//List out the relative position of each two particles, in a GridNumConsidered*GridNumConsidered box, peripheral particles are discarded
    vector<complex<double>> AnsVec;

    double SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1];
    double SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3];
    int N_X=BoxInfo[0].second[0]+1;
    int N_Y=BoxInfo[0].second[2]/N_X+1;
    double CharLength=BoxInfo[0].first[1]-BoxInfo[0].first[0];
    double MaxXYDist=GridNumConsidered/2*CharLength;

    auto IndexShift=[N_X,N_Y](int i,complex<int> Coor)->int{
      int
	n_X=i%N_X,
	n_Y=i/N_X,
	new_n_X=(n_X+Coor.real())%N_X,
	new_n_Y=(n_Y+Coor.imag())%N_Y;
      if(new_n_X<0) new_n_X+=N_X;
      if(new_n_Y<0) new_n_Y+=N_Y;
      return new_n_X+new_n_Y*N_X;
    };

    int cnt=0; //for controling the sample ratio
    
    for( int i_Box=0 ; i_Box<BoxInfo.size() ; i_Box++ ){
      
      //discard the periphery
      int n_X=i_Box%N_X, n_Y=i_Box/N_X;
      if(n_X-GridNumConsidered/2<0 or n_X+GridNumConsidered/2>=N_X or n_Y-GridNumConsidered/2<0 or n_Y+GridNumConsidered/2>=N_Y)
	continue;

      //collect all particles in (GridNumConsidered+1)*(GridNumConsidered+1)
      vector<vector<PolarParticle2D>> PortionOfMeshedPG;
      for( int i_Shift_X=-GridNumConsidered/2 ; i_Shift_X<=GridNumConsidered/2 ; i_Shift_X++ )
	for( int i_Shift_Y=-GridNumConsidered/2 ; i_Shift_Y<=GridNumConsidered/2 ; i_Shift_Y++ )
	  PortionOfMeshedPG.push_back(MeshedPG[IndexShift(i_Box,complex<int>(i_Shift_X,i_Shift_Y))]);
      auto Gathered_PortionOfMeshedPG=BoxMeshing_2D::Gathering(PortionOfMeshedPG);

      //calculate relative position
      for( auto p : MeshedPG[i_Box] ){
	if((cnt++)%SampleRatio!=0) continue;
	
	for( auto pp : Gathered_PortionOfMeshedPG ){
	  complex<double> RelativePos_=pp.Pos-p.Pos;
	  double Distance=abs(RelativePos_);
	  if(abs(RelativePos_.real())<MaxXYDist and abs(RelativePos_.imag())<MaxXYDist and Distance!=0)
	    AnsVec.push_back(RelativePos_);
	}
      }
    }

    return AnsVec;
  }
					  
  pair<vector<vector<double>>,vector<double>> AveragedPositionCorrelation(const vector<vector<PolarParticle2D>>& MeshedPG,
									  const vector<pair<vector<double>,vector<int>>>& BoxInfo,
									  int GridNumConsidered, /* an even number */
									  double BinSize, /* should be a factor of GridNumConsidered*CharLength */
									  int SampleRatio=1)
  {//Calculate the correlation function, in a GridNumConsidered*GridNumConsidered box, peripheral particles are discarded
   //return {Matrix of density, GridList},GridList is a vector with odd number of double

    double SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1];
    double SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3];
    int N_X=BoxInfo[0].second[0]+1;
    int N_Y=BoxInfo[0].second[2]/N_X+1;
    double CharLength=BoxInfo[0].first[1]-BoxInfo[0].first[0];
    double MaxXYDist=GridNumConsidered/2*CharLength;

    int N_Grids=GridNumConsidered*CharLength/BinSize;
    vector<vector<double>> AnsMat(N_Grids,vector<double>(N_Grids,0));

    //calculate GridList
    vector<double> GridList;
    for( int i=0 ; i<=N_Grids ; i++ )
      GridList.push_back(-GridNumConsidered/2*CharLength+i*BinSize);

    auto WhichGrid=[N_Grids,MaxXYDist,BinSize](complex<double> Coor)->pair<int,int>{
      double X=Coor.real()+MaxXYDist,Y=Coor.imag()+MaxXYDist;
      int nx=X/BinSize,ny=Y/BinSize;
      if(nx>=N_Grids)nx=nx-1;
      if(ny>=N_Grids)ny=ny-1;
      return {nx,ny};
    };

    auto IndexShift=[N_X,N_Y](int i,complex<int> Coor)->int{
      int
	n_X=i%N_X,
	n_Y=i/N_X,
	new_n_X=(n_X+Coor.real())%N_X,
	new_n_Y=(n_Y+Coor.imag())%N_Y;
      if(new_n_X<0) new_n_X+=N_X;
      if(new_n_Y<0) new_n_Y+=N_Y;
      return new_n_X+new_n_Y*N_X;
    };

    int cnt=0; //for controling the sample ratio
    int eff_cnt=0; //for count effective samples

    for( int i_Box=0 ; i_Box<BoxInfo.size() ; i_Box++ ){
      
      //discard the periphery
      int n_X=i_Box%N_X, n_Y=i_Box/N_X;
      if(n_X-GridNumConsidered/2<0 or n_X+GridNumConsidered/2>=N_X or n_Y-GridNumConsidered/2<0 or n_Y+GridNumConsidered/2>=N_Y)
	continue;

      //collect all particles in (GridNumConsidered+1)*(GridNumConsidered+1)
      vector<vector<PolarParticle2D>> PortionOfMeshedPG;
      for( int i_Shift_X=-GridNumConsidered/2 ; i_Shift_X<=GridNumConsidered/2 ; i_Shift_X++ )
	for( int i_Shift_Y=-GridNumConsidered/2 ; i_Shift_Y<=GridNumConsidered/2 ; i_Shift_Y++ )
	  PortionOfMeshedPG.push_back(MeshedPG[IndexShift(i_Box,complex<int>(i_Shift_X,i_Shift_Y))]);
      auto Gathered_PortionOfMeshedPG=BoxMeshing_2D::Gathering(PortionOfMeshedPG);

      //calculate relative position
      for( auto p : MeshedPG[i_Box] ){
	if((cnt++)%SampleRatio!=0) continue;
	eff_cnt++;
	
	for( auto pp : Gathered_PortionOfMeshedPG ){
	  complex<double> RelativePos_=pp.Pos-p.Pos;
	  double Distance=abs(RelativePos_);
	  if(abs(RelativePos_.real())<MaxXYDist and abs(RelativePos_.imag())<MaxXYDist and Distance!=0)
	    AnsMat[WhichGrid(RelativePos_).first][WhichGrid(RelativePos_).second]++;
	}
      }
    }

    //average 
    for( int i=0 ; i<N_Grids ; i++ )
      for( int j=0 ; j<N_Grids ; j++ )
	AnsMat[i][j]/=eff_cnt;
    
    return {AnsMat,GridList};
  }

  
  double DensityGradCompareAlongX(const PolarParticle2D_GT& PG,double SystemSize_X,double SystemSize_Y,double GridSize,
				  function<double(double)> ModulationFunc=[](double Diff)->double{return (Diff>0)?1:((Diff==0)?0:-1);})
  {//Calculate sum(ModulationFunc(Diff)), Diff=Density(x+Gridsize)-Density(x)
    double Ans=0;
    
    const auto MeshedPG=BoxMeshing_2D::Meshing(PG,SystemSize_X,SystemSize_Y,GridSize);
    const auto BoxInfo=BoxMeshing_2D::BoxInfo(SystemSize_X,SystemSize_Y,GridSize);
    
    for( int i_Box=0 ; i_Box<MeshedPG.size() ; i_Box++ ){
      double Density=MeshedPG[i_Box].size()/(GridSize*GridSize);
      int IndexNext=BoxInfo[i_Box].second[1];
      double DensityNext=MeshedPG[IndexNext].size()/(GridSize*GridSize);

      Ans+=ModulationFunc(DensityNext-Density);
    }

    return Ans;
  }
  
  //----------Not Recommended----------

  vector<pair<PolarParticle2D,PolarParticle2D>> ClosestNeighborList(vector<vector<PolarParticle2D>> MeshedPG,
								    vector<pair<vector<double>,vector<int>>> BoxInfo,
								    function<double(PolarParticle2D,PolarParticle2D)> DistanceOfTwoParticles)
  {//List out the closest particle of each particle: {{PolarParticle2D,Closest PolarParticle2D},...}
    vector<pair<PolarParticle2D,PolarParticle2D>> AnsVec;

    double CharLength=BoxInfo[0].first[1]-BoxInfo[0].first[0];
    int N_X=BoxInfo[0].second[0]+1;
    int N_Y=BoxInfo[0].second[2]/N_X+1;

    auto IndexShift=[N_X,N_Y](int i,complex<int> Coor)->int{
      int
	n_X=i%N_X,
	n_Y=i/N_X,
	new_n_X=(n_X+Coor.real())%N_X,
	new_n_Y=(n_Y+Coor.imag())%N_Y;
      if(new_n_X<0) new_n_X+=N_X;
      if(new_n_Y<0) new_n_Y+=N_Y;
      return new_n_X+new_n_Y*N_X;
    };
    
    for( int i_Box=0 ; i_Box<BoxInfo.size() ; i_Box++ )
      for( auto p : MeshedPG[i_Box] ){

	//find a near enough particle

	int cnt=0;
	auto Coor=_CoorPairList[cnt].first;
	while(MeshedPG[IndexShift(i_Box,Coor)].size()==0 or
	      (MeshedPG[IndexShift(i_Box,Coor)].size()==1 and DistanceOfTwoParticles(MeshedPG[IndexShift(i_Box,Coor)][0],p)==0)){
	  cnt++;
	  Coor=_CoorPairList[cnt].first;
	}

	PolarParticle2D NearEnoughP;
	if(DistanceOfTwoParticles(MeshedPG[IndexShift(i_Box,Coor)][0],p)!=0) NearEnoughP=MeshedPG[IndexShift(i_Box,Coor)][0];
	else NearEnoughP=MeshedPG[IndexShift(i_Box,Coor)][1];
	double dist=DistanceOfTwoParticles(NearEnoughP,p);
	
	//find a the closest particle

	for( int j_Box=cnt ; _CoorPairList[j_Box].second*CharLength<dist ; j_Box++ ){
	  auto Coor=_CoorPairList[j_Box].first;
	  
	  for( auto pp : MeshedPG[IndexShift(i_Box,Coor)] )
	    if(DistanceOfTwoParticles(pp,p)<dist and DistanceOfTwoParticles(pp,p)!=0){
	      dist=DistanceOfTwoParticles(pp,p);
	      NearEnoughP=pp;
	    }
	}

	AnsVec.push_back({p,NearEnoughP});
      }

    return AnsVec;
  }

};


#ifdef CheckingHeader
HeaderChecker(DataProcessing_hpp);
#endif //CheckingHeader

#endif //DataProcessing_hpp
