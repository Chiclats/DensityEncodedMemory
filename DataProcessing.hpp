/*--------------------
ver 250630
A collection of data processing functions.
--------------------*/

#ifndef DataProcessing_hpp
#define DataProcessing_hpp

#ifndef vs
#include<bits/stdc++.h>
#endif

#include"BasicTools.hpp"
#include"AgentRealization.hpp"
#include"StateProcessing.hpp"

#ifndef LocalProgram
#include<omp.h>
#endif//LocalProgram

using namespace std;

namespace DataProcessing
{//Tag processing: manipulations related to PolarParticle2D:.Tag

  vector<vector<PolarParticle2D>> ParticleGroupDivision(const vector<PolarParticle2D>& PG,int TagNum=-1)
  {//Divide a PG with particles of different tags into several PG, Tag should always start from 1 and be integers
   //-1 denotes to automatically detect TagNum

    if(TagNum==-1) // detect TagNum
      for( auto p : PG )
	TagNum=(p.Tag>TagNum)?p.Tag:TagNum;
    
    if(TagNum==1)
      return {PG};
    
    vector<vector<PolarParticle2D>> ansVecPG(TagNum);
    for( auto p : PG )
      ansVecPG[p.Tag-1].push_back(p);
      
    return ansVecPG;
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

};

namespace DataProcessing
{//Order parameter calculating

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

};

namespace DataProcessing
{//closest neighbors calculating
  
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

};

namespace DataProcessing::PositionCorrelationFunctions
{//two-poing correlation functions

  complex<double> RelativePos_Periodic(PolarParticle2D p1,PolarParticle2D p2,double SystemSize_X,double SystemSize_Y)
  {//Relative Position p2-p1, periodic boundary conditions
    double
      dx=p2.Pos.real()-p1.Pos.real(),
      dy=p2.Pos.imag()-p1.Pos.imag();

    if(abs(dx-SystemSize_X)<abs(dx) and abs(dx-SystemSize_X)<abs(dx+SystemSize_X))
      dx-=SystemSize_X;
    else if(abs(dx+SystemSize_X)<abs(dx))
      dx+=SystemSize_X;

    if(abs(dy-SystemSize_Y)<abs(dy) and abs(dy-SystemSize_Y)<abs(dy+SystemSize_Y))
      dy-=SystemSize_Y;
    else if(abs(dy+SystemSize_Y)<abs(dy))
      dy+=SystemSize_Y;

    return dx+ii*dy;
  }

  complex<double> RotateRelativePos(complex<double> RelativePos,double Theta)
  {//rotate the relative position counterclockwise
    return RelativePos*exp(ii*Theta);
  }

  int IndexShift(const int& index_From,const int& DX,const int& DY,const int& BoxNum_X,const int& BoxNum_Y)
  {//in the BoxInfo system, decide the box index from index_From + DX + DY * ii
   //ATTENTION: |DX|<BoxNum_X and |DY|<BoxNum_Y
    int
	n_X=index_From%BoxNum_X,
	n_Y=index_From/BoxNum_X,
	new_n_X=(n_X+DX+BoxNum_X)%BoxNum_X,
	new_n_Y=(n_Y+DY+BoxNum_Y)%BoxNum_Y;
      return new_n_X+new_n_Y*BoxNum_X;
  }
  
  vector<complex<double>> RelativePosListAroundAPosition(const complex<double>& Center,
							 const vector<vector<PolarParticle2D>>& MeshedPG,
							 const vector<pair<vector<double>,vector<int>>>& BoxInfo,
							 const tuple<double,double,double,int,int>& Parameters,
							 //{SystemSize_X,SystemSize_Y,CharLength,BoxNum_X,BoxNum_Y}
							 const int& GridNumConsidered /* an even number */)
  {//list out relative positions of neighboring particles around Center
    vector<complex<double>> AnsVec;

    const auto [SystemSize_X,SystemSize_Y,CharLength,BoxNum_X,BoxNum_Y]=Parameters;
    const PolarParticle2D p_Center(Center,0.); //an imaginary particle at Center
    const int index_From=BoxMeshing_2D::DecideBoxIndex(p_Center,SystemSize_X,SystemSize_Y,CharLength);
    const double MaxXYDist=GridNumConsidered/2*CharLength;

    //collect all particles in (GridNumConsidered+1)*(GridNumConsidered+1)
    vector<vector<PolarParticle2D>> PortionOfMeshedPG;
    for( int i_Shift_X=-GridNumConsidered/2 ; i_Shift_X<=GridNumConsidered/2 ; i_Shift_X++ )
      for( int i_Shift_Y=-GridNumConsidered/2 ; i_Shift_Y<=GridNumConsidered/2 ; i_Shift_Y++ )
	PortionOfMeshedPG.push_back(MeshedPG[IndexShift(index_From,i_Shift_X,i_Shift_Y,BoxNum_X,BoxNum_Y)]);
    auto Gathered_PortionOfMeshedPG=BoxMeshing_2D::Gathering(PortionOfMeshedPG);

    //calculate relative position	
    for( auto pp : Gathered_PortionOfMeshedPG ){
      complex<double> RelativePos_=RelativePos_Periodic(p_Center,pp,SystemSize_X,SystemSize_Y);
      double Distance=abs(RelativePos_);
      if(abs(RelativePos_.real())<MaxXYDist and abs(RelativePos_.imag())<MaxXYDist and Distance!=0)
	AnsVec.push_back(RelativePos_);
    }

    return AnsVec;
  }

  vector<vector<double>> PositionCorrelation_TwoPG(const PolarParticle2D_GT& PG_Center,
						   const PolarParticle2D_GT& PG_Target,
						   const tuple<double,double,double,double,double>& Parameters,
						   //GridSize RangeConsidered SystemSize_X SystemSize_Y CharLength(a factor of RangeConsidered)
						   const double& RotatingAngle=0,
						   const int& SamplePer=1 /*a sample per SamplePer particles*/)
  {//calculate p(r in PG_Target | r in PG_Center)
    
    const auto [GridSize,RangeConsidered,SystemSize_X,SystemSize_Y,CharLength]=Parameters;
    const auto BoxInfo=BoxMeshing_2D::BoxInfo(SystemSize_X,SystemSize_Y,CharLength);
    const auto Parameters_Inner=BoxMeshing_2D::UnzipBoxInfo(BoxInfo);
    const auto MeshedPG_Target=BoxMeshing_2D::Meshing(PG_Target,SystemSize_X,SystemSize_Y,CharLength);
    const int GridNumConsidered=RangeConsidered/CharLength;

    //runtime error checking
    if(RangeConsidered/2/GridSize-floor(RangeConsidered/2/GridSize)>1e-10){
      cerr<<"RangeConsidered: "<<RangeConsidered<<" GridSize: "<<GridSize<<endl;
      throw(runtime_error("GridSize is not a factor of RangeConsidered / 2"));
    }
    if(RangeConsidered/2/CharLength-floor(RangeConsidered/2/CharLength)>1e-10){
      cerr<<"RangeConsidered: "<<RangeConsidered<<" CharLength: "<<CharLength<<endl;
      throw(runtime_error("CharLength is not a factor of RangeConsidered / 2"));
    }

    //counting RelativePos
    vector<vector<double>> ansMat(int(RangeConsidered/GridSize),vector<double>(int(RangeConsidered/GridSize),0));
    
    for( int i_Center=0 ; i_Center<PG_Center.size() ; i_Center+=SamplePer ){
      auto RelativePosList=RelativePosListAroundAPosition(PG_Center[i_Center].Pos,MeshedPG_Target,BoxInfo,Parameters_Inner,GridNumConsidered);

      for( auto RelativePos : RelativePosList ){
	if(norm(RelativePos)>=RangeConsidered*RangeConsidered/4.) continue;

	RelativePos=RotateRelativePos(RelativePos,RotatingAngle);
	const int Index_X=(RelativePos.real()+RangeConsidered/2)/GridSize,Index_Y=(RelativePos.imag()+RangeConsidered/2)/GridSize;
	ansMat[Index_X][Index_Y]++;
      }
    }

    //rescaling to pdf
    const int NumSampledParticles=ceil(PG_Center.size()/SamplePer);
    for( int i_X=0 ; i_X<ansMat.size() ; i_X++ )
      for( int i_Y=0 ; i_Y<ansMat[0].size() ; i_Y++ )
	ansMat[i_X][i_Y]/=(NumSampledParticles*GridSize*GridSize);

    return ansMat;
  }
		   

};

namespace DataProcessing
{//Misc
  vector<complex<double>> RelativePosList(const vector<vector<PolarParticle2D>>& MeshedPG,
					  const vector<pair<vector<double>,vector<int>>>& BoxInfo,
					  const int& GridNumConsidered, /* an even number */
					  const int& SampleRatio /*a sample per SampleRatio particles*/) 
  {//List out the relative position of each two particles, in a GridNumConsidered*GridNumConsidered box, peripheral particles are discarded
    vector<complex<double>> AnsVec;

    auto [SystemSize_X,SystemSize_Y,CharLength,N_X,N_Y]=BoxMeshing_2D::UnzipBoxInfo(BoxInfo);
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

  vector<int> StatDirDistribution(const PolarParticle2D_GT& PG,const int& BinNum)
  {//Divide 2*pi into BinNum bins, and count the number of particles whose Dir is in a bin
    //Attention: 0 <= Dir <= 2*pi
    const double BinSize=2*pi/BinNum;
    vector<int> AnsVec(BinNum,0);
    int n;
    for( auto p : PG ){
      if(p.Dir<0 or 2*pi<p.Dir){
	cerr<<'('<<p.Pos.real()<<','<<p.Pos.imag()<<')'<<p.Dir<<endl;
	throw(runtime_error("In StatDirDistribution: Illegal value of Dir."));
      }
			  
      n=floor(p.Dir/BinSize);
      AnsVec[(n>=BinNum)?(BinNum-1):n]++;
    }

    return AnsVec;
  }

  vector<vector<int>> StatDirJumpingMat(const PolarParticle2D_GT& PG1,const PolarParticle2D_GT& PG2,const int& BinNum)
  {//Divide 2*pi into BinNum bins, and count the number of particles who in PG1 at Bin_i and in PG2 at Bin_j
    //Attention: 0 <= Dir <= 2*pi
    if(PG1.size()!=PG2.size())
      throw(runtime_error("In StatDirJumpingMat: PG1 and PG2 is of different length."));
  
    const double BinSize=2*pi/BinNum;
    vector<vector<int>> AnsMat(BinNum,vector<int>(BinNum,0));

    int i,j;
    for( int i_p=0 ; i_p<PG1.size() ; i_p++ ){
      if(PG1[i_p].Dir<0 or 2*pi<PG1[i_p].Dir or PG2[i_p].Dir<0 or 2*pi<PG2[i_p].Dir)
	throw(runtime_error("In StatDirJumpingMat: Illegal value of Dir."));

      i=floor(PG1[i_p].Dir/BinSize),j=floor(PG2[i_p].Dir/BinSize);
      AnsMat[(i>=BinNum)?(BinNum-1):i][(j>=BinNum)?(BinNum-1):j]++;
    }
  
    return AnsMat;
  }

};


#ifdef CheckingHeader
HeaderChecker(DataProcessing_hpp);
#endif //CheckingHeader

#endif //DataProcessing_hpp
