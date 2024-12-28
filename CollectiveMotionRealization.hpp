/*--------------------
ver 241222
--------------------*/

#ifndef CollectiveMotionRealization_hpp
#define CollectiveMotionRealization_hpp

#include<bits/stdc++.h>
#include<omp.h>
#include"BasicTools.hpp"
#include"AgentRealization.hpp"
#include"StateProcessing.hpp"
#include"DataProcessing.hpp"

using namespace std;

using CouplingT=function<PolarParticle2D(PolarParticle2D,PolarParticle2D)>;
using SelfEvolutionT=function<PolarParticle2D(PolarParticle2D,double)>;
using BoundaryT=function<PolarParticle2D(PolarParticle2D)>;
using DistanceT=function<double(PolarParticle2D,PolarParticle2D)>;

Meshed_PolarParticle2D_GT EvolveMeshedList_BoundaryCond_PolarParticle2D(const Meshed_PolarParticle2D_GT& list, const BoxInfoT& BoxInfo,
									const double CharLength,               const double TimeStep,
									const CouplingT& Interaction,          const SelfEvolutionT& SelfEvolution,
									const BoundaryT& BoundaryCondFunc,     const DistanceT& DistanceOfTwoParticles)
{/*To evaluate the configuration of a collective system after a timestep, The boundary condition is inputed as a function
   BoxInfo: result of BoxMeshing_2D_BoxInfo()
   Interaction: a function, with two agents of Type. Return the affect of the second agent on the first in the form of increase RATE.
   SelfEvolution: a function, with an agent of Type and a timestep. Return the self evolution of one agent in the form of increase VALUE.
   BoundaryCondition: a function to move particles out of the system back.
   Distanceoftwoparticles: a function to return the distance of two particles
  */
  double SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1];
  double SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3];

  vector<vector<PolarParticle2D>> AnsList(BoxInfo.size());
  
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ ){
    array<int,8> AdjacentBoxList={
      BoxInfo[BoxIndex].second[0],BoxInfo[BoxIndex].second[1],
      BoxInfo[BoxIndex].second[2],BoxInfo[BoxIndex].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[3],
    };
    //A list of index of adjacent boxes
    //[Warning] Assuming the box matrix has more than 3 rows and 3 collumns.

    for(int ip=0 ; ip<list[BoxIndex].size(); ip++){
      auto p=list[BoxIndex][ip];
      auto dp=SelfEvolution(p,TimeStep);

      PolarParticle2D dpInter;
      for(int AdjBoxIndex : AdjacentBoxList)
	for( auto p2 : list[AdjBoxIndex])
	  if(DistanceOfTwoParticles(p,p2)<=CharLength){
	    dpInter=Interaction(p,p2);
	    dp.Pos+=dpInter.Pos*TimeStep;
	    dp.Dir+=dpInter.Dir*TimeStep;
	  }

      PolarParticle2D p2;
      for( int i_sameBox=0 ; i_sameBox<list[BoxIndex].size() ; i_sameBox++ ){
	p2=list[BoxIndex][i_sameBox];
	if(i_sameBox!=ip and DistanceOfTwoParticles(p,p2)<=CharLength){
	  dpInter=Interaction(p,p2);
	  dp.Pos+=dpInter.Pos*TimeStep;
	  dp.Dir+=dpInter.Dir*TimeStep;
	}
      }
      
      dp.Pos+=p.Pos;
      dp.Dir+=p.Dir;
      dp.Tag=p.Tag;
      dp.DoubleIdentifier=p.DoubleIdentifier;

      auto dp_Processed=BoundaryCondFunc(dp);
      AnsList[BoxMeshing_2D::DecideBoxIndex(dp_Processed,SystemSize_X,SystemSize_Y,CharLength)].push_back(dp_Processed);
    }
  }

  return AnsList;
}

Meshed_PolarParticle2D_GT EvolveMeshedList_PolarParticle2D_VM(const Meshed_PolarParticle2D_GT& list,   const BoxInfoT& BoxInfo,
							      const double CharLength,                 const double Velocity,
							      const double Sigma,                      const BoundaryT& BoundaryCondFunc,
							      const DistanceT& DistanceOfTwoParticles, const double GuidingForceCoeff=0,
							      const int i_MP=0)
{/*For one species VM only!!!
   BoxInfo: result of BoxMeshing_2D_BoxInfo()
   BoundaryCondition: a function to move particles out of the system back.
   Distanceoftwoparticles: a function to return the distance of two particles
   GuidingForceCoeff: the intensity of guiding force in \pm y direction
 */
  double SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1];
  double SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3];

  vector<vector<PolarParticle2D>> AnsList(BoxInfo.size());
  
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ ){
    array<int,8> AdjacentBoxList={
      BoxInfo[BoxIndex].second[0],BoxInfo[BoxIndex].second[1],
      BoxInfo[BoxIndex].second[2],BoxInfo[BoxIndex].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[3],
    };
    //A list of index of adjacent boxes
    //[Warning] Assuming the box matrix has more than 3 rows and 3 collumns.

    for(int ip=0 ; ip<list[BoxIndex].size(); ip++){
      auto p=list[BoxIndex][ip];

      complex<double> SumDirComplex(0,0);
      for(int AdjBoxIndex : AdjacentBoxList)
	for( auto p2 : list[AdjBoxIndex])
	  if(DistanceOfTwoParticles(p,p2)<=CharLength)
	    SumDirComplex+=exp(ii*p2.Dir);

      PolarParticle2D p2;
      for( int i_sameBox=0 ; i_sameBox<list[BoxIndex].size() ; i_sameBox++ ){
	p2=list[BoxIndex][i_sameBox];
	if(DistanceOfTwoParticles(p,p2)<=CharLength)//the particle itself is also taken into account
	    SumDirComplex+=exp(ii*p2.Dir);
      }

      PolarParticle2D dp;
      dp.Dir=arg(SumDirComplex)+Sigma*(RandGen[i_MP].RandomDouble()*2*pi-pi);   
      dp.Dir+=GuidingForceCoeff*sin(2*dp.Dir);
      dp.Pos=p.Pos+exp(ii*dp.Dir)*Velocity;
      
      auto dp_Processed=BoundaryCondFunc(dp);
      AnsList[BoxMeshing_2D::DecideBoxIndex(dp_Processed,SystemSize_X,SystemSize_Y,CharLength)].push_back(dp_Processed);
    }
  }

  return AnsList;
}

Meshed_PolarParticle2D_GT EvolveMeshedList_PolarParticle2D_VM_FG(const Meshed_PolarParticle2D_GT& list,   const BoxInfoT& BoxInfo,
								 const double CharLength,                 const double Velocity,
								 const double Sigma,                      const BoundaryT& BoundaryCondFunc,
								 const DistanceT& DistanceOfTwoParticles, const double GuidingForceCoeff=0,
								 const int i_MP=0)
{/*For one species VM only!!! Guiding force will be calculated from the previous direction (matching VicsekModel.hpp)
   BoxInfo: result of BoxMeshing_2D_BoxInfo()
   BoundaryCondition: a function to move particles out of the system back.
   Distanceoftwoparticles: a function to return the distance of two particles
   GuidingForceCoeff: the intensity of guiding force in \pm y direction
 */
  double SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1];
  double SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3];

  vector<vector<PolarParticle2D>> AnsList(BoxInfo.size());
  
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ ){
    array<int,8> AdjacentBoxList={
      BoxInfo[BoxIndex].second[0],BoxInfo[BoxIndex].second[1],
      BoxInfo[BoxIndex].second[2],BoxInfo[BoxIndex].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[3],
    };
    //A list of index of adjacent boxes
    //[Warning] Assuming the box matrix has more than 3 rows and 3 collumns.

    for(int ip=0 ; ip<list[BoxIndex].size(); ip++){
      auto p=list[BoxIndex][ip];

      complex<double> SumDirComplex(0,0);
      for(int AdjBoxIndex : AdjacentBoxList)
	for( auto p2 : list[AdjBoxIndex])
	  if(DistanceOfTwoParticles(p,p2)<=CharLength)
	    SumDirComplex+=exp(ii*p2.Dir);

      PolarParticle2D p2;
      for( int i_sameBox=0 ; i_sameBox<list[BoxIndex].size() ; i_sameBox++ ){
	p2=list[BoxIndex][i_sameBox];
	if(DistanceOfTwoParticles(p,p2)<=CharLength)//the particle itself is also taken into account
	    SumDirComplex+=exp(ii*p2.Dir);
      }

      PolarParticle2D dp;
      dp.Dir=arg(SumDirComplex)+Sigma*(RandGen[i_MP].RandomDouble()*2*pi-pi);   
      dp.Dir+=GuidingForceCoeff*sin(2*p.Dir);//difference here
      dp.Pos=p.Pos+exp(ii*dp.Dir)*Velocity;
      
      auto dp_Processed=BoundaryCondFunc(dp);
      AnsList[BoxMeshing_2D::DecideBoxIndex(dp_Processed,SystemSize_X,SystemSize_Y,CharLength)].push_back(dp_Processed);
    }
  }

  return AnsList;
}

Meshed_PolarParticle2D_GT EvolveMeshedList_PolarParticle2D_VM_VecFluc(const Meshed_PolarParticle2D_GT& list,
								      const BoxInfoT& BoxInfo,
								      const double Velocity,
								      const double Sigma,
								      const BoundaryT& BoundaryCondFunc,
								      const DistanceT& DistanceOfTwoParticles,
								      const double GuidingForceCoeff=0,
								      const int i_MP=0)
{/*For one species VM with vectorial fluctuation only!!!
   BoxInfo: result of BoxMeshing_2D_BoxInfo()
   BoundaryCondition: a function to move particles out of the system back.
   Distanceoftwoparticles: a function to return the distance of two particles
   GuidingForceCoeff: the intensity of guiding force in \pm y direction
 */
  double
    SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1],
    SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3],
    CharLength=BoxInfo[0].first[1]-BoxInfo[0].first[0];

  Meshed_PolarParticle2D_GT AnsList(BoxInfo.size());

  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ ){
    array<int,8> AdjacentBoxList={
      BoxInfo[BoxIndex].second[0],BoxInfo[BoxIndex].second[1],
      BoxInfo[BoxIndex].second[2],BoxInfo[BoxIndex].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[3],
    };
    //A list of index of adjacent boxes
    //[Warning] Assuming the box matrix has more than 3 rows and 3 collumns.

    for(int ip=0 ; ip<list[BoxIndex].size(); ip++){
      auto p=list[BoxIndex][ip];

      complex<double> SumDirComplex(0,0);
      int SumDirCount=0;
      
      for(int AdjBoxIndex : AdjacentBoxList)
	for( auto p2 : list[AdjBoxIndex])
	  if(DistanceOfTwoParticles(p,p2)<=CharLength)
	    SumDirCount++,
	    SumDirComplex+=exp(ii*p2.Dir);

      PolarParticle2D p2;
      for( int i_sameBox=0 ; i_sameBox<list[BoxIndex].size() ; i_sameBox++ ){
	p2=list[BoxIndex][i_sameBox];
	if(DistanceOfTwoParticles(p,p2)<=CharLength)//the particle itself is also taken into account
	  SumDirCount++,
	  SumDirComplex+=exp(ii*p2.Dir);
      }

      PolarParticle2D dp;
      dp.Dir=arg(SumDirComplex+SumDirCount*Sigma*exp(ii*RandGen[i_MP].RandomDouble()*2.*pi));
      dp.Dir+=GuidingForceCoeff*sin(2*dp.Dir);
      dp.Pos=p.Pos+exp(ii*dp.Dir)*Velocity;
      
      auto dp_Processed=BoundaryCondFunc(dp);
      AnsList[BoxMeshing_2D::DecideBoxIndex(dp_Processed,SystemSize_X,SystemSize_Y,CharLength)].push_back(dp_Processed);
    }
  }

  return AnsList;
}

Meshed_PolarParticle2D_GT EvolveMeshedList_PolarParticle2D_RepulsiveVM(const Meshed_PolarParticle2D_GT& list,
								       const BoxInfoT& BoxInfo,     //CharLength is obtained from BoxInfo
								       const double AligningLength, //Radius of the aligning zone, less than CharLength in BoxInfo
								       const double Velocity,
								       const double Sigma,
								       const function<double(double)> RepulsiveForce, // only taken into account in CharLength
								       const BoundaryT& BoundaryCondFunc,
								       const DistanceT& DistanceOfTwoParticles,
								       const double GuidingForceCoeff=0,
								       const int& i_MP=0)
{/*For one species repulsive VM only!!!
   CharLength is obtained from BoxInfo, AligningLength must be less than CharLength
   RepulsiveForce: radius -> force, positive repulsive
   GuidingForceCoeff: the intensity of guiding force in \pm y direction
 */
  double SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1];
  double SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3];
  double CharLength=BoxInfo[0].first[1]-BoxInfo[0].first[0];

  Meshed_PolarParticle2D_GT AnsList(BoxInfo.size());

  auto EffRelativePos=[SystemSize_X,SystemSize_Y](complex<double> RelativePos)->complex<double>{
    double EffRelX=RelativePos.real(),EffRelY=RelativePos.imag();
    if(not(abs(EffRelX)<=abs(EffRelX+SystemSize_X) and abs(EffRelX)<=abs(EffRelX-SystemSize_X))){
      if(abs(EffRelX+SystemSize_X)<=abs(EffRelX-SystemSize_X)) EffRelX=EffRelX+SystemSize_X;
      else EffRelX=EffRelX-SystemSize_X;
    }
    if(not(abs(EffRelY)<=abs(EffRelY+SystemSize_Y) and abs(EffRelY)<=abs(EffRelY-SystemSize_Y))){
      if(abs(EffRelY+SystemSize_Y)<=abs(EffRelY-SystemSize_Y)) EffRelY=EffRelY+SystemSize_Y;
      else EffRelY=EffRelY-SystemSize_Y;
    }
    return EffRelX+ii*EffRelY;
  };// do deal with the relative position on the periphery
  
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ ){
    array<int,8> AdjacentBoxList={
      BoxInfo[BoxIndex].second[0],BoxInfo[BoxIndex].second[1],
      BoxInfo[BoxIndex].second[2],BoxInfo[BoxIndex].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[3],
    };
    //A list of index of adjacent boxes
    //[Warning] Assuming the box matrix has more than 3 rows and 3 collumns.

    for(int ip=0 ; ip<list[BoxIndex].size(); ip++){
      auto p=list[BoxIndex][ip];

      complex<double> SumDirComplex(0,0);
      double Distance;
      for(int AdjBoxIndex : AdjacentBoxList)
	for( auto p2 : list[AdjBoxIndex]){
	  Distance=DistanceOfTwoParticles(p,p2);
	  
	  if(Distance<=AligningLength)
	    SumDirComplex+=exp(ii*p2.Dir); //Alignment interaction

	  if(Distance<=CharLength)
	    SumDirComplex+=RepulsiveForce(Distance)*exp(ii*arg(EffRelativePos(p.Pos-p2.Pos))); //Repulsive interaction
	}

      PolarParticle2D p2;
      for( int i_sameBox=0 ; i_sameBox<list[BoxIndex].size() ; i_sameBox++ ){
	p2=list[BoxIndex][i_sameBox];
	Distance=DistanceOfTwoParticles(p,p2);
	
	if(Distance<=AligningLength)//the particle itself is also taken into account
	  SumDirComplex+=exp(ii*p2.Dir); //Alignment Interaction

	if(Distance<=CharLength and i_sameBox!=ip)
	  SumDirComplex+=RepulsiveForce(Distance)*exp(ii*arg(EffRelativePos(p.Pos-p2.Pos))); //Repulsive interaction
      }

      PolarParticle2D dp;
      dp.Dir=arg(SumDirComplex)+Sigma*(RandGen[i_MP].RandomDouble()*2*pi-pi);
      dp.Dir+=GuidingForceCoeff*sin(2*dp.Dir);
      dp.Pos=p.Pos+exp(ii*dp.Dir)*Velocity;
      
      auto dp_Processed=BoundaryCondFunc(dp);
      AnsList[BoxMeshing_2D::DecideBoxIndex(dp_Processed,SystemSize_X,SystemSize_Y,CharLength)].push_back(dp_Processed);
    }
  }

  return AnsList;
}


//----------MultiProcess codes----------

#ifndef LocalProgram //not compiled when in personal pc

Meshed_PolarParticle2D_GT EvolveMeshedList_PolarParticle2D_MP(const Meshed_PolarParticle2D_GT& list,
							      const BoxInfoT& BoxInfo,
							      const double CharLength,
							      const double TimeStep,
							      const CouplingT& Interaction,
							      const SelfEvolutionT& SelfEvolution,
							      const BoundaryT& BoundaryCondFunc,
							      const DistanceT& DistanceOfTwoParticles)
{/*To evaluate the configuration of a collective system after a timestep, The boundary condition is inputed as a function
   BoxInfo: result of BoxMeshing_2D_BoxInfo()
   Interaction: a function, with two agents of Type. Return the affect of the second agent on the first in the form of increase RATE.
   SelfEvolution: a function, with an agent of Type and a timestep. Return the self evolution of one agent in the form of increase VALUE.
   BoundaryCondition: a function to move particles out of the system back.
   Distanceoftwoparticles: a function to return the distance of two particles
 */
  double SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1];
  double SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3];

  vector<vector<PolarParticle2D>> AnsList(BoxInfo.size());

  #pragma omp parallel for
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ ){
    vector<int> AdjacentBoxList=BoxInfo[BoxIndex].second;
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[0]].second[2]);
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[0]].second[3]);
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[1]].second[2]);
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[1]].second[3]);
    //A list of index of adjacent boxes
    //[Warning] Assuming the box matrix has more than 3 rows and 3 collumns.

    for(int ip=0 ; ip<list[BoxIndex].size(); ip++){
      auto p=list[BoxIndex][ip];
      auto dp=SelfEvolution(p,TimeStep);

      PolarParticle2D dpInter;
      for(int AdjBoxIndex : AdjacentBoxList)
	for( auto p2 : list[AdjBoxIndex])
	  if(DistanceOfTwoParticles(p,p2)<=CharLength){
	    dpInter=Interaction(p,p2);
	    dp.Pos+=dpInter.Pos*TimeStep;
	    dp.Dir+=dpInter.Dir*TimeStep;
	  }

      PolarParticle2D p2;
      for( int i_sameBox=0 ; i_sameBox<list[BoxIndex].size() ; i_sameBox++ ){
	p2=list[BoxIndex][i_sameBox];
	if(i_sameBox!=ip and DistanceOfTwoParticles(p,p2)<=CharLength){
	  dpInter=Interaction(p,p2);
	  dp.Pos+=dpInter.Pos*TimeStep;
	  dp.Dir+=dpInter.Dir*TimeStep;
	}
      }
      
      dp.Pos+=p.Pos;
      dp.Dir+=p.Dir;
      dp.Tag=p.Tag;
      dp.DoubleIdentifier=p.DoubleIdentifier;

      auto dp_Processed=BoundaryCondFunc(dp);
      AnsList[BoxIndex].push_back(dp_Processed);
    }
  }

  vector<vector<PolarParticle2D>> FnlAnsList(BoxInfo.size());
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ )
    for( auto p : AnsList[BoxIndex] )
      FnlAnsList[BoxMeshing_2D::DecideBoxIndex(p,SystemSize_X,SystemSize_Y,CharLength)].push_back(p);
  
  return FnlAnsList;
}

Meshed_PolarParticle2D_GT EvolveMeshedList_PolarParticle2D_VM_MP(const Meshed_PolarParticle2D_GT& list,
								 const BoxInfoT& BoxInfo,
								 const double CharLength,
								 const double Velocity,
								 const double Sigma,
								 const BoundaryT& BoundaryCondFunc,
								 const DistanceT& DistanceOfTwoParticles,
								 const double GuidingForceCoeff=0)
{/*For one species VM only!!! MP ver
   BoxInfo: result of BoxMeshing_2D_BoxInfo()
   BoundaryCondition: a function to move particles out of the system back.
   Distanceoftwoparticles: a function to return the distance of two particles
   GuidingForceCoeff: the intensity of guiding force in \pm y direction
 */
  double SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1];
  double SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3];

  auto AnsList=list;

  //omp_set_num_threads(MPNum);
  #pragma omp parallel for
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ ){
    array<int,8> AdjacentBoxList={
      BoxInfo[BoxIndex].second[0],BoxInfo[BoxIndex].second[1],
      BoxInfo[BoxIndex].second[2],BoxInfo[BoxIndex].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[0]].second[3],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[2],
      BoxInfo[BoxInfo[BoxIndex].second[1]].second[3],
    };
    //A list of index of adjacent boxes
    //[Warning] Assuming the box matrix has more than 3 rows and 3 collumns.

    for(int ip=0 ; ip<list[BoxIndex].size(); ip++){
      auto p=list[BoxIndex][ip];

      complex<double> SumDirComplex(0,0);
      for(int AdjBoxIndex : AdjacentBoxList)
	for( auto p2 : list[AdjBoxIndex])
	  if(DistanceOfTwoParticles(p,p2)<=CharLength)
	    SumDirComplex+=exp(ii*p2.Dir);

      PolarParticle2D p2;
      for( int i_sameBox=0 ; i_sameBox<list[BoxIndex].size() ; i_sameBox++ ){
	p2=list[BoxIndex][i_sameBox];
	if(DistanceOfTwoParticles(p,p2)<=CharLength)//the particle itself is also taken into account
	    SumDirComplex+=exp(ii*p2.Dir);
      }

      PolarParticle2D dp;
      dp.Dir=arg(SumDirComplex)+Sigma*(RandGen[omp_get_thread_num()].RandomDouble()*2*pi-pi);
      //dp.Dir=arg(SumDirComplex)+Sigma*(RandGen[0].RandomDouble()*2*pi-pi);
      dp.Dir+=GuidingForceCoeff*sin(2*dp.Dir);
      dp.Pos=p.Pos+exp(ii*dp.Dir)*Velocity;
      
      auto dp_Processed=BoundaryCondFunc(dp);
      AnsList[BoxIndex][ip]=dp_Processed;
    }
  }

  vector<vector<PolarParticle2D>> FnlAnsList(BoxInfo.size());
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ )
    for( auto p : AnsList[BoxIndex] )
      FnlAnsList[BoxMeshing_2D::DecideBoxIndex(p,SystemSize_X,SystemSize_Y,CharLength)].push_back(p);

  return FnlAnsList;
}

Meshed_PolarParticle2D_GT EvolveMeshedList_PolarParticle2D_VM_VecFluc_MP(const Meshed_PolarParticle2D_GT& list,
									 const BoxInfoT& BoxInfo,
									 const double CharLength,
									 const double Velocity,
									 const double Sigma,
									 const BoundaryT& BoundaryCondFunc,
									 const DistanceT& DistanceOfTwoParticles,
									 const double GuidingForceCoeff=0)
{/*For one species VM with vectorial fluctuation only!!! MP ver
   BoxInfo: result of BoxMeshing_2D_BoxInfo()
   BoundaryCondition: a function to move particles out of the system back.
   Distanceoftwoparticles: a function to return the distance of two particles
   GuidingForceCoeff: the intensity of guiding force in \pm y direction
 */
  double SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1];
  double SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3];

  vector<vector<PolarParticle2D>> AnsList(BoxInfo.size());

  #pragma omp parallel for
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ ){
    vector<int> AdjacentBoxList=BoxInfo[BoxIndex].second;
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[0]].second[2]);
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[0]].second[3]);
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[1]].second[2]);
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[1]].second[3]);
    //A list of index of adjacent boxes
    //[Warning] Assuming the box matrix has more than 3 rows and 3 collumns.

    for(int ip=0 ; ip<list[BoxIndex].size(); ip++){
      auto p=list[BoxIndex][ip];

      complex<double> SumDirComplex(0,0);
      int SumDirCount=0;
      
      for(int AdjBoxIndex : AdjacentBoxList)
	for( auto p2 : list[AdjBoxIndex])
	  if(DistanceOfTwoParticles(p,p2)<=CharLength)
	    SumDirCount++,
	    SumDirComplex+=exp(ii*p2.Dir);

      PolarParticle2D p2;
      for( int i_sameBox=0 ; i_sameBox<list[BoxIndex].size() ; i_sameBox++ ){
	p2=list[BoxIndex][i_sameBox];
	if(DistanceOfTwoParticles(p,p2)<=CharLength)//the particle itself is also taken into account
	  SumDirCount++,
	  SumDirComplex+=exp(ii*p2.Dir);
      }

      PolarParticle2D dp;
      dp.Dir=arg(SumDirComplex+SumDirCount*Sigma*exp(ii*RandGen[omp_get_thread_num()].RandomDouble()*2.*pi));
      //dp.Dir=arg(SumDirComplex+SumDirCount*Sigma*exp(ii*RandGen[0].RandomDouble()*2.*pi));
      //vectorial fluctuation employed here
      dp.Dir+=GuidingForceCoeff*sin(2*dp.Dir);
      dp.Pos=p.Pos+exp(ii*dp.Dir)*Velocity;
      
      auto dp_Processed=BoundaryCondFunc(dp);
      AnsList[BoxIndex].push_back(dp_Processed);
    }
  }

  vector<vector<PolarParticle2D>> FnlAnsList(BoxInfo.size());
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ )
    for( auto p : AnsList[BoxIndex] )
      FnlAnsList[BoxMeshing_2D::DecideBoxIndex(p,SystemSize_X,SystemSize_Y,CharLength)].push_back(p);

  return FnlAnsList;
}

Meshed_PolarParticle2D_GT EvolveMeshedList_PolarParticle2D_RepulsiveVM_MP(const Meshed_PolarParticle2D_GT& list,
									  const BoxInfoT& BoxInfo,     //CharLength is obtained from BoxInfo
									  const double AligningLength, //Radius of the aligning zone, less than CharLength in BoxInfo
									  const double Velocity,
									  const double Sigma,
									  const function<double(double)> RepulsiveForce, // only taken into account in CharLength
									  const BoundaryT& BoundaryCondFunc,
									  const DistanceT& DistanceOfTwoParticles,
									  const double GuidingForceCoeff=0)
{/*For one species VM only!!! MP ver
   CharLength is obtained from BoxInfo, AligningLength must be less than CharLength
   RepulsiveForce: radius -> force, positive repulsive
   GuidingForceCoeff: the intensity of guiding force in \pm y direction
 */
  double SystemSize_X=BoxInfo[BoxInfo.size()-1].first[1];
  double SystemSize_Y=BoxInfo[BoxInfo.size()-1].first[3];
  double CharLength=BoxInfo[0].first[1]-BoxInfo[0].first[0];

  auto EffRelativePos=[SystemSize_X,SystemSize_Y](complex<double> RelativePos)->complex<double>{
    double EffRelX=RelativePos.real(),EffRelY=RelativePos.imag();
    if(not(abs(EffRelX)<=abs(EffRelX+SystemSize_X) and abs(EffRelX)<=abs(EffRelX-SystemSize_X))){
      if(abs(EffRelX+SystemSize_X)<=abs(EffRelX-SystemSize_X)) EffRelX=EffRelX+SystemSize_X;
      else EffRelX=EffRelX-SystemSize_X;
    }
    if(not(abs(EffRelY)<=abs(EffRelY+SystemSize_Y) and abs(EffRelY)<=abs(EffRelY-SystemSize_Y))){
      if(abs(EffRelY+SystemSize_Y)<=abs(EffRelY-SystemSize_Y)) EffRelY=EffRelY+SystemSize_Y;
      else EffRelY=EffRelY-SystemSize_Y;
    }
    return EffRelX+ii*EffRelY;
  };// do deal with the relative position on the periphery

  Meshed_PolarParticle2D_GT AnsList(BoxInfo.size());

  #pragma omp parallel for
  for( int BoxIndex=0 ; BoxIndex<BoxInfo.size() ; BoxIndex++ ){
    vector<int> AdjacentBoxList=BoxInfo[BoxIndex].second;
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[0]].second[2]);
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[0]].second[3]);
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[1]].second[2]);
    AdjacentBoxList.push_back(BoxInfo[AdjacentBoxList[1]].second[3]);
    //A list of index of adjacent boxes
    //[Warning] Assuming the box matrix has more than 3 rows and 3 collumns.

    for(int ip=0 ; ip<list[BoxIndex].size(); ip++){
      auto p=list[BoxIndex][ip];

      complex<double> SumDirComplex(0,0);
      double Distance;
      for(int AdjBoxIndex : AdjacentBoxList)
	for( auto p2 : list[AdjBoxIndex]){
	  Distance=DistanceOfTwoParticles(p,p2);
	  
	  if(Distance<=AligningLength)
	    SumDirComplex+=exp(ii*p2.Dir); //Alignment interaction

	  if(Distance<=CharLength)
	    SumDirComplex+=RepulsiveForce(Distance)*exp(ii*arg(EffRelativePos(p.Pos-p2.Pos))); //Repulsive interaction
	}

      PolarParticle2D p2;
      for( int i_sameBox=0 ; i_sameBox<list[BoxIndex].size() ; i_sameBox++ ){
	p2=list[BoxIndex][i_sameBox];
	Distance=DistanceOfTwoParticles(p,p2);
	
	if(Distance<=AligningLength)//the particle itself is also taken into account
	  SumDirComplex+=exp(ii*p2.Dir); //Alignment Interaction

	if(Distance<=CharLength and i_sameBox!=ip)
	  SumDirComplex+=RepulsiveForce(Distance)*exp(ii*arg(EffRelativePos(p.Pos-p2.Pos))); //Repulsive interaction
      }

      PolarParticle2D dp;
      dp.Dir=arg(SumDirComplex)+Sigma*(RandGen[omp_get_thread_num()].RandomDouble()*2*pi-pi);
      //dp.Dir=arg(SumDirComplex)+Sigma*(RandGen[0].RandomDouble()*2*pi-pi);
      dp.Dir+=GuidingForceCoeff*sin(2*dp.Dir);
      dp.Pos=p.Pos+exp(ii*dp.Dir)*Velocity;
      
      auto dp_Processed=BoundaryCondFunc(dp);
      AnsList[BoxIndex].push_back(dp_Processed);
    }
  }

  Meshed_PolarParticle2D_GT FnlAnsList(BoxInfo.size());
  for( auto Box : AnsList)
    for( auto p : Box)
      FnlAnsList[BoxMeshing_2D::DecideBoxIndex(p,SystemSize_X,SystemSize_Y,CharLength)].push_back(p);

  return FnlAnsList;
}

#endif //LocalProgram


#ifdef CheckingHeader
HeaderChecker(CollectiveMotionRealization_hpp);
#endif //CheckingHeader

#endif //Collectivemotionrealization_hpp
