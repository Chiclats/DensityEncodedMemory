/*--------------------
  ver. 250630
  --------------------*/

#ifndef AligningRepulsiveInteraction_hpp
#define AligningRepulsiveInteraction_hpp

#ifndef vs
#include<bits/stdc++.h>
#endif

#include"BasicTools.hpp"
#include"StateProcessing.hpp"
#include"ExtendedCollectiveMotionRealization.hpp"
#include"ParticleInField.hpp"

using namespace std;

namespace AligningRepulsiveInteraction{

  
  complex<double> _TriangleInterp_C(complex<double> p,
				    pair<complex<double>,complex<double>> p1_v1,
				    pair<complex<double>,complex<double>> p2_v2,
				    pair<complex<double>,complex<double>> p3_v3)
  {//bilinear/triangle interpolation
    double x1=p1_v1.first.real(),x2=p2_v2.first.real(),x3=p3_v3.first.real();
    double y1=p1_v1.first.imag(),y2=p2_v2.first.imag(),y3=p3_v3.first.imag();
    complex<double> z1=p1_v1.second,z2=p2_v2.second,z3=p3_v3.second;
    complex<double> a=-((y2*z1-y3*z1-y1*z2+y3*z2+y1*z3-y2*z3)/(x2*y1-x3*y1-x1*y2+x3*y2+x1*y3-x2*y3));
    complex<double> b=-((-x2*z1+x3*z1+x1*z2-x3*z2-x1*z3+x2*z3)/(x2*y1-x3*y1-x1*y2+x3*y2+x1*y3-x2*y3));
    complex<double> c=-((-x3*y2*z1+x2*y3*z1+x3*y1*z2-x1*y3*z2-x2*y1*z3+x1*y2*z3)/(x2*y1-x3*y1-x1*y2+x3*y2+x1*y3-x2*y3));
    double x=p.real(),y=p.imag();
    return a*x+b*y+c;
  }

  //---------------Evolution
  
  EvolutionT<ParticleInField> Evolution_AligningRepulsive(function<complex<double>(ParticleInField,ParticleInField)> PPRepulsion,//latter on former
							  function<complex<double>(ParticleInField)> WPRepulsion,
							  function<double(ParticleInField,ParticleInField_GT)> PDirUpdate,//DirChangeRate(p,Neighbors)
							  double PActiveVel,double HPCoeff,DistanceExtendedT<ParticleInField> Distance,double PNoise,
							  double PHCoeff,double PHRange,double HRelaxationTime,double TimeStep)
  {//Evolve the particle or field point in a TimeStep
   //Dir=PDirUpdate(...), Vp=PActiveVel*e(Dir)+PPRepulsion(...)+WPRepulsion(...)+HPCoeff*Vh+(noise)
   //\dot Vh=PHCoeff*\sum Vp*exp(-(rp-rh)^2/PHRange^2)-Vh/HRelaxationTime
    return [=](ParticleInField p,ParticleInField_GT Neighbors)->ParticleInField{
      if(p.IsField()){
	auto ap=p;
	ap.Vel-=TimeStep*ap.Vel/HRelaxationTime;

	for(auto np : Neighbors)
	  if(np.IsParticle())
	    ap.Vel+=TimeStep*PHCoeff*np.Vel*exp(-norm(ap.Pos-np.Pos)/PHRange/PHRange);

	return ap;
      }
      else{
	auto ap=p;
	ap.Dir+=TimeStep*PDirUpdate(p,Neighbors);
	
	ap.Vel=0;
	ap.Vel+=PActiveVel*Unit(ap.Dir);
	ap.Vel+=WPRepulsion(p);

	for(auto np : Neighbors)
	  if(np.IsParticle())
	    ap.Vel+=PPRepulsion(p,np);

	auto Neighbors_Field=RemoveIf(Neighbors,(function<bool(ParticleInField)>)ParticleInField::ParticleP);
	sort(Neighbors_Field.begin(),Neighbors_Field.end(),
	     [p,Distance](ParticleInField p1,ParticleInField p2)->bool{
	       return Distance(p1,p)<Distance(p2,p);});
	auto f1=Neighbors_Field[0],f2=Neighbors_Field[1],f3=Neighbors_Field[2];
	//[WARNING] Here field points must be homogeneously distributed.
	
	complex<double> Vh=_TriangleInterp_C(p.Pos,{f1.Pos,f1.Vel},{f2.Pos,f2.Vel},{f3.Pos,f3.Vel});
	ap.Vel+=HPCoeff*Vh;
	
	ap.Pos+=ap.Vel*TimeStep+PNoise*sqrt(TimeStep)*RandGen[0].GaussianRandomDouble();
	return ap;
      }
    };
  }

  //---------------WPRepulsion
  
  function<Complex(ParticleInField)> WallRepulsion_SquareUD(double SystemSize_X,double SystemSize_Y,double RepulsionRange,double Stiffness)
  {//Linear repulsion
    return [SystemSize_Y,RepulsionRange,Stiffness](ParticleInField p)->Complex{
      double y=p.Pos.imag();
      if(y<=RepulsionRange) return Stiffness*(RepulsionRange-y)*ii;
      else if(SystemSize_Y-RepulsionRange<y) return -Stiffness*(y-(SystemSize_Y-RepulsionRange))*ii;
      else return 0;
    };
  }

  //---------------PPRepulsion

  function<Complex(ParticleInField,ParticleInField)> PPRepulsion_Exp(double PPRange,double Stiffness)
  {//F = Stiffness * exp( - r / PPRange )
    return [PPRange,Stiffness](ParticleInField p,ParticleInField np)->Complex{
      Complex Dr=p.Pos-np.Pos;
      double r=abs(Dr);
      Complex UDr=Dr/r;
      return UDr*Stiffness*exp(-r/PPRange);
    };
  }

  //---------------PDirUpdate

  function<double(ParticleInField,ParticleInField_GT)> PDirUpdate_SineAlignment(double Gamma)
  {//sine alignment
    return [Gamma](ParticleInField p,ParticleInField_GT Neighbors)->double{
      double ans=0;
      for(auto np : Neighbors)
	if(np.IsParticle())
	  ans+=sin(np.Dir-p.Dir);
      return ans;
    };
  }

};

#ifdef CheckingHeader
HeaderChecker(AligningRepulsiveInteraction_hpp);
#endif //CheckingHeader

#endif//AligningRepulsiveInteraction_hpp
