/*--------------------
ver 250630
--------------------*/

#ifndef InteractionRealization_hpp
#define InteractionRealization_hpp

#ifndef vs
#include<bits/stdc++.h>
#endif

#include"BasicTools.hpp"
#include"AgentRealization.hpp"
#include"CollectiveMotionRealization.hpp"

using namespace std;

namespace Interaction{

  //----------Coupling---------------
  //WARNING: MP RandGen adaptation haven't been done

  CouplingT SineAngularCoupling_1Species(double Gamma)
  {//\Gamma\sin(\phi_j-\phi_i), for only one species
    return [Gamma](PolarParticle2D p1,PolarParticle2D p2)->PolarParticle2D{
	return Make_PolarParticle2D(0,0,Gamma*sin(p2.Dir-p1.Dir));
    };    
  }

  CouplingT SineAngularCoupling(double GammaA,double GammaB, double GammaAB)
  {//\Gamma_{ij}\sin(\phi_j-\phi_i)
    return [GammaA,GammaB,GammaAB](PolarParticle2D p1,PolarParticle2D p2)->PolarParticle2D{
      if(p1.Tag==1 and p2.Tag==1)
	return Make_PolarParticle2D(0,0,GammaA*sin(p2.Dir-p1.Dir));
      else if(p1.Tag==2 and p2.Tag==2)
	return Make_PolarParticle2D(0,0,GammaB*sin(p2.Dir-p1.Dir));
      else
	return Make_PolarParticle2D(0,0,GammaAB*sin(p2.Dir-p1.Dir));
    };
  }//checked 240502

  CouplingT Sine_LinearPert_AngularCoupling(double GammaA,double GammaB,double GammaAB,double c)
  {//\Gamma_{ij}(\sin(\phi_j-\phi_i)+c(1-(\phi_j-\phi_i)/\pi))
    return [GammaA,GammaB,GammaAB,c](PolarParticle2D p1,PolarParticle2D p2)->PolarParticle2D{
      double delta_Dir_reduced=p2.Dir-p1.Dir;
      while(delta_Dir_reduced>=2*pi) delta_Dir_reduced-=2*pi;
      while(delta_Dir_reduced<=0) delta_Dir_reduced+=2*pi;
    
      if(p1.Tag==1 and p2.Tag==1)
	return Make_PolarParticle2D(0,0,GammaA*(sin(delta_Dir_reduced)+c*(1-delta_Dir_reduced/pi)));
      else if(p1.Tag==2 and p2.Tag==2)
	return Make_PolarParticle2D(0,0,GammaB*(sin(delta_Dir_reduced)+c*(1-delta_Dir_reduced/pi)));
      else
	return Make_PolarParticle2D(0,0,GammaAB*(sin(delta_Dir_reduced)+c*(1-delta_Dir_reduced/pi)));
    };
  }//checked 240502

  CouplingT LinearAngularCoupling(double GammaA,double GammaB,double GammaAB)
  {//\Gamma_{ij}(1-(\phi_j-\phi_i)/\pi)
    return [GammaA,GammaB,GammaAB](PolarParticle2D p1,PolarParticle2D p2)->PolarParticle2D{
      double delta_Dir_reduced=p2.Dir-p1.Dir;
      while(delta_Dir_reduced>=2*pi) delta_Dir_reduced-=2*pi;
      while(delta_Dir_reduced<=0) delta_Dir_reduced+=2*pi;
    
      if(p1.Tag==1 and p2.Tag==1)
	return Make_PolarParticle2D(0,0,GammaA*(1-delta_Dir_reduced/pi));
      else if(p1.Tag==2 and p2.Tag==2)
	return Make_PolarParticle2D(0,0,GammaB*(1-delta_Dir_reduced/pi));
      else
	return Make_PolarParticle2D(0,0,GammaAB*(1-delta_Dir_reduced/pi));
    };
  }//checked 240502

  CouplingT Sin_CosPert_AngularCoupling(double GammaA,double GammaB,double GammaAB,double c)
  {//\Gamma_{ij}(\sin(\phi_j-\phi_i)+c\cos((\phi_j-\phi_i)/2))
    return [GammaA,GammaB,GammaAB,c](PolarParticle2D p1,PolarParticle2D p2)->PolarParticle2D{
      double delta_Dir_reduced=p2.Dir-p1.Dir;
      while(delta_Dir_reduced>=2*pi) delta_Dir_reduced-=2*pi;
      while(delta_Dir_reduced<=0) delta_Dir_reduced+=2*pi;
    
      if(p1.Tag==1 and p2.Tag==1)
	return Make_PolarParticle2D(0,0,GammaA*(sin(delta_Dir_reduced)+c*cos(delta_Dir_reduced/2)));
      else if(p1.Tag==2 and p2.Tag==2)
	return Make_PolarParticle2D(0,0,GammaB*(sin(delta_Dir_reduced)+c*cos(delta_Dir_reduced/2)));
      else
	return Make_PolarParticle2D(0,0,GammaAB*(sin(delta_Dir_reduced)+c*cos(delta_Dir_reduced/2)));
    };
  }//checked 240502

  CouplingT NFoldSineAngularCoupling(double GammaA,double GammaB, double GammaAB,int fold)
  {//\Gamma_{ij}(\sin(\phi_j-\phi_i)-\sin(n(\phi_j-\phi_i))/n)
    return [GammaA,GammaB,GammaAB,fold](PolarParticle2D p1,PolarParticle2D p2)->PolarParticle2D{
      if(p1.Tag==1 and p2.Tag==1)
	return Make_PolarParticle2D(0,0,GammaA*(sin(p2.Dir-p1.Dir)-sin(fold*(p2.Dir-p1.Dir))/fold));
      else if(p1.Tag==2 and p2.Tag==2)
	return Make_PolarParticle2D(0,0,GammaB*(sin(p2.Dir-p1.Dir)-sin(fold*(p2.Dir-p1.Dir))/fold));
      else
	return Make_PolarParticle2D(0,0,GammaAB*(sin(p2.Dir-p1.Dir)-sin(fold*(p2.Dir-p1.Dir))/fold));
    };
  }//checked 240502

  CouplingT SineNAngularCoupling(double GammaA,double GammaB,double GammaAB,int fold)
  {//\Gamma_{ij}\sin(n(\phi_j-\phi_i))
    return [GammaA,GammaB,GammaAB,fold](PolarParticle2D p1,PolarParticle2D p2)->PolarParticle2D{
      if(p1.Tag==1 and p2.Tag==1)
	return Make_PolarParticle2D(0,0,GammaA*sin(fold*(p2.Dir-p1.Dir)));
      else if(p1.Tag==2 and p2.Tag==2)
	return Make_PolarParticle2D(0,0,GammaB*sin(fold*(p2.Dir-p1.Dir)));
      else
	return Make_PolarParticle2D(0,0,GammaAB*sin(fold*(p2.Dir-p1.Dir)));
    };
  }

  //----------Self evolution---------------
  //WARNING: MP RandGen adaptation haven't been done

  SelfEvolutionT SelfDrivenEvolution(double v,double sigma,int i_MP=0)
  {//\dot x_i=v\cos\phi_i, \dot y_i=v\sin\phi_i, \phi_i+=\sigma\xi_i
    return [v,sigma,i_MP](PolarParticle2D p,double TimeStep)->PolarParticle2D{
      return Make_PolarParticle2D(cos(p.Dir)*v*TimeStep,sin(p.Dir)*v*TimeStep,sigma*sqrt(TimeStep)*RandGen[i_MP].GaussianRandomDouble());
    };
  }//checked 240502

  SelfEvolutionT NoMotionEvolution(double sigma,int i_MP=0)
  {//particles never moves, \dot\phi_i+=\sigma\xi_i
    return [sigma,i_MP](PolarParticle2D p,double TimeStep)->PolarParticle2D{
      return Make_PolarParticle2D(0,0,sigma*sqrt(TimeStep)*RandGen[i_MP].GaussianRandomDouble());
    };
  }//checked 240502

  SelfEvolutionT DirectedSelfDrivenEvolution(double v,double sigma,double epsilon,int i_MP=0)
  {//\dot x_i=v\cos\phi_i, \dot y_i=v\sin\phi_i, \phi_i+=\sigma\xi_i\dd t+\epsilon\sin 2\phi_i\dd t
    return [v,sigma,epsilon,i_MP](PolarParticle2D p,double TimeStep)->PolarParticle2D{
      return Make_PolarParticle2D(cos(p.Dir)*v*TimeStep,
				  sin(p.Dir)*v*TimeStep,
				  sigma*sqrt(TimeStep)*RandGen[i_MP].GaussianRandomDouble()+epsilon*sin(2*p.Dir)*TimeStep);
    };
  }

  //----------BoundaryCondFunc---------------
  // Attention here for BoundaryCondFunc, the "pull back" of Dir need also be considered!

  BoundaryT Boundary_Periodic(double SystemSize_X,double SystemSize_Y)
  {//4 periodic bdr
    return [SystemSize_X,SystemSize_Y](PolarParticle2D p)->PolarParticle2D{
      PolarParticle2D AnsP=p;

      if(AnsP.Pos.imag()>SystemSize_Y) AnsP.Pos-=SystemSize_Y*ii;
      if(AnsP.Pos.imag()<0) AnsP.Pos+=SystemSize_Y*ii;

      if(AnsP.Pos.real()<0) AnsP.Pos+=SystemSize_X;
      if(AnsP.Pos.real()>SystemSize_X) AnsP.Pos-=SystemSize_X;

      if(AnsP.Dir<0) AnsP.Dir+=2*pi;
      if(AnsP.Dir>2*pi) AnsP.Dir-=2*pi;
      
      return AnsP;
    };
  }

  BoundaryT Boundary_PRPR(double SystemSize_X,double SystemSize_Y)
  {//square system with 2 y=c boundaries periodic and 2 x=c boundaries reflective
    return [SystemSize_X,SystemSize_Y](PolarParticle2D p)->PolarParticle2D{
      PolarParticle2D AnsP=p;

      if(AnsP.Pos.imag()>SystemSize_Y) AnsP.Pos-=SystemSize_Y*ii;
      if(AnsP.Pos.imag()<0) AnsP.Pos+=SystemSize_Y*ii;

      if(AnsP.Pos.real()<0){
	AnsP.Pos=-AnsP.Pos.real()+AnsP.Pos.imag()*ii;
	AnsP.Dir=pi-AnsP.Dir;
	if(AnsP.Dir<0) AnsP.Dir+=2*pi;
      }
      if(AnsP.Pos.real()>SystemSize_X){
	AnsP.Pos=2*SystemSize_X-AnsP.Pos.real()+AnsP.Pos.imag()*ii;
	AnsP.Dir=pi-AnsP.Dir;
	if(AnsP.Dir<0) AnsP.Dir+=2*pi;
      }

      if(AnsP.Dir<0) AnsP.Dir+=2*pi;
      if(AnsP.Dir>2*pi) AnsP.Dir-=2*pi;
      
      return AnsP;
    };
  }

  BoundaryT Boundary_RPRP(double SystemSize_X,double SystemSize_Y)
  {//square system with 2 x=c boundaries periodic and 2 y=c boundaries reflective
    return [SystemSize_X,SystemSize_Y](PolarParticle2D p)->PolarParticle2D{
      PolarParticle2D AnsP=p;

      if(AnsP.Pos.real()>SystemSize_X) AnsP.Pos-=SystemSize_X;
      if(AnsP.Pos.real()<0) AnsP.Pos+=SystemSize_X;

      if(AnsP.Pos.imag()<0){
	AnsP.Pos=AnsP.Pos.real()-AnsP.Pos.imag()*ii;
	AnsP.Dir=-AnsP.Dir;
      }
      if(AnsP.Pos.imag()>SystemSize_Y){
	AnsP.Pos=AnsP.Pos.real()+(2*SystemSize_Y-AnsP.Pos.imag())*ii;
	AnsP.Dir=-AnsP.Dir;
      }

      if(AnsP.Dir<0) AnsP.Dir+=2*pi;
      if(AnsP.Dir>2*pi) AnsP.Dir-=2*pi;
      
      return AnsP;
    };
  }

  //----------Distance in system---------------

  DistanceT Distance_Periodic(double SystemSize_X,double SystemSize_Y)
  {//4 periodic bdr
    return [SystemSize_X,SystemSize_Y](PolarParticle2D p1,PolarParticle2D p2)->double{
      double
	d1=sqrt(norm(p2.Pos-p1.Pos)),
	d2=sqrt(norm(p2.Pos-p1.Pos+SystemSize_Y*ii)),
	d3=sqrt(norm(p2.Pos-p1.Pos-SystemSize_Y*ii)),
	d4=sqrt(norm(p2.Pos-p1.Pos+SystemSize_X)),
	d5=sqrt(norm(p2.Pos-p1.Pos-SystemSize_X)),
	d6=sqrt(norm(p2.Pos-p1.Pos+SystemSize_X+SystemSize_Y*ii)),
	d7=sqrt(norm(p2.Pos-p1.Pos+SystemSize_X-SystemSize_Y*ii)),
	d8=sqrt(norm(p2.Pos-p1.Pos-SystemSize_X+SystemSize_Y*ii)),
	d9=sqrt(norm(p2.Pos-p1.Pos-SystemSize_X-SystemSize_Y*ii));
      return MinVector(vector<double>{d1,d2,d3,d4,d5,d6,d7,d8,d9});
    };
  }

  DistanceT Distance_PRPR(double SystemSize_X,double SystemSize_Y)
  {//square system with 2 y=c boundaries periodic and 2 x=c boundaries reflective
    return [SystemSize_X,SystemSize_Y](PolarParticle2D p1,PolarParticle2D p2)->double{
      double
	d1=sqrt(norm(p2.Pos-p1.Pos)),
	d2=sqrt(norm(p2.Pos-p1.Pos+SystemSize_Y*ii)),
	d3=sqrt(norm(p2.Pos-p1.Pos-SystemSize_Y*ii));
      return min(min(d1,d2),d3);
    };
  }

  DistanceT Distance_RPRP(double SystemSize_X,double SystemSize_Y)
  {//square system with 2 x=c boundaries periodic and 2 y=c boundaries reflective
    return [SystemSize_X,SystemSize_Y](PolarParticle2D p1,PolarParticle2D p2)->double{
      double
	d1=sqrt(norm(p2.Pos-p1.Pos)),
	d2=sqrt(norm(p2.Pos-p1.Pos+SystemSize_X)),
	d3=sqrt(norm(p2.Pos-p1.Pos-SystemSize_X));
      return min(min(d1,d2),d3);
    };
  }

  //----------RelativePos---------------

  complex<double> RelativePos_Periodic(PolarParticle2D p1,PolarParticle2D p2,double SystemSize_X,double SystemSize_Y){
    vector<complex<double>> PosList={
      p2.Pos-p1.Pos,
      p2.Pos-p1.Pos+SystemSize_X,
      p2.Pos-p1.Pos-SystemSize_X,
      p2.Pos-p1.Pos+SystemSize_Y*ii,
      p2.Pos-p1.Pos-SystemSize_Y*ii,
      p2.Pos-p1.Pos+SystemSize_X+SystemSize_Y*ii,
      p2.Pos-p1.Pos+SystemSize_X-SystemSize_Y*ii,
      p2.Pos-p1.Pos-SystemSize_X+SystemSize_Y*ii,
      p2.Pos-p1.Pos-SystemSize_X-SystemSize_Y*ii
    };
    sort(PosList.begin(),PosList.end(),
	 [](complex<double> c1,complex<double> c2)->bool{
	   return abs(c1)<abs(c2);
	 });
    return PosList[0];
  }

  complex<double> RelativePos_RPRP(PolarParticle2D p1,PolarParticle2D p2,double SystemSize_X,double SystemSize_Y){
    vector<complex<double>> PosList={
      p2.Pos-p1.Pos,
      p2.Pos-p1.Pos+SystemSize_X,
      p2.Pos-p1.Pos-SystemSize_X
    };
    sort(PosList.begin(),PosList.end(),
	 [](complex<double> c1,complex<double> c2)->bool{
	   return abs(c1)<abs(c2);
	 });
    return PosList[0];
  }
  

//----------Not recommended---------------

  function<PolarParticle2D(PolarParticle2D,PolarParticle2D)> DoubleSineAngularCoupling(double GammaA,double GammaB, double GammaAB)
  {//\dot \phi_i+=\sum_{j\in\Omega_i} \Gamma_{ij}\sin(\phi_j-\phi_i)
    return [GammaA,GammaB,GammaAB](PolarParticle2D p1,PolarParticle2D p2)->PolarParticle2D{
      if(p1.Tag==1 and p2.Tag==1)
	return Make_PolarParticle2D(0,0,GammaA*(sin(p2.Dir-p1.Dir)-sin(2*(p2.Dir-p1.Dir))/2));
      else if(p1.Tag==2 and p2.Tag==2)
	return Make_PolarParticle2D(0,0,GammaB*(sin(p2.Dir-p1.Dir)-sin(2*(p2.Dir-p1.Dir))/2));
      else
	return Make_PolarParticle2D(0,0,GammaAB*(sin(p2.Dir-p1.Dir)-sin(2*(p2.Dir-p1.Dir))/2));
    };
  }  //Not recommended

  function<PolarParticle2D(PolarParticle2D,PolarParticle2D)> Sin_MixPert_AngularCoupling(double GammaA,double GammaB,double GammaAB,double c)
  {//\dot \phi_i+=\sum_{j\in\Omega_i} \Gamma_{ij}(\sin(\phi_j-\phi_i)+c\cos((\phi_j-\phi_i)/2))
    return [GammaA,GammaB,GammaAB,c](PolarParticle2D p1,PolarParticle2D p2)->PolarParticle2D{
      double delta_Dir_reduced=p2.Dir-p1.Dir;
      while(delta_Dir_reduced>=2*pi) delta_Dir_reduced-=2*pi;
      while(delta_Dir_reduced<=0) delta_Dir_reduced+=2*pi;
    
      if(p1.Tag==1 and p2.Tag==1)
	return Make_PolarParticle2D(0,0,GammaA*(sin(p2.Dir-p1.Dir)+c/2*cos(delta_Dir_reduced/2)+c/2*(1-delta_Dir_reduced/pi)));
      else if(p1.Tag==2 and p2.Tag==2)
	return Make_PolarParticle2D(0,0,GammaB*(sin(p2.Dir-p1.Dir)+c/2*cos(delta_Dir_reduced/2)+c/2*(1-delta_Dir_reduced/pi)));
      else
	return Make_PolarParticle2D(0,0,GammaAB*(sin(p2.Dir-p1.Dir)+c/2*cos(delta_Dir_reduced/2)+c/2*(1-delta_Dir_reduced/pi)));
    };
  }  //Not recommended
  
}

#ifdef CheckingHeader
HeaderChecker(InteractionRealization_hpp);
#endif //CheckingHeader

#endif //Interactionrealization_hpp
