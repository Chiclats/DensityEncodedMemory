/*--------------------
ver 241229
--------------------*/

#ifndef AgentRealization_hpp
#define AgentRealization_hpp

#include<bits/stdc++.h>
#include"BasicTools.hpp"

using namespace std;

//----------PolarParticle2D---------------

class PolarParticle2D
{
public:
  
  complex<double> Pos;
  double Dir;
  int Tag=1; //Tag default value 1
  vector<double> DoubleIdentifier; //For more complicated identification

  //constructors
  PolarParticle2D() {}
  PolarParticle2D(complex<double> Po_,double Di_,int Ta_,vector<double> DI_){
    Pos=Po_;
    Dir=Di_;
    Tag=Ta_;
    DoubleIdentifier=DI_;

    //Dir range [0,2pi)
    while(Dir>=2*pi) Dir-=2*pi;
    while(Dir<0) Dir+=2*pi;
  }
  PolarParticle2D(complex<double> Po_,double Di_,int Ta_) : PolarParticle2D(Po_,Di_,Ta_,{}) {}
  PolarParticle2D(complex<double> Po_,double Di_) : PolarParticle2D(Po_,Di_,1) {}
  PolarParticle2D(complex<double> Po_,int i_MP=0) : PolarParticle2D(Po_,RandGen[i_MP].RandomDouble()*2*pi) {}

  void SetPos(double x,double y)
  {
    Pos.real(x);
    Pos.imag(y);
    return;
  }//checked 240501

  void SetDir(double v)
  {
    Dir=v;
    return;
  }//checked 240501

  void SetTag(int tag)
  {
    Tag=tag;
    return ;
  }//checked 240501

  complex<double> DirComplex()
  {//return a unit vector (complex) of the velocity direction
    complex<double> ans(cos(Dir),sin(Dir));
    return ans;
  }//checked 240501

  void DirectlyOutput(bool ShowTag=true,bool ShowDoubleIdentifier=false)
  {
    cout<<Pos<<' '<<Dir;
    if(ShowTag)cout<<' '<<Tag;
    cout<<endl;
    if(ShowDoubleIdentifier)
      DirectlyOutput::Vec(DoubleIdentifier);
    return;
  }//checked 240501

};//checked 240501

using PolarParticle2D_GT=vector<PolarParticle2D>;
using Meshed_PolarParticle2D_GT=vector<vector<PolarParticle2D>>;

namespace DirectlyOutput{
  void a_PolarParticle2D(PolarParticle2D p,bool ShowTag=true,bool ShowDoubleIdentifier=false)
  {
    cout<<p.Pos<<' '<<p.Dir;
    if(ShowTag)cout<<' '<<p.Tag;
    cout<<endl;
    if(ShowDoubleIdentifier)
      DirectlyOutput::Vec(p.DoubleIdentifier);
    return;
  }//checked 240501

  void Vec(PolarParticle2D_GT PG,bool ShowTag=true,bool ShowDoubleIdentifier=false){
    for( auto p : PG )
      a_PolarParticle2D(p,ShowTag,ShowDoubleIdentifier);
    return;
  }

  void Vec(Meshed_PolarParticle2D_GT MeshedPG,bool ShowTag=true,bool ShowDoubleIdentifier=false){
    for( auto PG : MeshedPG)
      for( auto p : PG )
	a_PolarParticle2D(p,ShowTag,ShowDoubleIdentifier);
    return;
  }

  //----------not recommended----------

  void Vec_PolarParticle2D(PolarParticle2D_GT PG,bool ShowTag=true,bool ShowDoubleIdentifier=false)
  {//use Vec instead
    for( auto p : PG )
      a_PolarParticle2D(p,ShowTag,ShowDoubleIdentifier);
    return;
  }//checked 240501
};//checked 240501

namespace FileIO{
  
  vector<PolarParticle2D> InParticleVec(string Filename,bool ReadTag=true)
  {
    ifstream fin;
    fin.open(Filename);
    
    vector<PolarParticle2D> AnsVec;
    PolarParticle2D ptc;
    double datum;
    int index=0;
    
    while(fin>>datum){
      if(ReadTag){
	if(index==0)
	  ptc.Pos.real(datum);
	else if(index==1)
	  ptc.Pos.imag(datum);
	else if(index==2)
	  ptc.Dir=datum;
	else ptc.Tag=datum;
	
	index++;
	if(index>3){
	  AnsVec.push_back(ptc);
	  index=0;
	}
      }
      else{
	if(index==0)
	  ptc.Pos.real(datum);
	else if(index==1)
	  ptc.Pos.imag(datum);
	else if(index==2)
	  ptc.Dir=datum;
	
	index++;
	if(index>2){
	  AnsVec.push_back(ptc);
	  index=0;
	}
      }
    }
    
    fin.close();
    return AnsVec;
  }//checked 240501

  void OutParticleVec(string Filename,PolarParticle2D_GT PG,bool WriteTag=true)
  {
    ofstream fout;
    fout.open(Filename);
    for( int i=0 ; i<PG.size() ; i++){
      fout<<setprecision(DoublePrecession)<<real(PG[i].Pos)<<' '<<imag(PG[i].Pos)<<' '<<PG[i].Dir;
      if(WriteTag) fout<<' '<<PG[i].Tag;
      fout<<endl;
    }
    
    fout.close();
    return;
  }

  void OutParticleVec(string Filename,Meshed_PolarParticle2D_GT MeshedPG,bool WriteTag=true)
  {
    ofstream fout;
    fout.open(Filename);
    for( int j=0 ; j<MeshedPG.size() ; j++ )
      for( int i=0 ; i<MeshedPG[j].size() ; i++){
	fout<<setprecision(DoublePrecession)<<real(MeshedPG[j][i].Pos)<<' '<<imag(MeshedPG[j][i].Pos)<<' '<<MeshedPG[j][i].Dir;
	if(WriteTag) fout<<' '<<MeshedPG[j][i].Tag;
	fout<<endl;
      }
    
    fout.close();
    return;
  }
  
};

//----------not recommended----------

PolarParticle2D Make_PolarParticle2D(double x,double y,double dir,int Tag=1)//use the constructor PolarParticle2D
{// generate an instance of the type PolarParticle2D
    PolarParticle2D ans;
    ans.SetDir(dir);
    ans.SetPos(x,y);
    ans.SetTag(Tag);
    return ans;
}//checked 240501


#ifdef CheckingHeader
HeaderChecker(AgentRealization_hpp);
#endif //CheckingHeader

#endif //Agentrealization_hpp
