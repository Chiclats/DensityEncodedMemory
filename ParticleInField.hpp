/*--------------------
  ver. 250630
  --------------------*/

#ifndef ParticleInField_hpp
#define ParticleInField_hpp

#ifndef vs
#include<bits/stdc++.h>
#endif

#include"BasicTools.hpp"
#include"StateProcessing.hpp"

using namespace std;

struct ParticleInField{

  complex<double> Pos;
  complex<double> Vel;
  double Dir;

  int Tag;//Tag=0 particles, Tag=1 field points
  static const int Field=1;
  static const int Particle=0;

  ParticleInField(complex<double> P_,complex<double> V_,double D_,int T_=0){
    Pos=P_;
    Vel=V_;
    Dir=D_;
    Tag=T_;
  }
  ParticleInField(complex<double> P_,complex<double> V_) : ParticleInField(P_,V_,0,1) {}
  ParticleInField(complex<double> P_,double D_) : ParticleInField(P_,0,D_,0) {}
  ParticleInField() : ParticleInField(0,0,0,0) {}

  bool operator==(const ParticleInField& PF) const{
    return Pos==PF.Pos and Vel==PF.Vel and Dir==PF.Dir and Tag==PF.Tag;
  }
  bool operator!=(const ParticleInField& PF) const{
    return not(*this==PF);
  }

  bool IsField() const{
    return Tag==1;
  }
  bool IsParticle() const{
    return Tag==0;
  }
  static bool FieldP(ParticleInField p){
    return p.Tag==1;
  }
  static bool ParticleP(ParticleInField p){
    return p.Tag==0;
  }
};

using ParticleInField_GT=vector<ParticleInField>;
using ParticleInField_MGT=vector<vector<ParticleInField>>;

namespace DirectlyOutput{

  void a_ParticleInField(const ParticleInField& PF,bool ShowAll=false){
    if(ShowAll) cout<<PF.Pos<<"\t"<<PF.Vel<<"\t"<<PF.Dir<<"\t"<<PF.Tag<<endl;
    else if(PF.IsField()) cout<<PF.Pos<<"\t"<<PF.Vel<<"\tF\n";
    else cout<<PF.Pos<<"\t"<<PF.Vel<<"\t"<<PF.Dir<<"\tP\n";
  }

  void Vec(const vector<ParticleInField>& PFG,bool ShowAll=false){
    for(auto PF : PFG) a_ParticleInField(PF,ShowAll);
  }
  
};

namespace FileIO{

  template<>
  ParticleInField_GT InVec<ParticleInField>(const string& Filename){
    
    ifstream fin;
    fin.open(Filename);
    
    vector<ParticleInField> AnsVec;
    ParticleInField ptc;
    double datum;
    int index=0;
    
    while(fin>>datum){
      
      switch(index){
        case 0: ptc.Pos.real(datum); break;
        case 1: ptc.Pos.imag(datum); break;
        case 2: ptc.Vel.real(datum); break;
        case 3: ptc.Vel.imag(datum); break;
        case 4: ptc.Dir=datum; break;
        case 5: ptc.Tag=datum; break;
      }
	
      index++;
      if(index>5){
	AnsVec.push_back(ptc);
	index=0;
      }
    }
    
    fin.close();
    return AnsVec;
  }

  void OutParticleVec(const string& Filename,const ParticleInField_GT& PFG){
    
    ofstream fout;
    fout.open(Filename);
    
    for( int i=0 ; i<PFG.size() ; i++)
      fout<<setprecision(DoublePrecession)
	  <<real(PFG[i].Pos)<<' '<<imag(PFG[i].Pos)<<' '
	  <<real(PFG[i].Vel)<<' '<<imag(PFG[i].Vel)<<' '
	  <<PFG[i].Dir<<' '<<PFG[i].Tag<<endl;
    
    fout.close();
    return;
  }

  void OutParticleVec(const string& Filename,const ParticleInField_MGT& MPFG){
    
    ofstream fout;
    fout.open(Filename);
    
    for( int j=0 ; j<MPFG.size() ; j++ )
      for( int i=0 ; i<MPFG[j].size() ; i++)
	fout<<setprecision(DoublePrecession)
	    <<real(MPFG[j][i].Pos)<<' '<<imag(MPFG[j][i].Pos)<<' '
	    <<real(MPFG[j][i].Vel)<<' '<<imag(MPFG[j][i].Vel)<<' '
	    <<MPFG[j][i].Dir<<' '<<MPFG[j][i].Tag<<endl;
    
    fout.close();
    return;
  }
  
};

#ifdef CheckingHeader
HeaderChecker(ParticleInField_hpp);
#endif //CheckingHeader

#endif//ParticleInField_hpp
