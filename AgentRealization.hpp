/*--------------------
ver 250630
--------------------*/

#ifndef AgentRealization_hpp
#define AgentRealization_hpp

#ifndef vs
#include<bits/stdc++.h>
#endif

#include"BasicTools.hpp"

#include <cstring>
#include <cerrno>

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
    if(not fin.is_open()){
      cerr<<"Filename: "<<Filename<<endl;
      throw runtime_error("In FileIO::InParticleVec. Failed to open the file.");
    }
    
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
  }
  
  void OutParticleVec(string Filename,PolarParticle2D_GT PG,bool WriteTag=true)
  {
    ofstream fout;
    fout.open(Filename);
    if(not fout.is_open()){
      cerr<<"Filename: "<<Filename<<endl;
      cerr<<"In FileIO::OutParticleVec. Failed to open the file."<<endl;
      return;
    }

    fout<<fixed<<setprecision(DoublePrecession);
    
    for( int i=0 ; i<PG.size() ; i++){
      fout<<real(PG[i].Pos)<<' '<<imag(PG[i].Pos)<<' '<<PG[i].Dir;
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
    if(not fout.is_open()){
      cerr<<"Filename: "<<Filename<<endl;
      cerr<<"In FileIO::OutParticleVec. Failed to open the file."<<endl;
      return;
    }
    
    for( int j=0 ; j<MeshedPG.size() ; j++ )
      for( int i=0 ; i<MeshedPG[j].size() ; i++){
	fout<<setprecision(DoublePrecession)<<real(MeshedPG[j][i].Pos)<<' '<<imag(MeshedPG[j][i].Pos)<<' '<<MeshedPG[j][i].Dir;
	if(WriteTag) fout<<' '<<MeshedPG[j][i].Tag;
	fout<<endl;
      }
    
    fout.close();
    return;
  }

  void OutParticleVec_Bin(const string& Filename,const PolarParticle2D_GT& PG)
  {// output the particle group into a binary `.bpg' file. for each line {double,double,float}. no tag.
    #pragma pack(1)
    struct temp_P{
      double x,y;
      float Dir;
    };
    #pragma pack()

    vector<temp_P> temp_PG(PG.size());
    for( int i=0 ; i<temp_PG.size() ; i++ ){
      temp_PG[i].x=PG[i].Pos.real();
      temp_PG[i].y=PG[i].Pos.imag();
      temp_PG[i].Dir=PG[i].Dir;
    }

    std::ofstream outfile(Filename, std::ios::binary);
    if (!outfile.is_open()){
      cerr<<"Filename: "<<Filename<<endl;
      cerr<<"In FileIO::OutParticleVec_Bin. Failed to open the file."<<endl;
      return;
    }

    //first the line number of the file
    const size_t size = temp_PG.size();
    outfile.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
    //second all the data
    outfile.write(reinterpret_cast<const char*>(temp_PG.data()), size*sizeof(temp_P));

    if (!outfile.good()){
        std::cerr << "Error: Failed to write data to " << Filename << std::endl;
        throw std::runtime_error("In FileIO::OutParticleVec_Bin.");
    }

    outfile.close();
  }

  void OutParticleVec_Bin_Robust(const string& Filename,const PolarParticle2D_GT& PG)
  {// output the particle group into a binary `.bpg' file. for each line {double,double,float}. no tag.
    #pragma pack(1)
    struct temp_P{
      double x,y;
      float Dir;
    };
    #pragma pack()

    vector<temp_P> temp_PG(PG.size());
    for( int i=0 ; i<temp_PG.size() ; i++ ){
      temp_PG[i].x=PG[i].Pos.real();
      temp_PG[i].y=PG[i].Pos.imag();
      temp_PG[i].Dir=PG[i].Dir;
    }

    std::ofstream outfile(Filename, std::ios::binary);
    int Attempts=0;
    while(!outfile.is_open()){
      Attempts++;
      
      cerr<<"Filename: "<<Filename<<endl;
      cerr<<"In FileIO::OutParticleVec_Bin. Failed to open the file. Attempt="<<Attempts<<endl;

      //system error output
#ifndef vs
      std::cerr << "Error: " << strerror(errno) << std::endl;
#else
      char errMsg[256];
      strerror_s(errMsg, sizeof(errMsg), errno);
      std::cerr << "Error: " << errMsg << std::endl;
#endif // !vs

      if(Attempts==3){
	FileIO::OutParticleVec(Filename+".txt",PG);	
	return;
      }
	
      // play for a while before next attempt
      int i_Temp;
      for( int j_Temp=0 ; j_Temp<1e7 ; j_Temp++ ) i_Temp=j_Temp;
    }

    //first the line number of the file
    const size_t size = temp_PG.size();
    outfile.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
    //second all the data
    outfile.write(reinterpret_cast<const char*>(temp_PG.data()), size*sizeof(temp_P));

    if (!outfile.good()){
        std::cerr << "Error: Failed to write data to " << Filename << std::endl;
        throw std::runtime_error("In FileIO::OutParticleVec_Bin.");
    }

    outfile.close();
  }

  PolarParticle2D_GT InParticleVec_Bin(const string& Filename)
  {
    #pragma pack(1)
    struct temp_P{
      double x,y;
      float Dir;
    };
    #pragma pack()

    std::ifstream infile(Filename, std::ios::binary);
    if (!infile.is_open()) {
        std::cerr << "Error: Failed to open file " << Filename << std::endl;
        throw std::runtime_error("In FileIO::InParticleVec_Bin.");
    }

    size_t size;
    infile.read(reinterpret_cast<char*>(&size), sizeof(size_t));

    vector<temp_P> temp_PG(size);
    infile.read(reinterpret_cast<char*>(temp_PG.data()), size * sizeof(temp_P));

    PolarParticle2D_GT ansPG(size);
    for( int i=0 ; i<size ; i++ )
      ansPG[i]=PolarParticle2D(temp_PG[i].x+ii*temp_PG[i].y,(double) temp_PG[i].Dir);

    if (!infile.good()) {
      std::cerr << "Error: Failed to read data from " << Filename << std::endl;
      throw std::runtime_error("In FileIO::InParticleVec_Bin.");
    }
    infile.close();

    return ansPG;
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
