/*--------------------
ver 250222
--------------------*/

#ifndef BasicTools_hpp
#define BasicTools_hpp

#include<bits/stdc++.h>

#define ii complex<double>(0,1)//unit of imaginary numbers
#define pi 3.1415926535897932

#define HeaderChecker(HeaderName) \
  int _HeaderChecker##HeaderName (){ \
    cout<<#HeaderName<<" done\n";     \
    return 0; \
  } \
  int _Var_HeadChecker##HeaderName = _HeaderChecker##HeaderName();

using namespace std;

//---------------

//random number generator, each thread use one RandGen[]
struct RandomGenerator{
private:  
  random_device rd;
  mt19937 gen;

  uniform_real_distribution<> _dis_RandomDouble;
  normal_distribution<double> _dis_GaussianRandomDouble;

public:
  RandomGenerator(){
    gen=mt19937(rd());
    _dis_RandomDouble=uniform_real_distribution<>(0,1);
    _dis_GaussianRandomDouble=normal_distribution<double>(0,1);
  }

  double RandomDouble()
  {//generate a double datum, uniform dist in [0,1]
    return _dis_RandomDouble(gen);
  }

  double GaussianRandomDouble()
  {//generate a double datum, Gaussian dist N(0,1)
    return _dis_GaussianRandomDouble(gen);
  }

  int RandomInt(int a,int b)
  {//generate an int, uniformly dist in [a,b]
    uniform_int_distribution<int> dis(a, b);
    return dis(gen);
  }

  int RandomBinomialInt(int n,double p)
  {//generate an int, ~ Binomial(n,p)
    binomial_distribution<int> dis(n,p);
    return dis(gen);
  }
};
vector<RandomGenerator> RandGen(30);
//WARNING: max number of threads: 30

//---------------

complex<double> Unit(double rad){//angle to unit complex
  return cos(rad)+sin(rad)*ii;
}

complex<double> Pair2C(pair<double,double> r_i)
{//transform pair<double,double> to complex<double>
  return complex<double>(r_i.first,r_i.second);
}

pair<double,double> C2Pair(complex<double> c)
{//transform complex<double> to pair<double,double>
  return {c.real(),c.imag()};
}

//---------------

namespace DirectlyOutput
{// a collection of procedures to print data in the command-line

  template<typename Type>
  void Vec(vector<Type> vec){
    for(Type item : vec)
      cout<<item<<endl;
  }//checked 240501

  template<typename Type1, typename Type2>
  void Vec(vector<pair<Type1,Type2>> vec){
    for(auto item : vec)
      cout<<item.first<<"\t"<<item.second<<"\n";
  }

  template<typename Type>
  void Mat(vector<vector<Type>> mat){
    for(vector<Type> vec : mat){
      for(Type item : vec)
	cout<<item<<'\t';
      cout<<endl;
    }
  }//checked 240501
  
};//checked 240501

namespace FileIO
{// a collection of procedures to print data in the file

  short DoublePrecession=10;

  template<typename Type>
  vector<Type> InVec(const string& Filename)
  {//ATTENTION here the separator of Type=string is ' '
   //ATTENTION cannot use Type=int to input Type=double 
    ifstream fin;
    fin.open(Filename);
    if(not fin.is_open())
      throw runtime_error("In FileIO::InVec. Failed to open the file.");
    
    vector<Type> AnsVec;
    Type datum;
    while(fin>>datum)
      AnsVec.push_back(datum);
    
    fin.close();
    return AnsVec;
  }

  template<typename Type>
  vector<vector<Type>> InMat(const string& Filename,const int& cols)
  {
    vector<vector<Type>> AnsMat;
    vector<Type> TempVec(cols);
    int TempInt=0;
    
    ifstream fin;
    fin.open(Filename);
    if(not fin.is_open())
      throw runtime_error("In FileIO::InMat. Failed to open the file.");
    
    Type datum;
    while(fin>>datum){
      TempVec[TempInt]=datum;
      TempInt++;

      if(TempInt==cols){
	AnsMat.push_back(TempVec);
	TempInt=0;
      }
    }
    
    fin.close();
    
    return AnsMat;
  }

  template<typename Type>
  void OutVec(const string& Filename,const vector<Type>& List)
  {
    ofstream fout;
    fout.open(Filename);
    if(not fout.is_open()){
      //throw runtime_error("In FileIO::OutVec. Failed to open the file.");
      cerr<<"In FileIO::OutVec. Failed to open the file.\n";
      return;
    }
    
    for( int i=0 ; i<List.size() ; i++)
      fout<<setprecision(DoublePrecession)<<List[i]<<endl;
    fout.close();
    return;
  }

  template<typename Type>
  void OutVecPair(const string& Filename,const vector<pair<Type,Type>>& List)
  {
    ofstream fout;
    fout.open(Filename);
    if(not fout.is_open()){
      //throw runtime_error("In FileIO::OutVec. Failed to open the file.");
      cerr<<"In FileIO::OutVecPair. Failed to open the file.\n";
      return;
    }
    
    for( int i=0 ; i<List.size() ; i++)
      fout<<setprecision(DoublePrecession)<<List[i].first<<' '<<List[i].second<<endl;
    fout.close();
    return;
  }

  template<typename Type>
  void OutMatPair(const string& Filename,const vector<vector<pair<Type,Type>>>& List)
  {
    ofstream fout;
    fout.open(Filename);
    if(not fout.is_open()){
      //throw runtime_error("In FileIO::OutMatPair. Failed to open the file.");
      cerr<<"In FileIO::OutMatPair. Failed to open the file.\n";
      return;
    }
    
    for( auto v : List){
      for(auto p : v) fout<<setprecision(DoublePrecession)<<p.first<<' '<<p.second<<"\t";
      fout<<endl;
    }
    fout.close();
    return;
  }

  template<typename Type>
  void OutMat(const string& Filename,const vector<vector<Type>>& Mat)
  {
    ofstream fout;
    fout.open(Filename);
    if(not fout.is_open()){
      //throw runtime_error("In FileIO::OutMat. Failed to open the file.");
      cerr<<"In FileIO::OutMat. Failed to open the file.\n";
      return;
    }
    
    for( int i=0 ; i<Mat.size() ; i++){
      for(int j=0; j<Mat[0].size(); j++)
	fout<<setprecision(DoublePrecession)<<Mat[i][j]<<' ';
      fout<<endl;
    }
    fout.close();
    return;
  }

};

//---------------

template<typename Type>
vector<Type> SubVector(vector<Type> OriginVector,vector<int> IndexList)
{//input an IndexList, and return the SubVector
  vector<Type> AnsVec;
  for(auto i : IndexList)
    AnsVec.push_back(OriginVector[i]);
  return AnsVec;
}//checked 240501

template<typename Type>
Type MinVector(vector<Type> list,function<bool(Type,Type)> cf=[](Type i1,Type i2)->bool{return i1<i2;})
{//return the minimal of the list
  Type m=list[0];
  for( auto p : list)
    if(cf(p,m))
      m=p;
  return m;
} //Need further correction, can't be used when cf is set

double Mean(vector<double> list)
{//return the mean value of the list
  double sum=0;
  for( int i=0 ; i<list.size() ; i++ ) sum+=list[i];
  return sum/list.size();
}

double Var(vector<double> list)
{//return the variance value of the list
  double sum=0;
  for( int i=0 ; i<list.size() ; i++ ) sum+=list[i]*list[i];
  return -Mean(list)*Mean(list)+sum/list.size();
}

template<typename Type>
vector<Type> RemoveIf(vector<Type> OriginVector,const function<bool(Type)>& Predicate)
{//remove items in OriginVector satisfying Predicate(it)==true
  auto nit=remove_if(OriginVector.begin(),OriginVector.end(),Predicate);
  OriginVector.erase(nit,OriginVector.end());
  return OriginVector;
}

template<typename Type>
vector<Type> Gathering(vector<vector<Type>> VectorList)
{//gather all vectors into one vector

  vector<Type> AnsVec={};
  for( auto v : VectorList)
    AnsVec.insert(AnsVec.end(),v.begin(),v.end());

  return AnsVec;
}

//---------------

namespace Tictoc{

  chrono::time_point<chrono::high_resolution_clock> _StartTimePoint;

  void Tic(){
    _StartTimePoint=chrono::high_resolution_clock::now();
    return;
  }

  double Toc(string MarkerString="Time duration is %f ms.\n"){
    chrono::time_point<chrono::high_resolution_clock> EndTimePoint=chrono::high_resolution_clock::now();
    chrono::duration<double,milli> TimeDuration=EndTimePoint-_StartTimePoint;
    printf(MarkerString.c_str(),TimeDuration.count());
    return TimeDuration.count();
  }
  
};



#ifdef CheckingHeader
HeaderChecker(BasicTools_hpp);
#endif //CheckingHeader

#endif //BasicTools_hpp
