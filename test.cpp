#include"./ParticleSimulationHeaders.hpp"
#include"./ContinuousVicsekModel.hpp"

namespace VM=VicsekModel_2D;

int main(){

  const double
    SystemSize_X=10,
    SystemSize_Y=10,
    CharLength=1,
    Velocity=1,
    AligningCoeff=0,//1,
    Sigma=0,//0.1,
    GuidingCoeff=-1,
    TimeInterval=0.1,
    SqrtTimeInterval=sqrt(TimeInterval);

  const vector<double> Parameters={SystemSize_X,SystemSize_Y,CharLength,Velocity,AligningCoeff,Sigma,GuidingCoeff,TimeInterval,SqrtTimeInterval};

  const auto BoxInfo=BoxMeshing_2D::BoxInfo(SystemSize_X,SystemSize_Y,CharLength);
  const auto LDAdjBox=VM::LowerRightAdjacentBoxIndex(BoxInfo);
  const auto Distance=VM::Distance::Periodic(SystemSize_X,SystemSize_Y,CharLength);

  PolarParticle2D_GT PG={
    PolarParticle2D(5.5+5.5*ii,0.),
    PolarParticle2D(5.5+5.5*ii,0.),
    PolarParticle2D(5.5+5.5*ii,pi/3),
    PolarParticle2D(5.5+5.5*ii,-pi/2),
    PolarParticle2D(5.5+5.5*ii,-pi/4),
  };

  DirectlyOutput::Vec_PolarParticle2D(PG);
  std::cout<<"---------------\n";

  PG=VM::Evolve_CtnNml_VicsekModel_Periodic(PG,BoxInfo,Parameters,1);

  DirectlyOutput::Vec_PolarParticle2D(PG);
  std::cout<<"---------------\n";

}

//let me have a test!!
