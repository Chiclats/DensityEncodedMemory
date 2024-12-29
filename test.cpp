#include"./ParticleSimulationHeaders.hpp"
//#include"./ContinuousVicsekModel.hpp"
#include"./NematicAlgnPP.hpp"

namespace NAPP_2D=NematicAlgnPolarParticles_2D;
//namespace NAPP_2D=VicsekModel_2D;

int main(){

  const double
    SystemSize_X=10,
    SystemSize_Y=10,
    CharLength=1,
    Velocity=1,
    Sigma=0,//0.1,
    GuidingCoeff=-1;

  const vector<double> Parameters={SystemSize_X,SystemSize_Y,CharLength,Velocity,Sigma,GuidingCoeff};

  const auto BoxInfo=BoxMeshing_2D::BoxInfo(SystemSize_X,SystemSize_Y,CharLength);
  const auto LDAdjBox=NAPP_2D::LowerRightAdjacentBoxIndex(BoxInfo);

  PolarParticle2D_GT PG={
    // PolarParticle2D(5.5+5.5*ii,0.),
    // PolarParticle2D(5.5+5.5*ii,0.),
    // PolarParticle2D(5.5+5.5*ii,pi/3),
    // PolarParticle2D(7.5+5.5*ii,pi/2),
    // PolarParticle2D(7.5+5.5*ii,-pi/4),

    PolarParticle2D(5.5+5.5*ii,2*pi/3),
  };

  DirectlyOutput::Vec_PolarParticle2D(PG);
  std::cout<<"---------------\n";

   PG=NAPP_2D::Evolve_NAPP_Periodic(PG,BoxInfo,Parameters,1);
  //PG=NAPP_2D::Evolve_VicsekModel_Periodic(PG,BoxInfo,Parameters,1);

  DirectlyOutput::Vec_PolarParticle2D(PG);
  std::cout<<"---------------\n";

}

