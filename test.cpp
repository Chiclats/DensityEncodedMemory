//#define vs

#include"ParticleSimulationHeaders.hpp"

int main() {

	auto PG = StateProcessing::GenerateUniformRandomParticleGroup(1000, 1000, 1000 * 1000);
	FileIO::OutParticleVec("0.txt", PG);

	const vector<double> ParameterList = {
		1000,1000,1,0.5,0.15,-0.0001
	};
	const auto BoxInfo = BoxMeshing_2D::BoxInfo(1000, 1000, 1);

	PG = VicsekModel_2D::Evolve_VicsekModel_Periodic(
		PG,BoxInfo,ParameterList,10,0
	);

	FileIO::OutParticleVec("10.txt", PG);

}
