#pragma once

namespace PARAM{
	#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
	const short int Dim = 2;
	#else
	const short int Dim = 3;
	#endif
	const PS::F64 SMTH = 1.2;
	const PS::F64 C_CFL = 0.3;
	const PS::U64 NUMBER_OF_DENSITY_SMOOTHING_LENGTH_LOOP = 3;
	//Balsara (1995)'s switch
	const bool FLAG_B95  = true;
	//Rosswog et al. (2000)'s switch
	const bool FLAG_R00 = false;
	//Strength of AV
	const PS::F64 AV_STRENGTH = 1.0;
};
