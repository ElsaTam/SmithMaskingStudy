#pragma once

#include "NDFs/Microfacet.hpp"
#include "OptiX/OptixRenderer.h"
#include "tools/csvWriter.h"
#include "tools/statistics.h"


class Analyzer
{
private:
	typedef std::vector<scal> Vec_Sc;
	typedef std::vector<Vec_Sc> Vec_VSc;
	typedef std::vector<std::string> Vec_Str;

	struct DirectionData {
		scal phiStart;
		scal phiRange;
		scal thetaStart;
		scal thetaRange;
		int nPhi;
		int nTheta;
	};

	DirectionData D, DOut, DIn;

	TriangleMesh* mesh;
	OptixRenderer* optixRenderer;

	void prepareGPU();

	void logRenderingInfo() const;

	// utils
	//std::unique_ptr<MicrofacetDistribution> getTheoricalNDF() const;
	scal partial_error(scal ref, scal estimation) const;
	scal normalize_error(scal sum_E, int N) const;

	// tabulation
	void tabulateFunctions(std::vector<std::string> filenames, std::vector<scal(*)(void*, scal, scal, scal, scal)> T_functions, std::vector<void*> contexts, scal phiIn = 0, scal thetaIn = 0);


public:

	Analyzer(TriangleMesh* mesh = nullptr, bool useGPU = false);
	~Analyzer();

	void setGeo(TriangleMesh* mesh);

	/**
	 * @brief Compute masking with ray tracing (\f$G_1^{rc}\f$) and with Smith equation (\f$G_1\f$). 
	 * 
	 */
	void G1();

	/**
	 * @brief Compute GAF with ray tracing and with Smith equation (\f$G_1 * G_1\f$).
	 *
	 */
	void GAF();

	scal error();

	// Tabulation
	void writeTheta(csv::CSVWriter& writer, bool fromDistrib = false) const;
	void tabulate(bool D, bool G1_Ashikhmin, bool G1_RT);
	void tabulateDistrib();
	void tabulateG1_Ashikhmin();
	void tabulateG1_RT();
	void tabulateGAF_RT();

	// Ambient occlusion
	void ambientOcclusion();

	// Decoupe en sets
	void sets();

	// Statistics
	void statistics(bool computeError);

	void fullPipeline();
};

