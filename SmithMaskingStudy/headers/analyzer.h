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
	scal dSmooth(int N, int i) const;
	scal dLinear(int N, int i) const;
	std::unique_ptr<MicrofacetDistribution> getTheoricalNDF() const;
	std::string getFolder(const std::string& root) const;

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

	// Tabulation
	void writeTheta(csv::CSVWriter& writer, bool fromDistrib = false) const;
	void tabulateFunctions(std::vector<std::string> filenames, std::vector<scal (*)(void*, scal, scal, scal, scal)> T_functions, std::vector<void*> contexts, scal phiIn = 0, scal thetaIn = 0);
	void tabulate(bool D, bool G1_Ashikhmin, bool G1_RT, Discrete* NDF = nullptr);
	void tabulateDistrib(Discrete* NDF = nullptr);
	void tabulateG1_Ashikhmin(Discrete* NDF = nullptr);
	void tabulateG1_RT();
	void tabulateGAF_RT();
	// Heights histograms
	void tabulateHeights();

	// Ambient occlusion
	void ambientOcclusion();

	// Error
	ErrorStats error();
	scal error(scal theta, scal phi, const MicrofacetDistribution* NDF);

	// Decoupe en sets
	void sets();

	// Statistics (without GPU)
	void distrib() const;
	void statistics(csv::CSVWriter* writer = nullptr); // with GPU for error

	void fullPipeline();
};

