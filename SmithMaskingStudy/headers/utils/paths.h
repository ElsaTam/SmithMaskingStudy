#pragma once

#include <string>
#include "OptiX/LaunchParams.h"

namespace Path {

	// Inputs

	const std::string& objFolder();
	const std::string& hfFolder();

	std::string objFile(const std::string& surfName, int res);
	std::string hfFile(const std::string& surfName);

	const std::vector<std::string>& surfaceNames();
	const std::vector<int>& resolutions();

	// Outputs

	const std::string& outputRootFolder();

	std::string logs_Folder();
	std::string renders_Folder();
	std::string G1_Folder();
	std::string G1_2D_Folder();
	std::string G1_2D_Folder(int res);
	std::string G1_3D_Folder();
	std::string G1_3D_Folder(int res);
	std::string GAF_Folder();
	std::string GAF_2D_Folder();
	std::string GAF_2D_Folder(int res);
	std::string GAF_2D_Folder(int res, scal phiIn, scal thetaIn);
	std::string GAF_3D_Folder();
	std::string GAF_3D_Folder(int res);
	std::string GAF_3D_Folder(int res, scal phiIn, scal thetaIn);
	std::string tabulations_Folder();
	std::string tabulations_Folder(int res);
	std::string ambientOcclusion_Folder();
	std::string ambientOcclusion_Folder(int res);
	std::string statistics_Folder();
	std::string statistics_Folder(int res);

	const std::string& ptxFile();
	const std::string& gnuplotExe();

	std::string clusterImg(const std::string& surfName, int K);
	std::string renderImg(const std::string& surfName, const LaunchParams& launchParams, scal phi, scal theta);
	std::string ambientOcclusionImg(const std::string& surfName, int res = 0);

	std::string tabulationD(const std::string& surfName, int res = 0);
	std::string tabulationG1_rc(const std::string& surfName, int res = 0);
	std::string tabulationG1_smith(const std::string& surfName, int res = 0);
	std::string tabulationGAF_rc(const std::string& surfName, scal phiIn, scal thetaIn, int res = 0);
	std::string tabulationError(const std::string& surfName, int res = 0);
	std::string statisticsFile(int res = 0);


	// Utils

	std::string subdivisionFolderName(int res);
	std::string subdivisionFolderPath(const std::string& root, int res);

	bool exists(const std::string& path);
	bool isFile(const std::string& path);
	bool isFolder(const std::string& path);

	bool checkPaths(bool print = false);
}