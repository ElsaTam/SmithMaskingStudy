#include "utils/paths.h"
#include "utils/params.h"
#include "utils/console.h"

namespace Path {

	// Paths to folders

	const std::vector<int>& resolutions() {
		return Parameters::get()->currentParams()->pathParams.resolutions;
	}


	const std::string& objFolder() {
		return Parameters::get()->currentParams()->pathParams.objDir;
	}

	const std::string& hfFolder() {
		return Parameters::get()->currentParams()->pathParams.hfDir;
	}

	const std::string& outputRootFolder() {
		return Parameters::get()->currentParams()->pathParams.outputsDir;
	}

	std::string subdivisionFolderName(int res) {
		return std::to_string(res) + "_subdivisions/";
	}
	std::string subdivisionFolderPath(const std::string& root, int res) {
		return root + (res > 0 ? subdivisionFolderName(res) : "");
	}

	std::string logs_Folder() {
		return outputRootFolder() + "logs/";
	}
	std::string renders_Folder() {
		return outputRootFolder() + "renders/";
	}
	std::string G1_Folder() {
		return outputRootFolder() + "G1/";
	}
	std::string G1_2D_Folder() {
		return G1_Folder() + "2D/";
	}
	std::string G1_2D_Folder(int res) {
		return subdivisionFolderPath(G1_2D_Folder(), res);
	}
	std::string G1_3D_Folder() {
		return G1_Folder() + "3D/";
	}
	std::string G1_3D_Folder(int res) {
		return subdivisionFolderPath(G1_3D_Folder(), res);
	}
	std::string GAF_Folder() {
		return outputRootFolder() + "GAF/";
	}
	std::string GAF_2D_Folder() {
		return GAF_Folder() + "2D/";
	}
	std::string GAF_2D_Folder(int res) {
		return subdivisionFolderPath(GAF_2D_Folder(), res);
	}
	std::string GAF_2D_Folder(int res, scal phiIn, scal thetaIn) {
		return GAF_2D_Folder(res) + "phi_" + std::to_string(phiIn) + "_theta_" + std::to_string(thetaIn) + "/";
	}
	std::string GAF_3D_Folder() {
		return GAF_Folder() + "3D/";
	}
	std::string GAF_3D_Folder(int res) {
		return subdivisionFolderPath(GAF_3D_Folder(), res);
	}
	std::string GAF_3D_Folder(int res, scal phiIn, scal thetaIn) {
		return GAF_3D_Folder(res) + "phi_" + std::to_string(phiIn) + "_theta_" + std::to_string(thetaIn) + "/";
	}
	std::string tabulations_Folder() {
		return outputRootFolder() + "tabulations/";
	}
	std::string tabulations_Folder(int res) {
		return subdivisionFolderPath(tabulations_Folder(), res);
	}
	std::string ambientOcclusion_Folder() {
		return outputRootFolder() + "ambient_occlusion/";
	}
	std::string ambientOcclusion_Folder(int res) {
		return subdivisionFolderPath(ambientOcclusion_Folder(), res);
	}
	std::string features_Folder() {
		return outputRootFolder() + "features/";
	}
	std::string features_Folder(int res) {
		return subdivisionFolderPath(features_Folder(), res);
	}

	// Paths to files

	const std::vector<std::string>& surfaceNames() {
		return Parameters::get()->currentParams()->pathParams.surfNames;
	}

	std::string objFile(const std::string& surfName, int res) {
		return subdivisionFolderPath(objFolder(), res) + surfName + ".obj";
	}

	std::string hfFile(const std::string& surfName) {
		return hfFolder() + surfName + ".png";
	}

	std::string clusterImg(const std::string& surfName, int K) {
		return renders_Folder() + surfName + "_" + std::to_string(K) + "_clusters.png";
	}
	std::string renderImg(const std::string& surfName, const LaunchParams& launchParams, scal phi, scal theta) {
		std::string filename = renders_Folder();
		filename += surfName;
		filename += "_phi_" + std::to_string(phi);
		filename += "_theta_" + std::to_string(theta);
		filename += "_size_" + std::to_string(launchParams.frame.size[0]);
		filename += "_samples_" + std::to_string(launchParams.camera.nPixelSamples);
		filename += "_border_" + std::to_string(launchParams.sideEffect.borderPercentage);
		filename += ".png";
		return filename;
	}
	std::string ambientOcclusionImg(const std::string& surfName, int res) {
		return ambientOcclusion_Folder(res) + surfName + ".png";
	}

	std::string tabulationD(const std::string& surfName, int res) {
		return tabulations_Folder(res) + surfName + "_D.csv";
	}
	std::string tabulationG1_rc(const std::string& surfName, int res) {
		return tabulations_Folder(res) + surfName + "_G1-rc.csv";
	}
	std::string tabulationG1_smith(const std::string& surfName, int res) {
		return tabulations_Folder(res) + surfName + "_G1-smith.csv";
	}
	std::string tabulationGAF_rc(const std::string& surfName, scal phiIn, scal thetaIn, int res) {
		return tabulations_Folder(res) + surfName + "_wi-" + std::to_string(thetaIn) + "-" + std::to_string(phiIn) + "_GAF-rc.csv";
	}
	std::string tabulationError(const std::string& surfName, int res) {
		return tabulations_Folder(res) + surfName + "_E.csv";
	}
	std::string featuresFile(int res) {
		return features_Folder(res) + "features.csv";
	}

	const std::string& ptxFile() {
		return Parameters::get()->currentParams()->pathParams.ptxFile;
	}
	const std::string& gnuplotExe() {
		return Parameters::get()->currentParams()->pathParams.gnuplotPath;
	}

	// Utils

	bool exists(const std::string& path) {
		struct stat s;
		return (stat(path.c_str(), &s) == 0);
	}
	bool isFile(const std::string& path) {
		struct stat s;
		return (stat(path.c_str(), &s) == 0) && (s.st_mode & S_IFREG);
	}
	bool isFolder(const std::string& path) {
		struct stat s;
		return (stat(path.c_str(), &s) == 0) && (s.st_mode & S_IFDIR);
	}

	bool checkFolder(std::string path, bool print) {
		if (isFolder(path)) {
			if (print) {
				Console::out << path << " : ";
				Console::succ << "ok" << std::endl;
			}
			return true;
		}
		else {
			if (print) {
				Console::out << path << " : ";
				Console::err << "not found or not a directory" << std::endl;
			}
			return false;
		}
	}
	bool checkFile(std::string path, bool print) {
		if (isFile(path)) {
			if (print) {
				Console::out << path << " : ";
				Console::succ << "ok" << std::endl;
			}
			return true;
		}
		else {
			if (print) {
				Console::out << path << " : ";
				Console::err << "not found or not a file" << std::endl;
			}
			return false;
		}
	}
	bool checkPaths(bool print) {
		print = print && Parameters::get()->currentParams()->outLevel > OutLevel::NO_OUTPUT;
		bool allPathCorrects = true;
		if (print) Console::out << "----- Root folders -----" << std::endl;
		allPathCorrects = checkFolder(logs_Folder(), print)             && allPathCorrects;
		allPathCorrects = checkFolder(renders_Folder(), print)          && allPathCorrects;
		allPathCorrects = checkFolder(G1_Folder(), print)               && allPathCorrects;
		allPathCorrects = checkFolder(G1_2D_Folder(), print)            && allPathCorrects;
		allPathCorrects = checkFolder(G1_3D_Folder(), print)            && allPathCorrects;
		allPathCorrects = checkFolder(GAF_Folder(), print)              && allPathCorrects;
		allPathCorrects = checkFolder(GAF_2D_Folder(), print)           && allPathCorrects;
		allPathCorrects = checkFolder(GAF_3D_Folder(), print)           && allPathCorrects;
		allPathCorrects = checkFolder(tabulations_Folder(), print)      && allPathCorrects;
		allPathCorrects = checkFolder(ambientOcclusion_Folder(), print) && allPathCorrects;
		allPathCorrects = checkFolder(features_Folder(), print)       && allPathCorrects;

		if (print) Console::out << "----- Resolution folders -----" << std::endl;
		for (int res : resolutions()) {
			allPathCorrects = checkFolder(G1_2D_Folder(res), print)            && allPathCorrects;
			allPathCorrects = checkFolder(G1_3D_Folder(res), print)            && allPathCorrects;
			allPathCorrects = checkFolder(GAF_2D_Folder(res), print)           && allPathCorrects;
			allPathCorrects = checkFolder(GAF_3D_Folder(res), print)           && allPathCorrects;
			allPathCorrects = checkFolder(tabulations_Folder(res), print)      && allPathCorrects;
			allPathCorrects = checkFolder(ambientOcclusion_Folder(res), print) && allPathCorrects;
			allPathCorrects = checkFolder(features_Folder(res), print)       && allPathCorrects;
		}

		if (print) Console::out << "----- Param files -----" << std::endl;
		allPathCorrects = checkFile(ptxFile(), print) && allPathCorrects;
		//allPathCorrects = checkFile(gnuplotExe(), print) && allPathCorrects;

		return allPathCorrects;
	}
}