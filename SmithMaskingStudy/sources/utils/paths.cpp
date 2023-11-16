#include "utils/paths.h"
#include "utils/params.h"

namespace Path {
	std::string subdivisionFolder(int res) {
		return std::to_string(res) + "_subdivisions/";
	}

	std::string objPath(const std::string& name, int res) {
		return Parameters::userParams.pathParams.objFolder + subdivisionFolder(res) + name + ".obj";
	}

	std::string hfPath(const std::string& name) {
		return Parameters::userParams.pathParams.hfFolder + name + ".png";
	}

	bool fileExists(const std::string& path) {
		FILE* file;
		errno_t err = fopen_s(&file, path.c_str(), "r");
		if (err == 0) {
			fclose(file);
			return true;
		}
		else {
			return false;
		}
	}
}