#pragma once

#include <string>

namespace Path {

	std::string subdivisionFolder(int res);

	std::string objPath(const std::string& name, int res);

	std::string hfPath(const std::string& name);

	bool fileExists(const std::string& path);
}