#pragma once

#include "shapes/TriangleMesh.h"

namespace ObjWriter {

	void writeObj(const std::string& filename, const TriangleMesh* mesh);

};