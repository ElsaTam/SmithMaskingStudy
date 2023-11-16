#include "tools/objWriter.h"
#include "utils/console.h"

#include <fstream>

void ObjWriter::writeObj(const std::string& filename, const TriangleMesh* mesh)
{
	std::ofstream file;
	file.open(filename.c_str());

	file << "o " << mesh->name << std::endl;

	file << std::endl;

	for (const gdt::vec3sc& v : mesh->vertex)
	{
		file << "v " << v.x << " " << v.y << " " << v.z << std::endl;
	}

	file << std::endl;

	for (const gdt::vec3sc& n : mesh->vertex_normal)
	{
		file << "vn " << n.x << " " << n.y << " " << n.z << std::endl;
	}
	
	file << std::endl;

	for (const gdt::vec3i& idx : mesh->index)
	{
		file << "f " << (idx[0] + 1) << "//" << (idx[0] + 1) << " " << (idx[1] + 1) << "//" << (idx[1] + 1) << " " << (idx[2] + 1) << "//" << (idx[2] + 1) << std::endl;
	}

	file.close();
	Console::succ << Console::timeStamp << "OBJ written at: " << filename << std::endl;
}