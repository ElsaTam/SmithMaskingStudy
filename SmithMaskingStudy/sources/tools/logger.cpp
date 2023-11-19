#include "tools/logger.h"
#include "utils/params.h"
#include "utils/console.h"
#include "magic_enum.hpp"

#include <sstream>
#include <iomanip>
#include <fstream>
#include <boost/algorithm/string.hpp>


Logger& Logger::getInstance()
{
	static Logger instance;
	return instance;
}

void Logger::enable(bool _activated)
{
	this->activated = _activated;
}

void Logger::setFolder(const std::string& _folder)
{
	this->folder = _folder;
	this->generateNewId();
}

bool exist(const std::string& filename)
{
	struct stat buffer;
	return (stat(filename.c_str(), &buffer) == 0);
}

std::string findFilename(const std::string& folder, unsigned int id)
{
	std::ostringstream oss;
	oss << std::setfill('0') << std::setw(5) << id;
	std::string base = std::string(magic_enum::enum_name(Parameters::get()->currentParams()->methodParams.method));
	boost::to_lower(base);
	return folder + base + "_" + oss.str() + ".txt";
}

Logger::Logger()
{
	this->generateNewId();
}

void Logger::updateFilename()
{
	this->filename = findFilename(folder, currentId);
}

void Logger::generateNewId()
{
	this->currentId = 0;
	this->updateFilename();

	while (exist(this->filename))
	{
		this->currentId++;
		this->updateFilename();
	}
}

void Logger::changeCurrentLog(unsigned int newId)
{
	this->currentId = newId;
	this->updateFilename();
}

void Logger::removeLog(unsigned int id) const
{
	if (this->activated) {
		remove(findFilename(this->folder, id).c_str());
	}
}
void Logger::removeLog() const
{
	if (this->activated) {
		remove(this->filename.c_str());
	}
}

void Logger::clearLog() const
{
	if (this->activated) {
		std::ofstream myfile;
		myfile.open(this->filename);
		myfile.clear();
		myfile.close();
	}
}


void Logger::logMessage(const char* func, const std::string& message) const
{
	if (this->activated) {
		std::ofstream myfile;
		myfile.open(this->filename, std::ofstream::out | std::ofstream::app);

		myfile << "[" << func << "] " << message << std::endl;

		myfile.close();
	}
}

void Logger::logError(const char* file, const char* func, int line, const std::string& message) const
{
	if (this->activated) {
		std::ofstream myfile;
		myfile.open(this->filename, std::ofstream::out | std::ofstream::app);

		myfile << "ERROR: " << file << ":" << line << std::endl;
		myfile << "[" << func << "] " << message << std::endl;

		myfile.close();
	}
}

void Logger::logName(const std::string& folder, const std::string& name) const
{
	if (this->activated) {
		std::ofstream myfile;
		myfile.open(this->filename, std::ofstream::out | std::ofstream::app);

		myfile << std::endl << std::endl;
		myfile << folder << name << std::endl;
		myfile << "-----------------------" << std::endl;

		myfile.close();
	}
}

void Logger::logObj(const TriangleMesh* const mesh) const
{
	if (this->activated) {
		size_t numTriangles = mesh->index.size();;
		size_t numVertices = mesh->vertex.size();

		std::ofstream myfile;
		myfile.open(this->filename, std::ofstream::out | std::ofstream::app);
		myfile << *mesh << std::endl;
		myfile.close();
	}
}

void Logger::logStats(const StatisticsTool& stats) const
{
	if (this->activated) {
		std::ofstream myfile;
		myfile.open(this->filename, std::ofstream::out | std::ofstream::app);
		myfile << stats << std::endl;
		myfile.close();
	}
}

void Logger::logTime(const char* func, const std::string& objName,
					 time_point start, time_point end) const
{
	if (this->activated) {
		Duration duration(start, end);

		std::ofstream myfile;
		myfile.open(this->filename, std::ofstream::out | std::ofstream::app);

		myfile << "[" << func << "(" << objName << ")" << "] "
			   << duration.str()
			   << std::endl;

		myfile.close();
	}
}

void Logger::logTime(const char* func,
	time_point start, time_point end) const
{
	if (this->activated) {
		Duration duration(start, end);

		std::ofstream myfile;
		myfile.open(this->filename, std::ofstream::out | std::ofstream::app);

		myfile << "[" << func << "] "
			   << duration.str()
			   << std::endl;

		myfile.close();
	}
}

void Logger::hSeparation() const
{
	if (this->activated) {
		std::ofstream myfile;
		myfile.open(this->filename, std::ofstream::out | std::ofstream::app);

		myfile << std::endl;
		myfile << std::setfill('=') << std::setw(80) << "" << std::endl;
		myfile << std::setfill(' ') << std::setw(15) << "" << Console::timeStamp << std::endl;
		myfile << std::setfill('=') << std::setw(80) << "" << std::endl;

		myfile.close();
	}
}