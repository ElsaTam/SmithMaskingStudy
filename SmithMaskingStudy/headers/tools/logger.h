#pragma once

#include "shapes/TriangleMesh.h"
#include "tools/statistics.h"
#include "utils/duration.h"

#define LOG_MESSAGE( message )			Logger::getInstance().logMessage( __FUNCTION__, message)
#define LOG_ERROR( message )			Logger::getInstance().logError( __FILE__, __FUNCTION__, __LINE__, message)
#define LOG_TIME( start, end )			Logger::getInstance().logTime( __FUNCTION__, start, end)
#define LOG_TIME_OBJ( obj, start, end)	Logger::getInstance().logTime( __FUNCTION__, obj, start, end)
#define LOG_NAME( folder, name )		Logger::getInstance().logName( folder, name)
#define LOG_OBJ( obj )					Logger::getInstance().logObj( obj)
#define LOG_STATS( stats )				Logger::getInstance().logStats( stats)
#define LOG_HLINE						Logger::getInstance().hSeparation()
#define LOG_NEW							Logger::getInstance().generateNewId()
#define LOG_ID( id )					Logger::getInstance().changeCurrentLog(id);

class Logger
{
private:
	bool activated;
	unsigned int currentId;
	std::string filename;
	std::string folder;

	static Logger* singleton;
	Logger();
	void updateFilename();

public:
	Logger(const Logger&) = delete;
	Logger& operator=(const Logger&) = delete;
	static Logger& getInstance();

	void enable(bool _activated);
	void setFolder(const std::string& _folder);

	void generateNewId();
	void changeCurrentLog(unsigned int newId);

	void removeLog(unsigned int id) const;
	void removeLog() const;

	void clearLog() const;

	void logMessage(const char* func, const std::string& message) const;
	void logError(const char* file, const char* func, int line, const std::string& message) const;

	void logName(const std::string& folder, const std::string& name) const;
	void logObj(const TriangleMesh* const mesh) const;
	void logStats(const StatisticsTool& stats) const;

	void logTime(const char* func, const std::string& objName,
		time_point start, time_point end) const;

	void logTime(const char* func,
		time_point start, time_point end) const;

	void hSeparation() const;
};