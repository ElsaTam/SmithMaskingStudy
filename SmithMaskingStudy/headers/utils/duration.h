#pragma once

#include <chrono>

typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point;

class Duration
{
public:
	Duration(time_point start, time_point end);
	Duration(long milliseconds);

	std::string str() const;

	long milliseconds;
	double seconds;
	float minutes;
	float hours;
	float days;
};