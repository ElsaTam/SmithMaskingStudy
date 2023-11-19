#pragma once

#define NOMINMAX
#include <windows.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>

namespace Console {

	extern void printColor();

	struct OutStream {
		std::ostream* stream;
		int color;
		int backgroundColor;

		OutStream(std::ostream* dest, int _color = 15, int _backgroundColor = 0) {
			stream = dest;
			color = _color;
			backgroundColor = _backgroundColor;
		}

		template <typename T>
		OutStream& operator<<(const T& obj)
		{
			if (stream) {
				HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
				SetConsoleTextAttribute(hConsole, backgroundColor * 16 + color);
				*stream << obj;
				SetConsoleTextAttribute(hConsole, 15);
			}
			return *this;
		}

		OutStream& operator<<(std::ostream& (*f)(std::ostream&))
		{
			HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
			SetConsoleTextAttribute(hConsole, backgroundColor * 16 + color);
			*stream << *f;
			SetConsoleTextAttribute(hConsole, 15);
			return *this;
		}
	};

	extern OutStream out;   // default
	extern OutStream info;  // information
	extern OutStream prog;  // progress
	extern OutStream light; // light, discrete
	extern OutStream succ;  // success
	extern OutStream warn;  // warning
	extern OutStream err;   // error

	extern OutStream optixOut;
	extern OutStream optixInfo;
	extern OutStream optixSucc;
	extern OutStream optixWarn;
	extern OutStream optixErr;

	void print(OutLevel level, std::string message, const char* end = "\n");
	void printOptix(OutLevel level, std::string message, const char* end = "\n");

	struct TimeStamp {
		std::string str() {
			std::stringstream ss;
			ss << this;
			return ss.str();
		}
		friend std::ostream& operator<< (std::ostream& os, const TimeStamp& _timeStamp)
		{
			auto now = std::chrono::system_clock::now();
			std::time_t now_time = std::chrono::system_clock::to_time_t(now);
			std::tm ltime;
			localtime_s(&ltime, &now_time);
			os << "[" << std::setfill('0')
				<< std::setw(2) << ltime.tm_hour << ":"
				<< std::setw(2) << ltime.tm_min << ":"
				<< std::setw(2) << ltime.tm_sec << "] ";
			return os;
		}
	};

	extern TimeStamp timeStamp;
	extern std::string timePad;

	extern std::string line;
	extern std::string shortline;
	extern std::string indent;

	std::string fixed_size(float value, int nCharacters);
}