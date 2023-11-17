#include "utils/console.h"
#include <sstream>

namespace Console {
	int BLACK = 0;
	int DARK_BLUE = 1;
	int DARK_GREEN = 2;
	int GRAY_BLUE = 3;
	int RED = 4;
	int PURPLE = 5;
	int ORANGE = 6;
	int LIGHT_GRAY = 7;
	int GRAY = 8;
	int BLUE = 9;
	int GREEN = 10;
	int SKY_BLUE = 11;
	int SALMON = 12;
	int MAGENTA = 13;
	int IVORY = 14;
	int WHITE = 15;

	void printColor()
	{
		for (int b = 0; b < 16; ++b) {
			for (int t = 0; t < 16; ++t) {
				HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
				SetConsoleTextAttribute(hConsole, b * 16 + t);
				std::cout << "Hellow World";
				SetConsoleTextAttribute(hConsole, 15);
				std::cout << " [" << b << ":" << t << "]" << std::endl;
			}
		}
	}


	OutStream out{ &std::cout };			// default
	OutStream info{ &std::cout, SKY_BLUE }; // information
	OutStream prog{ &std::cout, IVORY };    // progress
	OutStream light{ &std::cout, GRAY };    // light, discrete
	OutStream succ{ &std::cout, GREEN };    // success
	OutStream warn{ &std::cout, ORANGE };   // warning
	OutStream err{ &std::cerr, RED };       // error

	OutStream optixOut{ &std::cout, WHITE, DARK_BLUE };
	OutStream optixInfo{ &std::cout, SKY_BLUE, DARK_BLUE };
	OutStream optixSucc{ &std::cout, GREEN, DARK_BLUE };
	OutStream optixWarn{ &std::cout, ORANGE, DARK_BLUE };
	OutStream optixErr{ &std::cerr, RED, DARK_BLUE };


	TimeStamp timeStamp{ };
	std::string timePad(std::string(11, ' '));

	std::string line(std::string(80, '-') + "\n");
	std::string shortline(std::string(20, '-') + "\n");
	std::string indent("|  ");

	std::string fixed_size(float value, int nCharacters)
	{
		// Remove one character if negative
		if (value < 0) nCharacters--;
		// Remove as many character as there is in the integer part
		int nDigitsInt = int(log10(std::max((int)std::floor(abs(value)), 1)) + 1);
		nCharacters -= nDigitsInt;
		// Remove one character for the point
		nCharacters--;

		std::stringstream stream;
		stream << std::fixed << std::setprecision(nCharacters) << value;
		return stream.str();
	}
}