#include <sstream>
#include <iomanip>

#include "utils/duration.h"

Duration::Duration(time_point start, time_point end) {
	auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	milliseconds = (long) int_ms.count();
	seconds = (double)(milliseconds / 1000.0);
	minutes = (float) seconds / 60.f;
	hours = minutes / 60.f;
	days = hours / 24.f;
}

Duration::Duration(long milliseconds) : milliseconds(milliseconds) {
	seconds = (double)(milliseconds / 1000.0);
	minutes = (float) seconds / 60.f;
	hours = minutes / 60.f;
	days = hours / 24.f;
}

std::string Duration::str() const {
	std::ostringstream ms; ms << std::setfill('0') << std::setw(3) << (int)( milliseconds % 1000);
	std::ostringstream s;   s << std::setfill('0') << std::setw(2) << (int)( (int)seconds % 60 );
	std::ostringstream m;   m << std::setfill('0') << std::setw(2) << (int)( ((int)seconds % 3600) / 60 );
	std::ostringstream h;   h << std::setfill('0') << std::setw(2) << (int)( seconds / 3600 );
	return h.str() + "h " + m.str() + "m " + s.str() + "s " + ms.str() + "ms";
}
