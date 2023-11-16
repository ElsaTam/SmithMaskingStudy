#include "tools/csvWriter.h"

#include <iomanip>
#include <sstream>
#include "utils/console.h"

std::string csv::elem::toString()
{
	switch (tag)
	{
	case elem::Tag::STRING:
		return std::get<std::string>(value);
	case elem::Tag::INT:
		return std::to_string(std::get<int>(value));
	case elem::Tag::SCAL_ARR:
	{
		const std::vector<scal>& arr = std::get<std::vector<scal>>(value);
		std::stringstream stream;
		stream << "[";
		stream << std::fixed << std::setprecision(std::numeric_limits<scal>::max_digits10) << arr[0];
		for (int i = 1; i < arr.size(); ++i) {
			stream << std::fixed << std::setprecision(std::numeric_limits<scal>::max_digits10) << ", " << arr[i];
		}
		stream << "]";
		return stream.str();
	}
	case elem::Tag::SCAL:
	{
		std::stringstream stream;
		stream << std::fixed << std::setprecision(std::numeric_limits<scal>::max_digits10) << std::get<scal>(value);
		return stream.str();
	}
	case elem::Tag::FLOAT:
	{
		std::stringstream stream;
		stream << std::fixed << std::setprecision(std::numeric_limits<float>::max_digits10) << std::get<float>(value);
		return stream.str();
	}
	case elem::Tag::DOUBLE:
	{
		std::stringstream stream;
		stream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << std::get<double>(value);
		return stream.str();
	}
	case elem::Tag::LONG_DOUBLE:
	{
		std::stringstream stream;
		stream << std::fixed << std::setprecision(std::numeric_limits<long double>::max_digits10) << std::get<long double>(value);
		return stream.str();
	}
	default:
		return "";
	}
}

csv::CSVWriter::CSVWriter() { }

csv::CSVWriter::CSVWriter(const std::string& _filename, std::ios_base::openmode mode)
	: filename(_filename)
{
	file.open(filename, mode);

	lines = 0;

	std::ifstream in(filename);
	std::string unused;
	while (std::getline(in, unused))
		++lines;
	in.close();
}

csv::CSVWriter::~CSVWriter()
{
	close();
}

void csv::CSVWriter::writeRow(std::vector<elem> datas)
{
	if (!file.is_open()) return;

	for (int i = 0; i < datas.size() - 1; ++i) {
		file << datas[i].toString() << separator;
	}
	file << datas.back().toString() << std::endl;
	++lines;
}

void csv::CSVWriter::close()
{
	if (!file.is_open()) return;
	file.close();

	Console::succ << Console::timePad << "CSV saved to " << filename << " ... done." << std::endl;
}

int csv::CSVWriter::numberOfLines() const
{
	return lines;
}