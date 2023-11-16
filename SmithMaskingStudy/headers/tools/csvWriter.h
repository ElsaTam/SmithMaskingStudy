#pragma once

#include <fstream>
#include <variant>
#include <vector>
#include "gdt/math/vec.h"

namespace csv {
    struct elem
    {
        enum class Tag { STRING, SCAL_ARR, INT, SCAL, FLOAT, DOUBLE, LONG_DOUBLE } tag;
        std::variant<std::string, std::vector<scal>, int, float, double, long double> value;

        std::string toString();
    };

    class CSVWriter {
    private:
        char separator = ',';
        std::string filename = "";
        std::ofstream file;
        int lines = 0;

    public:
        CSVWriter();
        CSVWriter(const std::string& filename, std::ios_base::openmode mode = std::ios_base::trunc);
        ~CSVWriter();
        void writeRow(std::vector<elem> datas);
        void close();
        int numberOfLines() const;
        const std::string& getFilename() const { return filename; }
    };
};