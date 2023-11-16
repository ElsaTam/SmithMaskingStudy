#pragma once

#include <string>
#include <vector>
#include "tools/csvWriter.h"
#include "shapes/TriangleMesh.h"

typedef long double ld;
typedef gdt::vec_t<ld, 2> vec2_ld;
typedef gdt::vec_t<ld, 3> vec3_ld;
typedef gdt::box_t<vec2_ld> box2_ld;

struct ErrorStats
{
    ld integral = 0;
    ld mean = 0;
    ld std = 0;
    ld max = 0;
    ld min = 1;
};

/**
 * @brief Class used to compute and store statistics of datas.
 * 
 */
class StatisticsTool
{
private:
    friend std::ostream& operator<<(std::ostream& output, const StatisticsTool& stats);
    const TriangleMesh* mesh;

    // error
    ErrorStats m_Error{ };

    // scores
    gdt::vec2i m_heightSums{ 0 };
    ld m_slopeScore{ 0 };
    ld m_anisotropyScore{ 0 };

    // normal cartesian coordinates
    vec3_ld m_nXYZMean{ 0 };
    vec3_ld m_nXYZStd{ 0 };
    box2_ld m_nXYRanges;
    ld m_nXYCorr{ 0 };

    // normal spherical coordinates
    vec2_ld m_nTPMean{ 0 };
    vec2_ld m_nTPStd{ 0 };
    box2_ld m_nTPRanges;
    ld m_nTPCorr{ 0 };

    // thetas
    ld m_TCV{ 0 };     // Coefficient of Variation (or Relative Standard Deviation)
    ld m_TMAD{ 0 };    // Median Absolute Deviation
    ld m_TQ1{ 0 };     // First quartile (25%)
    ld m_TMedian{ 0 }; // Second quartile (50%)
    ld m_TQ3{ 0 };     // Third quartile (75%)
    ld m_TIQR{ 0 };    // Interquartile range
    ld m_TQCD{ 0 };    // Quartile coefficient of dispersion

    // slopes
    vec2_ld m_slopeMean{ 0 };
    vec2_ld m_slopeStd{ 0 };
    box2_ld m_slopeRanges;
    ld m_slopeCorr{ 0 };

    // areas
    ld m_totalA{ 0 };
    ld m_areaMean{ 0 };
    ld m_areaStd{ 0 };
    gdt::interval<ld> m_areaRange;

    // cluster
    ld m_lowH{ 0 };

    // heights
    ld m_heightMean{ 0 };
    ld m_heightStd{ 0 };
    ld m_maxDeltaHeight{ 0 };
    gdt::interval<ld> m_heightRange;
    ld m_heightCV{ 0 };     // Coefficient of Variation (or Relative Standard Deviation)
    ld m_heightMAD{ 0 };    // Median Absolute Deviation
    ld m_heightQ1{ 0 };     // First quartile (25%)
    ld m_heightMedian{ 0 }; // Second quartile (50%)
    ld m_heightQ3{ 0 };     // Third quartile (75%)
    ld m_heightIQR{ 0 };    // Interquartile range
    ld m_heightQCD{ 0 };    // Quartile coefficient of dispersion

    // shape parameters
    ld m_skewness{ 0 };
    ld m_kurtosis{ 0 };

    // Correlations
    ld m_ATCorr{ 0 }; // Correlation between theta and area
    ld m_THCorr{ 0 }; // Correlation between theta and height

    static int wHead;
    static int wCell;
    static int nColumns;

    scal computeAnisotropy() const;
    
public:
    StatisticsTool(const TriangleMesh* _mesh, ErrorStats error = { }, int n_features = -1);
    ~StatisticsTool();

    // set variables
    void compute5Statistics();
    void compute25Statistics();
    void computeStatistics();

    // Export data
    static void CSVHeader(csv::CSVWriter* writer);
    void toCSV(csv::CSVWriter* writer);

    void print();
};