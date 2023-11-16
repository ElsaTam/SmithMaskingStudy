#pragma once

#include <string>
#include <vector>
#include "tools/csvWriter.h"
#include "shapes/TriangleMesh.h"


/**
 * @brief Class used to compute and store statistics of datas.
 * 
 */
class StatisticsTool
{
    csv::elem::Tag DECIMAL_TAG = csv::elem::Tag::LONG_DOUBLE;
    typedef long double ld;
    typedef gdt::vec_t<ld, 2> vec2_ld;
    typedef gdt::vec_t<ld, 3> vec3_ld;
    typedef gdt::box_t<vec2_ld> box2_ld;

private:
    friend std::ostream& operator<<(std::ostream& output, const StatisticsTool& stats);
    const TriangleMesh* mesh;

    // error
    scal m_error{ 0 };

    // phi
    ld m_phiAnisotropy{ 0 };

    // thetas
    ld m_thetaMean{ 0 };
    ld m_thetaStd{ 0 };
    gdt::interval<ld> m_thetaRange;
    ld m_thetaCV{ 0 };     // Coefficient of Variation (or Relative Standard Deviation)
    ld m_thetaMAD{ 0 };    // Median Absolute Deviation
    ld m_thetaQ1{ 0 };     // First quartile (25%)
    ld m_thetaMedian{ 0 }; // Second quartile (50%)
    ld m_thetaQ3{ 0 };     // Third quartile (75%)
    ld m_thetaIQR{ 0 };    // Interquartile range
    ld m_thetaQCD{ 0 };    // Quartile coefficient of dispersion
    ld m_thetaSkewness{ 0 };
    ld m_thetaKurtosis{ 0 };

    // areas
    ld m_areaTotal{ 0 };
    ld m_areaMean{ 0 };
    ld m_areaStd{ 0 };
    gdt::interval<ld> m_areaRange;
    ld m_areaCV{ 0 };     // Coefficient of Variation (or Relative Standard Deviation)
    ld m_areaMAD{ 0 };    // Median Absolute Deviation
    ld m_areaQ1{ 0 };     // First quartile (25%)
    ld m_areaMedian{ 0 }; // Second quartile (50%)
    ld m_areaQ3{ 0 };     // Third quartile (75%)
    ld m_areaIQR{ 0 };    // Interquartile range
    ld m_areaQCD{ 0 };    // Quartile coefficient of dispersion

    // heights
    ld m_heightMean{ 0 };
    ld m_heightStd{ 0 };
    gdt::interval<ld> m_heightRange;
    ld m_heightCV{ 0 };     // Coefficient of Variation (or Relative Standard Deviation)
    ld m_heightMAD{ 0 };    // Median Absolute Deviation
    ld m_heightQ1{ 0 };     // First quartile (25%)
    ld m_heightMedian{ 0 }; // Second quartile (50%)
    ld m_heightQ3{ 0 };     // Third quartile (75%)
    ld m_heightIQR{ 0 };    // Interquartile range
    ld m_heightQCD{ 0 };    // Quartile coefficient of dispersion
    ld m_heightCb{ 0 };     // Lowest cluster proportion

    // Correlations
    ld m_ATCorr{ 0 }; // Correlation between theta and area
    ld m_THCorr{ 0 }; // Correlation between theta and height
    ld m_AHCorr{ 0 }; // Correlation between areas and height

    static int wHead;
    static int wCell;
    static int nColumns;

    scal computeAnisotropy() const;
    
public:
    StatisticsTool(const TriangleMesh* _mesh, scal error = 0, int n_features = -1);
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