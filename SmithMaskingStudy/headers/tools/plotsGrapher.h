#pragma once

#include "gdt/math/vec.h"
#include <vector>
#include <string>

/**
 * @brief Namespace containing methods used to plot graphs, such as the
 * normal distribution graphs.
 * 
 * Look for available palettes at https://github.com/Gnuplotting/gnuplot-palettes
 * 
 */
namespace PlotsGrapher
{
    // -------------------------- CLOUD POINTS --------------------------------
    // scatter 2D points cloud
    void scatter2D(const std::string& filename, const std::vector<scal>& X, const std::vector<scal>& Y);
    // scatter 3D points cloud
    void scatter3D(const std::string& filename, const std::vector<scal>& theta, const std::vector<scal>& phi);


    // ---------------------------- 2D CURVES ---------------------------------
    // plot 2D curves on one plot
    void plot(const std::string& filename,
        const std::vector<scal>& X,
        const std::vector<std::vector<scal>>& Y,
        const std::vector<std::string>& titles);
    void plotG1(const std::string& filename,
        const std::vector<scal>& X,
        const std::vector<std::vector<scal>>& Y,
        const std::vector<std::string>& titles);
    void plotDistrib(const std::string& filename,
        const std::vector<scal>& X,
        const std::vector<scal>& Y1,
        const std::vector<scal>& Y2 = std::vector<scal>());

    // --------------------------- 3D SURFACES --------------------------------
    // plot one 3D surface
    void splot(const std::string& filename,
        const std::vector<std::vector<gdt::vec3sc>>& points);
    // plot three 3D surface to visualize a difference
    void splotDiff(const std::string& filename,
        const std::vector<std::vector<gdt::vec3sc>>& points);

    // -------------------------- PROJ COLORMAP -------------------------------
    // plot colormap for a visibility term (in 0:1)
    void cmplot(const std::string& filename, const std::vector<gdt::vec3sc>& col);
    // plot colormap for a diff between two visibility terms (in 0:diffMax)
    void cmplotDiff(const std::string& filename, const std::vector<gdt::vec3sc>& col, std::string palette, scal diffMax = -1);
};