#include "tools/plotsGrapher.h"

#include "gnuplot-iostream.h"
#include "utils/console.h"
#include "utils/params.h"
#include "utils/paths.h"




// Tells gnuplot-iostream how to print objects of class vec3sc.
namespace gnuplotio {
    template<typename T>
    struct BinfmtSender<gdt::vec_t<T, 3>> {
        static void send(std::ostream& stream) {
            BinfmtSender<T>::send(stream);
            BinfmtSender<T>::send(stream);
            BinfmtSender<T>::send(stream);
        }
    };

    template <typename T>
    struct BinarySender<gdt::vec_t<T, 3>> {
        static void send(std::ostream& stream, const gdt::vec_t<T, 3>& v) {
            BinarySender<T>::send(stream, v.x);
            BinarySender<T>::send(stream, v.y);
            BinarySender<T>::send(stream, v.z);
        }
    };

    // We don't use text mode in this demo.  This is just here to show how it would go.
    template<typename T>
    struct TextSender<gdt::vec_t<T, 3>> {
        static void send(std::ostream& stream, const gdt::vec_t<T, 3>& v) {
            TextSender<T>::send(stream, v.x);
            stream << " ";
            TextSender<T>::send(stream, v.y);
            stream << " ";
            TextSender<T>::send(stream, v.z);
        }
    };
} // namespace gnuplotio

// ---------- CHECK IF A CORRECT PATH TO GNUPLOT HAS BEN PROVIDED --------------

bool gnuplot_found() {
    if (Path::exists(Path::gnuplotExe())) {
        return true;
    }
    else {
        Console::print(OutLevel::ERR, "Could not find gnuplot executable at " + Path::gnuplotExe() + ". Plots can't be created.");
        return false;
    }
}

// ----------------------------- OUTPUT TO CONSOLE ----------------------------

void outputSuccess(const std::string& filename)
{
    Console::print(OutLevel::SUCCESS, "Plot rendered, and saved to " + filename  + " ... done.");
}


// ------------------------------- CLOUD POINTS -------------------------------


void PlotsGrapher::scatter2D(const std::string& filename, const std::vector<scal>& X, const std::vector<scal>& Y)
{
    if (gnuplot_found()) {
        {
            Gnuplot gp(Path::gnuplotExe());
            gp << "set terminal png truecolor enhanced size 1000,1000\n";
            gp << "set output \"" << filename << "\"\n";
            gp << "set xtics font \"Verdana,15\"\n";
            gp << "set xlabel font \"Verdana,15\"\n";
            gp << "set ytics font \"Verdana,15\"\n";
            gp << "set ylabel font \"Verdana,15\"\n";
            //gp << "set size ratio -1\n";
            gp << "set yrange [-1.4:1.4]\n";
            gp << "set pointsize 0.5\n";
            gp << "set xlabel \"x\"\n";
            gp << "set ylabel \"h\"\n";
            gp << "plot '-' with points pt 7 lc rgb \"#fd000000\" notitle\n";
            gp.send1d(std::tuple(X, Y));
        }
        outputSuccess(filename);
    }

}

void PlotsGrapher::scatter3D(const std::string& filename, const std::vector<scal>& theta, const std::vector<scal>& phi)
{
    if (gnuplot_found()) {
        {
            Gnuplot gp(Path::gnuplotExe());
            gp << "set terminal png truecolor enhanced size 1000,500\n";
            gp << "set output \"" << filename << "\"\n";
            gp << "set mapping spherical\n";
            gp << "set xtics font \"Verdana,15\"\n";
            gp << "set xlabel font \"Verdana,15\"\n";
            gp << "set ytics font \"Verdana,15\"\n";
            gp << "set ylabel font \"Verdana,15\"\n";
            gp << "set pointsize 0.3\n";
            gp << "set xlabel \"x\"\n";
            gp << "set ylabel \"y\"\n";
            gp << "splot '-' with points pt 7 lc rgb \"#fc000000\" notitle\n";
            gp.send1d(std::tuple(theta, phi));
        }

        outputSuccess(filename);
    }
}


// -------------------------------- 2D CURVES ---------------------------------


void PlotsGrapher::plot(const std::string& filename,
    const std::vector<scal>& X,
    const std::vector<std::vector<scal>>& Y,
    const std::vector<std::string>& titles)
{
    if (gnuplot_found()) {
        scal yMin = std::numeric_limits<scal>::max();
        scal yMax = std::numeric_limits<scal>::min();
        for (const std::vector<scal>& v : Y) {
            std::pair<std::vector<scal>::const_iterator, std::vector<scal>::const_iterator> _minmax = std::minmax_element(v.begin(), v.end());
            yMin = std::min(yMin, *(_minmax.first));
            yMax = std::max(yMax, *(_minmax.second));
        }

        {
            Gnuplot gp(Path::gnuplotExe());
            gp << "set encoding utf8\n";
            gp << "set terminal png size 1500,1100 enhanced font 'Verdana,50'\n";
            //gp << "set terminal png size 1300,1100\n";
            gp << "set output \"" << filename << "\"\n";
            gp << "set xrange [" << *std::min_element(X.begin(), X.end()) << ":" << *std::max_element(X.begin(), X.end()) << "]\n";
            gp << "set xtic (-1.5, 0, 1.5)\n";
            gp << "set yrange [" << std::min((scal)0, yMin) << ":" << std::max((scal)1, yMax) << "]\n";
            gp << "set ytics (0, 1)\n";
            //gp << "set key outside\n";
            gp << "set key right top\n";
            gp << "plot";
            for (int i = 0; i < Y.size(); ++i) {
                if (i < titles.size()) {
                    /*if (i == Y.size() - 1)
                        gp << gp.file1d(std::tuple(X, Y[i])) << "with lines lw 5 lc \"red\" title \"" << titles[i] << "\"";
                    else*/
                    gp << gp.file1d(std::tuple(X, Y[i])) << "with lines lw 5 title \"" << titles[i] << "\"";
                }
                else
                    gp << gp.file1d(std::tuple(X, Y[i])) << " with lines lw 5 notitle ";
                if (i < Y.size() - 1)
                    gp << ", ";
            }
        }

        outputSuccess(filename);
    }
}

void PlotsGrapher::plotG1(const std::string& filename,
    const std::vector<scal>& X,
    const std::vector<std::vector<scal>>& Y,
    const std::vector<std::string>& titles)
{
    if (gnuplot_found()) {
        scal yMin = std::numeric_limits<scal>::max();
        scal yMax = std::numeric_limits<scal>::min();
        for (const std::vector<scal>& v : Y) {
            std::pair<std::vector<scal>::const_iterator, std::vector<scal>::const_iterator> _minmax = std::minmax_element(v.begin(), v.end());
            yMin = std::min(yMin, *(_minmax.first));
            yMax = std::max(yMax, *(_minmax.second));
        }

        {
            Gnuplot gp(Path::gnuplotExe());
            gp << "set encoding utf8\n";
            gp << "set terminal png size 1500,1100 enhanced font 'Verdana,50'\n";
            //gp << "set terminal png size 1300,1100\n";
            gp << "set output \"" << filename << "\"\n";
            gp << "set xrange [" << *std::min_element(X.begin(), X.end()) << ":" << *std::max_element(X.begin(), X.end()) << "]\n";
            gp << "set xtic (-1.5, 0, 1.5)\n";
            gp << "set yrange [" << std::min((scal)0, yMin) << ":" << std::max((scal)1, yMax) << "]\n";
            //gp << "set ytics (0, 1)\n";
            //gp << "set key outside\n";
            gp << "set key center bottom\n";
            gp << "plot";
            for (int i = 0; i < Y.size(); ++i) {
                if (i < titles.size()) {
                    /*if (i == Y.size() - 1)
                        gp << gp.file1d(std::tuple(X, Y[i])) << "with lines lw 5 lc \"red\" title \"" << titles[i] << "\"";
                    else*/
                    gp << gp.file1d(std::tuple(X, Y[i])) << "with lines lw 5 title \"" << titles[i] << "\"";
                }
                else
                    gp << gp.file1d(std::tuple(X, Y[i])) << " with lines lw 5 notitle ";
                if (i < Y.size() - 1)
                    gp << ", ";
            }
        }

        outputSuccess(filename);
    }
}

void PlotsGrapher::plotDistrib(const std::string& filename,
    const std::vector<scal>& X,
    const std::vector<scal>& Y1,
    const std::vector<scal>& Y2)
{
    if (gnuplot_found()) {
        scal yMin = std::numeric_limits<scal>::max();
        scal yMax = std::numeric_limits<scal>::min();
        std::pair<std::vector<scal>::const_iterator, std::vector<scal>::const_iterator> _minmax = std::minmax_element(Y1.begin(), Y1.end());
        yMin = std::min(yMin, *(_minmax.first));
        yMax = std::max(yMax, *(_minmax.second));
        if (Y2.size() > 0) {
            _minmax = std::minmax_element(Y2.begin(), Y2.end());
            yMin = std::min(yMin, *(_minmax.first));
            yMax = std::max(yMax, *(_minmax.second));
        }

        {
            Gnuplot gp(Path::gnuplotExe());
            gp << "set encoding utf8\n";
            gp << "set terminal png size 1300,1100 enhanced font 'Verdana,30'\n";
            gp << "set output \"" << filename << "\"\n";
            gp << "unset ytics\n";
            //gp << "set ytics nomirror tc lt 1\n";
            //gp << "set ylabel 'Surface' tc lt 1\n";
            if (Y2.size() > 0) {
                gp << "set y2tics nomirror tc lt 2\n";
                gp << "set y2label 'Theory' tc lt 2\n";
            }
            gp << "plot";
            gp << gp.file1d(std::tuple(X, Y1)) << "axis x1y1 with lines lt 1 lw 3 notitle ,";
            if (Y2.size() > 0) {
                gp << gp.file1d(std::tuple(X, Y2)) << "axis x1y2 with lines lt 2 lw 3 notitle\n";
            }
        }

        outputSuccess(filename);
    }
}


// ------------------------------ 3D SURFACES ---------------------------------


void PlotsGrapher::splot(const std::string& filename, const std::vector<std::vector<gdt::vec3sc>>& points)
{
    if (gnuplot_found()) {
        {
            Gnuplot gp(Path::gnuplotExe());
            gp << "set terminal png size 1000, 1000\n";
            gp << "set output \"" << filename << "\"\n";
            gp << "set hidden3d\n";
            gp << "splot" << gp.file2d(points) << "with lines lt rgb \"blue\" notitle\n";
        }

        outputSuccess(filename);
    }
}

void PlotsGrapher::splotDiff(const std::string& filename, const std::vector<std::vector<gdt::vec3sc>>& points)
{
    if (gnuplot_found()) {
        {
            Gnuplot gpTop(Path::gnuplotExe());
            gpTop << "set terminal png size 1000, 1000 enhanced font 'Verdana,30'\n";
            gpTop << "set output \"" << filename << "_top.png\"\n";
            gpTop << "set hidden3d\n";
            gpTop << "set view equal xyz\n";
            gpTop << "set view 0, 0\n";
            gpTop << "set xtics offset 0, 0, 0\n";
            gpTop << "set ytics offset 2, 0, 0\n";
            gpTop << "unset ztics\n";
            gpTop << "set xzeroaxis\n";
            gpTop << "set yzeroaxis\n";
            gpTop << "splot" << gpTop.file2d(points) << "with lines lt rgb \"red\" notitle\n";
        }

        outputSuccess(filename + "_top.png");

        {
            Gnuplot gpFront(Path::gnuplotExe());
            gpFront << "set terminal png size 1000, 1000 enhanced font 'Verdana,30'\n";
            gpFront << "set output \"" << filename << "_front.png\"\n";
            gpFront << "set hidden3d\n";
            gpFront << "set view equal xyz\n";
            gpFront << "set view 90, 0\n";
            gpFront << "set xtics offset 0,-1,0\n";
            gpFront << "unset ytics\n";
            gpFront << "set format z \"%.2f\"; set ztics 0,0.01,1\n";
            gpFront << "splot" << gpFront.file2d(points) << "with lines lt rgb \"red\" notitle\n";
        }

        outputSuccess(filename + "_front.png");

        {
            Gnuplot gpRight(Path::gnuplotExe());
            gpRight << "set terminal png size 1000, 1000 enhanced font 'Verdana,30'\n";
            gpRight << "set output \"" << filename << "_right.png\"\n";
            gpRight << "set hidden3d\n";
            gpRight << "set view equal xyz\n";
            gpRight << "set view 90, 90.1\n";
            gpRight << "set xrange [-1:1]\n";
            gpRight << "set xtics offset 0,-1,0\n";
            gpRight << "unset xtics\n";
            gpRight << "set format z \"%.2f\"; set ztics 0,0.01,1\n";
            gpRight << "splot" << gpRight.file2d(points) << "with lines lt rgb \"red\" notitle\n";
        }

        outputSuccess(filename + "_right.png");
    }
}


// ----------------------------- PROJ COLORMAP --------------------------------


void PlotsGrapher::cmplot(const std::string& filename, const std::vector<gdt::vec3sc>& col)
{
    if (gnuplot_found()) {
        {
            Gnuplot gp(Path::gnuplotExe());
            gp << "set terminal pngcairo size 2000,2000\n";
            gp << "set output '" << filename << "'\n";
            gp << "set size square\n";

            gp << "set dgrid3d 200,200,2\n";
            gp << "set table $dat\n";
            gp << "splot" << gp.file1d(col) << "\n";
            gp << "unset table\n";
            gp << "unset dgrid3d\n";

            gp << "set pm3d map\n";
            gp << "set cbrange[0:1.2]\n";
            gp << "set palette model RGB\n";
            gp << "set palette model RGB defined(0 \"black\", 0.33 \"blue\", 0.66 \"dark-green\", 1 \"yellow\", 1 \"red\", 1.2 \"dark-red\")\n";
            gp << "circle(x, y, z) = x**2 + y**2 > 1 ? 0 : z\n";
            gp << "splot $dat u 1:2:(circle($1, $2, $3))\n";
        }

        outputSuccess(filename);
    }
}

void PlotsGrapher::cmplotDiff(const std::string& filename, const std::vector<gdt::vec3sc>& col, std::string palette, scal diffMax)
{
    if (gnuplot_found()) {
        {
            Gnuplot gp(Path::gnuplotExe());
            gp << "set loadpath \"Z:/gnuplot/palettes\"\n";
            //gp << "load \"" << palette << ".pal\"\n";
            gp << "set palette gray\n";

            gp << "set terminal pngcairo size 2100,2000 enhanced font 'Verdana,50'\n";
            gp << "set output '" << filename << "'\n";
            gp << "set size square\n";

            gp << "set xtics (-1, 0, 1)\n";
            gp << "set ytics (-1, 0, 1)\n";

            if (diffMax > 0) {
                gp << "set cbrange[0:" << diffMax << "]\n";
                gp << "set cbtics (0," << diffMax << ")\n";
            }

            gp << "set dgrid3d 500,500,2\n";
            //gp << "set dgrid3d 200,200,2\n";
            gp << "set table $dat\n";
            gp << "splot" << gp.file1d(col) << "\n";
            gp << "unset table\n";
            gp << "unset dgrid3d\n";

            gp << "set pm3d map\n";
            gp << "circle(x, y, z) = x**2 + y**2 > 1 ? 0 : z\n";
            gp << "splot $dat u 1:2:(circle($1, $2, $3))\n";
        }

        outputSuccess(filename);
    }
}
