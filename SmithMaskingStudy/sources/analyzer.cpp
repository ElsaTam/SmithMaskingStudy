#include "analyzer.h"

#include <sys/stat.h>

#include "utils/console.h"
#include "utils/math/math.h"
#include "utils/math/random.h"
#include "utils/params.h"
#include "tools/logger.h"
#include "tools/plotsGrapher.h"
#include "HDFs/HeightsDiscrete.hpp"


#define ENTER                                                                 \
    if (Parameters::userParams.outLevel >= OutLevel::TRACE) { PING }          \
    const auto start = std::chrono::high_resolution_clock::now();

#define GPU_ENTER                                                             \
    ENTER                                                                     \
    prepareGPU();                                                             \
    logRenderingInfo();                                                       \

#define EXIT                                                                  \
    LOG_TIME(start, std::chrono::high_resolution_clock::now());



void INFO_PHI(scal phi, int i, size_t nPhi)
{
    if (Parameters::userParams.outLevel >= OutLevel::INFO) {
        Console::light << Console::timePad << Console::indent;
        Console::light << "(" << std::setfill(' ')
            << std::setw(std::to_string(nPhi).length()) << i << " / " << nPhi << ") ";
        Console::light << "phi = " << phi << std::endl;
    }
}


Analyzer::Analyzer(TriangleMesh* _mesh, bool useGPU)
{
    DIn = {};
    DIn.phiStart = Parameters::userParams.directionInParams.phiStart;
    DIn.phiRange = Parameters::userParams.directionInParams.phiEnd - DIn.phiStart;
    DIn.thetaStart = Parameters::userParams.directionInParams.thetaStart;
    DIn.thetaRange = Parameters::userParams.directionInParams.thetaEnd - DIn.thetaStart;
    DIn.nPhi = Parameters::userParams.directionInParams.nAzimuthSamples;
    DIn.nTheta = Parameters::userParams.directionInParams.nElevationSamples;

    DOut = {};
    DOut.phiStart = Parameters::userParams.directionOutParams.phiStart;
    DOut.phiRange = Parameters::userParams.directionOutParams.phiEnd - DOut.phiStart;
    DOut.thetaStart = Parameters::userParams.directionOutParams.thetaStart;
    DOut.thetaRange = Parameters::userParams.directionOutParams.thetaEnd - DOut.thetaStart;
    DOut.nPhi = Parameters::userParams.directionOutParams.nAzimuthSamples;
    DOut.nTheta = Parameters::userParams.directionOutParams.nElevationSamples;

    D = DOut;

    if (useGPU) {
        time_point tmp_start = std::chrono::high_resolution_clock::now();
        optixRenderer = new OptixRenderer(mesh);
        time_point tmp_end = std::chrono::high_resolution_clock::now();
        Duration tmp_duration = Duration(tmp_start, tmp_end);
        LOG_MESSAGE("OptixRenderer built in " + tmp_duration.str());
    }
    setGeo(_mesh);
}

Analyzer::~Analyzer()
{
    if (optixRenderer != nullptr) {
        delete optixRenderer;
    }
}


void Analyzer::setGeo(TriangleMesh* _mesh)
{
    mesh = _mesh;
    if (mesh && optixRenderer) {
        time_point tmp_start = std::chrono::high_resolution_clock::now();
        optixRenderer->setGeo(mesh);
        optixRenderer->resize(Parameters::userParams.renderingParams.renderSize);
        time_point tmp_end = std::chrono::high_resolution_clock::now();
        Duration tmp_duration = Duration(tmp_start, tmp_end);
        LOG_MESSAGE("Analyzer geo prepared in " + tmp_duration.str());
    }
    Console::out << std::endl;
}


// ------------------------------------------------------------------------- //
// ---------------------------      SETUP      ----------------------------- //
// ------------------------------------------------------------------------- //

void Analyzer::prepareGPU()
{
    if (optixRenderer == nullptr) {
        if (Parameters::userParams.outLevel >= OutLevel::TRACE) { PING }
        time_point tmp_start = std::chrono::high_resolution_clock::now();
        optixRenderer = new OptixRenderer(mesh);
        time_point tmp_end = std::chrono::high_resolution_clock::now();
        Duration tmp_duration = Duration(tmp_start, tmp_end);
        LOG_MESSAGE("OptixRenderer built in " + tmp_duration.str());
    }
    optixRenderer->resize(Parameters::userParams.renderingParams.renderSize);
}

// ------------------------------------------------------------------------- //
// ----------------------------      LOGS      ----------------------------- //
// ------------------------------------------------------------------------- //

void Analyzer::logRenderingInfo() const
{
    switch (Parameters::userParams.method)
    {
    case Method::G1:
    case Method::FULL_PIPELINE:
        LOG_MESSAGE("nPhi=" + std::to_string(Parameters::userParams.directionOutParams.nAzimuthSamples)
            + ", phiStart=" + std::to_string(Parameters::userParams.directionOutParams.phiStart)
            + ", phiEnd=" + std::to_string(Parameters::userParams.directionOutParams.phiEnd));
        LOG_MESSAGE("nTheta=" + std::to_string(Parameters::userParams.directionOutParams.nElevationSamples)
            + ", thetaStart=" + std::to_string(Parameters::userParams.directionOutParams.thetaStart)
            + ", thetaEnd=" + std::to_string(Parameters::userParams.directionOutParams.thetaEnd));
        break;

    case Method::GAF:
    case Method::TABULATION:
    case Method::STATISTICS:
        LOG_MESSAGE("IN_nPhi=" + std::to_string(Parameters::userParams.directionInParams.nAzimuthSamples)
            + ", IN_phiStart=" + std::to_string(Parameters::userParams.directionInParams.phiStart)
            + ", IN_phiEnd=" + std::to_string(Parameters::userParams.directionInParams.phiEnd));
        LOG_MESSAGE("IN_nTheta=" + std::to_string(Parameters::userParams.directionInParams.nElevationSamples)
            + ", IN_thetaStart=" + std::to_string(Parameters::userParams.directionInParams.thetaStart)
            + ", IN_thetaEnd=" + std::to_string(Parameters::userParams.directionInParams.thetaEnd));
        LOG_MESSAGE("OUT_nPhi=" + std::to_string(Parameters::userParams.directionOutParams.nAzimuthSamples)
            + ", OUT_phiStart=" + std::to_string(Parameters::userParams.directionOutParams.phiStart)
            + ", OUT_phiEnd=" + std::to_string(Parameters::userParams.directionOutParams.phiEnd));
        LOG_MESSAGE("OUT_nTheta=" + std::to_string(Parameters::userParams.directionOutParams.nElevationSamples)
            + ", OUT_thetaStart=" + std::to_string(Parameters::userParams.directionOutParams.thetaStart)
            + ", OUT_thetaEnd=" + std::to_string(Parameters::userParams.directionOutParams.thetaEnd));
        break;

    case Method::AMBIENT_OCCLUSION:
        break;

    default:
        break;
    }

    LOG_MESSAGE("NUM_PIXEL_SAMPLES=" + std::to_string(Parameters::userParams.renderingParams.nPixelSamples)
        + ", Frame=[" + std::to_string(Parameters::userParams.renderingParams.renderSize.x) + "x" + std::to_string(Parameters::userParams.renderingParams.renderSize.y) + "]"
        + ", N_SAMPLES=" + std::to_string(Parameters::userParams.renderingParams.nPixelSamples * Parameters::userParams.renderingParams.renderSize.x * Parameters::userParams.renderingParams.renderSize.y));
}

// ------------------------------------------------------------------------- //
// ----------------------------      UTILS      ---------------------------- //
// ------------------------------------------------------------------------- //

scal Analyzer::dSmooth(int N, int i) const
{
    if (N <= 1) return 0.5;
    float m = 1;
    float a = 0.5 / (atan(m * (scal)(N - 1) * 0.5));
    float y = a * atan(m * ((scal)i - (scal)(N - 1) * 0.5)) + 0.5;
    return y;
}

scal Analyzer::dLinear(int N, int i) const
{
    if (N <= 1) return 0.f;
    return (scal)i / (scal)(N - 1);
}

std::unique_ptr<MicrofacetDistribution> Analyzer::getTheoricalNDF() const
{
    vec2sc alpha = Parameters::userParams.ndfParams.alpha;
    MicrofacetProfil profil = Parameters::userParams.ndfParams.profil;
    if (Parameters::userParams.ndfParams.findAlphaFromFileName) {
        std::string fname = mesh->name;

        std::vector<std::string> tokens;
        std::istringstream f(mesh->name);
        std::string s;
        while (std::getline(f, s, '-')) {
            tokens.push_back(s);
        }

        if (tokens.size() != 3 && tokens.size() != 5) {
            Console::err << Console::timePad
                << "Name of the .obj file does not follow the format such as beckmann-0-50, no analytic NDF could be found. Beckmann with roughness=0.50 is used." << std::endl;
            alpha.x = alpha.y = 0.5f;
            profil = MicrofacetProfil::BECKMANN;
        }
        else {
            if (tokens[0].compare("beckmann") == 0) {
                profil = MicrofacetProfil::BECKMANN;
            }
            else if (tokens[0].compare("ggx") == 0) {
                profil = MicrofacetProfil::GGX;
            }

            alpha.x = std::stoi(tokens[1]) + std::pow(10, -(int)(tokens[2].length())) * std::stoi(tokens[2]);

            if (tokens.size() == 5)
                alpha.y = std::stoi(tokens[3]) + std::pow(10, -(int)(tokens[4].length())) * std::stoi(tokens[4]);
            else
                alpha.y = alpha.x;
        }
    }

    switch (profil) {
    case MicrofacetProfil::BECKMANN:
    {
        return std::make_unique<Beckmann>(alpha.x, alpha.y);
    }
    case MicrofacetProfil::GGX:
    {
        return std::make_unique<GGX>(alpha.x, alpha.y);
    }
    default:
        return nullptr;
    }
}

std::string Analyzer::getFolder(const std::string& root) const
{
    std::string folder = Parameters::userParams.pathParams.outputsFolder + root + "subd" + std::to_string(mesh->subdivisions()) + "/";
    struct stat sb;
    if (stat(folder.c_str(), &sb) == 0)
        return folder;
    else {
        Console::err << "Folder " << folder << " does not exist." << std::endl;
        exit(1);
    }
}


// ------------------------------------------------------------------------- //
// -----------------------------      G1      ------------------------------ //
// ------------------------------------------------------------------------- //


void Analyzer::G1()
{
    GPU_ENTER

    // Create the discrete NDF
    std::unique_ptr<Discrete> discrete_ndf(new Discrete(*mesh, nullptr, Parameters::userParams.sideEffectParams.borderPercentage));

    // Init csv writer to store results
    csv::CSVWriter* writer_rc    = new csv::CSVWriter(getFolder("tabulations/") + mesh->name + "_G1-rc.csv");
    csv::CSVWriter* writer_smith = new csv::CSVWriter(getFolder("tabulations/") + mesh->name + "_G1-smith.csv");
    // Write first row (theta)
    writeTheta(*writer_rc, false);
    writeTheta(*writer_smith, false);

    std::vector<std::vector<gdt::vec3sc>> pts_rc(D.nPhi), pts_smith(D.nPhi), pts_diff(D.nPhi);
    std::vector<gdt::vec3sc> col_rc, col_smith, col_diff;

    for (int i = 0; i < D.nPhi; ++i) {
        scal phi = D.phiStart + i * D.phiRange / (scal)D.nPhi; // [ -pi, ..., ..., ... ], pi
        INFO_PHI(phi, i, D.nPhi);

        // Vectors for the 2D plots
        Vec_Sc pts_theta, pts_G1_rc, pts_G1_smith;

        // One vector for each row
        std::vector<csv::elem> vector_g1_rc;
        std::vector<csv::elem> vector_g1_smith;
        // Push back the first value: phi (at the start of the cell)
        vector_g1_rc.push_back({ csv::elem::Tag::SCAL, phi });
        vector_g1_smith.push_back({ csv::elem::Tag::SCAL, phi });

        // Move phi at the center of the cell
        phi = D.phiStart + (i + 0.5f) * D.phiRange / (scal)D.nPhi; // [ -pi+d, ...+d, ...+d, ...+d ], pi+d

        for (int j = 0; j < D.nTheta; ++j) {
            scal theta = D.thetaStart + (j + 0.5f) * D.thetaRange / (scal)D.nTheta;
            // Compute raytracing (reference) and theoretical (Smith) masking
            scal G1_rc = optixRenderer->G1(phi, theta);
            scal G1_smith = discrete_ndf->G1(Conversion::polar_to_cartesian(theta, phi), { 0, 0, 1 });
            // for final 3D plot
            pts_rc[i].push_back(Conversion::polar_to_cartesian(theta, phi, G1_rc));
            col_rc.push_back(Conversion::polar_to_colormap(theta, phi, G1_rc));
            pts_smith[i].push_back(Conversion::polar_to_cartesian(theta, phi, G1_smith));
            col_smith.push_back(Conversion::polar_to_colormap(theta, phi, G1_smith));
            pts_diff[i].push_back(Conversion::polar_to_cartesian(theta, phi, G1_rc - G1_smith));
            col_diff.push_back(Conversion::polar_to_colormap(theta, phi, G1_rc - G1_smith));
            // for 2D plots
            pts_theta.push_back(theta);
            pts_G1_rc.push_back(G1_rc);
            pts_G1_smith.push_back(G1_smith);
            // for csv files
            vector_g1_smith.push_back({ csv::elem::Tag::SCAL, G1_smith });
            vector_g1_rc.push_back({ csv::elem::Tag::SCAL, G1_rc });
        }

        // Write rows
        writer_rc->writeRow(vector_g1_rc);
        writer_smith->writeRow(vector_g1_smith);

        // Plot 2D graph
        PlotsGrapher::plot(getFolder("G1/2D/") + mesh->name + "_" + std::to_string(phi) + ".png", pts_theta, { pts_G1_rc, pts_G1_smith }, { "Ray casting", "Smith"});
    }

    // Close csv writers
    writer_rc->close();
    writer_smith->close();
    delete writer_rc;
    delete writer_smith;

    // Plot 3D graph
    PlotsGrapher::splot     (getFolder("G1/3D/") + mesh->name + "_RC.png" ,                pts_rc);
    PlotsGrapher::cmplot    (getFolder("G1/3D/") + mesh->name + "_RC_colormap.png",        col_rc);
    PlotsGrapher::splot     (getFolder("G1/3D/") + mesh->name + "_Smith.png",              pts_smith);
    PlotsGrapher::cmplot    (getFolder("G1/3D/") + mesh->name + "_Smith_colormap.png",     col_smith);
    PlotsGrapher::splotDiff (getFolder("G1/3D/") + mesh->name + "_ashikhmin_diff",         pts_diff);
    PlotsGrapher::cmplotDiff(getFolder("G1/3D/") + mesh->name + "_ashikhmin_colormap.png", col_diff, "magma");

    EXIT
}


// ------------------------------------------------------------------------- //
// -----------------------------      GAF      ----------------------------- //
// ------------------------------------------------------------------------- //

void Analyzer::GAF()
{
    GPU_ENTER

    // Create the discrete NDF
    std::unique_ptr<Discrete> discrete_ndf(new Discrete(*mesh, nullptr, Parameters::userParams.sideEffectParams.borderPercentage));


    for (int i1 = 0; i1 < DIn.nPhi; ++i1) {
        scal phiIn = DIn.phiStart + i1 * DIn.phiRange / (scal)DIn.nPhi;

        for (int j1 = 0; j1 < DIn.nTheta; ++j1) {
            scal thetaIn = DIn.thetaStart + j1 * DIn.thetaRange / (scal)DIn.nTheta;

            // Init csv writer to store results
            csv::CSVWriter* writer_rc = new csv::CSVWriter(getFolder("tabulations/") + mesh->name + "_G1-rc.csv");
            csv::CSVWriter* writer_smith = new csv::CSVWriter(getFolder("tabulations/") + mesh->name + "_G1-smith.csv");
            // Write first row (theta)
            writeTheta(*writer_rc, false);
            writeTheta(*writer_smith, false);

            std::vector<std::vector<gdt::vec3sc>> pts_rc(D.nPhi), pts_smith(D.nPhi), pts_diff(D.nPhi);
            std::vector<gdt::vec3sc> col_rc, col_smith, col_diff;

            // first loop for the outgoing direction an ingoing direction is already set
            for (int i2 = 0; i2 < DOut.nPhi; ++i2) {
                scal phiOut = DOut.phiStart + i2 * DOut.phiRange / (scal)DOut.nPhi; // [ -pi, ..., ..., ... ], pi
                INFO_PHI(phiOut, i2, DOut.nPhi);

                // Vectors for the 2D plots
                Vec_Sc pts_theta, pts_G1_rc, pts_G1_smith;

                // One vector for each row
                std::vector<csv::elem> vector_g1_rc;
                std::vector<csv::elem> vector_g1_smith;
                // Push back the first value: phi (at the start of the cell)
                vector_g1_rc.push_back({ csv::elem::Tag::SCAL, phiOut });
                vector_g1_smith.push_back({ csv::elem::Tag::SCAL, phiOut });

                // Move phi at the center of the cell
                phiOut = DOut.phiStart + (i2 + 0.5f) * DOut.phiRange / (scal)DOut.nPhi; // [ -pi+d, ...+d, ...+d, ...+d ], pi+d

                // last loop : compute the gaf
                for (int j2 = 0; j2 < DOut.nTheta; ++j2) {
                    scal thetaOut = DOut.thetaStart + (j2 + 0.5f) * DOut.thetaRange / (scal)DOut.nTheta;
                    // Compute raytracing (reference) and theoretical (Smith) masking
                    scal GAF_rc = optixRenderer->GAF(phiIn, thetaIn, phiOut, thetaOut);
                    scal G1_i_smith = discrete_ndf->G1(Conversion::polar_to_cartesian(thetaIn, phiIn), { 0, 0, 1 });
                    scal G1_o_smith = discrete_ndf->G1(Conversion::polar_to_cartesian(thetaOut, phiOut), { 0, 0, 1 });
                    scal GAF_smith = G1_i_smith * G1_o_smith;
                    // for final 3D plot
                    pts_rc[i2].push_back(Conversion::polar_to_cartesian(thetaOut, phiOut, GAF_rc));
                    col_rc.push_back(Conversion::polar_to_colormap(thetaOut, phiOut, GAF_rc));
                    pts_smith[i2].push_back(Conversion::polar_to_cartesian(thetaOut, phiOut, GAF_smith));
                    col_smith.push_back(Conversion::polar_to_colormap(thetaOut, phiOut, GAF_smith));
                    pts_diff[i2].push_back(Conversion::polar_to_cartesian(thetaOut, phiOut, GAF_rc - GAF_smith));
                    col_diff.push_back(Conversion::polar_to_colormap(thetaOut, phiOut, GAF_rc - GAF_smith));
                    // for 2D plots
                    pts_theta.push_back(thetaOut);
                    pts_G1_rc.push_back(GAF_rc);
                    pts_G1_smith.push_back(GAF_smith);
                    // for csv files
                    vector_g1_smith.push_back({ csv::elem::Tag::SCAL, GAF_smith });
                    vector_g1_rc.push_back({ csv::elem::Tag::SCAL, GAF_rc });
                }

                // Write rows
                writer_rc->writeRow(vector_g1_rc);
                writer_smith->writeRow(vector_g1_smith);

                // plot a graph for an ingoing direction, an outgoing phi and N outgoing theta
                std::string folderDir = getFolder("GAF/2D/") + "phi_" + std::to_string(phiIn) + "_theta_" + std::to_string(thetaIn) + "/";
                CreateDirectoryA(folderDir.c_str(), NULL);
                PlotsGrapher::plot(folderDir + mesh->name + "_" + std::to_string(phiOut) + ".png", pts_theta, { pts_G1_rc, pts_G1_smith }, { });
            }

            // Close csv writers
            writer_rc->close();
            writer_smith->close();
            delete writer_rc;
            delete writer_smith;

            // Plot 3D graph
            std::string folderDir = getFolder("GAF/3D/") + "phi_" + std::to_string(phiIn) + "_theta_" + std::to_string(thetaIn) + "/";
            CreateDirectoryA(folderDir.c_str(), NULL);
            PlotsGrapher::splot     (folderDir + mesh->name + "_RC.png",                 pts_rc);
            PlotsGrapher::cmplot    (folderDir + mesh->name + "_RC_colormap.png",        col_rc);
            PlotsGrapher::splot     (folderDir + mesh->name + "_Smith.png",              pts_smith);
            PlotsGrapher::cmplot    (folderDir + mesh->name + "_Smith_colormap.png",     col_smith);
            PlotsGrapher::splotDiff (folderDir + mesh->name + "_ashikhmin_diff",         pts_diff);
            PlotsGrapher::cmplotDiff(folderDir + mesh->name + "_ashikhmin_colormap.png", col_diff, "magma");
        }
    }

    EXIT
}



// ------------------------------------------------------------------------- //
// -------------------------      TABULATION      -------------------------- //
// ------------------------------------------------------------------------- //

/*
                 THETA
           0                  (pi/2)
           ^                   ^
    - pi <  -------------------
           |   |   |   |   |   |
            -------------------
           |   |   |   |   |   |
  P         -------------------
  H        |   |   |   |   |   |
  I         -------------------
           |   |   |   |   |   |
            -------------------
           |   |   |   |   |   |
  (+ pi) <  -------------------

*/

void Analyzer::writeTheta(csv::CSVWriter& writer, bool fromDistrib) const
{
    scal thetaStart, thetaEnd;
    size_t nTheta;
    if (fromDistrib) {
        thetaStart = Discrete::thetaStart();
        thetaEnd = Discrete::thetaEnd();
        nTheta = Discrete::thetaSize();
    }
    else {
        thetaStart = D.thetaStart;
        thetaEnd = D.thetaStart + D.thetaRange;
        nTheta = D.nTheta;
    }
    const scal thetaStep = (thetaEnd - thetaStart) / (scal)nTheta;

    std::vector<csv::elem> thetas;
    thetas.push_back({ csv::elem::Tag::STRING, "" });
    for (int j = 0; j < nTheta; ++j) {
        scal theta = thetaStart + j * thetaStep;
        thetas.push_back({ csv::elem::Tag::SCAL, theta });
    }
    writer.writeRow(thetas);
}

scal forwarderRendererG1(void* renderer, scal phiIn, scal thetaIn, scal phiOut, scal thetaOut) {
    return static_cast<OptixRenderer*>(renderer)->G1(phiOut, thetaOut);
}
scal forwarderRendererGAF(void* renderer, scal phiIn, scal thetaIn, scal phiOut, scal thetaOut) {
    return static_cast<OptixRenderer*>(renderer)->GAF(phiIn, thetaIn, phiOut, thetaOut);
}
scal forwarderDistribG1(void* distrib, scal phiIn, scal thetaIn, scal phiOut, scal thetaOut) {
    vec3sc dir = Conversion::polar_to_cartesian(thetaOut, phiOut);
    return static_cast<Discrete*>(distrib)->G1(dir, {0, 0, 1});
}
scal forwarderDistribD(void* distrib, scal phiIn, scal thetaIn, scal phiOut, scal thetaOut) {
    vec3sc dir = Conversion::polar_to_cartesian(thetaOut, phiOut);
    return static_cast<Discrete*>(distrib)->D(dir);
}

void Analyzer::tabulateFunctions(std::vector<std::string> filenames, std::vector<scal(*)(void*, scal, scal, scal, scal)> T_functions, std::vector<void*> contexts, scal phiIn, scal thetaIn)
{
    size_t N_functions = T_functions.size();

    std::vector<csv::CSVWriter*> writers;
    for (int f = 0; f < N_functions; ++f) {
        writers.push_back(new csv::CSVWriter(filenames[f]));
        writeTheta(*writers[f]);
    }

    const scal phiStart = Discrete::phiStart(), phiEnd = Discrete::phiEnd();
    const scal thetaStart = Discrete::thetaStart(), thetaEnd = Discrete::thetaEnd();
    const size_t nPhi = Discrete::phiSize(), nTheta = Discrete::thetaSize();
    const scal phiStep = (phiEnd - phiStart) / (scal)nPhi;
    const scal thetaStep = (thetaEnd - thetaStart) / (scal)nTheta;

    for (int i = 0; i < nPhi; ++i) {
        scal phiOut = phiStart + i * phiStep; // [ -pi, ..., ..., ... ], pi
        INFO_PHI(phiOut, i, nPhi);

        std::vector<std::vector<csv::elem>> g1s(N_functions);
        for (int f = 0; f < N_functions; ++f) { g1s[f].push_back({ csv::elem::Tag::SCAL, phiOut }); }
        phiOut += 0.5 * phiStep; // [ -pi+d, ...+d, ...+d, ...+d ], pi+d
        for (int j = 0; j < nTheta; ++j) {
            scal thetaOut = thetaStart + j * thetaStep + 0.5 * thetaStep; // [0+d, ...+d, ...+d, ...+d ], pi/2+d
            for (int f = 0; f < N_functions; ++f)
                g1s[f].push_back({ csv::elem::Tag::SCAL, T_functions[f](contexts[f], phiIn, thetaIn, phiOut, thetaOut) });
        }
        for (int f = 0; f < N_functions; ++f) { writers[f]->writeRow(g1s[f]); }
    }
    for (int f = 0; f < N_functions; ++f) { delete writers[f]; }
}


void Analyzer::tabulate(bool D, bool G1_Ashikhmin, bool G1_RT, Discrete* NDF)
{
    ENTER
    if (G1_RT) {
        GPU_ENTER
        std::map<std::string, optionalParam> params;
        params["programType"] = ProgramType::G1;
        optixRenderer->setLaunchParams(params);
    }

    if (!NDF && (D || G1_Ashikhmin)) NDF = new Discrete(*mesh);

    std::vector<std::string> filenames;
    std::vector<scal(*)(void*, scal, scal, scal, scal)> functions;
    std::vector<void*> contexts;
    if (D) {
        filenames.push_back(getFolder("tabulations/") + mesh->name + "_D.csv");
        functions.push_back(&forwarderDistribD);
        contexts.push_back(NDF);
    }
    if (G1_Ashikhmin) {
        filenames.push_back(getFolder("tabulations/") + mesh->name + "_G1-ashikhmin.csv");
        functions.push_back(&forwarderDistribG1);
        contexts.push_back(NDF);
    }
    if (G1_RT) {
        filenames.push_back(getFolder("tabulations/") + mesh->name + "_G1-rt.csv");
        functions.push_back(&forwarderRendererG1);
        contexts.push_back(optixRenderer);
    }

    tabulateFunctions(filenames, functions, contexts);

    EXIT
}

void Analyzer::tabulateDistrib(Discrete* NDF)
{
    tabulate(true, false, false, NDF);
}

void Analyzer::tabulateG1_Ashikhmin(Discrete* NDF)
{
    tabulate(false, true, false, NDF);
}

void Analyzer::tabulateG1_RT()
{
    tabulate(false, false, true);
}

void Analyzer::tabulateGAF_RT()
{
    GPU_ENTER

    std::map<std::string, optionalParam> params;
    params["programType"] = ProgramType::GAF;
    optixRenderer->setLaunchParams(params);

    const scal phiStart = -m_pi, phiEnd = m_pi;   // phiEnd excluded
    const scal thetaStart = 0, thetaEnd = m_pi_2 - 0.001; // thetaEnd included
    const size_t nPhi = 16, nTheta = 30;
    const scal phiStep = (phiEnd - phiStart) / (scal)nPhi;
    const scal thetaStep = (thetaEnd - thetaStart) / ((scal)nTheta - 1);

    for (int i = 0; i < nPhi; ++i) {
        scal phiIn = phiStart + i * phiStep; // [ -pi, ..., ..., ... ], pi
        INFO_PHI(phiIn, i, nPhi);
        for (int j = 0; j < nTheta; ++j) {
            scal thetaIn = thetaStart + j * thetaStep; // [0, ..., ..., ... ], pi/2
            tabulateFunctions(
                { getFolder("tabulations/") + mesh->name + "_wi-" + std::to_string(thetaIn) + "-" + std::to_string(phiIn) + "_GAF.csv" },
                { &forwarderRendererGAF },
                { optixRenderer },
                phiIn, thetaIn);
        }
    }

    EXIT
}


void Analyzer::tabulateHeights()
{
    ENTER

    // Create the height NDF
    Console::out << Console::timeStamp << "Analyzer computing Histograms..." << std::endl;
    std::unique_ptr<HeightsDiscrete> hdf(new HeightsDiscrete(*mesh, Parameters::userParams.sideEffectParams.borderPercentage));

    // Create CSV writers
    csv::CSVWriter* writer = new csv::CSVWriter(getFolder("tabulations/") + mesh->name + "_H.csv");

    // Write first row (theta)
    writeTheta(*writer);

    // Prepare iteration values
    const scal phiStart = Discrete::phiStart(), phiEnd = Discrete::phiEnd();
    const scal thetaStart = Discrete::thetaStart(), thetaEnd = Discrete::thetaEnd();
    const size_t nPhi = Discrete::phiSize(), nTheta = Discrete::thetaSize();
    const scal phiStep = (phiEnd - phiStart) / (scal)nPhi;
    const scal thetaStep = (thetaEnd - thetaStart) / (scal)nTheta;

    // Start double loop

    for (int i = 0; i < nPhi; ++i) {
        scal phi = phiStart + i * phiStep; // [ -pi, ..., ..., ... ], pi

        // One vector for each row
        std::vector<csv::elem> vector_h;

        // Push back the first value: phi (at the start of the cell)
        vector_h.push_back({ csv::elem::Tag::SCAL, phi });

        // Move phi at the center of the cell
        phi += 0.5 * phiStep; // [ -pi+d, ...+d, ...+d, ...+d ], pi+d

        for (int j = 0; j < nTheta; ++j) {
            scal theta = thetaStart + j * thetaStep + 0.5 * thetaStep; // [0+d, ...+d, ...+d, ...+d ], pi/2+d
            vector_h.push_back({ csv::elem::Tag::SCAL_ARR, hdf->hist(Conversion::polar_to_cartesian(theta, phi)) });
        }

        // Write rows
        writer->writeRow(vector_h);

        if (Parameters::userParams.outLevel >= OutLevel::INFO) {
            Console::light << Console::timePad << Console::indent;
            Console::light << "(" << std::setfill(' ') << std::setw(std::to_string(nPhi).length()) << i << " / " << nPhi << ") ";
            Console::light << "     ###     phi = " << std::to_string(phi) << std::endl;
        }
    }

    writer->close();
    delete writer;

    EXIT
}


// ------------------------------------------------------------------------- //
// ---------------------      AMBIENT OCCLUSION      ----------------------- //
// ------------------------------------------------------------------------- //

void Analyzer::ambientOcclusion()
{
    GPU_ENTER

        optixRenderer->render(getFolder("ambient_occlusion/") + mesh->name + ".png");

    EXIT
}



// ------------------------------------------------------------------------- //
// ----------------------------      ERROR     ----------------------------- //
// ------------------------------------------------------------------------- //

scal Analyzer::error(scal theta, scal phi, const MicrofacetDistribution* NDF)
{
    scal measured = optixRenderer->G1(phi, theta);
    scal computed = NDF->G1(Conversion::polar_to_cartesian(theta, phi), { 0, 0, 1 });
    return abs(measured - computed);
}

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
using namespace boost::accumulators;

ErrorStats Analyzer::error()
{
    GPU_ENTER

    std::unique_ptr<Discrete> discrete_ndf(new Discrete(*mesh, nullptr, Parameters::userParams.sideEffectParams.borderPercentage));

    Console::out << Console::timeStamp << "Analyzer computing error..." << std::endl;

    ErrorStats E;
    E.integral = 0;

    scal dPhi = D.phiRange / (scal)(D.nPhi);
    scal dTheta = D.thetaRange / (scal)(D.nTheta);

    std::unique_ptr<csv::CSVWriter> writer(new csv::CSVWriter(getFolder("tabulations/") + mesh->name + "_Error.csv"));
    accumulator_set<scal, stats<tag::min, tag::max, tag::immediate_mean, tag::variance> > accError;

    std::vector<csv::elem> thetas;
    thetas.push_back({ csv::elem::Tag::STRING, "" });
    for (int j = 0; j < D.nTheta; ++j) {
        scal theta = D.thetaStart + j * dTheta;
        thetas.push_back({ csv::elem::Tag::SCAL, theta });
    }
    writer->writeRow(thetas);

    for (int i = 0; i < D.nPhi; ++i) {
        scal phi = D.phiStart + i * dPhi;

        std::vector<csv::elem> errors;
        errors.push_back({ csv::elem::Tag::SCAL, phi });

        for (int j = 0; j < D.nTheta; ++j) {
            scal theta = D.thetaStart + j * dTheta;
            scal e = error(theta, phi, discrete_ndf.get());
            errors.push_back({ csv::elem::Tag::SCAL, e });
            accError(e);
            E.integral += e * abs(sin(theta));
        }

        writer->writeRow(errors);

        if (Parameters::userParams.outLevel >= OutLevel::INFO) {
            Console::light << Console::timePad << Console::indent;
            Console::light << "(" << std::setfill(' ') << std::setw(std::to_string(D.nPhi).length()) << i << " / " << D.nPhi << ") ";
            Console::light << "Current error = " << std::to_string(E.integral);
            Console::light << "     ###     phi = " << std::to_string(phi) << std::endl;
        }
    }
    E.integral *= dPhi * dTheta;
    E.mean = extract::mean(accError);
    E.std = sqrt(extract::variance(accError));
    E.min = extract::min(accError);
    E.max = extract::max(accError);

    Console::prog << Console::timePad << "Error = " << E.integral << std::endl;

    EXIT
    return E;
}



// ------------------------------------------------------------------------- //
// -----------------------      DECOUPE EN SETS     ------------------------ //
// ------------------------------------------------------------------------- //

//input: ratio is between 0.0 to 1.0
//output: rgb color
vec3sc rgb(scal ratio)
{
    vec3sc first = { 1, 0, 0 }; // red
    vec3sc middle = { 1, 1, 0 }; // yellow
    vec3sc last = { 0, 1, 0 }; // green
    if (ratio < 0.5)
        return first + ratio * 2 * (middle - first);
    else if (ratio > 0.5)
        return middle + (ratio * 2 - 1) * (last - middle);
    else
        return middle;
}

void Analyzer::sets()
{
    GPU_ENTER

        for (int K = 2; K < 5; ++K) {
            Vec_Sc centroids(K);
            std::vector<std::set<int>> sets = mesh->heightSeparation(K, centroids); // each set contains the face's ids

            std::map<scal, int> mapCentroidToSet;
            for (int k = 0; k < K; ++k) {
                mapCentroidToSet[centroids[k]] = k;
            }
            std::vector<int> orderedSetsIds;
            for (std::map<scal, int>::iterator it = mapCentroidToSet.begin(); it != mapCentroidToSet.end(); ++it) {
                orderedSetsIds.push_back(it->second);
            }

            size_t flatFaces = 0;

            std::vector<int> colorsId(mesh->index.size());
            std::vector<vec3sc> colorBySet(K + 1);
            for (int i = 0; i < K; ++i) {
                int k = orderedSetsIds[i];
                colorBySet[k] = rgb((scal)i / (scal)(K - 1));
            }
            colorBySet[K] = { 0 };

            for (int faceID = 0; faceID < mesh->index.size(); ++faceID) {
                bool inSet = false;
                for (int i = 0; i < K && !inSet; ++i) {
                    int k = orderedSetsIds[i];
                    if (sets[k].find(faceID) != sets[k].end()) {
                        colorsId[faceID] = k;
                        flatFaces++;
                        inSet = true;
                    }
                }
                if (!inSet) { colorsId[faceID] = K; }
            }

            std::map<std::string, optionalParam> params;
            params["setsColors"] = colorBySet;
            params["colorsId"] = colorsId;
            optixRenderer->setLaunchParams(params);

            optixRenderer->render(Parameters::userParams.pathParams.outputsFolder + mesh->name + "_" + std::to_string(K) + "_clusters.png");
        }

    EXIT
}



// ------------------------------------------------------------------------- //
// ----------------------------      DATA      ----------------------------- //
// ------------------------------------------------------------------------- //


void Analyzer::distrib() const
{
    ENTER

    // Create distributions
    std::unique_ptr<MicrofacetDistribution> analytic_ndf = std::make_unique<GGX>(0.5, 0.5);

    Discrete mesh_1024_ndf(*mesh, nullptr, 0);

    TriangleMesh* mesh_1024CC = createMesh("12_subdivisions/VirtualGonio/ggx-0-50_1024_CatmullClark.obj");
    Discrete mesh_1024CC_ndf(*mesh_1024CC, nullptr, 0);
    delete mesh_1024CC;

    TriangleMesh* mesh_4096 = createMesh("12_subdivisions/VirtualGonio/ggx-0-50_4096.obj");
    Discrete mesh_4096_ndf(*mesh_4096, nullptr, 0);
    delete mesh_4096;

    const scal thetaStep = (mesh_1024_ndf.thetaEnd() - mesh_1024_ndf.thetaStart()) / ((scal)mesh_1024_ndf.thetaSize());
    const scal phiStep = (mesh_1024_ndf.phiEnd() - mesh_1024_ndf.phiStart()) / ((scal)mesh_1024_ndf.phiSize());

    bool D = true;
    std::vector<gdt::vec3sc> colAnalyticD, colMesh1024, colMesh1024CC, colMesh4096;
    scal analyticMaxD = 0, mesh1024MaxD = 0, mesh1024CCMaxD = 0, mesh4096MaxD = 0;

    for (int i = 0; i < mesh_1024_ndf.phiSize(); ++i) {
        for (int j = 0; j < mesh_1024_ndf.thetaSize(); ++j) {
            const scal phi = mesh_1024_ndf.phiStart() + phiStep * (i + 0.5);
            const scal theta = mesh_1024_ndf.thetaStart() + thetaStep * (j + 0.5);
            const vec3sc wn = Conversion::polar_to_cartesian(theta, phi);

            if (D) {
                scal D_analytic = analytic_ndf->D(wn);
                scal D_1024 = mesh_1024_ndf.D(wn);
                scal D_1024CC = mesh_1024CC_ndf.D(wn);
                scal D_4096 = mesh_4096_ndf.D(wn);

                colAnalyticD.push_back(Conversion::polar_to_colormap(theta, phi, D_analytic));
                if (analyticMaxD < D_analytic) analyticMaxD = D_analytic;

                colMesh1024.push_back(Conversion::polar_to_colormap(theta, phi, D_1024));
                if (mesh1024MaxD < D_1024) mesh1024MaxD = D_1024;

                colMesh1024CC.push_back(Conversion::polar_to_colormap(theta, phi, D_1024CC));
                if (mesh1024CCMaxD < D_1024CC) mesh1024CCMaxD = D_1024CC;

                colMesh4096.push_back(Conversion::polar_to_colormap(theta, phi, D_4096));
                if (mesh4096MaxD < D_4096) mesh4096MaxD = D_4096;
            }
        }
    }

    std::stringstream file;
    if (D) {
        scal maxValue = std::max(std::max(std::max(analyticMaxD, mesh1024MaxD), mesh1024CCMaxD), mesh4096MaxD);

        file.str(""); file.clear();
        file << Parameters::userParams.pathParams.outputsFolder << "distribs/D_" << mesh->name << "_analytic.png";
        PlotsGrapher::cmplotDiff(file.str(), colAnalyticD, "plasma", maxValue);

        file.str(""); file.clear();
        file << Parameters::userParams.pathParams.outputsFolder << "distribs/D_" << mesh->name << "_mesh_1024.png";
        //PlotsGrapher::cmplotDiff(file.str(), colMesh1024, "plasma", maxValue);

        file.str(""); file.clear();
        file << Parameters::userParams.pathParams.outputsFolder << "distribs/D_" << mesh->name << "_mesh_1024_CC.png";
        PlotsGrapher::cmplotDiff(file.str(), colMesh1024CC, "plasma", maxValue);

        file.str(""); file.clear();
        file << Parameters::userParams.pathParams.outputsFolder << "distribs/D_" << mesh->name << "_mesh_4096.png";
        PlotsGrapher::cmplotDiff(file.str(), colMesh4096, "plasma", maxValue);
    }

    EXIT
}

void Analyzer::statistics(csv::CSVWriter* writer)
{
    //ErrorStats E = error();
    ErrorStats E;

    ENTER

    Console::out << Console::timeStamp << "Computing statistics..." << std::endl;
    StatisticsTool stats(mesh, E);

    Console::out << std::endl;
    stats.print();
    Console::out << std::endl;

    if (writer != nullptr) {
        Console::out << Console::timePad << "Writing statistics..." << std::endl;
        if (writer->numberOfLines() == 0) stats.CSVHeader(writer);
        stats.toCSV(writer);

        Console::succ << Console::timePad
            << "Statistics append in " << writer->getFilename() << std::endl;
    }

    EXIT
}




// ------------------------------------------------------------------------- //
// ------------------------      FULL PIPELINE      ------------------------ //
// ------------------------------------------------------------------------- //

void Analyzer::fullPipeline()
{
    const bool render = true;
    // Prepare GPU
    GPU_ENTER
    if (render) {
        std::map<std::string, optionalParam> params;
        params["programType"] = ProgramType::G1;
        optixRenderer->setLaunchParams(params);
    }

    // Create the discrete NDF
    std::unique_ptr<Discrete> discrete_ndf(new Discrete(*mesh, nullptr, Parameters::userParams.sideEffectParams.borderPercentage));
    
    // Create CSV writers
    csv::CSVWriter* distributionWriter = new csv::CSVWriter(getFolder("tabulations/") + mesh->name + "_D.csv");
    csv::CSVWriter* G1AshikhminWriter = new csv::CSVWriter(getFolder("tabulations/") + mesh->name + "_G1-ashikhmin.csv");
    csv::CSVWriter* G1RTWriter = render ? new csv::CSVWriter(getFolder("tabulations/") + mesh->name + "_G1-rt.csv") : new csv::CSVWriter();
    csv::CSVWriter* errorWriter = render ? new csv::CSVWriter(getFolder("tabulations/") + mesh->name + "_Error.csv") : new csv::CSVWriter();
    csv::CSVWriter* statWriter = new csv::CSVWriter(Parameters::userParams.pathParams.outputsFolder + "statistics/data_HR.csv", std::ios_base::app);

    // Write first row (theta)
    writeTheta(*distributionWriter, true);
    writeTheta(*G1AshikhminWriter, false);
    writeTheta(*G1RTWriter, false);
    writeTheta(*errorWriter, false);

    Console::out << Console::timeStamp << "Analyzer computing G1s..." << std::endl;

    // Prepare error
    ErrorStats E;
    E.integral = 0;
    accumulator_set<scal, stats<tag::min, tag::max, tag::immediate_mean, tag::variance> > accError;

    // Prepare iteration values
    const scal phiStart = Discrete::phiStart(), phiEnd = Discrete::phiEnd();
    const scal thetaStart = Discrete::thetaStart(), thetaEnd = Discrete::thetaEnd();
    const size_t nPhi = Discrete::phiSize(), nTheta = Discrete::thetaSize();
    const scal phiStep = (phiEnd - phiStart) / (scal)nPhi;
    const scal thetaStep = (thetaEnd - thetaStart) / (scal)nTheta;

    // Start double loop

    for (int i = 0; i < nPhi; ++i) {
        scal phi = phiStart + i * phiStep; // [ -pi, ..., ..., ... ], pi

        // One vector for each row
        std::vector<csv::elem> vector_d;
        std::vector<csv::elem> vector_g1_ashikhmin;
        std::vector<csv::elem> vector_g1_rt;
        std::vector<csv::elem> vector_error;

        // Push back the first value: phi (at the start of the cell)
        vector_d.push_back({            csv::elem::Tag::SCAL, phi });
        vector_g1_ashikhmin.push_back({ csv::elem::Tag::SCAL, phi });
        vector_g1_rt.push_back({        csv::elem::Tag::SCAL, phi });
        vector_error.push_back({        csv::elem::Tag::SCAL, phi });
        
        // Move phi at the center of the cell
        phi += 0.5 * phiStep; // [ -pi+d, ...+d, ...+d, ...+d ], pi+d

        for (int j = 0; j < nTheta; ++j) {
            scal theta = thetaStart + j * thetaStep + 0.5 * thetaStep; // [0+d, ...+d, ...+d, ...+d ], pi/2+d

            vec3sc dir = Conversion::polar_to_cartesian(theta, phi);

            // Compute each value
            scal D            = discrete_ndf->D(dir);
            scal G1_Ashikhmin = discrete_ndf->G1(dir, { 0, 0, 1 });
            scal G1_RT        = render ? optixRenderer->G1(phi, theta) : 0;
            scal e            = render ? abs(G1_RT - G1_Ashikhmin) : 0;

            // Add each value to the row vector
            vector_d.push_back({            csv::elem::Tag::SCAL, D });
            vector_g1_ashikhmin.push_back({ csv::elem::Tag::SCAL, G1_Ashikhmin });
            vector_g1_rt.push_back({        csv::elem::Tag::SCAL, G1_RT });
            vector_error.push_back({        csv::elem::Tag::SCAL, e });

            // Increase error
            accError(e);
            E.integral += e * abs(sin(theta));
        }

        // Write rows
        distributionWriter->writeRow(vector_d);
        G1AshikhminWriter->writeRow(vector_g1_ashikhmin);
        G1RTWriter->writeRow(vector_g1_rt);
        errorWriter->writeRow(vector_error);

        if (Parameters::userParams.outLevel >= OutLevel::INFO) {
            Console::light << Console::timePad << Console::indent;
            Console::light << "(" << std::setfill(' ') << std::setw(std::to_string(nPhi).length()) << i << " / " << nPhi << ") ";
            if (render)
                Console::light << "Current error = " << std::to_string(E.integral);
            Console::light << "     ###     phi = " << std::to_string(phi) << std::endl;
        }
    }

    // End error statistics
    if (render) {
        E.integral *= phiStep * thetaStep;
        E.mean = extract::mean(accError);
        E.std = sqrt(extract::variance(accError));
        E.min = extract::min(accError);
        E.max = extract::max(accError);
        Console::prog << Console::timePad << "Error = " << E.integral << std::endl;
    }

    Console::out << Console::timeStamp << "Computing statistics..." << std::endl;
    StatisticsTool stats(mesh, E);

    Console::out << Console::timePad << "Writing statistics..." << std::endl;
    if (statWriter->numberOfLines() == 0) stats.CSVHeader(statWriter);
    stats.toCSV(statWriter);

    Console::succ << Console::timePad << "Statistics append in " << statWriter->getFilename() << std::endl;

    distributionWriter->close();
    G1AshikhminWriter->close();
    G1RTWriter->close();
    errorWriter->close();
    statWriter->close();

    delete distributionWriter;
    delete G1AshikhminWriter;
    delete G1RTWriter;
    delete errorWriter;
    delete statWriter;

    EXIT
}

