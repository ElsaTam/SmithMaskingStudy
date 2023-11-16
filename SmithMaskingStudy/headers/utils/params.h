#pragma once

#include "tools/jasonParser.h"
#include "gdt/math/vec.h"

// choice of method
enum class Method {
    G1_2D,
    G1_3D,
    ASHIKHMIN_DIFF_3D,
    COMPARE_SMITH,
    COMPARE_ASHIKHMIN,
    COMPARE_EDB,
    COMPARE_SAMPLES,
    GAF_2D,
    GAF_3D,
    TABULATION,
    HISTOGRAMS,
    AMBIENT_OCCLUSION,
    NORMAL_DISTRIB,
    SLOPE_DISTRIB,
    PLOT_DISTRIB,
    STATISTICS,
    GENERATE_MICROFLAKES,
    FULL_PIPELINE,
    BENCHMARK
};

enum class OutLevel {
    NO_OUTPUT = 0,
    ERR,
    WARNING,
    INFO,
    TRACE
};

enum class MicrofacetProfil {
    BECKMANN = 0,
    GGX
};

// path (input and output files and foldes) parameters
struct PathParams {
    bool resetOutput                    { false };
    std::string inputsFolder            { "Z:/obj/" };
    std::string outputsFolder           { "Z:/outputs/SurfAnalyzer/" };
    std::string gnuplotPath             { "" };
    std::vector<std::string> objNames   { "10 subdivisions/PerTex/001.obj" };
    std::string ptxFile                 { "./sources/cuda/devicePrograms.cu.ptx" };
};
bool operator==(const PathParams& lhs, const PathParams& rhs);
bool operator!=(const PathParams& lhs, const PathParams& rhs);

// direction parameters
struct DirectionParams {
    // azimuth
    scal phiStart                       { 0 };
    scal phiEnd                         { m_pi };
    int nAzimuthSamples                 { 1 };
    // elevation
    scal thetaStart                     { -m_pi_2+0.01f };
    scal thetaEnd                       { m_pi_2-0.01f };
    int nElevationSamples               { 50 };
};
bool operator==(const DirectionParams& lhs, const DirectionParams& rhs);
bool operator!=(const DirectionParams& lhs, const DirectionParams& rhs);

// side effect parameters
struct SideEffectParams {
    float borderPercentage              { (scal)0.2 }; // [0; 1]
    bool directional                    { false };
};
bool operator==(const SideEffectParams& lhs, const SideEffectParams& rhs);
bool operator!=(const SideEffectParams& lhs, const SideEffectParams& rhs);

// rendering parameters
struct RenderingParams {
    gdt::vec2i renderSize               { 1024, 1024 };
    int nPixelSamples                   { 8 };
    bool createPicture                  { false };
    scal tMax                           { (scal)std::numeric_limits<scal>::max() };
    bool useSmooth                      { false };
};
bool operator==(const RenderingParams& lhs, const RenderingParams& rhs);
bool operator!=(const RenderingParams& lhs, const RenderingParams& rhs);

// NDFs parameters
struct NDFParams {
    MicrofacetProfil profil             { MicrofacetProfil::BECKMANN };
    bool findAlphaFromFileName          { true };
    gdt::vec2sc alpha                   { gdt::vec2sc(0.5, 0.5) };
};
bool operator==(const NDFParams& lhs, const NDFParams& rhs);
bool operator!=(const NDFParams& lhs, const NDFParams& rhs);

struct BinaryFilesParams {
    bool readFromBinary                 { false };
    bool createBinary                   { false };
    bool overwriteBinary                { true };
};
bool operator==(const BinaryFilesParams& lhs, const BinaryFilesParams& rhs);
bool operator!=(const BinaryFilesParams& lhs, const BinaryFilesParams& rhs);

struct UserParams {
    Method method                       { Method::COMPARE_ASHIKHMIN };
    OutLevel outLevel                   { OutLevel::TRACE };
    bool log                            { false };
    PathParams pathParams               { };
    DirectionParams directionInParams   { };
    DirectionParams directionOutParams  { }; // this is the one used when we need only one direction
    SideEffectParams sideEffectParams   { };
    RenderingParams renderingParams     { };
    NDFParams ndfParams                 { };
    BinaryFilesParams binaryFilesParams { };
};
bool operator==(const UserParams& lhs, const UserParams& rhs);
bool operator!=(const UserParams& lhs, const UserParams& rhs);


class Parameters {
public:
    static UserParams userParams;

    Parameters();
    Parameters(const std::vector<std::string>& paths);

    void addDefault();
    size_t getNumberOfLaunchs() const;
    const UserParams& getParamsForLaunch(int i) const;

private:
    void parse(const std::string& path);
    UserParams createUserParams(jParser::jValue jValue) const;
    PathParams createPathParamsForOBJ(jParser::jValue jValue) const;
    PathParams createPathParamsForHF(jParser::jValue jValue) const;
    DirectionParams createDirectionParams(jParser::jValue jValue) const;
    SideEffectParams createSideEffectParams(jParser::jValue jValue) const;
    RenderingParams createRenderingParams(jParser::jValue jValue) const;
    NDFParams createNDFParams(jParser::jValue jValue) const;
    BinaryFilesParams createBinaryFilesParams(jParser::jValue jValue) const;

    std::vector<UserParams> userParamsVector;
};

// precision
constexpr auto Epsilon = 1e-3f;