#pragma once

#include "tools/jasonParser.h"
#include "gdt/math/vec.h"

// choice of method
enum class Method {
    G1,
    GAF,
    D_TABULATION,
    FEATURES,
    FULL_PIPELINE,
    AMBIENT_OCCLUSION,
    GENERATE_MICROFLAKES
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

struct MethodParams {
    Method method                       { Method::G1 };
    bool computeError                   { false };
    bool jiggleFlakes                   { false };
};
bool operator==(const MethodParams& lhs, const MethodParams& rhs);
bool operator!=(const MethodParams& lhs, const MethodParams& rhs);

// path (input and output files and foldes) parameters
struct PathParams {
    std::string objDir               { "Z:/obj/" };
    std::string hfDir                { "Z:/hf/" };
    std::string outputsDir           { "Z:/outputs/SurfAnalyzer/" };
    std::string gnuplotPath             { "gnuplot" };
    std::vector<std::string> surfNames  { };
    std::vector<int> resolutions        { };
    std::string ptxFile                 { "./sources/cuda/devicePrograms.cu.ptx" };
};
bool operator==(const PathParams& lhs, const PathParams& rhs);
bool operator!=(const PathParams& lhs, const PathParams& rhs);

// direction parameters
struct DirectionParams {
    // azimuth
    scal phiStart                       { -m_pi };
    scal phiEnd                         { m_pi };
    int nPhiSamples                     { 1 };
    // elevation
    scal thetaStart                     { 0 };
    scal thetaEnd                       { m_pi_2 };
    int nThetaSamples                   { 90 };
};
bool operator==(const DirectionParams& lhs, const DirectionParams& rhs);
bool operator!=(const DirectionParams& lhs, const DirectionParams& rhs);

// NDF parameters
struct NDFParams {
    int nPhi{ 400 };
    int nTheta{ 100 };
};
bool operator==(const NDFParams& lhs, const NDFParams& rhs);
bool operator!=(const NDFParams& lhs, const NDFParams& rhs);

// side effect parameters
struct SideEffectParams {
    float borderPercentage              { (scal)0.2 }; // [0; 1]
    bool BBox                           { false };
};
bool operator==(const SideEffectParams& lhs, const SideEffectParams& rhs);
bool operator!=(const SideEffectParams& lhs, const SideEffectParams& rhs);

// rendering parameters
struct RenderingParams {
    gdt::vec2i renderSize               { 1024, 1024 };
    int nPixelSamples                   { 4 };
    bool createPicture                  { false };
    bool useSmooth                      { false };
};
bool operator==(const RenderingParams& lhs, const RenderingParams& rhs);
bool operator!=(const RenderingParams& lhs, const RenderingParams& rhs);

struct UserParams {
    OutLevel outLevel                   { OutLevel::WARNING };
    bool log                            { false };
    MethodParams methodParams           { };
    PathParams pathParams               { };
    DirectionParams directionInParams   { };
    DirectionParams directionOutParams  { }; // this is the one used when we need only one direction
    NDFParams ndfParams                 { };
    SideEffectParams sideEffectParams   { };
    RenderingParams renderingParams     { };
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
    bool checkRequired(jParser::jValue jValue) const;
    UserParams createUserParams(jParser::jValue jValue) const;
    MethodParams createMethodParams(jParser::jValue jValue) const;
    PathParams createPathParamsForOBJ(jParser::jValue jValue) const;
    PathParams createPathParamsForHF(jParser::jValue jValue) const;
    DirectionParams createDirectionParams(jParser::jValue jValue) const;
    NDFParams createNDFParams(jParser::jValue jValue) const;
    SideEffectParams createSideEffectParams(jParser::jValue jValue) const;
    RenderingParams createRenderingParams(jParser::jValue jValue) const;

    std::vector<UserParams> userParamsVector;
};

// precision
constexpr auto Epsilon = 1e-3f;