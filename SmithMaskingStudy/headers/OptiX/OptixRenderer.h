#pragma once

#include "OptiX/SBT.h"
#include "OptiX/Context.h"
#include <variant>

using optionalParam = std::variant<
    scal,
    int,
    std::vector<size_t>,
    std::vector<vec3sc>,
    std::vector<int>,
    bool,
    ProgramType
>;

class OptixRenderer
{
private:
    const TriangleMesh* mesh;

    Context context;
    Pipeline* pipeline;
    SBT*	sbt;
    AccelerationStructure* as;
    CUstream stream;

    /*! @{ our launch parameters, on the host, and the buffer to store
        them on the device */
    LaunchParams launchParams;
    CUDABuffer   launchParamsBuffer;
    /*! @} */

    /*! @{ input buffers to pass info to the device */
    CUDABuffer targetedTrianglesBuffer;
    CUDABuffer setsColorBuffer; // one color per set
    CUDABuffer colorsIdBuffer; // one color id by face
    /*! @} */

    /*! @{ output buffers to store computed info from the device */
    CUDABuffer colorBuffer;
    CUDABuffer visibilityBuffer;
    /*! @} */

public:
	OptixRenderer(const TriangleMesh* mesh = nullptr);
    ~OptixRenderer();

    void setGeo(const TriangleMesh* mesh);

    void resize(const vec2i& newSize);
    void setUserParams();
    void setCamera();

    void render(const std::string& fileName = "");
    scal G1(scal phi, scal theta);
    scal GAF(scal phiIn, scal thetaIn, scal phiOut, scal thetaOut);

    void setLaunchParams(std::map<std::string, optionalParam> params);

protected:
    void render(scal phi, scal theta);
    void render(scal phiIn, scal thetaIn, scal phiOut, scal thetaOut);
    bool checkReady() const;
    scal downloadVisibility();
};