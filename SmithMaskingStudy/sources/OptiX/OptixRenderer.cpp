#include "OptiX/OptixRenderer.h"
#include "utils/console.h"
#include "utils/math/math.h"
#include "utils/params.h"
#include "utils/paths.h"

// this include may only appear in a single source file:
#include <optix_function_table_definition.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

OptixRenderer::OptixRenderer(const TriangleMesh* _mesh) :
    context(Context()),
    pipeline(nullptr),
    sbt(nullptr),
    as(nullptr),
    mesh(nullptr)
{
    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Building OptixRenderer...");

	const char* ptx_path = Path::ptxFile().c_str();

	// Initialize optix and create the context (static)
	context = Context();
    Console::print(OutLevel::TRACE, Console::timePad + "Context created");

	// CUDA stream that optix pipeline will run on
	CUDA_CHECK(StreamCreate(&stream));

	// Build the pipeline (create the module and program groups)
	pipeline = new Pipeline(ptx_path);
    Console::print(OutLevel::TRACE, Console::timePad + "Pipeline created");

	// Build the SBT (raygen and miss records)
	sbt = new SBT(pipeline);
    Console::print(OutLevel::TRACE, Console::timePad + "SBT created");

    if (_mesh)
        setGeo(_mesh);
}

OptixRenderer::~OptixRenderer()
{
    delete as;
    delete sbt;
    delete pipeline;
}

void OptixRenderer::setGeo(const TriangleMesh* _mesh)
{
    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Setting the geometry in the renderer...");

    mesh = _mesh;

    if (as) delete as;

	// Build the acceleration structure (Bottom-Level Acceleration Structure)
	as = new AccelerationStructure(mesh);
    Console::print(OutLevel::TRACE, Console::timePad + "Acceleration Structure created");

	// Update the SBT (hitgroup records)
	sbt->updateAS(*as);
    Console::print(OutLevel::TRACE, Console::timePad + "SBT updated with AS");

	// Initialize the launch parameters
	launchParams.traversable = as->getHandle();
	launchParamsBuffer.alloc(sizeof(launchParams));

    if (checkReady()) {
        Console::print(OutLevel::SUCCESS, "[OptixRenderer] Fully built and ready to render!");
    }
}

void OptixRenderer::resize(const vec2i& newSize)
{
    // resize our cuda frame buffer
    visibilityBuffer.resize(newSize.x * newSize.y * sizeof(scal));
    colorBuffer.resize(newSize.x * newSize.y * sizeof(uint32_t));

    // update the launch parameters that we'll pass to the optix launch:
    launchParams.frame.size = newSize;
    launchParams.frame.visibilityBuffer = (scal*)visibilityBuffer.d_pointer();
    launchParams.frame.colorBuffer = (uint32_t*)colorBuffer.d_pointer();

    // and re-set the camera, since aspect may have changed
    setCamera();
}

void OptixRenderer::setUserParams()
{
    // Visibility
    launchParams.visibility.useSmooth = Parameters::get()->currentParams()->renderingParams.useSmooth;

    // Program type
    switch (Parameters::get()->currentParams()->methodParams.method) {
    case Method::G1:
    case Method::GENERATE_MICROFLAKES:
    case Method::FEATURES:
    case Method::FULL_PIPELINE:
        launchParams.camera.programType = ProgramType::G1;
        break;
    case Method::GAF:
        launchParams.camera.programType = ProgramType::GAF;
        break;
    case Method::AMBIENT_OCCLUSION:
        launchParams.camera.programType = ProgramType::AMBIENT_OCCLUSION;
        break;
    default:
        break;
    }
}

void OptixRenderer::setCamera()
{
    if (!mesh) {
        Console::print(OutLevel::WARNING, Console::timePad + "Setting the camera, but mesh is null. Camera won't be centered.");
        return;
    }
    vec3sc from = vec3sc( mesh->bounds.center().x, mesh->bounds.center().y, mesh->bounds.upper.z + 5. );
    vec3sc at{ mesh->bounds.center() };
    vec3sc up{ 0., 1., 0. };
    scal cosFovy = max(mesh->bounds.span().x, mesh->bounds.span().y);

    launchParams.camera.nPixelSamples = Parameters::get()->currentParams()->renderingParams.nPixelSamples;
    launchParams.camera.position = from;
    launchParams.camera.direction = normalize(at - from);
    launchParams.camera.cosFovy = cosFovy;
    const float aspect = launchParams.frame.size.x / float(launchParams.frame.size.y);
    launchParams.camera.horizontal
        = (scal)(cosFovy * aspect) * normalize(cross(launchParams.camera.direction, up));
    launchParams.camera.vertical
        = (scal)(cosFovy) * normalize(cross(launchParams.camera.horizontal, launchParams.camera.direction));
}

void OptixRenderer::render(const std::string& fileName)
{
    if (!checkReady()) return;

    setUserParams();
    setCamera();
    launchParamsBuffer.upload(&launchParams, 1);
    OPTIX_CHECK(optixLaunch(
        /*! pipeline we're launching launch: */
        pipeline->getPipeline(),
        stream,
        /*! parameters and SBT */
        launchParamsBuffer.d_pointer(),
        launchParamsBuffer.sizeInBytes,
        &sbt->getSBT(),
        /*! dimensions of the launch: */
        launchParams.frame.size.x,
        launchParams.frame.size.y,
        1
    ));
    CUDA_SYNC_CHECK();

    // create PNG
    if (Parameters::get()->currentParams()->renderingParams.createPicture)
    {
        std::vector<uint32_t> pixels;
        const size_t npixels = launchParams.frame.size.x * launchParams.frame.size.y;
        pixels.resize(npixels);
        colorBuffer.download(pixels.data(), npixels);

        stbi_write_png(fileName.c_str(), launchParams.frame.size.x, launchParams.frame.size.y, 4, pixels.data(), launchParams.frame.size.x * sizeof(uint32_t));
        Console::print(OutLevel::SUCCESS, Console::timePad + "Image rendered, and saved to " + fileName);
    }
}

void OptixRenderer::render(scal phi, scal theta)
{
    if (!checkReady()) return;

    // Set outgoing direction
    vec3sc dir(-sin(theta) * cos(phi), -sin(theta) * sin(phi), cos(theta));
    dir = normalize(dir);
    dir = Geometry::rotateAlongNormal(dir, mesh->meso_normal);
    launchParams.visibility.directionOut = dir;

    std::string fileName = Parameters::get()->currentParams()->renderingParams.createPicture ? Path::renderImg(mesh->name, launchParams, phi, theta) : "";
    render(fileName);
}

void OptixRenderer::render(scal phiIn, scal thetaIn, scal phiOut, scal thetaOut)
{
    if (!checkReady()) return;

    // Set incoming direction
    vec3sc dirIn(-sin(thetaIn) * cos(phiIn), -sin(thetaIn) * sin(phiIn), cos(thetaIn));
    dirIn = normalize(dirIn);
    dirIn = Geometry::rotateAlongNormal(dirIn, mesh->meso_normal);
    launchParams.visibility.directionIn = normalize(dirIn);

    render(phiOut, thetaOut);
}

scal OptixRenderer::G1(scal phi, scal theta)
{
    render(phi, theta);
    return downloadVisibility();
}

scal OptixRenderer::GAF(scal phiIn, scal thetaIn, scal phiOut, scal thetaOut)
{
    render(phiIn, thetaIn, phiOut, thetaOut);
    return downloadVisibility();
}

scal OptixRenderer::downloadVisibility()
{
    if (!checkReady()) return 0.0;

    std::vector<scal> visibilities;
    visibilities.resize(launchParams.frame.size.x * launchParams.frame.size.y);
    visibilityBuffer.download(visibilities.data(),
        launchParams.frame.size.x * launchParams.frame.size.y);
    int concernedPixels = 0;
    scal visibility = 0.f;
    for (scal v : visibilities) {
        if (v >= 0.f) {
            concernedPixels++;
            visibility += v;
        }
    }
    if (concernedPixels > 0)
        visibility /= concernedPixels;
    return visibility;
}

void OptixRenderer::setLaunchParams(std::map<std::string, optionalParam> params)
{
    for (auto it = params.begin(); it != params.end(); ++it)
    {
        if (it->first == "borderPercentage") {
            launchParams.sideEffect.borderPercentage = std::get<scal>(it->second);
        }
        else if (it->first == "BBox") {
            launchParams.sideEffect.BBox = std::get<bool>(it->second);
        }
        else if (it->first == "targetedTriangles") {
            if (launchParams.visibility.numberOfTargetedTriangles > 0) {
                targetedTrianglesBuffer.freeBuffer();
            }
            std::vector<size_t> triangleIDs = std::get<std::vector<size_t>>(it->second);
            launchParams.visibility.numberOfTargetedTriangles = triangleIDs.size();
            if (launchParams.visibility.numberOfTargetedTriangles > 0) {
                targetedTrianglesBuffer.alloc_and_upload(triangleIDs);
                launchParams.visibility.targetedTriangles = (size_t*)targetedTrianglesBuffer.d_pointer();
            }
        }
        else if (it->first == "setsColors") {
            if (launchParams.sets.numberSets > 0) {
                setsColorBuffer.freeBuffer();
            }
            std::vector<vec3sc> colors = std::get<std::vector<vec3sc>>(it->second);
            launchParams.sets.numberSets = colors.size();
            if (launchParams.sets.numberSets > 0) {
                setsColorBuffer.alloc_and_upload(colors);
                launchParams.sets.colorsBySet = (vec3sc*)setsColorBuffer.d_pointer();
            }
        }
        else if (it->first == "colorsId") {
            if (colorsIdBuffer.d_ptr != nullptr) {
                colorsIdBuffer.freeBuffer();
            }
            std::vector<int> colorsId = std::get<std::vector<int>>(it->second);
            if (colorsId.size() > 0) {
                colorsIdBuffer.alloc_and_upload(colorsId);
                launchParams.sets.colorsId = (int*)colorsIdBuffer.d_pointer();
            }
        }
        else if (it->first == "useSmooth") {
            launchParams.visibility.useSmooth = std::get<bool>(it->second);
        }
        else if (it->first == "programType") {
            launchParams.camera.programType = std::get<ProgramType>(it->second);
        }
    }
}



bool OptixRenderer::checkReady() const
{
    bool ready = true;
    if (!as) {
        Console::print(OutLevel::ERR, "[OptixRenderer] Acceleration Structure is not built.");
        ready = false;
    }
    if (!mesh) {
        Console::print(OutLevel::ERR, "[OptixRenderer] No mesh given.");
        ready = false;
    }
    if (!sbt) {
        Console::print(OutLevel::ERR, "[OptixRenderer] SBT is not built.");
        ready = false;
    }
    if (!pipeline) {
        Console::print(OutLevel::ERR, "[OptixRenderer] Pipeline is not built.");
        ready = false;
    }
    if (!Context::optixContext) {
        Console::print(OutLevel::ERR, "[OptixRenderer] Context is not built.");
        ready = false;
    }

    if (!ready) {
        exit(EXIT_FAILURE);
    }
    return ready;
}