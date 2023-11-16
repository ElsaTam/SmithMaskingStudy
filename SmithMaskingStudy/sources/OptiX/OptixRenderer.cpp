#include "OptiX/OptixRenderer.h"
#include "utils/console.h"
#include "utils/math/math.h"
#include "utils/params.h"

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
    Console::out << Console::timeStamp
                 << "Building OptixRenderer..." << std::endl;

	const char* ptx_path = Parameters::userParams.pathParams.ptxFile.c_str();

	// Initialize optix and create the context (static)
	context = Context();
    Console::light << Console::timePad << "Context created" << std::endl;

	// CUDA stream that optix pipeline will run on
	CUDA_CHECK(StreamCreate(&stream));

	// Build the pipeline (create the module and program groups)
	pipeline = new Pipeline(ptx_path);
    Console::light << Console::timePad << "Pipeline created" << std::endl;

	// Build the SBT (raygen and miss records)
	sbt = new SBT(pipeline);
    Console::light << Console::timePad << "SBT created" << std::endl;

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
    Console::out << Console::timeStamp
        << "Setting the geometry in the renderer..." << std::endl;

    mesh = _mesh;

    if (as) delete as;

	// Build the acceleration structure (Bottom-Level Acceleration Structure)
	as = new AccelerationStructure(mesh);
    Console::light << Console::timePad << "Acceleration Structure created" << std::endl;

	// Update the SBT (hitgroup records)
	sbt->updateAS(*as);
    Console::light << Console::timePad << "SBT updated with AS" << std::endl;

	// Initialize the launch parameters
	launchParams.traversable = as->getHandle();
	launchParamsBuffer.alloc(sizeof(launchParams));

    if (checkReady()) {
        Console::succ << Console::timePad << "SUCCESS: [OptixRenderer] Fully built and ready to render!" << std::endl;
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
    launchParams.visibility.tMax = Parameters::userParams.renderingParams.tMax;
    launchParams.visibility.useSmooth = Parameters::userParams.renderingParams.useSmooth;

    // Side effect
    if (Parameters::userParams.method != Method::COMPARE_EDB) {
        launchParams.sideEffect.borderPercentage = Parameters::userParams.sideEffectParams.borderPercentage;
        launchParams.sideEffect.directional = Parameters::userParams.sideEffectParams.directional;
    }

    // Program type
    switch (Parameters::userParams.method) {
    case Method::G1:
    case Method::COMPARE_EDB:
    case Method::COMPARE_SAMPLES:
    case Method::ASHIKHMIN_DIFF_3D:
    case Method::GENERATE_MICROFLAKES:
    case Method::STATISTICS:
    case Method::FULL_PIPELINE:
    case Method::BENCHMARK:
        launchParams.camera.programType = ProgramType::G1;
        break;
    case Method::GAF_2D:
    case Method::GAF_3D:
        launchParams.camera.programType = ProgramType::GAF;
        break;
    case Method::AMBIENT_OCCLUSION:
        launchParams.camera.programType = ProgramType::AMBIENT_OCCLUSION;
        break;
    case Method::TABULATION:
        if (Parameters::userParams.outLevel >= OutLevel::TRACE) {
            //Console::out << "userParams.method == Method::TABULATION. The programType is set to ";
            switch (launchParams.camera.programType) {
            case ProgramType::G1:
                //Console::out << "ProgramType::G1" << std::endl;
                break;
            case ProgramType::GAF:
                //Console::out << "ProgramType::GAF" << std::endl;
                break;
            case ProgramType::AMBIENT_OCCLUSION:
                //Console::out << "ProgramType::AMBIENT_OCCLUSION" << std::endl;
                break;
            default:
                //Console::warn << "unknown" << std::endl;
                break;
            }
        }
        break;
    default:
        break;
    }
}

void OptixRenderer::setCamera()
{
    if (!mesh) {
        Console::warn << "Setting the camera, but mesh is null. "
            << "Camera won't be centered." << std::endl;
        return;
    }
    vec3sc from = vec3sc( mesh->bounds.center().x, mesh->bounds.center().y, mesh->bounds.upper.z + 5. );
    vec3sc at{ mesh->bounds.center() };
    vec3sc up{ 0., 1., 0. };
    scal cosFovy = max(mesh->bounds.span().x, mesh->bounds.span().y);

    launchParams.camera.nPixelSamples = Parameters::userParams.renderingParams.nPixelSamples;
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
    if (Parameters::userParams.renderingParams.createPicture)
    {
        std::vector<uint32_t> pixels;
        const size_t npixels = launchParams.frame.size.x * launchParams.frame.size.y;
        pixels.resize(npixels);
        colorBuffer.download(pixels.data(), npixels);

        stbi_write_png(fileName.c_str(), launchParams.frame.size.x, launchParams.frame.size.y, 4, pixels.data(), launchParams.frame.size.x * sizeof(uint32_t));
        Console::succ << std::endl
            << "Image rendered, and saved to " << fileName << " ... done."
            << std::endl << std::endl;
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

    std::string fileName = "";
    if (Parameters::userParams.renderingParams.createPicture)
    {
        fileName = Parameters::userParams.pathParams.outputsFolder;
        fileName += mesh->name;
        fileName += "_phi_" + std::to_string(phi);
        fileName += "_theta_" + std::to_string(theta);
        fileName += "_size_" + std::to_string(launchParams.frame.size[0]); 
        fileName += "_samples_" + std::to_string(launchParams.camera.nPixelSamples);
        fileName += "_border_" + std::to_string(launchParams.sideEffect.borderPercentage);
        fileName += ".png";
    }
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
        else if (it->first == "directional") {
            launchParams.sideEffect.directional = std::get<bool>(it->second);
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
        else if (it->first == "nPixelSamples") {
            launchParams.camera.nPixelSamples = std::get<int>(it->second);
            Parameters::userParams.renderingParams.nPixelSamples = std::get<int>(it->second);
        }
    }
}



bool OptixRenderer::checkReady() const
{
    bool ready = true;
    if (!as) {
        Console::err << "ERROR: [OptixRenderer] Acceleration Structure is not built." << std::endl;
        ready = false;
    }
    if (!mesh) {
        Console::err << "ERROR: [OptixRenderer] No mesh given." << std::endl;
        ready = false;
    }
    if (!sbt) {
        Console::err << "ERROR: [OptixRenderer] SBT is not built." << std::endl;
        ready = false;
    }
    if (!pipeline) {
        Console::err << "ERROR: [OptixRenderer] Pipeline is not built." << std::endl;
        ready = false;
    }
    if (!Context::optixContext) {
        Console::err << "ERROR: [OptixRenderer] Context is not built." << std::endl;
        ready = false;
    }

    if (!ready) {
        exit(EXIT_FAILURE);
    }
    return ready;
}