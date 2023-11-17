#include "analyzer.h"
#include "shapes/TriangleMesh.h"
#include "tools/microflakesGenerator.h"
#include "tools/csvWriter.h"
#include "tools/logger.h"
#include "tools/objWriter.h"
#include "utils/console.h"
#include "utils/params.h"
#include "utils/paths.h"
#include <filesystem>

void help() {
    std::cout << "Usage: SurfAnalyzer [OPTION] ... [FILES]" << std::endl;
    std::cout << "--help           (-h)" << std::endl;
    std::cout << "        Print this help information on this tool." << std::endl;
    std::cout << "--example        (-e)" << std::endl;
    std::cout << "        Print the example of a .json file for the default parameters." << std::endl;
}

void example() {
    std::cout << "{" << std::endl;
    std::cout << "    \"userParams\": {" << std::endl;
    std::cout << "        \"method\": \"GLOBAL_VISIBILITY\"," << std::endl;
    std::cout << "        \"outLevel\" : \"TRACE\"," << std::endl;
    std::cout << "        \"log\" : false," << std::endl;
    std::cout << "        \"pathParams\": {" << std::endl;
    std::cout << "            \"objNames\": [\"PerTex/001.obj\", \"PerTex/002.obj\"]," << std::endl;
    std::cout << "            \"objMin\" : 3," << std::endl;
    std::cout << "            \"objMax\" : 10," << std::endl;
    std::cout << "            \"inputsFolder\" : \"../../inputs/\"," << std::endl;
    std::cout << "            \"outputsFolder\" : \"../../outputs/\"" << std::endl;
    std::cout << "        }," << std::endl;
    std::cout << "        \"directionParams\" : {" << std::endl;
    std::cout << "            \"phiStart\": 0," << std::endl;
    std::cout << "            \"phiEnd\" : 3.141593," << std::endl;
    std::cout << "            \"nAzimuthSamples\" : 1," << std::endl;
    std::cout << "            \"thetaStart\" : -1.570796," << std::endl;
    std::cout << "            \"thetaEnd\" : 1.570796," << std::endl;
    std::cout << "            \"nElevationSamples\" : 10" << std::endl;
    std::cout << "        }," << std::endl;
    std::cout << "        \"sideEffectParams\": {" << std::endl;
    std::cout << "            \"borderPercentage\": 0.5," << std::endl;
    std::cout << "            \"directional\" : true" << std::endl;
    std::cout << "        }," << std::endl;
    std::cout << "        \"renderingParams\" : {" << std::endl;
    std::cout << "            \"renderSize\": [1000, 1000] ," << std::endl;
    std::cout << "            \"nPixelSamples\" : 16," << std::endl;
    std::cout << "            \"createPicture\" : false" << std::endl;
    std::cout << "        }" << std::endl;
    std::cout << "    }" << std::endl;
    std::cout << "}" << std::endl;
}

void moveFiles(const std::string& src, const std::string& target)
{
    try
    {
        std::filesystem::copy(src, target, std::filesystem::copy_options::overwrite_existing | std::filesystem::copy_options::recursive);
    }
    catch (std::exception& e)
    {
        Console::err << "[ERROR] " << e.what() << std::endl;
    }
}

void cleanDirectory(const std::string& path)
{
    try
    {
        for (const auto& entry : std::filesystem::directory_iterator(path))
            std::filesystem::remove_all(entry.path());
    }
    catch (std::exception& e) {
        Console::err << "[ERROR] " << e.what() << std::endl;
    }
}

void run(const UserParams& params) {

    bool useGPU = (params.method == Method::G1)
        || (params.method == Method::GAF)
        //|| (params.method == Method::TABULATION)
        //|| (params.method == Method::STATISTICS)
        || (params.method == Method::AMBIENT_OCCLUSION)
        || (params.method == Method::GENERATE_MICROFLAKES);

    Analyzer analyzer(nullptr, useGPU);

    // for each surface
    for (const std::string& surfName : Path::surfaceNames())
    {
        // for each resolution
        for (int res : Path::resolutions())
        {
            Console::info << Console::line << Console::line;
            Console::info << surfName << ", " << res << " subdivisions" << std::endl;
            Console::info << Console::line << std::endl;

            try {
                TriangleMesh* mesh;
                if (params.method == Method::GENERATE_MICROFLAKES) {
                    LOG_NAME(Path::hfFolder(), surfName);
                    MicroflakesGenerator generator(Path::hfFile(surfName));
                    int size = pow(2, res);
                    mesh = generator.createModel({ size, size });
                }
                else {
                    LOG_NAME(Path::objFolder(), surfName);
                    mesh = createMesh(Path::objFile(surfName, res));
                }

                LOG_OBJ(mesh);

                analyzer.setGeo(mesh);

                switch (params.method) {
                case Method::G1:
                    analyzer.G1();
                    break;
                case Method::GAF:
                    analyzer.GAF();
                    break;
                case Method::D_TABULATION:
                    analyzer.tabulateDistrib();
                    break;
                case Method::AMBIENT_OCCLUSION:
                    analyzer.ambientOcclusion();
                    break;
                case Method::STATISTICS:
                    analyzer.statistics(true);
                    break;
                case Method::GENERATE_MICROFLAKES:
                {
                    ObjWriter::writeObj(Path::objFile(mesh->name, res), mesh);
                    //analyzer.normals();
                    //analyzer.statistics();
                    break;
                }
                case Method::FULL_PIPELINE:
                    analyzer.fullPipeline();
                    break;
                default:
                    break;
                }
                delete mesh;

                Console::succ << Console::timeStamp << "Done with " << surfName << std::endl;
            }
            catch (const std::exception& e) {
                if (params.outLevel >= OutLevel::ERR)
                    Console::err << e.what() << std::endl;
                LOG_ERROR(e.what());
            }
        }
    }
}

void createSubdvisionFolders(const std::string& root) {
    for (int i : Path::resolutions()) {
        std::string dir = root + Path::subdivisionFolderName(i);
        CreateDirectoryA(dir.c_str(), NULL);
    }
}

void createFolder(const std::string& path, bool withSubdivisions) {
    CreateDirectoryA(path.c_str(), NULL);
    if (withSubdivisions) createSubdvisionFolders(path);
}

bool createOutputFolders() {
    Console::out << "Create folders : " << Path::outputRootFolder() << std::endl;
    // check if root exists :
    if (! Path::exists(Path::outputRootFolder())) {
        Console::err << "Incorrect output path: " << Path::outputRootFolder() << std::endl;
        return false;
    }
    if (! Path::isFolder(Path::outputRootFolder())) {
        Console::err << "Output path is not a directory: " << Path::outputRootFolder() << std::endl;
        return false;
    }

    // logs
    createFolder(Path::logs_Folder(), false);

    // renders
    createFolder(Path::renders_Folder(), false);

    // G1
    createFolder(Path::G1_Folder(), false);
    createFolder(Path::G1_2D_Folder(), true);
    createFolder(Path::G1_3D_Folder(), true);

    // GAF
    createFolder(Path::GAF_Folder(), false);
    createFolder(Path::GAF_2D_Folder(), true);
    createFolder(Path::GAF_3D_Folder(), true);

    // tabulations
    createFolder(Path::tabulations_Folder(), true);

    // ambient occlusion
    createFolder(Path::ambientOcclusion_Folder(), true);

    // statistics
    createFolder(Path::statistics_Folder(), true);

    if (!Path::checkPaths()) {
        Path::checkPaths(true);
        Console::err << "There were issues while creating folders. Please correct the paths manually before restarting the program." << std::endl;
        return false;
    }

    return true;
}

int main(int argc, char* argv[]) {

#ifdef _DEBUG
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

    if (argc == 1) {
        Console::err << "No input file has been given." << std::endl;
        return 0;
    }

    srand(time(0));


    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
        if (files.empty() && (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0)) {
            help();
            return 0;
        }
        else if (files.empty() && (std::strcmp(argv[i], "-e") == 0 || std::strcmp(argv[i], "--example") == 0)) {
            example();
            return 0;
        }
        else {
            files.push_back(argv[i]);
        }
    }

    Parameters parameters(files);

    for (int launch = 0; launch < parameters.getNumberOfLaunchs(); ++launch)
    {
        const UserParams& userParams = parameters.getParamsForLaunch(launch);

        if (createOutputFolders()) {
            Logger::getInstance().setFolder(Path::logs_Folder());
            Logger::getInstance().enable(userParams.log);
            run(userParams);
            Console::out << std::endl << std::endl;
        }
    }

    return 0;
}