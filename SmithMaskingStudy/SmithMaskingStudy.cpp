#include "analyzer.h"
#include "shapes/TriangleMesh.h"
#include "tools/microflakesGenerator.h"
#include "tools/csvWriter.h"
#include "tools/logger.h"
#include "tools/objWriter.h"
#include "utils/console.h"
#include "utils/params.h"
#include "utils/paths.h"
#include "magic_enum.hpp"
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
    std::cout << "\t'userParams': [" << std::endl;
    std::cout << "\t\t{" << std::endl;
    std::cout << "\t\t\t'methodParams': {" << std::endl;
    std::cout << "\t\t\t\t'name': 'G1'" << std::endl;
    std::cout << "\t\t\t}," << std::endl;
    std::cout << "\t\t\t'pathParams': {" << std::endl;
    std::cout << "\t\t\t\t'objDir'     : 'path/to/obj/directory/'," << std::endl;
    std::cout << "\t\t\t\t'outputsDir' : 'path/to/outputs/directory/'," << std::endl;
    std::cout << "\t\t\t\t'surfNames'  : ['subpath/to/surfA.obj', 'subpath/to/surfB.obj']," << std::endl;
    std::cout << "\t\t\t\t'surfMin'    : 1," << std::endl;
    std::cout << "\t\t\t\t'surfMax'    : 10," << std::endl;
    std::cout << "\t\t\t\t'resolutions': [8, 9, 10]," << std::endl;
    std::cout << "\t\t\t\t'ptxFile'    : './sources/cuda/devicePrograms.cu.ptx'," << std::endl;
    std::cout << "\t\t\t}," << std::endl;
    std::cout << "\t\t\t'directionParams': {" << std::endl;
    std::cout << "\t\t\t\t'phiStart'     : 0," << std::endl;
    std::cout << "\t\t\t\t'phiEnd'       : 6.283185," << std::endl;
    std::cout << "\t\t\t\t'nPhiSamples'  : 400," << std::endl;
    std::cout << "\t\t\t\t'thetaStart'   : 0," << std::endl;
    std::cout << "\t\t\t\t'thetaEnd'     : 1.570796," << std::endl;
    std::cout << "\t\t\t\t'nThetaSamples': 100" << std::endl;
    std::cout << "\t\t\t}," << std::endl;
    std::cout << "\t\t\t'sideEffectParams': {" << std::endl;
    std::cout << "\t\t\t\t'borderPercentage': 0.2," << std::endl;
    std::cout << "\t\t\t\t'BBox'            : false" << std::endl;
    std::cout << "\t\t\t}," << std::endl;
    std::cout << "\t\t\t'renderingParams': {" << std::endl;
    std::cout << "\t\t\t\t'renderSize'   : [1024, 1024]," << std::endl;
    std::cout << "\t\t\t\t'nPixelSamples': 4," << std::endl;
    std::cout << "\t\t\t\t'createPicture': false," << std::endl;;
    std::cout << "\t\t\t\t'useSmooth'    : false" << std::endl;
    std::cout << "\t\t\t}," << std::endl;
    std::cout << "\t\t\t'outLevel': 'INFO'," << std::endl;
    std::cout << "\t\t\t'log'     : false" << std::endl;
    std::cout << "\t\t}" << std::endl;
    std::cout << "\t]" << std::endl;
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

void run() {

    Console::info << Console::line << Console::line;
    Console::info << "Method " << std::string(magic_enum::enum_name(Parameters::get()->currentParams()->methodParams.method)) << std::endl;
    Console::info << Console::line << std::endl;

    bool useGPU = (Parameters::get()->currentParams()->methodParams.method == Method::G1)
        || (Parameters::get()->currentParams()->methodParams.method == Method::GAF)
        || (Parameters::get()->currentParams()->methodParams.method == Method::FEATURES && Parameters::get()->currentParams()->methodParams.computeError)
        || (Parameters::get()->currentParams()->methodParams.method == Method::FULL_PIPELINE)
        || (Parameters::get()->currentParams()->methodParams.method == Method::AMBIENT_OCCLUSION);

    Analyzer analyzer(nullptr, useGPU);

    // for each surface
    for (const std::string& surfName : Path::surfaceNames())
    {
        // for each resolution
        for (int res : Path::resolutions())
        {
            Console::out << Console::shortline;
            Console::out << surfName << ", " << res << " subdivisions" << std::endl;
            Console::out << Console::shortline << std::endl;

            try {
                TriangleMesh* mesh;
                if (Parameters::get()->currentParams()->methodParams.method == Method::GENERATE_MICROFLAKES) {
                    LOG_NAME(Path::hfFolder(), surfName);
                    const std::string& hfPath = Path::hfFile(surfName);
                    if (!Path::isFile(hfPath)) {
                        throw std::runtime_error("File not found: " + hfPath);
                    }
                    MicroflakesGenerator generator(hfPath);
                    int size = pow(2, res);
                    mesh = generator.createModel({ size, size });
                }
                else {
                    LOG_NAME(Path::objFolder(), surfName);
                    const std::string& objPath = Path::objFile(surfName, res);
                    if (!Path::isFile(objPath)) {
                        throw std::runtime_error("File not found: " + objPath);
                    }
                    mesh = createMesh(objPath);
                }

                LOG_OBJ(mesh);

                analyzer.setGeo(mesh);

                switch (Parameters::get()->currentParams()->methodParams.method) {
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
                case Method::FEATURES:
                    analyzer.features();
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
                if (Parameters::get()->currentParams()->outLevel >= OutLevel::ERR)
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
    createFolder(Path::features_Folder(), true);

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


    std::string file = "";
    for (int i = 1; i < argc; ++i) {
        if (file.empty() && (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0)) {
            help();
            return 0;
        }
        else if (file.empty() && (std::strcmp(argv[i], "-e") == 0 || std::strcmp(argv[i], "--example") == 0)) {
            example();
            return 0;
        }
        else if (file.empty()) {
            file = argv[i];
        }
        else {
            Console::err << "Unknown argument " << argv[i] << " (parameters file has already been provided: " << file << ")" << std::endl;
        }
    }

    Parameters* parameters = Parameters::get();
    parameters->parse(file);

    for (int launch = 0; launch < parameters->size(); ++launch)
    {
        parameters->setIndex(launch);

        if (createOutputFolders()) {
            Logger::getInstance().setFolder(Path::logs_Folder());
            Logger::getInstance().enable(Parameters::get()->currentParams()->log);
            run();
            Console::out << std::endl << std::endl;
        }
    }

    return 0;
}