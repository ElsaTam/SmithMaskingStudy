#include "analyzer.h"
#include "shapes/TriangleMesh.h"
#include "tools/csvWriter.h"
#include "tools/logger.h"
#include "tools/objWriter.h"
#include "utils/console.h"
#include "utils/params.h"
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

    bool useGPU = (params.method == Method::G1_2D)
        || (params.method == Method::G1_3D)
        || (params.method == Method::GAF_2D)
        || (params.method == Method::GAF_3D)
        //|| (params.method == Method::TABULATION)
        //|| (params.method == Method::STATISTICS)
        || (params.method == Method::COMPARE_SMITH)
        || (params.method == Method::COMPARE_ASHIKHMIN)
        || (params.method == Method::COMPARE_EDB)
        || (params.method == Method::COMPARE_SAMPLES)
        || (params.method == Method::AMBIENT_OCCLUSION)
        || (params.method == Method::GENERATE_MICROFLAKES);

    Analyzer analyzer(nullptr, useGPU);

    // for each model
    for (const std::string& objName : params.pathParams.objNames)
    {
        Console::info << Console::line << Console::line;
        Console::info << objName << std::endl;
        LOG_NAME(params.pathParams.inputsFolder, objName);

        try {
            TriangleMesh* mesh = createMesh(objName);

            LOG_OBJ(mesh);

            analyzer.setGeo(mesh);

            switch (params.method) {
            case Method::G1_2D:
                analyzer.G1_2D();
                break;
            case Method::G1_3D:
                analyzer.G1_3D();
                break;
            case Method::ASHIKHMIN_DIFF_3D:
                analyzer.ashikhminDiff_3D();
                break;
            case Method::GAF_2D:
                analyzer.GAF_2D();
                break;
            case Method::GAF_3D:
                analyzer.GAF_3D();
                break;
            case Method::TABULATION:
            {
                bool D = true;
                bool G1_Ashikhmin = true;
                bool G1_RT = true;
                bool GAF = false;
                Discrete* NDF = G1_Ashikhmin ? new Discrete(*mesh, nullptr, 0) : nullptr;
                if (D || G1_Ashikhmin || G1_RT) analyzer.tabulate(D, G1_Ashikhmin, G1_RT, NDF);
                if (NDF) delete NDF;
                if (GAF) analyzer.tabulateGAF_RT();
                break;
            }
            case Method::HISTOGRAMS:
                analyzer.tabulateHeights();
                break;
            case Method::COMPARE_SMITH:
                analyzer.compareToSmith();
                break;
            case Method::COMPARE_ASHIKHMIN:
                analyzer.compareToAshikhmin();
                //analyzer.compareToAshikhminEquations();
                //analyzer.compareAshikhminToHeitz();
                break;
            case Method::COMPARE_EDB:
                //analyzer.compareEdB_border();
                analyzer.compareEdB_directional();
                break;
            case Method::COMPARE_SAMPLES:
                analyzer.compareRaycastingSamples();
                break;
            case Method::AMBIENT_OCCLUSION:
                analyzer.ambientOcclusion();
                break;
            case Method::NORMAL_DISTRIB:
                analyzer.normals();
                break;
            case Method::SLOPE_DISTRIB:
                analyzer.slopes();
                break;
            case Method::PLOT_DISTRIB:
                //analyzer.distribCuts();
                analyzer.distrib();
                //analyzer.compareDistribTmp();
                break;
            case Method::STATISTICS:
            {
                csv::CSVWriter* csvWriter = new csv::CSVWriter(params.pathParams.outputsFolder + "statistics/heights_and_thetas.csv", std::ios_base::app);
                analyzer.statistics(csvWriter);
                csvWriter->close();
                delete csvWriter;

                break;
            }
            case Method::GENERATE_MICROFLAKES:
                //analyzer.normals();
                //analyzer.statistics();
                analyzer.compareToAshikhmin();
                break;
            case Method::FULL_PIPELINE:
                analyzer.fullPipeline();
                break;
            case Method::BENCHMARK:
                //analyzer.benchmark();
                analyzer.benchmark_n_features();
                break;
            default:
                break;
            }
            delete mesh;

            Console::succ << Console::timeStamp << "Done with " << objName << std::endl;
        }
        catch (const std::exception& e) {
            if (params.outLevel >= OutLevel::ERR)
                Console::err << e.what() << std::endl;
            LOG_ERROR(e.what());
        }
    }
}

void createSubdvisionFolders(const std::string& root) {
    for (int i = 6; i < 13; ++i) {
        std::string dir = root + "subd" + std::to_string(i) + "/";
        CreateDirectoryA(dir.c_str(), NULL);
    }
}

bool createOutputFolders(const std::string& root) {
    Console::out << "Create folders : " << root << std::endl;
    // check if root exists :
    DWORD ftyp = GetFileAttributesA(root.c_str());
    if (ftyp == INVALID_FILE_ATTRIBUTES) {
        Console::err << "Incorrect output path: " << root << std::endl;
        return false;  //something is wrong with your path!
    }
    if (!(ftyp & FILE_ATTRIBUTE_DIRECTORY)) {
        Console::err << "Output path is not a directory: " << root << std::endl;
        return false;   // this is not a directory!
    }

    // logs
    std::string currentDir = root + "logs/";
    CreateDirectoryA(currentDir.c_str(), NULL);

    // visibility G1
    currentDir = root + "G1/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    currentDir = root + "G1/" + "2D/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);

    // G1
    currentDir = root + "G1/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    currentDir = root + "G1/" + "2D/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);
    currentDir = root + "G1/" + "3D/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);
    currentDir = root + "G1/" + "Smith/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);
    currentDir = root + "G1/" + "Ashikhmin/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);
    currentDir = root + "G1/" + "derivative/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);
    currentDir = root + "G1/" + "derivative2/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);

    // GAF
    currentDir = root + "GAF/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    currentDir = root + "GAF/" + "2D/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);
    currentDir = root + "GAF/" + "3D/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);

    // tabulations
    currentDir = root + "tabulations/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);

    // ambient occlusion
    currentDir = root + "ambient_occlusion/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);

    // statistics
    currentDir = root + "statistics/";
    CreateDirectoryA(currentDir.c_str(), NULL);

    // normals
    currentDir = root + "normals/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);

    // slopes
    currentDir = root + "slopes/";
    CreateDirectoryA(currentDir.c_str(), NULL);
    createSubdvisionFolders(currentDir);

    // distribs
    currentDir = root + "distribs/";
    CreateDirectoryA(currentDir.c_str(), NULL);

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

        if (createOutputFolders(userParams.pathParams.outputsFolder)) {
            Logger::getInstance().setFolder(userParams.pathParams.outputsFolder + "logs/");
            Logger::getInstance().enable(userParams.log || (userParams.method == Method::BENCHMARK));
            run(userParams);
            Console::out << std::endl << std::endl;
        }
    }

    return 0;
}