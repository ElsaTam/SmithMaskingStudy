#include "utils/params.h"

#include <iomanip>
#include <fstream>
#include <sstream>

#include "utils/paths.h"


bool operator==(const PathParams& lhs, const PathParams& rhs) {
	if (lhs.objFolder.compare(rhs.objFolder) != 0) return false;
	if (lhs.hfFolder.compare(rhs.hfFolder) != 0) return false;
	if (lhs.outputsFolder.compare(rhs.outputsFolder) != 0) return false;
	if (lhs.surfNames != rhs.surfNames) return false;
	return true;
}

bool operator!=(const PathParams& lhs, const PathParams& rhs) {
	return !(lhs == rhs);
}

bool operator==(const DirectionParams& lhs, const DirectionParams& rhs)
{
	if (lhs.phiStart != rhs.phiStart) return false;
	if (lhs.phiEnd != rhs.phiEnd) return false;
	if (lhs.nAzimuthSamples != rhs.nAzimuthSamples) return false;
	if (lhs.thetaStart != rhs.thetaStart) return false;
	if (lhs.thetaEnd != rhs.thetaEnd) return false;
	if (lhs.nElevationSamples != rhs.nElevationSamples) return false;
	return true;
}

bool operator!=(const DirectionParams& lhs, const DirectionParams& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const SideEffectParams& lhs, const SideEffectParams& rhs)
{
	if (lhs.borderPercentage != rhs.borderPercentage) return false;
	if (lhs.directional != rhs.directional) return false;
	return true;
}

bool operator!=(const SideEffectParams& lhs, const SideEffectParams& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const RenderingParams& lhs, const RenderingParams& rhs)
{
	if (lhs.renderSize != rhs.renderSize) return false;
	if (lhs.nPixelSamples != rhs.nPixelSamples) return false;
	if (lhs.createPicture != rhs.createPicture) return false;
	if (lhs.tMax != rhs.tMax) return false;
	if (lhs.useSmooth != rhs.useSmooth) return false;
	return true;
}

bool operator!=(const RenderingParams& lhs, const RenderingParams& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const UserParams& lhs, const UserParams& rhs)
{
	if (lhs.method != rhs.method) return false;
	if (lhs.outLevel != rhs.outLevel) return false;
	if (lhs.log != rhs.log) return false;
	if (lhs.pathParams != rhs.pathParams) return false;
	if (lhs.directionInParams != rhs.directionInParams) return false;
	if (lhs.directionOutParams != rhs.directionOutParams) return false;
	if (lhs.sideEffectParams != rhs.sideEffectParams) return false;
	if (lhs.renderingParams != rhs.renderingParams) return false;
	return true;
}

bool operator!=(const UserParams& lhs, const UserParams& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const NDFParams& lhs, const NDFParams& rhs)
{
	if (lhs.profil != rhs.profil) return false;
	if (lhs.alpha != rhs.alpha) return false;
	return true;
}

bool operator!=(const NDFParams& lhs, const NDFParams& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const BinaryFilesParams& lhs, const BinaryFilesParams& rhs)
{
	if (lhs.readFromBinary != rhs.readFromBinary) return false;
	if (lhs.createBinary != rhs.createBinary) return false;
	if (lhs.overwriteBinary != rhs.overwriteBinary) return false;
	return true;
}

bool operator!=(const BinaryFilesParams& lhs, const BinaryFilesParams& rhs)
{
	return !(lhs == rhs);
}



UserParams Parameters::userParams = UserParams();


Parameters::Parameters()
{
	userParamsVector.push_back(userParams);
}

Parameters::Parameters(const std::vector<std::string>& paths)
{
	userParamsVector = std::vector<UserParams>(0);
	for (const std::string& path : paths)
	{
		parse(path);
	}
	userParams = userParamsVector[0];
}

void Parameters::addDefault()
{
	userParamsVector.push_back(UserParams());
}



void Parameters::parse(const std::string& path)
{
	std::ifstream in(path.c_str());
	if (!in.is_open()) return;

	std::string str = "";
	std::string tmp;
	while (std::getline(in, tmp)) str += tmp; // we stock the whole file in a single string
	in.close();

	jParser::jValue allUserParams = jParser::parser::parse(str);
	// {
	//     "userParams" : {
	// 	       ...
	//     },
	//     "userParams" : {
	// 	       ...
	//     }
	// }
	for (int i = 0; i < allUserParams.size(); ++i)
	{
		userParamsVector.push_back(createUserParams(allUserParams[i]));
	}
}

// {
//     "method" : "GLOBAL_VISIBILITY",
//     "outLevel" : "TRACE",
//     "log" : false,
//     "pathParams": {
//      	"objNames": ["PerTex/001.obj"] ,
//     		"objMin" : 1,
//     		"objMax" : 10,
//     		"inputsFolder" : "../../inputs/",
//     		"outputsFolder" : "../../outputs/"
//     },
//     "directionParams" : {
//     	   ...
//     },
//     "sideEffectParams" : {
//         ...
//     },
//     "renderingParams" : {
//         ...
//     },
//	   "ndfParams" : {
//		   ...
//     }
// }
UserParams Parameters::createUserParams(jParser::jValue jValue) const
{
	// Method
	if (jValue.contains("method")) {
		const std::string methodStr = jValue["method"].as_string();
		if (methodStr.compare("GAF") == 0) userParams.method = Method::GAF;
		else if (methodStr.compare("G1") == 0) userParams.method = Method::G1;
		else if (methodStr.compare("D_TABULATION") == 0) userParams.method = Method::D_TABULATION;
		else if (methodStr.compare("AMBIENT_OCCLUSION") == 0) userParams.method = Method::AMBIENT_OCCLUSION;
		else if (methodStr.compare("STATISTICS") == 0) userParams.method = Method::STATISTICS;
		else if (methodStr.compare("GENERATE_MICROFLAKES") == 0) userParams.method = Method::GENERATE_MICROFLAKES;
		else if (methodStr.compare("FULL_PIPELINE") == 0) userParams.method = Method::FULL_PIPELINE;
	}

	// Out level
	if (jValue.contains("outLevel")) {
		const std::string outLevelStr = jValue["outLevel"].as_string();
		if (outLevelStr.compare("NO_OUTPUT") == 0) userParams.outLevel = OutLevel::NO_OUTPUT;
		else if (outLevelStr.compare("ERR") == 0) userParams.outLevel = OutLevel::ERR;
		else if (outLevelStr.compare("WARNING") == 0) userParams.outLevel = OutLevel::WARNING;
		else if (outLevelStr.compare("TRACE") == 0) userParams.outLevel = OutLevel::TRACE;
		else if (outLevelStr.compare("INFO") == 0) userParams.outLevel = OutLevel::INFO;
	}

	// Log boolean
	if (jValue.contains("log")) {
		userParams.log = jValue["log"].as_bool();
	}

	// path parameters
	if (jValue.contains("pathParams")) {
		if (userParams.method == Method::GENERATE_MICROFLAKES)
			userParams.pathParams = createPathParamsForHF(jValue["pathParams"]);
		else
			userParams.pathParams = createPathParamsForOBJ(jValue["pathParams"]);
	}

	// Direction Parameters
	if (jValue.contains("directionInParams")) {
		userParams.directionInParams = createDirectionParams(jValue["directionInParams"]);
	}
	if (jValue.contains("directionParams")) {
		userParams.directionOutParams = createDirectionParams(jValue["directionParams"]);
	}
	if (jValue.contains("directionOutParams")) {
		userParams.directionOutParams = createDirectionParams(jValue["directionOutParams"]);
	}

	// Side Effect Parameters
	if (jValue.contains("sideEffectParams")) {
		userParams.sideEffectParams = createSideEffectParams(jValue["sideEffectParams"]);
	}

	// Rendering Parameters
	if (jValue.contains("renderingParams")) {
		userParams.renderingParams = createRenderingParams(jValue["renderingParams"]);
	}

	// NDF Parameters
	if (jValue.contains("ndfParams")) {
		userParams.ndfParams = createNDFParams(jValue["ndfParams"]);
	}

	// Binary Files Parameters
	if (jValue.contains("binaryFilesParams")) {
		userParams.binaryFilesParams = createBinaryFilesParams(jValue["binaryFilesParams"]);
	}


	return userParams;
}

//  {
//  	"objNames": ["PerTex/001.obj"] ,
//      "objMin" : 1,
//  	"objMax" : 10,
//  	"inputsFolder" : "../../inputs/",
//  	"outputsFolder" : "../../outputs/"
//  }
PathParams Parameters::createPathParamsForOBJ(jParser::jValue jValue) const
{
	PathParams pathParams;

	if (jValue.contains("objNames") || (jValue.contains("objMin") && jValue.contains("objMax")))
		pathParams.surfNames.clear();

	// Mesh resolution
	std::vector<std::string> resolutionFolders;
	if (jValue.contains("resMin") && jValue.contains("resMax")) {
		for (int i = jValue["resMin"].as_int(); i <= jValue["resMax"].as_int(); ++i) {
			pathParams.resolutions.push_back(i);
		}
	}
	else if (jValue.contains("res")) {
		pathParams.resolutions.push_back(jValue["res"].as_int());
	}

	// obj names
	if (jValue.contains("objNames")) {
		for (int i = 0; i < jValue["objNames"].size(); ++i) {
			pathParams.surfNames.push_back(jValue["objNames"][i].as_string());
		}
	}

	// PerTex range
	if (jValue.contains("objMin") && jValue.contains("objMax")) {
		for (int i = jValue["objMin"].as_int(); i <= jValue["objMax"].as_int(); ++i) {
		    std::ostringstream oss;
		    oss << std::setfill('0') << std::setw(3) << i;
			pathParams.surfNames.push_back("PerTex/" + oss.str());
		}
	}

	// Folders
	if (jValue.contains("resetOutput")) {
		pathParams.resetOutput = jValue["resetOutput"].as_bool();
	}
	if (jValue.contains("objFolder")) {
		pathParams.objFolder = jValue["objFolder"].as_string();
	}
	if (jValue.contains("outputsFolder")) {
		pathParams.outputsFolder = jValue["outputsFolder"].as_string();
	}

	// Gnuplot executable
	if (jValue.contains("gnuplotPath") && jValue["gnuplotPath"].as_string().size() > 0) {
		pathParams.gnuplotPath = jValue["gnuplotPath"].as_string();
	}
	else {
		pathParams.gnuplotPath = "gnuplot";
	}

	// Cuda source
	if (jValue.contains("ptxFile")) {
		pathParams.ptxFile = jValue["ptxFile"].as_string();
	}

	return pathParams;
}


PathParams Parameters::createPathParamsForHF(jParser::jValue jValue) const
{
	PathParams pathParams;

	if (jValue.contains("objNames") || (jValue.contains("objMin") && jValue.contains("objMax")))
		pathParams.surfNames.clear();

	// hf names
	if (jValue.contains("objNames")) {
		for (int i = 0; i < jValue["objNames"].size(); ++i) {
			pathParams.surfNames.push_back(jValue["objNames"][i].as_string());
		}
	}

	// PerTex range
	if (jValue.contains("objMin") && jValue.contains("objMax")) {
		for (int i = jValue["objMin"].as_int(); i <= jValue["objMax"].as_int(); ++i) {
			std::ostringstream oss;
			oss << std::setfill('0') << std::setw(3) << i;
			pathParams.surfNames.push_back(oss.str());
		}
	}
	
	// Output mesh resolutions
	if ((jValue.contains("resMin") && jValue.contains("resMax")) || jValue.contains("res")) {
		pathParams.resolutions.clear();
	}
	if (jValue.contains("resMin") && jValue.contains("resMax")) {
		for (int i = jValue["resMin"].as_int(); i <= jValue["resMax"].as_int(); ++i) {
			pathParams.resolutions.push_back(i);
		}
	}
	else if (jValue.contains("res")) {
		pathParams.resolutions.push_back(jValue["res"].as_int());
	}

	// Folders
	if (jValue.contains("hfFolder")) {
		pathParams.hfFolder = jValue["hfFolder"].as_string();
	}
	if (jValue.contains("objFolder")) {
		pathParams.objFolder = jValue["objFolder"].as_string();
	}
	if (jValue.contains("outputsFolder")) {
		pathParams.outputsFolder = jValue["outputsFolder"].as_string();
	}

	return pathParams;
}

// {
//     "phiStart" : 0,
//     "phiEnd" : 3.141593,
//     "nAzimuthSamples" : 1,
//     "thetaStart" : -1.570796,
//     "thetaEnd" : 1.570796,
//     "nElevationSamples" : 10
// }
DirectionParams Parameters::createDirectionParams(jParser::jValue jValue) const
{
	DirectionParams directionParams;

	if (jValue.contains("phiStart")) {
		directionParams.phiStart = jValue["phiStart"].as_float();
	}
	if (jValue.contains("phiEnd")) {
		directionParams.phiEnd = jValue["phiEnd"].as_float();
	}
	if (jValue.contains("nAzimuthSamples")) {
		directionParams.nAzimuthSamples = jValue["nAzimuthSamples"].as_int();
	}

	if (jValue.contains("thetaStart")) {
		directionParams.thetaStart = jValue["thetaStart"].as_float();
	}
	if (jValue.contains("thetaEnd")) {
		directionParams.thetaEnd = jValue["thetaEnd"].as_float();
	}
	if (jValue.contains("nElevationSamples")) {
		directionParams.nElevationSamples = jValue["nElevationSamples"].as_int();
	}

	return directionParams;
}

// {
// 	    "borderPercentage" : 0.5,
// 		"directional" : true
// }
SideEffectParams Parameters::createSideEffectParams(jParser::jValue jValue) const
{
	SideEffectParams sideEffectParams;

	if (jValue.contains("borderPercentage")) {
		sideEffectParams.borderPercentage = jValue["borderPercentage"].as_float();
	}
	if (jValue.contains("directional")) {
		sideEffectParams.directional = jValue["directional"].as_bool();
	}

	return sideEffectParams;
}

// {
// 	   "renderSize" : [1000, 1000],
// 	   "nPixelSamples" : 16,
// 	   "createPicture" : false
// }
RenderingParams Parameters::createRenderingParams(jParser::jValue jValue) const
{
	RenderingParams renderingParams;

	if (jValue.contains("renderSize")) {
		if (jValue["renderSize"].size() == 2) {
			renderingParams.renderSize[0] = jValue["renderSize"][0].as_int();
			renderingParams.renderSize[1] = jValue["renderSize"][1].as_int();
		}
	}
	if (jValue.contains("nPixelSamples")) {
		renderingParams.nPixelSamples = jValue["nPixelSamples"].as_int();
	}
	if (jValue.contains("createPicture")) {
		renderingParams.createPicture = jValue["createPicture"].as_bool();
	}
	if (jValue.contains("useSmooth")) {
		renderingParams.useSmooth = jValue["useSmooth"].as_bool();
	}

	return renderingParams;
}

// {
// 	   "profil" : "BECKMANN",
//	   "findAlphaFromFileName" : true,
// 	   "alpha" : [0.5, 0.5]
// }
NDFParams Parameters::createNDFParams(jParser::jValue jValue) const
{
	NDFParams ndfParams;

	if (jValue.contains("profil")) {
		const std::string profilStr = jValue["profil"].as_string();
		if (profilStr.compare("BECKMANN") == 0) ndfParams.profil = MicrofacetProfil::BECKMANN;
		else if (profilStr.compare("GGX") == 0) ndfParams.profil = MicrofacetProfil::GGX;
	}
	if (jValue.contains("findAlphaFromFileName")) {
		ndfParams.findAlphaFromFileName = jValue["findAlphaFromFileName"].as_bool();
	}
	if (jValue.contains("alpha")) {
		if (jValue["alpha"].size() == 2) {
			ndfParams.alpha.x = jValue["alpha"][0].as_float();
			ndfParams.alpha.y = jValue["alpha"][1].as_float();
		}
	}

	return ndfParams;
}

// {
// 	   "readFromBinary" : true,
//	   "createBinary" : true,
// 	   "overwrite_binary" : true
// }
BinaryFilesParams Parameters::createBinaryFilesParams(jParser::jValue jValue) const
{
	BinaryFilesParams binaryFilesParams;

	if (jValue.contains("readFromBinary")) {
		binaryFilesParams.readFromBinary = jValue["readFromBinary"].as_bool();
	}
	if (jValue.contains("createBinary")) {
		binaryFilesParams.createBinary = jValue["createBinary"].as_bool();
	}
	if (jValue.contains("overwriteBinary")) {
		binaryFilesParams.overwriteBinary = jValue["overwriteBinary"].as_bool();
	}

	return binaryFilesParams;
}






size_t Parameters::getNumberOfLaunchs() const
{
	return userParamsVector.size();
}

const UserParams& Parameters::getParamsForLaunch(int i) const
{
	if (i < userParamsVector.size())
		userParams = userParamsVector[i];
	else
		userParams = UserParams();
	return userParams;
}