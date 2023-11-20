#include "utils/params.h"

#include <iomanip>
#include <fstream>
#include <sstream>

#include "utils/paths.h"
#include "utils/console.h"

Parameters* Parameters::instancePtr = nullptr;


bool operator==(const MethodParams& lhs, const MethodParams& rhs) {
	if (lhs.method != rhs.method) return false;
	if (lhs.computeError != rhs.computeError) return false;
	if (lhs.jiggleFlakes != rhs.jiggleFlakes) return false;
	return true;
}

bool operator!=(const MethodParams& lhs, const MethodParams& rhs) {
	return !(lhs == rhs);
}

bool operator==(const PathParams& lhs, const PathParams& rhs) {
	if (lhs.objDir.compare(rhs.objDir) != 0) return false;
	if (lhs.hfDir.compare(rhs.hfDir) != 0) return false;
	if (lhs.outputsDir.compare(rhs.outputsDir) != 0) return false;
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
	if (lhs.nPhiSamples != rhs.nPhiSamples) return false;
	if (lhs.thetaStart != rhs.thetaStart) return false;
	if (lhs.thetaEnd != rhs.thetaEnd) return false;
	if (lhs.nThetaSamples != rhs.nThetaSamples) return false;
	return true;
}

bool operator!=(const DirectionParams& lhs, const DirectionParams& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const NDFParams& lhs, const NDFParams& rhs)
{
	if (lhs.nTheta != rhs.nTheta) return false;
	if (lhs.nPhi != rhs.nPhi) return false;
	return true;
}

bool operator!=(const NDFParams& lhs, const NDFParams& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const SideEffectParams& lhs, const SideEffectParams& rhs)
{
	if (lhs.borderPercentage != rhs.borderPercentage) return false;
	if (lhs.BBox != rhs.BBox) return false;
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
	if (lhs.useSmooth != rhs.useSmooth) return false;
	return true;
}

bool operator!=(const RenderingParams& lhs, const RenderingParams& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const UserParams& lhs, const UserParams& rhs)
{
	if (lhs.outLevel != rhs.outLevel) return false;
	if (lhs.log != rhs.log) return false;
	if (lhs.methodParams != rhs.methodParams) return false;
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

Parameters* Parameters::get()
{
	if (instancePtr == nullptr) {
		instancePtr = new Parameters();
		return instancePtr;
	}
	else {
		return instancePtr;
	}
}

void Parameters::parse(const std::string& path)
{
	std::ifstream in(path.c_str());
	if (!in.is_open()) return;

	userParamsVector.clear();

	std::string str = "";
	std::string tmp;
	while (std::getline(in, tmp))
		str += tmp; // we stock the whole file in a single string

	jParser::jValue allUserParams = jParser::parser::parse(str)["userParams"];
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
		if (checkRequired(allUserParams[i])) {
			createUserParams(allUserParams[i]);
		}
	}

	userParams = &userParamsVector[0];
}

// Required :
// - methodParams.name
// - pathParams.objDir
// - pathParams.hfDir (if methodParams.name == "GENERATE_MICROFLAKES")
// - pathParams.outputsDir
// - pathParams.surfNames || (pathParams.surfMin && pathParams.surfMax)
bool Parameters::checkRequired(jParser::jValue jValue) const {
	bool requirements = true;
	std::string methodName = "";
	if (!jValue.contains("methodParams") || !jValue["methodParams"].contains("name")) {
		Console::print(OutLevel::ERR, "methodParams.name is required.");
		requirements = false;
	}
	else {
		methodName = jValue["methodParams"]["name"].as_string();
	}
	if (jValue.contains("pathParams")) {
		if (!jValue["pathParams"].contains("objDir")) {
			Console::print(OutLevel::ERR, "pathParams.objDir is required.");
			requirements = false;
		}
		if ((methodName.compare("GENERATE_MICROFLAKES") == 0) && (!jValue["pathParams"].contains("hfDir"))) {
			Console::print(OutLevel::ERR, "pathParams.hfDir is required with method GENERATE_MICROFLAKES.");
			requirements = false;
		}
		if (!jValue["pathParams"].contains("outputsDir")) {
			Console::print(OutLevel::ERR, "pathParams.outputsDir is required.");
			requirements = false;
		}
		if (jValue["pathParams"].contains("surfMin") || jValue["pathParams"].contains("surfMax")) {
			if (jValue["pathParams"]["surfMin"].as_int() > jValue["pathParams"]["surfMax"].as_int() && !jValue["pathParams"].contains("surfNames")) {
				Console::print(OutLevel::ERR, "At least on input surface is required (currently, pathParams.surfMin > pathParams.surfMax, and pathParams.surfNames is undefined. No surface names can be selected.");
				requirements = false;
			}
		}
		else if (!jValue["pathParams"].contains("surfNames") && !(jValue["pathParams"].contains("surfMin") || jValue["pathParams"].contains("surfMax"))) {
			Console::print(OutLevel::ERR, "At least on input surface is required (can be defined with pathParams.surfNames, or pathParams.surfMin and pathParams.surfMax).");
			requirements = false;
		}
		return requirements;
	}
}

// {
//     "method" : "GLOBAL_VISIBILITY",
//     "outLevel" : "TRACE",
//     "log" : false,
//     "pathParams": {
//      	"surfNames": ["PerTex/001.obj"] ,
//     		"surfMin" : 1,
//     		"surfMax" : 10,
//     		"inputsFolder" : "../../inputs/",
//     		"outputsDir" : "../../outputs/"
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
void Parameters::createUserParams(jParser::jValue jValue)
{
	userParamsVector.push_back(UserParams());
	userParams = &userParamsVector.back();
	// Out level
	if (jValue.contains("outLevel")) {
		const std::string outLevelStr = jValue["outLevel"].as_string();
		if (outLevelStr.compare("NO_OUTPUT") == 0) userParams->outLevel = OutLevel::NO_OUTPUT;
		else if (outLevelStr.compare("ERR") == 0) userParams->outLevel = OutLevel::ERR;
		else if (outLevelStr.compare("WARNING") == 0) userParams->outLevel = OutLevel::WARNING;
		else if (outLevelStr.compare("TRACE") == 0) userParams->outLevel = OutLevel::TRACE;
		else if (outLevelStr.compare("INFO") == 0) userParams->outLevel = OutLevel::INFO;
	}

	// Log boolean
	if (jValue.contains("log")) {
		userParams->log = jValue["log"].as_bool();
	}

	// Method Parameters
	if (jValue.contains("methodParams")) {
		createMethodParams(jValue["methodParams"]);
	}

	// path parameters
	if (jValue.contains("pathParams")) {
		createPathParams(jValue["pathParams"]);
	}

	// Direction Parameters
	if (jValue.contains("directionInParams")) {
		createDirectionParams(jValue["directionInParams"], userParams->directionInParams);
	}
	if (jValue.contains("directionParams")) {
		createDirectionParams(jValue["directionParams"], userParams->directionOutParams);
	}
	if (jValue.contains("directionOutParams")) {
		createDirectionParams(jValue["directionOutParams"], userParams->directionOutParams);
	}

	// NDF Parameters
	if (jValue.contains("ndfParams")) {
		createNDFParams(jValue["ndfParams"]);
	}

	// Side Effect Parameters
	if (jValue.contains("sideEffectParams")) {
		createSideEffectParams(jValue["sideEffectParams"]);
	}

	// Rendering Parameters
	if (jValue.contains("renderingParams")) {
		createRenderingParams(jValue["renderingParams"]);
	}
}

// "methodParams": {
//     "name": "G1" | "GAF" | "D_TABULATION" | "FEATURES" | "FULL_PIPELINE" | "AMBIENT_OCCLUSION" | "GENERATE_MICROFLAKES",
//     "computeError" : bool,
//     "jiggleFlakes" : bool
// }
void Parameters::createMethodParams(jParser::jValue jValue) {
	// Method
	if (jValue.contains("name")) {
		const std::string methodStr = jValue["name"].as_string();
		if (methodStr.compare("GAF") == 0) userParams->methodParams.method = Method::GAF;
		else if (methodStr.compare("G1") == 0) userParams->methodParams.method = Method::G1;
		else if (methodStr.compare("D_TABULATION") == 0) userParams->methodParams.method = Method::D_TABULATION;
		else if (methodStr.compare("FEATURES") == 0) userParams->methodParams.method = Method::FEATURES;
		else if (methodStr.compare("FULL_PIPELINE") == 0) userParams->methodParams.method = Method::FULL_PIPELINE;
		else if (methodStr.compare("AMBIENT_OCCLUSION") == 0) userParams->methodParams.method = Method::AMBIENT_OCCLUSION;
		else if (methodStr.compare("GENERATE_MICROFLAKES") == 0) userParams->methodParams.method = Method::GENERATE_MICROFLAKES;
	}
	if (jValue.contains("computeError")) {
		userParams->methodParams.computeError = jValue["computeError"].as_bool();
	}
	if (jValue.contains("jiggleFlakes")) {
		userParams->methodParams.jiggleFlakes = jValue["jiggleFlakes"].as_bool();
	}
}

//  {
//  	"surfNames": ["PerTex/001.obj"] ,
//      "surfMin" : 1,
//  	"surfMax" : 10,
//  	"inputsFolder" : "../../inputs/",
//  	"outputsDir" : "../../outputs/"
//  }
void Parameters::createPathParams(jParser::jValue jValue)
{
	if (jValue.contains("surfNames") || (jValue.contains("surfMin") && jValue.contains("surfMax")))
		userParams->pathParams.surfNames.clear();

	// surface names
	if (jValue.contains("surfNames")) {
		for (int i = 0; i < jValue["surfNames"].size(); ++i) {
			userParams->pathParams.surfNames.push_back(jValue["surfNames"][i].as_string());
		}
	}

	// PerTex range
	if (jValue.contains("surfMin") && jValue.contains("surfMax")) {
		for (int i = jValue["surfMin"].as_int(); i <= jValue["surfMax"].as_int(); ++i) {
		    std::ostringstream oss;
		    oss << std::setfill('0') << std::setw(3) << i;
			if (userParams->methodParams.method == Method::GENERATE_MICROFLAKES) {
				userParams->pathParams.surfNames.push_back(oss.str());
			}
			else {
				userParams->pathParams.surfNames.push_back("PerTex/" + oss.str());
			}
		}
	}

	// Input mesh resolutions
	if (jValue.contains("resolutions")) {
		userParams->pathParams.resolutions.clear();
		for (int i = 0; i < jValue["resolutions"].size(); ++i) {
			userParams->pathParams.resolutions.push_back(jValue["resolutions"][i].as_int());
		}
	}

	// Folders
	if (jValue.contains("objDir")) {
		userParams->pathParams.objDir = jValue["objDir"].as_string();
	}
	if (jValue.contains("hfDir")) {
		userParams->pathParams.hfDir = jValue["hfDir"].as_string();
	}
	if (jValue.contains("outputsDir")) {
		userParams->pathParams.outputsDir = jValue["outputsDir"].as_string();
	}

	// Gnuplot executable
	if (jValue.contains("gnuplotPath") && jValue["gnuplotPath"].as_string().size() > 0) {
		userParams->pathParams.gnuplotPath = jValue["gnuplotPath"].as_string();
	}

	// Cuda source
	if (jValue.contains("ptxFile")) {
		userParams->pathParams.ptxFile = jValue["ptxFile"].as_string();
	}
}

// {
//     "phiStart" : 0,
//     "phiEnd" : 3.141593,
//     "nPhiSamples" : 1,
//     "thetaStart" : -1.570796,
//     "thetaEnd" : 1.570796,
//     "nThetaSamples" : 10
// }
void Parameters::createDirectionParams(jParser::jValue jValue, DirectionParams& directionParams) 
{
	if (jValue.contains("phiStart")) {
		directionParams.phiStart = jValue["phiStart"].as_float();
	}
	if (jValue.contains("phiEnd")) {
		directionParams.phiEnd = jValue["phiEnd"].as_float();
	}
	if (jValue.contains("nPhiSamples")) {
		directionParams.nPhiSamples = jValue["nPhiSamples"].as_int();
	}

	if (jValue.contains("thetaStart")) {
		directionParams.thetaStart = jValue["thetaStart"].as_float();
	}
	if (jValue.contains("thetaEnd")) {
		directionParams.thetaEnd = jValue["thetaEnd"].as_float();
	}
	if (jValue.contains("nThetaSamples")) {
		directionParams.nThetaSamples = jValue["nThetaSamples"].as_int();
	}
}

// {
// 	    "nTheta" : 100,
// 		"nPhi" : 400
// }
void Parameters::createNDFParams(jParser::jValue jValue)
{
	if (jValue.contains("nTheta")) {
		userParams->ndfParams.nTheta = jValue["nTheta"].as_int();
	}
	if (jValue.contains("nPhi")) {
		userParams->ndfParams.nPhi = jValue["nPhi"].as_int();
	}
}

// {
// 	    "borderPercentage" : 0.2,
// 		"BBox" : false
// }
void Parameters::createSideEffectParams(jParser::jValue jValue)
{
	if (jValue.contains("borderPercentage")) {
		userParams->sideEffectParams.borderPercentage = jValue["borderPercentage"].as_float();
	}
	if (jValue.contains("BBox")) {
		userParams->sideEffectParams.BBox = jValue["BBox"].as_bool();
	}
}

// {
// 	   "renderSize" : [1000, 1000],
// 	   "nPixelSamples" : 16,
// 	   "createPicture" : false
// }
void Parameters::createRenderingParams(jParser::jValue jValue)
{
	if (jValue.contains("renderSize")) {
		if (jValue["renderSize"].size() == 2) {
			userParams->renderingParams.renderSize[0] = jValue["renderSize"][0].as_int();
			userParams->renderingParams.renderSize[1] = jValue["renderSize"][1].as_int();
		}
	}
	if (jValue.contains("nPixelSamples")) {
		userParams->renderingParams.nPixelSamples = jValue["nPixelSamples"].as_int();
	}
	if (jValue.contains("createPicture")) {
		userParams->renderingParams.createPicture = jValue["createPicture"].as_bool();
	}
	if (jValue.contains("useSmooth")) {
		userParams->renderingParams.useSmooth = jValue["useSmooth"].as_bool();
	}
}



void Parameters::setIndex(int i)
{
	userParams = &userParamsVector[i];
}