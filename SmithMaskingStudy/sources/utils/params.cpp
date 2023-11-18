#include "utils/params.h"

#include <iomanip>
#include <fstream>
#include <sstream>

#include "utils/paths.h"
#include "utils/console.h"


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
	while (std::getline(in, tmp))
		str += tmp; // we stock the whole file in a single string
	in.close();

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
			userParamsVector.push_back(createUserParams(allUserParams[i]));
		}
	}
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
		Console::err << "methodParams.name is required." << std::endl;
		requirements = false;
	}
	else {
		methodName = jValue["methodParams"]["name"].as_string();
	}
	if (jValue.contains("pathParams")) {
		if (!jValue["pathParams"].contains("objDir")) {
			Console::err << "pathParams.objDir is required." << std::endl;
			requirements = false;
		}
		if ((methodName.compare("GENERATE_MICROFLAKES") == 0) && (!jValue["pathParams"].contains("hfDir"))) {
			Console::err << "pathParams.hfDir is required with method GENERATE_MICROFLAKES." << std::endl;
			requirements = false;
		}
		if (!jValue["pathParams"].contains("outputsDir")) {
			Console::err << "pathParams.outputsDir is required." << std::endl;
			requirements = false;
		}
		if (jValue["pathParams"].contains("surfMin") || jValue["pathParams"].contains("surfMax")) {
			if (jValue["pathParams"]["surfMin"].as_int() > jValue["pathParams"]["surfMax"].as_int() && !jValue["pathParams"].contains("surfNames")) {
				Console::err << "At least on input surface is required (currently, pathParams.surfMin > pathParams.surfMax, and pathParams.surfNames is undefined. No surface names can be selected." << std::endl;
				requirements = false;
			}
		}
		else if (!jValue["pathParams"].contains("surfNames") && !(jValue["pathParams"].contains("surfMin") || jValue["pathParams"].contains("surfMax"))) {
			Console::err << "At least on input surface is required (can be defined with pathParams.surfNames, or pathParams.surfMin and pathParams.surfMax)." << std::endl;
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
UserParams Parameters::createUserParams(jParser::jValue jValue) const
{
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

	// Method Parameters
	if (jValue.contains("methodParams")) {
		userParams.methodParams = createMethodParams(jValue["methodParams"]);
	}

	// path parameters
	if (jValue.contains("pathParams")) {
		if (userParams.methodParams.method == Method::GENERATE_MICROFLAKES)
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

	// NDF Parameters
	if (jValue.contains("ndfParams")) {
		userParams.ndfParams = createNDFParams(jValue["ndfParams"]);
	}

	// Side Effect Parameters
	if (jValue.contains("sideEffectParams")) {
		userParams.sideEffectParams = createSideEffectParams(jValue["sideEffectParams"]);
	}

	// Rendering Parameters
	if (jValue.contains("renderingParams")) {
		userParams.renderingParams = createRenderingParams(jValue["renderingParams"]);
	}

	return userParams;
}

// "methodParams": {
//     "name": "G1" | "GAF" | "D_TABULATION" | "FEATURES" | "FULL_PIPELINE" | "AMBIENT_OCCLUSION" | "GENERATE_MICROFLAKES",
//     "computeError" : bool,
//     "jiggleFlakes" : bool
// }
MethodParams Parameters::createMethodParams(jParser::jValue jValue) const {
	MethodParams methodParams;

	// Method
	if (jValue.contains("name")) {
		const std::string methodStr = jValue["name"].as_string();
		if (methodStr.compare("GAF") == 0) methodParams.method = Method::GAF;
		else if (methodStr.compare("G1") == 0) methodParams.method = Method::G1;
		else if (methodStr.compare("D_TABULATION") == 0) methodParams.method = Method::D_TABULATION;
		else if (methodStr.compare("FEATURES") == 0) methodParams.method = Method::FEATURES;
		else if (methodStr.compare("FULL_PIPELINE") == 0) methodParams.method = Method::FULL_PIPELINE;
		else if (methodStr.compare("AMBIENT_OCCLUSION") == 0) methodParams.method = Method::AMBIENT_OCCLUSION;
		else if (methodStr.compare("GENERATE_MICROFLAKES") == 0) methodParams.method = Method::GENERATE_MICROFLAKES;
	}
	if (jValue.contains("computeError")) {
		methodParams.computeError = jValue["computeError"].as_bool();
	}
	if (jValue.contains("jiggleFlakes")) {
		methodParams.jiggleFlakes = jValue["jiggleFlakes"].as_bool();
	}

	return methodParams;
}

//  {
//  	"surfNames": ["PerTex/001.obj"] ,
//      "surfMin" : 1,
//  	"surfMax" : 10,
//  	"inputsFolder" : "../../inputs/",
//  	"outputsDir" : "../../outputs/"
//  }
PathParams Parameters::createPathParamsForOBJ(jParser::jValue jValue) const
{
	PathParams pathParams;

	if (jValue.contains("surfNames") || (jValue.contains("surfMin") && jValue.contains("surfMax")))
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

	// surface names
	if (jValue.contains("surfNames")) {
		for (int i = 0; i < jValue["surfNames"].size(); ++i) {
			pathParams.surfNames.push_back(jValue["surfNames"][i].as_string());
		}
	}

	// PerTex range
	if (jValue.contains("surfMin") && jValue.contains("surfMax")) {
		for (int i = jValue["surfMin"].as_int(); i <= jValue["surfMax"].as_int(); ++i) {
		    std::ostringstream oss;
		    oss << std::setfill('0') << std::setw(3) << i;
			pathParams.surfNames.push_back("PerTex/" + oss.str());
		}
	}

	// Folders
	if (jValue.contains("objDir")) {
		pathParams.objDir = jValue["objDir"].as_string();
	}
	if (jValue.contains("outputsDir")) {
		pathParams.outputsDir = jValue["outputsDir"].as_string();
	}

	// Gnuplot executable
	if (jValue.contains("gnuplotPath") && jValue["gnuplotPath"].as_string().size() > 0) {
		pathParams.gnuplotPath = jValue["gnuplotPath"].as_string();
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

	if (jValue.contains("surfNames") || (jValue.contains("surfMin") && jValue.contains("surfMax")))
		pathParams.surfNames.clear();

	// hf names
	if (jValue.contains("surfNames")) {
		for (int i = 0; i < jValue["surfNames"].size(); ++i) {
			pathParams.surfNames.push_back(jValue["surfNames"][i].as_string());
		}
	}

	// PerTex range
	if (jValue.contains("surfMin") && jValue.contains("surfMax")) {
		for (int i = jValue["surfMin"].as_int(); i <= jValue["surfMax"].as_int(); ++i) {
			std::ostringstream oss;
			oss << std::setfill('0') << std::setw(3) << i;
			pathParams.surfNames.push_back(oss.str());
		}
	}
	
	// Output mesh resolutions
	if (jValue.contains("resolutions")) {
		pathParams.resolutions.clear();
		for (int i = 0; i < jValue["resolutions"].size(); ++i) {
			pathParams.resolutions.push_back(jValue["resolutions"][i].as_int());
		}
	}

	// Folders
	if (jValue.contains("hfDir")) {
		pathParams.hfDir = jValue["hfDir"].as_string();
	}
	if (jValue.contains("objDir")) {
		pathParams.objDir = jValue["objDir"].as_string();
	}
	if (jValue.contains("outputsDir")) {
		pathParams.outputsDir = jValue["outputsDir"].as_string();
	}

	// Gnuplot executable
	if (jValue.contains("gnuplotPath") && jValue["gnuplotPath"].as_string().size() > 0) {
		pathParams.gnuplotPath = jValue["gnuplotPath"].as_string();
	}

	// Cuda source
	if (jValue.contains("ptxFile")) {
		pathParams.ptxFile = jValue["ptxFile"].as_string();
	}

	return pathParams;
}

// {
//     "phiStart" : 0,
//     "phiEnd" : 3.141593,
//     "nPhiSamples" : 1,
//     "thetaStart" : -1.570796,
//     "thetaEnd" : 1.570796,
//     "nThetaSamples" : 10
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

	return directionParams;
}

// {
// 	    "nTheta" : 100,
// 		"nPhi" : 400
// }
NDFParams Parameters::createNDFParams(jParser::jValue jValue) const
{
	NDFParams ndfParams;

	if (jValue.contains("nTheta")) {
		ndfParams.nTheta = jValue["nTheta"].as_int();
	}
	if (jValue.contains("nPhi")) {
		ndfParams.nPhi = jValue["nPhi"].as_int();
	}

	return ndfParams;
}

// {
// 	    "borderPercentage" : 0.2,
// 		"BBox" : false
// }
SideEffectParams Parameters::createSideEffectParams(jParser::jValue jValue) const
{
	SideEffectParams sideEffectParams;

	if (jValue.contains("borderPercentage")) {
		sideEffectParams.borderPercentage = jValue["borderPercentage"].as_float();
	}
	if (jValue.contains("BBox")) {
		sideEffectParams.BBox = jValue["BBox"].as_bool();
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