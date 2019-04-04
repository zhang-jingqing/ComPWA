/*! program : Y -> pi0 pi0 D0 D0bar 
* @file Pi0Pi0D0D0bar.cpp
*/

// Standard header files go here
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <utility>

#include <boost/program_options.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

// Boost header files
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>


// Core header files go here
#include "Core/Event.hpp"
//#include "Core/Particle.hpp"
//#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
//#include "Core/FunctionTree.hpp"
//#include "Core/TableFormater.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"

// ComPWA header files go here
#include "Data/DataSet.hpp"
#include "Data/RootIO/RootDataIO.hpp"

#include "Physics/ParticleList.hpp"
//#include "Physics/HelicityFormalism.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/CoherentIntensity.hpp"
#include "Physics/IntensityBuilderXML.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"

#include "Tools/RootGenerator.hpp"
#include "Tools/Generate.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/ParameterTools.hpp"
#include "Tools/ComponentIntensity.hpp"
#include "Tools/Plotting/ROOT/RootPlotData.hpp"

using namespace std;
using namespace ComPWA;
using namespace ComPWA::Data;
using ComPWA::Data::DataSet;
using ComPWA::Data::RootDataIO;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;
using namespace ComPWA::Physics;
using namespace ComPWA::Physics::HelicityFormalism;

// Enable serialization of MinuitResult. For some reason has to be outside
// any namespaces.
BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)

/*****************************************************************************
 *
 * The main function.
 *
 *****************************************************************************/

int main(int argc, char **argv) {

  // read configuration for program setting/arguments
  namespace po = boost::program_options; 

  std::string config_file;
  po::options_description generic("Generic options");
  generic.add_options()("help,h", "produce help message");
  generic.add_options()(
      "config,c", po::value<string>(&config_file)->default_value("config.cfg"),
      "name of a file of a configuration.");

  std::string logLevel;
  bool saveLogFile;

  std::string logFileName;
  std::string plotDataName;
  std::string plotFitName;

  std::string dataFile;
  std::string dataTree;
  std::string mcFile;
  std::string mcTree;
  std::string mcTrueFile;
  std::string mcTrueTree;

  std::string execMode;
  std::string modelFile;

  int nPhspMC, nData;
  int randomSeed;
  bool useMinos, useHesse, useRandomStartValues, usePhspRangeLimit;
  int noErrorSampling = -1;
  
  po::options_description config("General settings");
  config.add_options()(
      "logLevel", po::value<std::string>(&logLevel)->default_value("trace"),
      "set log level: error|warning|info|DEBUG");
  config.add_options()("saveLogFile",
                       po::value<bool>(&saveLogFile)->default_value(0),
                       "write log file? 0/1");
  config.add_options()(
      "logFileName", po::value<std::string>(&logFileName)
      ->default_value("out.log"), "The files logFileName are created");
  config.add_options()(
      "plotDataName", po::value<std::string>(&plotDataName)
      ->default_value("plot_data.root"), "root for data plot");
  config.add_options()(
      "plotFitName", po::value<std::string>(&plotFitName)
      ->default_value("plot_fit.root"), "root for fit plot");

  config.add_options()(
      "dataFile", po::value<std::string>(&dataFile)->default_value(""),
      "set input data file.");
  config.add_options()(
      "dataTree", po::value<std::string>(&dataTree)->default_value("tree"),
      "set tree name of input data.");
  config.add_options()(
      "mcFile", po::value<std::string>(&mcFile)->default_value(""),
      "set input phsp file.");
  config.add_options()(
      "mcTree", po::value<std::string>(&mcTree)->default_value("tree"),
      "set tree name of input mc.");
  config.add_options()(
      "mcTrueFile", po::value<std::string>(&mcTrueFile)->default_value(""),
      "set input phsp file.");
  config.add_options()(
      "mcTrueTree", po::value<std::string>(&mcTrueTree)->default_value("tree"),
      "set tree name of input mc.");
  
  config.add_options()(
      "execMode", po::value<std::string>(&execMode)->default_value("Fit"),
      "set execute mode; default is Fit.");

  config.add_options()(
     "modelFile", po::value<std::string>(&modelFile)
     ->default_value("model.xml"), "xml model for fit/generation");

  config.add_options()(
      "nPhspMC", po::value<int>(&nPhspMC)->default_value(10000),
      "Number of input MC; default is 10000");
  config.add_options()(
      "nData", po::value<int>(&nData)->default_value(1000),
      "Number of Data/toyData");

  config.add_options()(
      "randomSeed", po::value<int>(&randomSeed)->default_value(-1),
      "set random number randomSeed");

  config.add_options()(
      "noErrorSampling", po::value<int>(&noErrorSampling)->default_value(-1),
      "Number of Sampling when calculate statistical errors");

  config.add_options()(
      "useRandomStartValues", po::value<bool>(&useRandomStartValues)
      ->default_value(false), "use random start values or not in the fit");

  config.add_options()(
      "usePhspRangeLimit", po::value<bool>(&usePhspRangeLimit)
      ->default_value(true), "if remove the events beyond the phsp range");

  config.add_options()(
      "useMinos", po::value<bool>(&useMinos)->default_value(false),
      "use minos or not in fit");
  config.add_options()(
      "useHesse", po::value<bool>(&useHesse)->default_value(true),
      "use hesse or not in fit");

  po::options_description cmdline_options;
  cmdline_options.add(generic);
  cmdline_options.add(config);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm); 

  ifstream fin(config_file.c_str());
  if (fin) {
    store(parse_config_file(fin, config), vm);
    notify(vm);
  }
  fin.close();

  //if (!saveLogFile) logFileName = "";
  Logging log(logFileName, logLevel); // initialize logging

  LOG(INFO) << "PWA program starts";
  //LOG(INFO) << "minuitstragety.xml (optional) is configuration for minuit2 ";
  LOG(INFO) << "execMode = " << execMode << std::endl;

  if (!(execMode == "GenPhspMC" || execMode == "GenToyData" || execMode == "Fit"
      || execMode == "Calc" || execMode == "Draw" || execMode == "FitCalc" 
      || execMode == "FitDraw" || execMode == "FitCalcDraw")) {
    LOG(INFO) << "Unknown execMode! Possible choices are:" << std::endl;
    LOG(INFO) << "GenPhspMC, GenToyData, " << std::endl;
    LOG(INFO) << "Fit, FitDraw, FitCalc, FitCalcDraw, " << std::endl;
    LOG(INFO) << "Calc, Draw" << std::endl;
    LOG(INFO) << "Calc means calculate fit fractions/yields and their errors" 
        << std::endl;
    return 0;
  }
  
  ComPWA::Physics::IntensityBuilderXML Builder;

  boost::property_tree::ptree modelTree;
  boost::property_tree::xml_parser::read_xml(modelFile, modelTree);
  std::vector<pid> initialState, finalState;
  auto partList = std::make_shared<ComPWA::PartList>();
  ReadParticles(partList, modelTree);
  auto heliKins = Builder.createHelicityKinematics(partList,
      modelTree.get_child("HelicityKinematics"));
  //global/total intensity, StrengthIntensity
  auto intensity = Builder.createIntensity(partList, heliKins,
      modelTree.get_child("Intensity"));
  auto gen = std::make_shared<ComPWA::Tools::RootGenerator>(
      heliKins->getParticleStateTransitionKinematicsInfo(), randomSeed);

  std::shared_ptr<ComPWA::Data::DataSet> dataSample;
  std::shared_ptr<ComPWA::Data::DataSet> phspSample;
  std::shared_ptr<ComPWA::Data::DataSet> phspTrueSample;
  RootDataIO phspIO(mcTree, nPhspMC);
  RootDataIO phspTrueIO(mcTrueTree, -1);
  RootDataIO dataIO(dataTree, nData);

  //read MC/gen phsp mc
  if (mcFile.empty() || execMode == "GenPhspMC") {
    if (execMode == "GenPhspMC" && nPhspMC <= 0) {
      std::cout << " nPhspMC <= 0, no MC sample is generated" << std::endl;
      return 0;
    }
    phspSample = ComPWA::Tools::generatePhsp(nPhspMC, gen);
    if (execMode == "GenPhspMC") {
      if (mcFile.empty()) mcFile == "genPhspMC.root";
      phspIO.writeData(phspSample, mcFile);
      std::cout << "FINISHED" << std::endl;
      return 0;
    }
  } else {
    phspSample = phspIO.readData(mcFile);
  }
  if (usePhspRangeLimit) phspSample->reduceToPhaseSpace(heliKins);

  //read true MC
  if (execMode == "FitCalc" || execMode == "FitCalcDraw" 
      || execMode == "CalcFitFraction") {
    if (mcTrueFile.empty()) {
      LOG(INFO) << " no mcTrueFile, use mcFile instead.";
      LOG(INFO) << " But result is not precision. " << std::endl;
      phspTrueSample = phspSample;
    } else {
      phspTrueSample = phspTrueIO.readData(mcTrueFile);
    }
  }
 
  if (execMode == "GenToyData") {
    if (nData <= 0) {
      std::cout << " nData <= 0, no toy data sample is generated" << std::endl;
      return 0;
    }
    dataSample = ComPWA::Tools::generate(nData, heliKins, gen, intensity,
        phspSample);
    if (dataFile.empty()) dataFile == "genToyData.root";
    dataIO.writeData(dataSample, dataFile);
    LOG(INFO) << "FINISHED" << std::endl;
    return 0;
  }
  //read data
  dataSample = dataIO.readData(dataFile);
  if (usePhspRangeLimit) dataSample->reduceToPhaseSpace(heliKins);

  phspSample->convertEventsToParameterList(heliKins);
  dataSample->convertEventsToParameterList(heliKins);

  LOG(DEBUG) << " Data/MC input ok " << std::endl;

  auto estimator = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
      intensity, dataSample, phspSample);   
  LOG(DEBUG) << " createMinLogLH ok" << std::endl;
  LOG(INFO) << estimator->print(25);

  std::shared_ptr<FitResult> result;

  ParameterList fitPars, finalPars;
  intensity->addUniqueParametersTo(fitPars);
  LOG(DEBUG) << fitPars.to_str();

  LOG(INFO) << " Initial LH = " << estimator->evaluate() << ". " << std::endl;
    
  if (execMode == "Fit" || execMode == "FitDraw" || execMode == "FitCalc" 
      || execMode == "FitCalcDraw") {
    auto minuitIF = new Optimizer::Minuit2::MinuitIF(estimator, fitPars);

    LOG(DEBUG) << " minuitIF OK " << std::endl;
    minuitIF->setUseHesse(useHesse);
    minuitIF->setUseMinos(useMinos);
    //std::cout << fitPars << std::endl;
    LOG(DEBUG) << fitPars << std::endl;
    LOG(DEBUG) << " minuitIF OK";
//    setErrorOnParameterList(fitPars, 0.05, useMinos);
//      LOG(DEBUG) << " setErrorOnParameterList OK " << std::endl;
    LOG(INFO) << "useRandomStarValues = " << useRandomStartValues << std::endl;
    if (useRandomStartValues) {
      randomStartValues(fitPars);
    }
    result = minuitIF->exec(fitPars);
    LOG(DEBUG) << " minuitIF->exec(fitPars)";
    LOG(DEBUG) << fitPars << std::endl;
//    std::cout << fitPars << std::endl;
    finalPars = result->finalParameters();
    LOG(DEBUG) << " result->finalParameters() OK" << std::endl;
  }

  if (finalPars.numParameters() == 0) {
    finalPars = fitPars;
  }

  LOG(INFO) << "intensity->updateParametersFrom(fianlPars) OK" << std::endl;
  intensity->updateParametersFrom(finalPars);

  // SetComponentPaterns
  void SetComponentPaterns(std::vector<std::string> &componentNames,
      std::vector<std::vector<std::vector<std::string>>> &decayPaterns,
      std::vector<std::vector<std::pair<int, int>>> &decayLRanges,
      std::vector<std::vector<std::pair<int, int>>> &decaySRanges);

  std::vector<std::string> componentNames;
  std::vector<std::vector<std::vector<std::string>>> decayPaterns;
  std::vector<std::vector<std::pair<int, int>>> decayLRanges;
  std::vector<std::vector<std::pair<int, int>>> decaySRanges;

  SetComponentPaterns(componentNames, decayPaterns, decayLRanges, decaySRanges);
  LOG(INFO) << "SetComponentPaterns OK" << std::endl;

  // Get IncoherentIntensity to extract components
  // If we want to calc fit errors, Strength of StrengthIncoherentIntensity
  // must be fixed (which is the default setting)
  auto incoherentIntensity = Builder.createIncoherentIntensity(
      partList, heliKins, modelTree.get_child("Intensity.Intensity"));
  LOG(INFO) << "Get IncoherentIntensity OK" << std::endl;
  // use same parameters as intensity
  incoherentIntensity->updateParametersFrom(finalPars);
  LOG(INFO) << "incoherentIntensity->updateParametersFrom(finalPars) OK"
      << std::endl;
  
  //get component's intensity
  std::vector<std::shared_ptr<ComPWA::Intensity>> components(
      componentNames.size(), std::shared_ptr<ComPWA::Intensity>());
  for (std::size_t icomp = 0; icomp < componentNames.size(); ++icomp) {
    LOG(INFO) << "Getting component " << componentNames.at(icomp) << std::endl;
    auto icomponent = 
        ComPWA::Tools::getComponentIntensityFromIncoherentIntensity(
        std::dynamic_pointer_cast<ComPWA::Physics::IncoherentIntensity>(
        incoherentIntensity), componentNames.at(icomp), decayPaterns.at(icomp),
        decayLRanges.at(icomp), decaySRanges.at(icomp)); 
    components.at(icomp) = icomponent;
  }
  LOG(INFO) << "Get Component Intensity OK" << std::endl;

  if (execMode == "Calc" || execMode == "FitCalc" || execMode 
      == "FitCalcDraw") {
//121 ComPWA::ParameterList calculateFitFractions(
//122     std::shared_ptr<const ComPWA::Physics::CoherentIntensity> intensity,
//123     std::shared_ptr<ComPWA::Data::DataSet> sample,
//124     const std::vector<std::string> &components = {}) {
//    std::vector<std::string> components;
// do not know if it is right or not
    
  }

  if (execMode == "Draw" || execMode == "FitDraw" || execMode 
      == "FitCalcDraw") {
    ComPWA::Tools::Plotting::RootPlotData 
        plotData(heliKins->getParticleStateTransitionKinematicsInfo(), 
        plotDataName, "RECREATE");
    plotData.writeData(*dataSample);
    LOG(INFO) << "Plot Data OK" << std::endl;

    //global fit
    ComPWA::Tools::Plotting::RootPlotData
        plotWeightedMC(heliKins->getParticleStateTransitionKinematicsInfo(), 
        plotFitName, "RECREATE");
    plotWeightedMC.writeIntensityWeightedPhspSample(*phspSample, intensity);
    LOG(INFO) << "Plot Global Fit OK" << std::endl;
    //fit components
    for (std::size_t icomp = 0; icomp < componentNames.size(); ++icomp) {
      LOG(INFO) << "Plotting component " << componentNames.at(icomp)
          << std::endl;
      plotWeightedMC.writeIntensityWeightedPhspSample(*phspSample,
          incoherentIntensity, std::map<std::string,
          std::shared_ptr<const ComPWA::Intensity>>{{componentNames.at(icomp),
          components.at(icomp)}});
    }
    LOG(INFO) << "Plot components OK" << std::endl;
  }

  LOG(INFO) << "FINISHED" << std::endl;
  return 0;
}

void SetComponentPaterns(std::vector<std::string> &componentNames,
    std::vector<std::vector<std::vector<std::string>>> &decayPaterns,
    std::vector<std::vector<std::pair<int, int>>> &decayLRanges,
    std::vector<std::vector<std::pair<int, int>>> &decaySRanges) {
  //a component will be set by 
  //a std::vector<std::vector<std::string>>,
  //and two std::vector<std::pair<int, int>>

  componentNames = std::vector<std::string>({"a10", "a11", "a12", "a32"});

  //a10
  std::vector<std::string> decay1({"EpEm", "D*(2007)0", "D*(2007)0bar"});
  //if there are more requirement in the decay chain,
  //we need more decay/LSRange in patern1/LSRanges
  std::vector<std::vector<std::string>> patern1(1, decay1);
  decayPaterns.push_back(patern1);
  decayLRanges.push_back(
      std::vector<std::pair<int, int>>(1, std::pair<int, int>(1, 1)));
  decaySRanges.push_back(
      std::vector<std::pair<int, int>>(1, std::pair<int, int>(0, 0)));

  //a_11
  decayPaterns.push_back(patern1);
  decayLRanges.push_back(
      std::vector<std::pair<int, int>>(1, std::pair<int, int>(1, 1)));
  decaySRanges.push_back(
      std::vector<std::pair<int, int>>(1, std::pair<int, int>(1, 1)));

  //a_12
  decayPaterns.push_back(patern1);
  decayLRanges.push_back(
      std::vector<std::pair<int, int>>(1, std::pair<int, int>(1, 1)));
  decaySRanges.push_back(
      std::vector<std::pair<int, int>>(1, std::pair<int, int>(2, 2)));

  //a_32
  decayPaterns.push_back(patern1);
  decayLRanges.push_back(
      std::vector<std::pair<int, int>>(1, std::pair<int, int>(3, 3)));
  decaySRanges.push_back(
      std::vector<std::pair<int, int>>(1, std::pair<int, int>(2, 2)));
}
