#include "FitApp.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <limits>
#include <iomanip>
#include <fstream>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"
#include "Data/Root/RootDataIO.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Physics/ParticleList.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Plotting/RootPlotData.hpp"
#include "Tools/UpdatePTreeParameter.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

#include "TRandom3.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

using namespace ComPWA;
using ComPWA::Optimizer::Minuit2::MinuitResult;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

int main(int argc, char **argv) {
  //./FitApp path_and_name_to_log path_and_name_to_config
  std::string ConfPath, LogPath;
  if (argc >= 2) {
    LogPath = std::string(argv[1]);
  } else {
    LogPath = "./fit.log";
  }
  if (argc == 3) {
    ConfPath = std::string(argv[2]);
  } else {
    ConfPath = "./config.cfg";
  }

  //initialize logging
  ComPWA::Logging Log("info", "./log"); //LogPath);
  
  LOG(INFO) << std::setprecision(std::numeric_limits<double>::max_digits10);

  LOG(INFO) << "Start program!";

  LOG(INFO) << "Read Configuration file...";
  //read configuration
  std::map<std::string, std::string> MapConf = readConfig(ConfPath);
  std::map<std::string, std::string>::iterator MapIter;
  MapIter = MapConf.find("LogLevel");
  if (MapIter != MapConf.end()) {
    Log.setLogLevel(MapIter->second);   
  }
  
  std::string FitResultXml("./output/fit_result.xml");
  getValue(MapConf, "FitResultXml", FitResultXml);
  std::string FitResultLog("./output/fit_result.log");
  getValue(MapConf, "FitResultLog", FitResultLog);
  std::fstream Fout(FitResultLog, std::ios::out);

  std::string ModelFile("./model.xml");
  getValue(MapConf, "ModelFile", ModelFile);

  int Seed = 1;
  MapIter = MapConf.find("Seed");
  getValue(MapConf, "Seed", Seed);
  int NData = 1000;
  getValue(MapConf, "NData", NData);
  int NMC = 10000;
  getValue(MapConf, "NMC", NMC);
  int NMCTrue = 10000;
  getValue(MapConf, "NMCTrue", NMCTrue);

  std::string DataFile("./io_data.root");
  getValue(MapConf, "DataFile", DataFile);
  std::string DataTree("t");
  getValue(MapConf, "DataTree", DataTree);
  std::string MCFile("./io_phsp.root");
  getValue(MapConf, "MCFile", MCFile);
  std::string MCTree("t");
  getValue(MapConf, "MCTree", MCTree);
  std::string MCTrueFile("./io_phsp.root");
  getValue(MapConf, "MCTrueFile", MCTrueFile);
  std::string MCTrueTree("t");
  getValue(MapConf, "MCTrueTree", MCTrueTree);

  bool LoadFitResult(false);
  getValue(MapConf, "LoadFitResult", LoadFitResult);

  bool GenPhspMC(false), GenPhspTrueMC(false), GenToyData(false);
  getValue(MapConf, "GenPhspMC", GenPhspMC);
  getValue(MapConf, "GenPhspTrueMC", GenPhspTrueMC);
  getValue(MapConf, "GenToyData", GenToyData);
  bool GenToyDataUseInputPhsp(false);
  getValue(MapConf, "GenToyDataUseInputPhsp", GenToyDataUseInputPhsp);
  if (GenToyDataUseInputPhsp && !GenPhspMC)
    GenPhspMC = true;
  bool FitData(false), DrawPlot(false);
  getValue(MapConf, "FitData", FitData);
  getValue(MapConf, "DrawPlot", DrawPlot);
  bool CalcFF(false), CalcFFError(false);
  getValue(MapConf, "CalcFF", CalcFF);
  getValue(MapConf, "CalcFFError", CalcFFError);
  int NParSet = 0;
  getValue(MapConf, "NParSet", NParSet);

  int TotalFit = 1000, ValidFit = 100, CovergedFit = 10;
  getValue(MapConf, "TotalFit", TotalFit);
  getValue(MapConf, "ValidFit", ValidFit);
  getValue(MapConf, "CovergedFit", CovergedFit);

  bool FixA10(false), FixA11(false), FixA12(false), FixA32(false);
  getValue(MapConf, "FixA10", FixA10);
  getValue(MapConf, "FixA11", FixA11);
  getValue(MapConf, "FixA12", FixA12);
  getValue(MapConf, "FixA32", FixA32);

  getValue(MapConf, "A10Mag", A10Mag);
  getValue(MapConf, "A10Phase", A10Phase);
  getValue(MapConf, "A11Mag", A11Mag);
  getValue(MapConf, "A11Phase", A11Phase);
  getValue(MapConf, "A12Mag", A12Mag);
  getValue(MapConf, "A12Phase", A12Phase);
  getValue(MapConf, "A32Mag", A32Mag);
  getValue(MapConf, "A32Phase", A32Phase);

  double A10MagInit = 1.0, A10PhaseInit = 0.0;
  getValue(MapConf, "A10MagInit", A10MagInit);
  getValue(MapConf, "A10PhaseInit", A10PhaseInit);
  double A11MagInit = 1.0, A11PhaseInit = 0.0;
  getValue(MapConf, "A11MagInit", A11MagInit);
  getValue(MapConf, "A11PhaseInit", A11PhaseInit);
  double A12MagInit = 1.0, A12PhaseInit = 0.0;
  getValue(MapConf, "A12MagInit", A12MagInit);
  getValue(MapConf, "A12PhaseInit", A12PhaseInit);
  double A32MagInit = 1.0, A32PhaseInit = 0.0;
  getValue(MapConf, "A32MagInit", A32MagInit);
  getValue(MapConf, "A32PhaseInit", A32PhaseInit);

  bool UseHesse(true), UseMinos(false);
  getValue(MapConf, "UseHesse", UseHesse);
  getValue(MapConf, "UseMinos", UseMinos);

  LOG(INFO) << "Read configuration OK";
  MapConf.clear();
  
  LOG(INFO) << "Read physics model configuration...";
  //read model configuration/xml
  boost::property_tree::ptree ModelTree;
  boost::property_tree::xml_parser::read_xml(ModelFile, ModelTree);
  boost::property_tree::ptree IntensityTree = ModelTree.get_child("Intensity");
  ComPWA::Tools::updateParameterRangeByType(IntensityTree, "Magnitude", 0, 10);
  ComPWA::Tools::updateParameterRangeByType(IntensityTree,
      "Phase", -3.14159, 3.14159);

  ComPWA::Tools::updateParameterValue(IntensityTree, A10Mag,
      A10MagInit);
  ComPWA::Tools::updateParameterValue(IntensityTree, A10Phase,
      A10PhaseInit);
  ComPWA::Tools::updateParameterValue(IntensityTree, A11Mag,
      A11MagInit);
  ComPWA::Tools::updateParameterValue(IntensityTree, A11Phase,
      A11PhaseInit);
  ComPWA::Tools::updateParameterValue(IntensityTree, A12Mag,
      A12MagInit);
  ComPWA::Tools::updateParameterValue(IntensityTree, A12Phase,
      A12PhaseInit);
  ComPWA::Tools::updateParameterValue(IntensityTree, A32Mag,
      A32MagInit);
  ComPWA::Tools::updateParameterValue(IntensityTree, A32Phase,
      A32PhaseInit);
  if (FixA10) {
    ComPWA::Tools::fixParameter(IntensityTree, A10Mag, 
        A10MagInit);
    ComPWA::Tools::fixParameter(IntensityTree, A10Phase, 
        A10PhaseInit);
  }
  if (FixA11) {
    ComPWA::Tools::fixParameter(IntensityTree, A11Mag, 
        A11MagInit);
    ComPWA::Tools::fixParameter(IntensityTree, A11Phase, 
        A11PhaseInit);
  }
  if (FixA12) {
    ComPWA::Tools::fixParameter(IntensityTree, A12Mag, 
        A12MagInit);
    ComPWA::Tools::fixParameter(IntensityTree, A12Phase, 
        A12PhaseInit);
  }
  if (FixA32) {
    ComPWA::Tools::fixParameter(IntensityTree, A32Mag, 
        A32MagInit);
    ComPWA::Tools::fixParameter(IntensityTree, A32Phase, 
        A32PhaseInit);
  }

  boost::property_tree::ptree UnnormalizedIntensityTree 
      = IntensityTree.get_child("Intensity.Intensity");

  //Read particle information
  auto PartList = std::make_shared<ComPWA::PartList>();
  ReadParticles(PartList, ComPWA::Physics::defaultParticleList);
  ReadParticles(PartList, ModelTree); //myParticles particleList_xml

  //ComPWA::Physics::IntensityBuilderXML BuilderPhsp;
  ComPWA::Physics::IntensityBuilderXML Builder;

  LOG(INFO) << "Build Helicity Kinematics...";
  //Create kinematics object
  auto HeliKine = 
      Builder.createHelicityKinematics(PartList,
      ModelTree.get_child("HelicityKinematics"));

  LOG(INFO) << "Create generator and RootDataIO...";
  //Generate phase space sample (toy MC, if needed)
  ComPWA::Data::Root::RootGenerator MyGenerator(
      HeliKine.getParticleStateTransitionKinematicsInfo());
  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(Seed);
  ComPWA::Data::Root::RootDataIO PhspSampleIO(MCTree, NMC);
  ComPWA::Data::Root::RootDataIO PhspTrueSampleIO(MCTrueTree, NMCTrue);
  
  LOG(INFO) << "PhspSample and PhspSampleIO...";
  std::vector<ComPWA::Event> PhspSample, PhspTrueSample;
  if (GenPhspMC) {
    LOG(INFO) << "GenPhspMC...";
    PhspSample = ComPWA::Data::generatePhsp(NMC, MyGenerator, RandomGenerator);
    PhspSampleIO.writeData(PhspSample, MCFile);
  } else {
    LOG(INFO) << "ReadPhspMC...";
    if (NMC > 0) PhspSample = PhspSampleIO.readData(MCFile);
  }
  
  if (GenPhspTrueMC) {
    LOG(INFO) << "GenPhspTrueMC...";
    PhspTrueSample = ComPWA::Data::generatePhsp(NMCTrue, MyGenerator,
        RandomGenerator);
    PhspTrueSampleIO.writeData(PhspTrueSample, MCTrueFile);
  } else {
    LOG(INFO) << "ReadPhspTrueMC...";
    if (NMCTrue > 0) PhspTrueSample = PhspTrueSampleIO.readData(MCTrueFile);
    else {
      LOG(INFO) << "No PhspTrueSample is assigned!!";
      LOG(INFO) << "Use PhspSample as PhspTrueMC!!!";
      PhspTrueSample = PhspSample;
    }
  }
  auto BuilderPhsp = ComPWA::Physics::IntensityBuilderXML(PhspSample);

  LOG(INFO) << "Create model intensity...";
  //Create Intensity from model configuration file
  auto ModelIntensity = BuilderPhsp.createIntensity(PartList, HeliKine, 
      IntensityTree);
  LOG(INFO) << "Create model intensity over";

  ComPWA::FitResult BestResult;

  if (LoadFitResult) {
    LOG(INFO) << "LoadFitResult...";
    Fout << "LoadFitResult" << std::endl;
    std::ifstream Ifs(FitResultXml);
    boost::archive::xml_iarchive Ia(Ifs);
    Ia >> BOOST_SERIALIZATION_NVP(BestResult);
    Ifs.close();
  }
//
  LOG(INFO) << "DataSample and DataSampleIO...";
  ComPWA::Data::Root::RootDataIO DataSampleIO(DataTree, NData);
  std::vector<ComPWA::Event> DataSample;
  //Generate a data sample given intensity and kinematics
  if (GenToyData) {
    LOG(INFO) << "GenToyData...";
    if (!GenToyDataUseInputPhsp) {
      DataSample = ComPWA::Data::generate(NData, HeliKine, MyGenerator,
          ModelIntensity, RandomGenerator);
    } else {
      DataSample = ComPWA::Data::generate(NData, HeliKine, RandomGenerator,
          ModelIntensity, PhspSample, PhspSample);
    }
    DataSampleIO.writeData(DataSample, DataFile);
  } else {
    LOG(INFO) << "ReadData...";
    if (NData > 0) DataSample = DataSampleIO.readData(DataFile);
  }

  LOG(INFO) << "Convert events to data set...";
  //Fit to the data and print the result
  auto PhspSampleSet = ComPWA::Data::convertEventsToDataSet(PhspSample, 
      HeliKine);
  auto PhspTrueSampleSet = ComPWA::Data::convertEventsToDataSet(PhspTrueSample,
      HeliKine);
  auto DataSampleSet = Data::convertEventsToDataSet(DataSample, HeliKine);
  LOG(INFO) << "DataSet.size = " << DataSampleSet.Weights.size()
      << " PhspSampleSet.size = " << PhspSampleSet.Weights.size();
  std::vector<ComPWA::Event> TinySample;
  TinySample.push_back(DataSample[0]);
  auto TinySampleSet = Data::convertEventsToDataSet(TinySample, HeliKine);

  LOG(INFO) << "Create estimator and Minuit2...";
  auto MyEstimator = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
      ModelIntensity, DataSampleSet); //the parameters in intensity then is 
                                      //put to MyEstimator (FitParList, with fix or not info),
                                      //together with Tree
                                      //BUT the only connection between parameters 
                                      //in Estimator.first(FunTreeIntensity), and Estimator.second(FitParList)
                                      //now is the same name of a Fit Parameter, no address/pointer connection
  auto MyMinuitIF = Optimizer::Minuit2::MinuitIF(UseHesse, UseMinos);
  
  LOG(INFO) << "Start minimization...";
  //Start minimization

  if (FitData) {
    LOG(INFO) << "FitData...";
    ComPWA::FitParameterList RandomFitParList;
    int ValidFitCounter = 0, CovergedFitCounter = 0;
    double MinFinalNLL(9999.99), CovergedNLL(9999.99);
    for (int IFit = 0; IFit < TotalFit; ++IFit) {
      if (CovergedFit <= 0 && ValidFitCounter > ValidFit) break;
      if (CovergedFitCounter > CovergedFit) break;
      //Generate Random initial Values
      RandomFitParList = genRandomPars(MyEstimator.second, RandomGenerator);
      //if (notRandom) RandomFitParList = MyEstimator.second;
      //auto Result = MyMinuitIF.optimize(MyEstimator.first, MyEstimator.second);
      auto Result = MyMinuitIF.optimize(MyEstimator.first, RandomFitParList);
      if (!Result.IsValid) continue;
      if (!Result.HasValidCov) continue;
      if (!Result.HasAccCov) continue;
      ++ValidFitCounter;
      if (CovergedFitCounter == 0) {
        CovergedNLL = Result.FinalEstimatorValue;
      } else {
        double NLL = Result.FinalEstimatorValue;
        if (abs(CovergedNLL - NLL) > 0.5) {
          if (CovergedNLL > NLL) {
            CovergedNLL = NLL;
            CovergedFitCounter = 0;
          }
          continue;
        }
      }
      ++CovergedFitCounter;
      if (MinFinalNLL <= Result.FinalEstimatorValue) continue;
      BestResult = Result;
      MinFinalNLL = Result.FinalEstimatorValue;
      LOG(INFO) << "Arxiv...";
      //Save fit result
      std::ofstream Ofs(FitResultXml, std::ios::out);
      boost::archive::xml_oarchive Oa(Ofs);
      Oa << BOOST_SERIALIZATION_NVP(Result);
      Ofs.close();
      LOG(INFO) << "Arxiv ok";
    }
    Fout << "MinNLL: " << MinFinalNLL << std::endl;
    Fout << "--------------------------" << std::endl;
    Fout << "Fit Result: " << std::endl;
    Fout << "Name\t\t" << "InitValue\t" << "FinalValue\t" << "+- "
        << "Error" << std::endl;
    for (std::size_t i = 0; i < BestResult.FinalParameters.size(); ++i) {
      double ErrLow = BestResult.FinalParameters[i].Error.first;
      double ErrUp = BestResult.FinalParameters[i].Error.second;
      Fout << BestResult.FinalParameters[i].Name << "\t"
          << BestResult.InitialParameters[i].Value << "\t"
          << BestResult.FinalParameters[i].Value << "\t";
      if (fabs(ErrLow - ErrUp) < 1e-6)
          Fout << " +- " << ErrLow;
      else 
          Fout << " (- " << ErrLow << ", + " << ErrUp << ")";
      Fout << std::endl;
    }
    Fout << "--------------------------" << std::endl;
  }

  if (!CalcFF && !CalcFFError) {
    BestResult.print(std::cout);
  }
  if (!CalcFF && !CalcFFError && !DrawPlot) {
    LOG(INFO) << "DONE";
    return 0;
  }

  LOG(INFO) << "Filter components";
  //Intensity Component a10, a11, a12, a32
  //Need IncoherentIntensity, and NO NormalizedIntensity
  auto AtotIntensity = Builder.createIntensity(PartList, HeliKine, 
      UnnormalizedIntensityTree); 
  auto A10Intensity = Builder.createIntensity(PartList, HeliKine, 
      UnnormalizedIntensityTree); 
  auto A11Intensity = Builder.createIntensity(PartList, HeliKine, 
      UnnormalizedIntensityTree); 
  auto A12Intensity = Builder.createIntensity(PartList, HeliKine, 
      UnnormalizedIntensityTree); 
  auto A32Intensity = Builder.createIntensity(PartList, HeliKine, 
      UnnormalizedIntensityTree); 
  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> TempTree;
  ComPWA::FunctionTree::ParameterList TempParList;
  //to get the ParameterList, and the function tree;
  //the TinySampleSet.Data is not used in below calc ff or draw plot
  std::tie(TempTree, TempParList) = AtotIntensity.bind(TinySampleSet.Data);
  //get the FitPars corresponding to the UnnormalziedItensity
  ComPWA::FitParameterList ComponentFitPars = ComPWA::FunctionTree
      ::createFitParameterList(TempParList);

  std::map<std::string, std::shared_ptr<ComPWA::Intensity>> MapComponent;
  //(sample, intensity, map<string, intensity> component)
  //once AmpPtr created by make_shared, it is not related with AmpIntnesity any more
  auto AtotPtr = std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensity>(
      AtotIntensity);
  auto A10Ptr = std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensity>(
      A10Intensity);
  auto A11Ptr = std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensity>(
      A11Intensity);
  auto A12Ptr = std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensity>(
      A12Intensity);
  auto A32Ptr = std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensity>(
      A32Intensity);
  MapComponent[A10] = A10Ptr;
  MapComponent[A11] = A11Ptr;
  MapComponent[A12] = A12Ptr;
  MapComponent[A32] = A32Ptr;

  //to update the FitPar values to the fit result
  if (LoadFitResult || FitData) {
    for (std::size_t i = 0; i < ComponentFitPars.size(); ++i) {
      std::string Name = ComponentFitPars[i].Name;
      bool Updated(false);
      for (std::size_t j = 0; j < BestResult.FinalParameters.size(); ++j) {
        auto FitPar = BestResult.FinalParameters[j];
        if (Name == FitPar.Name) {
          ComponentFitPars[i].Value = FitPar.Value;
          Updated = true;
          break;
        }
      }
      if (!Updated) {
        LOG(INFO) << "Paramter " << Name
            << " In ComponentIntensity is not updated! Exit!!!" << std::endl;
        return 0;
      }
    }
  }

  std::vector<double> DoubleParList(ComponentFitPars.size(), 0);
  filterIntensity("Atot", ComponentFitPars, DoubleParList);
  AtotPtr->forceUpdateParametersFrom(DoubleParList);
  filterIntensity("A10", ComponentFitPars, DoubleParList);
  A10Ptr->forceUpdateParametersFrom(DoubleParList);
  filterIntensity("A11", ComponentFitPars, DoubleParList);
  A11Ptr->forceUpdateParametersFrom(DoubleParList);
  filterIntensity("A12", ComponentFitPars, DoubleParList);
  A12Ptr->forceUpdateParametersFrom(DoubleParList);
  filterIntensity("A32", ComponentFitPars, DoubleParList);
  A32Ptr->forceUpdateParametersFrom(DoubleParList);

  //Plot data sample and intensity
  if (DrawPlot) {
    LOG(INFO) << "DrawPlot...";
    ComPWA::Tools::Plotting::RootPlotData
        PlotData(HeliKine.getParticleStateTransitionKinematicsInfo(),
        "plot_data.root", "RECREATE");
    PlotData.writeData(DataSampleSet);
    ComPWA::Tools::Plotting::RootPlotData
        PlotFit(HeliKine.getParticleStateTransitionKinematicsInfo(),
        "plot_fit.root", "RECREATE");
    ////the Intensity of writeIntensityWeightedPhspSample is not same to those in Physics
    PlotFit.writeIntensityWeightedPhspSample(PhspSampleSet,
        *AtotPtr, "fit", MapComponent);
  }

  ComPWA::Tools::FitFractionList FFList;
  if (CalcFF) {
    LOG(INFO) << "CalcFF...";
    FFList = ComPWA::Tools::calculateFitFractions(MapComponent,
        AtotPtr, PhspTrueSampleSet);
  }
  if (CalcFFError) {
    LOG(INFO) << "CalcFFError...";
    //std::vector<std::vector<double>> NXVec; //for debug
    //sampled parameters given CovarianceMatrix and Mean:
    // X = \mu + L Z, with CovM = L * L^T,
    // Z is n-dimension indpenent normal distributed variables
    //covariance matrix
    std::vector<std::vector<double>> CovMatrix = BestResult.CovarianceMatrix;
    std::size_t NFreePars = CovMatrix[0].size();

    //Fout << "Cov i\\j\t";
    //for (std::size_t i = 0; i < NFreePars; ++i)
    //  Fout << i << "\t";
    //Fout << std::endl;
    //for (std::size_t i = 0; i < NFreePars; ++i) {
    //  Fout << i << "\t";
    //  for (std::size_t j = 0; j < NFreePars; ++j) {
    //    Fout << CovMatrix[i][j] << "\t";
    //  }
    //  Fout << std::endl;
    //}
    std::vector<std::vector<double>> LMatrix(NFreePars,
        std::vector<double>(NFreePars, 0));
    
    //Convert to GSL_matrix
    gsl_matrix *GslCovMatrix = gsl_matrix_alloc(NFreePars, NFreePars);
    for (std::size_t row = 0; row < NFreePars; ++row) {
      for (std::size_t col = 0; col < NFreePars; ++col) {
        gsl_matrix_set(GslCovMatrix, row, col, CovMatrix[row][col]);
      }
    }
    //Cholesky decomposition to obtain L, CovM = L L^T
    //on output, lower triangle contains L, upper triangle unmodified;
    // lower triangle, M[i][j], i >= j 21 22 23 upper triangle, M[i][j], i < j
    int Status = gsl_linalg_cholesky_decomp(GslCovMatrix);
    if (Status == GSL_EDOM) {
      LOG(ERROR) << "Decomposition has failed!";
    }
    //Get LMatrix
    for (std::size_t row = 0; row < NFreePars; ++row) {
      for (std::size_t col = row; col < NFreePars; ++col) {
        LMatrix[col][row] = gsl_matrix_get(GslCovMatrix, col, row);
        if (col != row) LMatrix[row][col] = 0;
      }
    }
    //debug
    //Fout << "GSL_M after decomp" << std::endl;
    //Fout << "col\\row\t";
    //for (std::size_t row = 0; row < NFreePars; ++row)
    //  Fout << row << "\t";
    //Fout << std::endl;
    //for (std::size_t row = 0; row < NFreePars; ++row) {
    //    Fout << row << "\t";
    //  for (std::size_t col = 0; col < NFreePars; ++col) {
    //    Fout << gsl_matrix_get(GslCovMatrix, row, col) << "\t";
    //  }
    //  Fout << std::endl;
    //}
    //Fout << std::endl;
    //Fout << "LMatrix: " << std::endl;
    //Fout << "col\\row\t";
    //for (std::size_t row = 0; row < NFreePars; ++row)
    //  Fout << row << "\t";
    //Fout << std::endl;
    //for (std::size_t row = 0; row < NFreePars; ++row) {
    //  Fout << row << "\t";
    //  for (std::size_t col = 0; col < NFreePars; ++col) {
    //    Fout << LMatrix[row][col] << "\t";
    //  }
    //  Fout << std::endl;
    //}
    //Fout << std::endl;

    gsl_matrix_free(GslCovMatrix);
    TRandom3 Rand3(Seed);
    //NFree independent normal distributed random variables
    std::vector<double> ZVec(NFreePars, 0), XVec(NFreePars, 0);
    std::vector<double> MuVec(NFreePars, 0);
    std::vector<ComPWA::Tools::FitFractionList> FFListVec;
    ComPWA::Tools::FitFractionList FFListTemp;
    LOG(INFO) << "gen random pars and calc ffs...";
    std::vector<std::pair<double, double>> BoundsVec;
    for (std::size_t i = 0; i < ComponentFitPars.size(); ++i) {
      if (ComponentFitPars[i].IsFixed) continue;
      BoundsVec.push_back(ComponentFitPars[i].Bounds);
    }
    for (int iSet = 0; iSet < NParSet; ++iSet) {
      LOG(INFO) << "iSet = " << iSet;
      //Generate parameters with given covariance matrix
      for (std::size_t i = 0; i < NFreePars; ++i) {
        ZVec[i] = Rand3.Gaus(0, 1);//Normal distribution 
        XVec[i] = 0;
      }
      for (std::size_t i = 0, j = 0; i < ComponentFitPars.size(); ++i) {
        if (ComponentFitPars[i].IsFixed) continue;
        MuVec[j] = ComponentFitPars[i].Value;
        ++j; 
      }
      // X = L * Z + mu
      bool InBounds(true);
      for (std::size_t i = 0; i < NFreePars; ++i) {
        XVec[i] += MuVec[i];
        for (std::size_t j = 0; j < NFreePars; ++j) {
          XVec[i] += LMatrix[i][j] * ZVec[j]; 
        }
        //check if x in bounds
        if (XVec[i] <= BoundsVec[i].first || XVec[i] >= BoundsVec[i].second) {
          InBounds = false;
          break;
        }
      }
      if (!InBounds) {
        iSet--;
        continue;
      }
      //NXVec.push_back(XVec);
      //update parameters to intensity wrapper
      //auto ComponentNewFitPars = ComponentFitPars;
      ComPWA::FitParameterList ComponentNewFitPars(ComponentFitPars);
      for (std::size_t i = 0, j = 0; i < ComponentNewFitPars.size(); ++i) {
        if (ComponentNewFitPars[i].IsFixed) {
          DoubleParList[i] = ComponentFitPars[i].Value;
        } else {
          DoubleParList[i] = XVec[j];
          ComponentNewFitPars[i].Value = XVec[j];
          ++j;
        }
      }
      LOG(INFO) << "Random pars ok";
      AtotPtr->forceUpdateParametersFrom(DoubleParList);
      filterIntensity("A10", ComponentNewFitPars, DoubleParList);
      A10Ptr->forceUpdateParametersFrom(DoubleParList);
      filterIntensity("A11", ComponentNewFitPars, DoubleParList);
      A11Ptr->forceUpdateParametersFrom(DoubleParList);
      filterIntensity("A12", ComponentNewFitPars, DoubleParList);
      A12Ptr->forceUpdateParametersFrom(DoubleParList);
      filterIntensity("A32", ComponentNewFitPars, DoubleParList);
      A32Ptr->forceUpdateParametersFrom(DoubleParList);
      
      FFListTemp = ComPWA::Tools::calculateFitFractions(MapComponent,
          AtotPtr, PhspTrueSampleSet);
      LOG(INFO) << "iParSet = " << iSet << " this ff ok";
      FFListVec.push_back(FFListTemp);
    }
    //calculate standard deviation of the NParSet FF samples
    //I think maybe should divide NParSet - 1, not NParSet
    for (std::size_t jFF = 0; jFF < FFList.size(); ++jFF) {
      double Mean = 0, SqSum = 0, StdDev = 0;
      for (std::size_t iSet = 0; iSet < FFListVec.size(); ++iSet) {
        double Temp = FFListVec[iSet][jFF].Value;
        Mean += Temp;
      }
      Mean /= NParSet;
      for (std::size_t iSet = 0; iSet < FFListVec.size(); ++iSet) {
        double Temp = FFListVec[iSet][jFF].Value;
        SqSum += (Temp - Mean) * (Temp - Mean);
      }
      StdDev = std::sqrt(SqSum / NParSet);
      FFList[jFF].Error = StdDev;
    }
    LOG(INFO) << "FFError Setted" << std::endl;
    // debug, check sigma of sampled pars,
    //std::vector<double> MeanVec(NXVec[0].size(), 0);
    //std::vector<double> SigmaVec(NXVec[0].size(), 0);
    //std::vector<double> SumDiff2(NXVec[0].size(), 0);
    //for (std::size_t i = 0; i < NXVec.size(); ++i) {
    //  for (std::size_t j = 0; j < NXVec[i].size(); ++j) {
    //    MeanVec[j] += NXVec[i][j];  
    //  }
    //}
    //for (std::size_t j = 0; j < MeanVec.size(); ++j) 
    //  MeanVec[j] /= NParSet;
    //for (std::size_t i = 0; i < NXVec.size(); ++i) {
    //  for (std::size_t j = 0; j < NXVec[i].size(); ++j) {
    //    SumDiff2[j] += (NXVec[i][j] - MeanVec[j]) * (NXVec[i][j] - MeanVec[j]); 
    //  }
    //}
    //for (std::size_t j = 0; j < SigmaVec.size(); ++j) {
    //  SigmaVec[j] = std::sqrt(SumDiff2[j]/NParSet);
    //}
    //Fout << "Mean and Sigma of sampled pars: " << std::endl;
    //Fout << "No.\tMean\tSigma" << std::endl;
    //for (std::size_t j = 0; j < SigmaVec.size(); ++j) {
    //  Fout << j << "\t" << MeanVec[j] << "\t" << SigmaVec[j] << std::endl;
    //}
  }

  Fout << "===========================" << std::endl;
  Fout << "FitFraction:" << std::endl;
  Fout << "Component:\t\t" << "Value\t\t" << "Error" << std::endl;
  for (auto & FF : FFList) {
    Fout << FF.Name << "\t\t" << FF.Value << "\t\t" 
        << "+- " << FF.Error << std::endl;;
  }
  Fout << "===========================" << std::endl;


  LOG(INFO) << "DONE";
  return 0;
}

std::map<std::string, std::string> readConfig(const std::string fname) {
  std::map<std::string, std::string> MapConf;
  std::ifstream fin;
  fin.open(fname, std::ios::in);
  if (fin.fail()) {
    std::cout << "WARNING: Configure file not exists!" << std::endl;
    exit(1);
  }
  std::string Line;
  std::size_t IEq, ICm;
  std::string Key, Val, temp;
  while (getline(fin, Line)) {
    ICm = Line.find("#");
    IEq = Line.find("=");
    if (IEq == std::string::npos) continue; //no = in this line
    if (ICm != std::string::npos && IEq > ICm) continue; //e.g., #key = value
    std::istringstream iss(Line.substr(0, ICm));
    iss >> Key >> temp >> Val;
    if (Key.empty()) continue; //no key
    if (Val.empty()) continue; //no value
    MapConf[Key] = Val;
  }
  fin.close();
  return MapConf;
}
bool boolStr(const std::string val) {
  bool res(false);
  if (val == "1" || val == "Yes" || val == "YES"
      || val == "TRUE" || val == "True" || val == "true") {
    res = true;
  } else if (val == "0" || val == "NO" || val == "No" || val == "no"
      || val == "FALSE" || val == "False" || val == "false") {
    res = false;
  } else {
    std::cout << "WARNING: wrong bool flag" << std::endl;
    exit(1);
  }
  return res;
}
ComPWA::FitParameterList genRandomPars(ComPWA::FitParameterList &FitParList,
    ComPWA::Data::Root::RootUniformRealGenerator & RandomGen) {
  ComPWA::FitParameterList RandomFitParList(FitParList);
  for (std::size_t i = 0; i < FitParList.size(); ++i) {
    auto FitPar = FitParList[i];
    if (FitPar.IsFixed) {
      RandomFitParList[i].Value = FitPar.Value;
      continue;
    }
    double RandomVal = RandomGen(); //uniform in (0, 1)
    if (!FitPar.HasBounds) {
      RandomFitParList[i].Value = RandomVal;
      continue;
    }
    double Lower = FitPar.Bounds.first;
    double Upper = FitPar.Bounds.second;
    RandomFitParList[i].Value = Lower + RandomVal * (Upper - Lower);
    continue;
  }
  return RandomFitParList;
}
void filterIntensity(const std::string Amp,
    const ComPWA::FitParameterList & ParList,
    std::vector<double> &OutDoublePars) {
    //Amp = A10, A11, A12, A32
  for (std::size_t i = 0; i < ParList.size(); ++i) {
    auto Par = ParList[i];
    OutDoublePars[i] = Par.Value;
    std::string Name = Par.Name;
    //if fit parameter is not Magnituder or Phase, then continue
    if (Name.find("Magnitude") == std::string::npos
        && Name.find("Phase") == std::string::npos) continue;
    if (Amp == "Atot") continue;
    if (Amp == "A10" && 
        (Name.find(A10Mag) != std::string::npos 
        || Name.find(A10Phase) != std::string::npos))
      continue;
    if (Amp == "A11" && 
        (Name.find(A11Mag) != std::string::npos
        || Name.find(A11Phase) != std::string::npos))
      continue;
    if (Amp == "A12" && 
        (Name.find(A12Mag) != std::string::npos
        || Name.find(A12Phase) != std::string::npos))
      continue;
    if (Amp == "A32" && 
        (Name.find(A32Mag) != std::string::npos
        || Name.find(A32Phase) != std::string::npos))
      continue;
    //std::cout << "Set " << Amp << " " << Name << " To zero " << std::endl;
    OutDoublePars[i] = 0.0;
  }
}
void getValue(std::map<std::string, std::string> &MapConf,
    std::string Key, std::string &Out) {
  std::map<std::string, std::string>::iterator MapIter = MapConf.find(Key);
  if (MapIter != MapConf.end())
    Out = MapIter->second;
}
void getValue(std::map<std::string, std::string> &MapConf,
    std::string Key, int &Out) {
  std::map<std::string, std::string>::iterator 
      MapIter = MapConf.find(Key);
  if (MapIter != MapConf.end())
    Out = std::stoi(MapIter->second);
}
void getValue(std::map<std::string, std::string> &MapConf,
    std::string Key, double &Out) {
  std::map<std::string, std::string>::iterator 
      MapIter = MapConf.find(Key);
  if (MapIter != MapConf.end())
    Out = std::stof(MapIter->second);
}
void getValue(std::map<std::string, std::string> &MapConf,
    std::string Key, bool &Out) {
  std::map<std::string, std::string>::iterator 
      MapIter = MapConf.find(Key);
  if (MapIter != MapConf.end())
    Out = boolStr(MapIter->second);
}
