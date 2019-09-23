#ifndef FIT_APP_HPP
#define FIT_APP_HPP

#include <string>
#include <vector>

#include "Core/FunctionTree/FitParameter.hpp"
#include "Core/FunctionTree/ParameterList.hpp"
#include "Data/Root/RootGenerator.hpp"

const std::string A10("A10");
const std::string A11("A11");
const std::string A12("A12");
const std::string A32("A32");
std::string A10Mag;
std::string A10Phase;
std::string A11Mag;
std::string A11Phase;
std::string A12Mag;
std::string A12Phase;
std::string A32Mag;
std::string A32Phase;

std::map<std::string, std::string> readConfig(const std::string fname);
void getValue(std::map<std::string, std::string> &MapConf,
    std::string Key, std::string &Out);
void getValue(std::map<std::string, std::string> &MapConf,
    std::string Key, int &Out);
void getValue(std::map<std::string, std::string> &MapConf,
    std::string Key, double &Out);
void getValue(std::map<std::string, std::string> &MapConf,
    std::string Key, bool &Out);
bool boolStr(const std::string val);
ComPWA::FitParameterList genRandomPars(ComPWA::FitParameterList &FitParList,
    ComPWA::Data::Root::RootUniformRealGenerator & RandomGen);
void filterIntensity(const std::string Amp,
    const ComPWA::FitParameterList & ParList,
    std::vector<double> &OutDoublePars);
//sampling to calc fit errors
//https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
//https://stats.stackexchange.com/questions/120179/generating-data-with-a-given-sample-covariance-matrix/120227#120227
//https://stats.stackexchange.com/questions/160054/how-to-use-the-cholesky-decomposition-or-an-alternative-for-correlated-data-si?noredirect=1&lq=1
#endif
