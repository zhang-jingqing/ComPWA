#include <algorithm>

#include "Core/Exceptions.hpp"
#include "Core/Generator.hpp"
#include "Core/Intensity.hpp"
#include "Core/Kinematics.hpp"
#include "Core/ProgressBar.hpp"
#include "Data/DataSet.hpp"
#include "Tools/Generate.hpp"

#include "ThirdParty/parallelstl/include/pstl/algorithm"
#include "ThirdParty/parallelstl/include/pstl/execution"

namespace ComPWA {
namespace Tools {

std::shared_ptr<ComPWA::Data::DataSet>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::Intensity> Intensity) {
  std::vector<ComPWA::Event> events;
  if (NumberOfEvents <= 0)
    return std::make_shared<ComPWA::Data::DataSet>(events);
  // initialize generator output vector
  unsigned int EventBunchSize(5000);
  events.reserve(NumberOfEvents);

  double SafetyMargin(0.05);
  double generationMaxValue(0.0);
  unsigned int initialSeed = Generator->getSeed();

  std::vector<ComPWA::Event> tmp_events(EventBunchSize);
  std::vector<double> Intensities(tmp_events.size());
  std::vector<double> RandomNumbers(Intensities.size());

  LOG(INFO) << "Generating hit-and-miss sample: [" << NumberOfEvents
            << " events] ";
  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    // generate events
    std::generate(
        tmp_events.begin(), tmp_events.end(),
        [Generator]() -> ComPWA::Event { return Generator->generate(); });

    // evaluate function
    // Note: some event generators create events outside of the phase space
    // boundary (due to numerical instability and precision). These events have
    // to be ignored!
    std::transform(pstl::execution::par_unseq, tmp_events.begin(),
                   tmp_events.end(), Intensities.begin(),
                   [Kinematics, Intensity](const ComPWA::Event &evt) -> double {
                     ComPWA::DataPoint point(Kinematics->convert(evt));
                     if (!Kinematics->isWithinPhaseSpace(point))
                       return 0.0;
                     return evt.Weight * Intensity->evaluate(point);
                   });
    // determine maximum
    double BunchMax(*std::max_element(pstl::execution::par_unseq,
                                      Intensities.begin(), Intensities.end()));
    // restart generation if we got above the current maximum
    if (BunchMax > generationMaxValue) {
      generationMaxValue = (1.0 + SafetyMargin) * BunchMax;
      if (events.size() > 0) {
        events.clear();
        Generator->setSeed(initialSeed);
        bar = ComPWA::ProgressBar(NumberOfEvents);
        LOG(INFO) << "Tools::generate() | Error in HitMiss "
                     "procedure: Maximum value of random number generation "
                     "smaller then amplitude maximum! We raise the maximum "
                     "to "
                  << generationMaxValue << " value and restart generation!";
        continue;
      }
    }
    // do hit and miss
    // first generate random numbers (no multithreading here, to ensure
    // deterministic behavior independent on the number of threads)
    std::generate(RandomNumbers.begin(), RandomNumbers.end(),
                  [Generator, generationMaxValue]() -> double {
                    return Generator->uniform(0, generationMaxValue);
                  });

    for (unsigned int i = 0; i < tmp_events.size(); ++i) {
      if (RandomNumbers[i] < Intensities[i]) {
        events.push_back(tmp_events[i]);
        events.back().Weight = 1.0;
        bar.next();
        if (events.size() == NumberOfEvents)
          break;
      }
    }
    if (events.size() == NumberOfEvents)
      break;
  }
  auto DataSample = std::make_shared<ComPWA::Data::DataSet>(events);
  DataSample->convertEventsToDataPoints(Kinematics);
  return DataSample;
}

std::shared_ptr<ComPWA::Data::DataSet>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::Intensity> Intensity,
         std::shared_ptr<ComPWA::Data::DataSet> phsp,
         std::shared_ptr<ComPWA::Data::DataSet> phspTrue) {

  // Doing some checks
  if (NumberOfEvents <= 0)
    throw std::runtime_error("Tools::generate() negative number of events: " +
                             std::to_string(NumberOfEvents));
  if (!Intensity)
    throw std::runtime_error("Tools::generate() | Amplitude not valid");
  if (!Generator)
    throw std::runtime_error("Tools::generate() | Generator not valid");
  if (!phsp)
    throw std::runtime_error("Tools::generate() | No phase space sample given");
  if (phspTrue &&
      phspTrue->getEventList().size() != phsp->getEventList().size())
    throw std::runtime_error(
        "Tools::generate() | We have a sample of true "
        "phsp events, but the sample size doesn't match that one of "
        "the phsp sample!");

  std::vector<ComPWA::Event> events;
  events.reserve(NumberOfEvents);

  double SafetyMargin(0.05);

  double maxSampleWeight(ComPWA::getMaximumSampleWeight(phsp->getEventList()));
  if (phspTrue) {
    double temp_maxweight(
        ComPWA::getMaximumSampleWeight(phspTrue->getEventList()));
    if (temp_maxweight > maxSampleWeight)
      maxSampleWeight = temp_maxweight;
  }

  if (maxSampleWeight <= 0.0)
    throw std::runtime_error("Tools::generate() Sample maximum value is zero!");
  double generationMaxValue(maxSampleWeight * (1.0 + SafetyMargin));
  unsigned int initialSeed = Generator->getSeed();

  LOG(INFO) << "Tools::generate() | Using " << generationMaxValue
            << " as maximum value of the intensity.";

  auto const &PhspEvents = phsp->getEventList();
  unsigned int limit(PhspEvents.size());

  unsigned int EventBunchSize(5000);
  if (PhspEvents.size() < EventBunchSize)
    EventBunchSize = PhspEvents.size();
  std::vector<ComPWA::Event> TrueEventsBunch(EventBunchSize);
  std::vector<double> Intensities(TrueEventsBunch.size());
  std::vector<double> RandomNumbers(Intensities.size());

  auto CurrentStartIterator = PhspEvents.begin();
  auto CurrentTrueStartIterator = PhspEvents.begin();
  if (phspTrue)
    CurrentTrueStartIterator = phspTrue->getEventList().begin();
  unsigned int CurrentStartIndex(0);
  LOG(INFO) << "Generating hit-and-miss sample: [" << NumberOfEvents
            << " events] ";
  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    if (CurrentStartIndex + EventBunchSize > limit)
      EventBunchSize = limit - CurrentStartIndex;

    // evaluate function
    // Note: some event generators create events outside of the phase space
    // boundary (due to numerical instability and precision). These events have
    // to be ignored!
    std::transform(pstl::execution::par_unseq, CurrentTrueStartIterator,
                   CurrentTrueStartIterator + EventBunchSize,
                   Intensities.begin(),
                   [Kinematics, Intensity](const ComPWA::Event &evt) -> double {
                     ComPWA::DataPoint point(Kinematics->convert(evt));
                     if (!Kinematics->isWithinPhaseSpace(point))
                       return 0.0;
                     return Intensity->evaluate(point);
                   });

    // determine maximum
    double BunchMax(*std::max_element(pstl::execution::par_unseq,
                                      Intensities.begin(), Intensities.end()));
    // restart generation if we got above the current maximum
    if (maxSampleWeight * BunchMax > generationMaxValue) {
      generationMaxValue = maxSampleWeight * (1.0 + SafetyMargin) * BunchMax;
      LOG(INFO) << "We raise the maximum to " << generationMaxValue;
      if (events.size() > 0) {
        events.clear();
        Generator->setSeed(initialSeed);
        CurrentStartIterator = PhspEvents.begin();
        CurrentTrueStartIterator = PhspEvents.begin();
        if (phspTrue)
          CurrentTrueStartIterator = phspTrue->getEventList().begin();
        CurrentStartIndex = 0;
        bar = ComPWA::ProgressBar(NumberOfEvents);
        LOG(INFO) << "Tools::generate() | Error in HitMiss "
                     "procedure: Maximum value of random number generation "
                     "smaller then amplitude maximum! Restarting generation!";
      }
      continue;
    }

    // do hit and miss
    // first generate random numbers (no multithreading here, to ensure
    // deterministic behavior independent on the number of threads)
    std::generate(RandomNumbers.begin(), RandomNumbers.end(),
                  [Generator, generationMaxValue]() -> double {
                    return Generator->uniform(0, generationMaxValue);
                  });

    for (unsigned int i = 0; i < TrueEventsBunch.size(); ++i) {
      if (RandomNumbers[i] < CurrentStartIterator->Weight * Intensities[i]) {
        events.push_back(*CurrentStartIterator);
        events.back().Weight = 1.0;
        bar.next();
        if (events.size() == NumberOfEvents)
          break;
      }
      ++CurrentStartIterator;
    }

    // increment true iterator
    std::advance(CurrentTrueStartIterator, EventBunchSize);
    CurrentStartIndex += EventBunchSize;

    if (events.size() == NumberOfEvents)
      break;

    if (CurrentStartIndex >= limit)
      break;
  }
  double gen_eff = (double)events.size() / NumberOfEvents;
  if (CurrentStartIndex > NumberOfEvents) {
    gen_eff = (double)events.size() / CurrentStartIndex;
  }
  LOG(INFO) << "Efficiency of toy MC generation: " << gen_eff;

  auto DataSample = std::make_shared<ComPWA::Data::DataSet>(events);
  DataSample->convertEventsToDataPoints(Kinematics);
  return DataSample;
}

std::shared_ptr<ComPWA::Data::DataSet>
generatePhsp(unsigned int nEvents, std::shared_ptr<ComPWA::Generator> gen) {
  std::vector<ComPWA::Event> sample;

  LOG(INFO) << "Generating phase-space MC: [" << nEvents << " events] ";

  ComPWA::ProgressBar bar(nEvents);
  for (unsigned int i = 0; i < nEvents; ++i) {
    ComPWA::Event tmp = gen->generate();
    double ampRnd = gen->uniform(0, 1);
    if (ampRnd > tmp.Weight) {
      --i;
      continue;
    }

    // Reset weights: weights are taken into account by hit&miss. The
    // resulting sample is therefore unweighted
    tmp.Weight = 1.0;

    sample.push_back(tmp);
    bar.next();
  }
  return std::make_shared<ComPWA::Data::DataSet>(sample);
}

std::shared_ptr<ComPWA::Data::DataSet>
generateImportanceSampledPhsp(unsigned int NumberOfEvents,
                              std::shared_ptr<ComPWA::Kinematics> Kinematics,
                              std::shared_ptr<ComPWA::Generator> Generator,
                              std::shared_ptr<ComPWA::Intensity> Intensity) {
  std::vector<ComPWA::Event> events;
  if (NumberOfEvents <= 0)
    return std::make_shared<ComPWA::Data::DataSet>(events);
  // initialize generator output vector
  unsigned int EventBunchSize(5000);
  events.reserve(NumberOfEvents);

  double SafetyMargin(0.05);
  double generationMaxValue(0.0);
  unsigned int initialSeed = Generator->getSeed();
  double WeightSum(0.0);

  std::vector<ComPWA::Event> tmp_events(EventBunchSize);
  std::vector<double> Intensities(tmp_events.size());
  std::vector<double> RandomNumbers(Intensities.size());
  LOG(INFO)
      << "Generating phase space sample (hit-and-miss importance sampled): ["
      << NumberOfEvents << " events] ";
  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    // generate events
    std::generate(
        tmp_events.begin(), tmp_events.end(),
        [Generator]() -> ComPWA::Event { return Generator->generate(); });

    // evaluate function
    // Note: some event generators create events outside of the phase space
    // boundary (due to numerical instability and precision). These events have
    // to be ignored!
    std::transform(pstl::execution::par_unseq, tmp_events.begin(),
                   tmp_events.end(), Intensities.begin(),
                   [Kinematics, Intensity](const ComPWA::Event &evt) -> double {
                     ComPWA::DataPoint point(Kinematics->convert(evt));
                     if (!Kinematics->isWithinPhaseSpace(point))
                       return 0.0;
                     return evt.Weight * Intensity->evaluate(point);
                   });
    // determine maximum
    double BunchMax(*std::max_element(pstl::execution::par_unseq,
                                      Intensities.begin(), Intensities.end()));
    // restart generation if we got above the current maximum
    if (BunchMax > generationMaxValue) {
      generationMaxValue = (1.0 + SafetyMargin) * BunchMax;
      if (events.size() > 0) {
        events.clear();
        WeightSum = 0.0;
        Generator->setSeed(initialSeed);
        bar = ComPWA::ProgressBar(NumberOfEvents);
        LOG(INFO)
            << "Tools::generateImportanceSampledPhsp() | Error in HitMiss "
               "procedure: Maximum value of random number generation "
               "smaller then amplitude maximum! We raise the maximum "
               "to "
            << generationMaxValue << " value and restart generation!";
        continue;
      }
    }
    // do hit and miss
    // first generate random numbers (no multithreading here, to ensure
    // deterministic behavior independent on the number of threads)
    std::generate(RandomNumbers.begin(), RandomNumbers.end(),
                  [Generator, generationMaxValue]() -> double {
                    return Generator->uniform(0, generationMaxValue);
                  });

    for (unsigned int i = 0; i < tmp_events.size(); ++i) {
      if (RandomNumbers[i] < Intensities[i]) {
        events.push_back(tmp_events[i]);
        double weight(tmp_events[i].Weight / Intensities[i]);
        events.back().Weight = weight;
        WeightSum += weight;
        bar.next();
        if (events.size() == NumberOfEvents)
          break;
      }
    }
    if (events.size() == NumberOfEvents)
      break;
  }
  // now just rescale the event weights so that sum(event weights) = # events
  double rescale_factor(NumberOfEvents / WeightSum);
  for (auto &evt : events) {
    evt.Weight = evt.Weight * rescale_factor;
  }

  auto DataSample = std::make_shared<ComPWA::Data::DataSet>(events);
  DataSample->convertEventsToDataPoints(Kinematics);
  return DataSample;
}

} // namespace Tools
} // namespace ComPWA
