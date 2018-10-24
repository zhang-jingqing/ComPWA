#include <algorithm>

#include "Core/Exceptions.hpp"
#include "Tools/Generate.hpp"

#include "ThirdParty/parallelstl/include/pstl/algorithm"
#include "ThirdParty/parallelstl/include/pstl/execution"

namespace ComPWA {
namespace Tools {

std::shared_ptr<ComPWA::Data::Data>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::AmpIntensity> Intensity) {
  std::shared_ptr<ComPWA::Data::Data> data(new ComPWA::Data::Data);
  if (NumberOfEvents <= 0)
    return data;
  // initialize generator output vector
  unsigned int EventBunchSize(5000);
  std::vector<ComPWA::Event> events;
  events.reserve(NumberOfEvents);

  double SafetyMargin(0.05);
  double generationMaxValue(0.0);
  unsigned int initialSeed = Generator->getSeed();

  std::vector<ComPWA::Event> tmp_events(EventBunchSize);
  std::vector<double> Intensities(tmp_events.size());
  std::vector<double> RandomNumbers(Intensities.size());

  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    // generate events
    std::generate(
        tmp_events.begin(), tmp_events.end(),
        [Generator]() -> ComPWA::Event { return Generator->generate(); });

    // evaluate function
    std::transform(
        pstl::execution::par_unseq, tmp_events.begin(), tmp_events.end(),
        Intensities.begin(),
        [Kinematics, Intensity](const ComPWA::Event &evt) -> double {
          ComPWA::DataPoint point;
          try {
            Kinematics->convert(evt, point);
          } catch (ComPWA::BeyondPhsp &ex) { // event outside phase, remove
            LOG(TRACE) << ex.what();
          }
          return evt.weight() * Intensity->intensity(point);
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
        events.back().setWeight(1.);
        bar.next();
        if (events.size() == NumberOfEvents)
          break;
      }
    }
    if (events.size() == NumberOfEvents)
      break;
  }
  data->add(events);
  return data;
}

std::shared_ptr<ComPWA::Data::Data>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::AmpIntensity> Intensity,
         std::shared_ptr<ComPWA::Data::Data> phsp,
         std::shared_ptr<ComPWA::Data::Data> phspTrue) {

  // Doing some checks
  if (NumberOfEvents <= 0)
    throw std::runtime_error("Tools::generate() negative number of events: " +
                             std::to_string(NumberOfEvents));
  if (!phsp)
    throw std::runtime_error("Tools::generate() No phsp sample given!");
  if (!Intensity)
    throw std::runtime_error("Tools::generate() | Amplitude not valid");
  if (!Generator)
    throw std::runtime_error("Tools::generate() | Generator not valid");
  if (phspTrue && !phsp)
    throw std::runtime_error("Tools::generate() | We have a sample of true"
                             " phsp events, but no phsp sample!");
  if (!phspTrue)
    phspTrue = phsp;
  if (phspTrue && phspTrue->numEvents() != phsp->numEvents())
    throw std::runtime_error(
        "Tools::generate() | We have a sample of true "
        "phsp events, but the sample size doesn't match that one of "
        "the phsp sample!");

  std::shared_ptr<ComPWA::Data::Data> data(new ComPWA::Data::Data);
  std::vector<ComPWA::Event> events;
  events.reserve(NumberOfEvents);

  double SafetyMargin(0.05);
  double maxSampleWeight(phsp->maximumWeight());
  if (phspTrue->maximumWeight() > maxSampleWeight)
    maxSampleWeight = phspTrue->maximumWeight();
  if (maxSampleWeight <= 0.0)
    throw std::runtime_error("Tools::generate() Sample maximum value is zero!");
  double generationMaxValue(maxSampleWeight * (1.0 + SafetyMargin));
  unsigned int initialSeed = Generator->getSeed();

  LOG(INFO) << "Tools::generate() | Using " << generationMaxValue
             << " as maximum value of the intensity.";

  unsigned int limit(phsp->numEvents());

  unsigned int EventBunchSize(5000);
  if (phsp->numEvents() < EventBunchSize)
    EventBunchSize = phsp->numEvents();
  std::vector<ComPWA::Event> TrueEventsBunch(EventBunchSize);
  std::vector<double> Intensities(TrueEventsBunch.size());
  std::vector<double> RandomNumbers(Intensities.size());

  auto CurrentStartIterator = phsp->events().begin();
  auto CurrentTrueStartIterator = phspTrue->events().begin();
  unsigned int CurrentStartIndex(0);
  LOG(INFO) << "Generating hit-and-miss sample: [" << NumberOfEvents
            << " events] ";
  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    if (CurrentStartIndex + EventBunchSize > limit)
      EventBunchSize = limit - CurrentStartIndex;

    // evaluate function
    std::transform(
        pstl::execution::par_unseq, CurrentTrueStartIterator,
        CurrentTrueStartIterator + EventBunchSize, Intensities.begin(),
        [Kinematics, Intensity](const ComPWA::Event &evt) -> double {
          ComPWA::DataPoint point;
          try {
            Kinematics->convert(evt, point);
          } catch (ComPWA::BeyondPhsp &ex) { // event outside phase, remove
            LOG(TRACE) << ex.what();
          }
          return Intensity->intensity(point);
        });

    // determine maximum
    double BunchMax(*std::max_element(pstl::execution::par_unseq,
                                      Intensities.begin(), Intensities.end()));
    // restart generation if we got above the current maximum
    if (BunchMax > generationMaxValue) {
      generationMaxValue = maxSampleWeight * (1.0 + SafetyMargin) * BunchMax;
      LOG(INFO) << "We raise the maximum to " << generationMaxValue;
      if (events.size() > 0) {
        events.clear();
        Generator->setSeed(initialSeed);
        CurrentStartIterator = phsp->events().begin();
        CurrentTrueStartIterator = phspTrue->events().begin();
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
      if (RandomNumbers[i] < CurrentStartIterator->weight() * Intensities[i]) {
        events.push_back(*CurrentStartIterator);
        events.back().setWeight(1.);
        events.back().setEfficiency(1.);
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
  data->add(events);
  return data;
}

} // namespace Tools
} // namespace ComPWA