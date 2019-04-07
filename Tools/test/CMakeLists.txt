
# Check if all requirements are found
if( ${ROOT_FOUND} AND TARGET Tools
)

add_executable(RootGeneratorTest RootGeneratorTest.cpp)
target_link_libraries(RootGeneratorTest
    Tools Boost::unit_test_framework
)
set_target_properties(RootGeneratorTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)
add_test(NAME RootGeneratorTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/RootGeneratorTest
)
endif()

# Check if all requirements are found
if(TARGET Tools)

add_executable(EvtGenGeneratorTest EvtGenGeneratorTest.cpp)
target_link_libraries(EvtGenGeneratorTest
    Tools
    Boost::unit_test_framework
)
set_target_properties(EvtGenGeneratorTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)

add_test(NAME EvtGenGeneratorTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/EvtGenGeneratorTest
)
endif()

if(TARGET Tools)
  # --------------------- ComponentIntensity ---------------------- #
  add_executable(ComponentIntensityTests ComponentIntensityTest.cpp)

  target_link_libraries(ComponentIntensityTests
    HelicityFormalism
    Tools
    Boost::unit_test_framework
  )

  target_include_directories(ComponentIntensityTests
    PUBLIC ${ROOT_INCLUDE_DIR} ${BOOST_INCLUDE_DIR}
  )

  # Move testing binaries into a testBin directory
  set_target_properties(ComponentIntensityTests
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  )

  add_test(NAME ComponentIntensityTests
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/ComponentIntensityTests)
endif()

if(TARGET Tools)
  # --------------------- FitFractions ---------------------- #
  add_executable(FitFractionsTests FitFractionsTest.cpp)

  target_link_libraries(FitFractionsTests
    HelicityFormalism
    Tools
    Integration
    Boost::unit_test_framework
  )

  target_include_directories(FitFractionsTests
    PUBLIC ${ROOT_INCLUDE_DIR} ${BOOST_INCLUDE_DIR}
  )

  # Move testing binaries into a testBin directory
  set_target_properties(FitFractionsTests
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  )

  add_test(NAME FitFractionsTests
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/FitFractionsTests)
endif()
