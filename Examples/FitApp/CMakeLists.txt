###########################################################
# Create executable for Dalitz plot fit of J/psi -> ypipi #
###########################################################

#
# FindBoost
#
find_package(Boost COMPONENTS
    program_options
    REQUIRED QUIET
)

if( TARGET Minuit2IF
    AND TARGET MinLogLH
    AND TARGET Tools
    AND TARGET Integration
    AND TARGET Plotting
    AND TARGET HelicityFormalism
    AND TARGET RootData
)

file( GLOB lib_srcs FitApp.cpp )
file( GLOB lib_headers ../../Physics/ParticleList.hpp FitApp.hpp )

add_executable ( FitApp
    FitApp.cpp ${lib_headers} )

target_link_libraries( FitApp
    Core
    Minuit2IF
    MinLogLH
    RootData
    HelicityFormalism
    Tools
    Integration
    Plotting
    ROOT::RIO
    Boost::program_options
    ${ROOT_LIBRARIES} ${Boot_LIBRARIES}
    #Threads::Threads
    GSL::gsl
    #GSL::gslcblas
)

target_include_directories( FitApp
  PUBLIC 
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

install(TARGETS FitApp
    RUNTIME DESTINATION bin
)

else ()
  message( WARNING "Required targets not found! Not building\
                    FitApp executable!")
endif()
