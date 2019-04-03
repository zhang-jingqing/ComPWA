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

if( TARGET Core AND TARGET RootDataIO AND TARGET Minuit2IF AND TARGET MinLogLH
    AND TARGET Tools AND TARGET Integration AND TARGET HelicityFormalism
  )

file( GLOB lib_srcs Dst0Dst0.cpp )
#file( GLOB lib_headers ../../Physics/ParticleList.hpp )

add_executable ( Dst0Dst0
    Dst0Dst0.cpp ${lib_headers} )

target_link_libraries( Dst0Dst0
    Core 
    RootDataIO
    Minuit2IF
    MinLogLH
    Plotting
    Integration
    Tools
    HelicityFormalism
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
)

target_include_directories( Dst0Dst0
  PUBLIC 
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

install(TARGETS Dst0Dst0
    RUNTIME DESTINATION bin
)

else ()
  message( WARNING "Required targets not found! Not building\
                    Dst0Dst0 executable!")
endif()
