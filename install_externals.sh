#!/bin/bash

#set to the install path of your choice and remove the return statement below
#make sure that you have write access in that directory
install_path=$1

source ./checker.sh
source ./install_dependencies.sh

# wget, tar, cmake 2.8+ and a newer gcc (c++11 support) are required for this installation script
# these dependencies are checked here... of course the dependency tree is not checked any further...
check_install_dir
check_dependencies
check_compiler

#the temporary install dir is just the current path
cd ${install_path}
mkdir temp
temp_workdir=${install_path}/temp

# determine number of processors (use 1 less for building)
if [ ! $numproc ]; then
  numproc=`less /proc/cpuinfo | grep processor | wc -l`
  if [ "$numproc" -gt 1 ]; then
    numproc=$(($numproc-1));
  fi
fi

install_boost
install_geneva
install_root
install_qftpp
install_minuit2

cd ${install_path}
rm -rf temp

#export paths (this should be done in a config script of compwa or the users bashrc etc)
echo "Installation is finished! Please add the following path to your environment variables:"
echo "export PATH=${boost_root}/bin:$ROOTSYS/bin:\$PATH"
echo "export LD_LIBRARY_PATH=${geneva_root}/lib:${boost_root}/lib:$ROOTSYS/lib:\$LD_LIBRARY_PATH"