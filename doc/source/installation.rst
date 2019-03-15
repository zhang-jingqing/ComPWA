Installation
============

.. note::
   A detailed guide can be found below. The installation basically boils down to:

   * `git clone https://github.com/ComPWA/ComPWA.git <COMPWA_SOURCE_PATH>`
   * `cd <COMPWA_SOURCE_PATH> && git submodule init && git submodule update`
   * `mkdir build && cd build`
   * `cmake ../<COMPWA_SOURCE_PATH>`
   * `make`


Prerequisites
-------------

ComPWA is supposed to run on most modern unix systems (including MacOS). The following packages are mandatory:

* git (optional, for easier updates and if you want to contribute)
* cmake ( > 3.3 )
* gcc (> 5.1) or clang
* `Boost <http://www.boost.org/users/download/>`_\ , version >= 1.54

.. tip::
   For a more feature rich installation, the following packages are highly recommended:

   * python3 + virtualenv (for ComPWA expert system and python interface as well as a python plotting module)
   * `GSL <https://www.gnu.org/software/gsl/>`_
   * `ROOT <http://root.cern.ch/drupal/content/downloading-root>`_\ , version 5.34, 6.08
   * `Minuit2 <http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/>`_\ , version 5.34, 6.08
   * `Geneva <https://launchpad.net/geneva/+download>`_\ , version 1.8.0

In case that some dependencies are not met on your system use your package manager to install them or use the manual procedures described below. For MacOS you can use e.g. `MacPorts <https://www.macports.org>`_ as package manager.
You can also try to use different (newer) versions than the ones stated above, however those are not tested.

In order to install all dependencies it is probably also useful to have a look
on the `installation instructions file <https://github.com/ComPWA/ComPWA/blob/master/.travis.yml>`_ for TravisCI.


Manual installation of dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Boost**: to install Boost follow 
  `these <http://www.boost.org/doc/libs/1_54_0/more/getting_started/unix-variants.html#easy-build-and-install>`_ 
  instructions. The ``--prefix=path/to/installation/prefix`` option is useful
  as you can use the same path to point ComPWA to this specific Boost
  installation.

* **ROOT**: to install Root follow
  `these <http://root.cern.ch/drupal/content/installing-root-source>`_
  instructions.

* **Minuit2** is included in most ROOT installations. In case you want to
  install a stand-alone version follow
  `these <http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/gettingStarted/autoconf.html>`_
  instructions. In addition you should use
  ``./configure --prefix=/your/minuit2/installation/path`` followed by
  ``make install`` to specify an installation directory which ComPWA later
  needs to find the library and to install all needed files to this location.

* **Geneva**: to install Geneva follow 
  `these <http://www.gemfony.eu/index.php?id=genevainstallation>`_ 
  instructions. The ``DCMAKE_INSTALL_PREFIX="/where/the/library/goes"`` option
  is useful as you can use the same path to point ComPWA to this specific 
  Geneva installation:

  * Extract Geneva in a folder of your choice, lets call it "YOUR_GENEVA_PATH"
  * Go to ``YOUR_GENEVA_PATH/build``
  * Execute ``mkdir install``
  * Execute ``cmake ../ -DCMAKE_INSTALL_PREFIX="YOUR_GENEVA_PATH/build/install"``
  * Execute ``make``
  * Execute ``make install``
  * Execute ``cp install/CMakeModules/FindGeneva.cmake YOUR_COMPWA_PATH/cmake/Modules/``
  * Before compiling ComPWA, execute ``export GENEVA_ROOT=YOUR_GENEVA_PATH/build/install``
  * Note for Fedora 25: The Geneva tests are build by default but might have trouble finding the boost test libraries of the Fedora boost package. A workaround is to disable them within ``YOUR_GENEVA_PATH/CMakeModules/CommonGenevaBuild.cmake, line 55`` (replace the line with ``SET( GENEVA_BUILD_TESTS FALSE )``.
  * Alternatively you can follow the instructions from the Geneva `manual <http://www.gemfony.eu/fileadmin/documentation/geneva-manual.pdf>`_\ :

    * Go to ``/your/geneva/source/path/build``
    * ``cp ../scripts/genevaConfig.gcfg``
    * Modify ``./genevaConfig.gcfg`` to fit your needs.
    * ``../scripts/prepareBuild.sh ./genevaConfig.gcfg``

ComPWA installation
-------------------

Getting ComPWA
^^^^^^^^^^^^^^

To get the most recent version of the ComPWA framework clone its GitHub repository:

.. code-block:: shell

   git clone git@github.com:ComPWA/ComPWA <COMPWA_SOURCE_PATH>

this will clone the repository to the subfolder ``<COMPWA_SOURCE_PATH>`` within the current directory.
For multithreading ComPWA uses the parallel stl algorithms of c++17. Unfortunately the current compilers do not have any implementations for this. Here ComPWA currently relies on `TBB <https://github.com/01org/tbb>`_ and `parallelstl <https://github.com/intel/parallelstl>`_\ , which are included in ComPWA as git submodules. After the clone it is necessary to obtain these software modules via the commands

.. code-block:: shell

   cd <COMPWA_SOURCE_PATH>
   git submodule init && git submodule update

.. _setup-venv-label:

Setting up a python virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to use the python interface to ComPWA and/or use the python modules of ComPWA, setting up a virtual environment (venv) is highly recommended. Below are the setup instructions. Simply replace **$PATH_OF_YOUR_VENV** with the path where the venv should be installed.

.. code-block:: shell

   virtualenv -p python3 <PATH_OF_YOUR_VENV>
   source <PATH_OF_YOUR_VENV>/bin/activate
   pip install virtualenvwrapper

Now the virtual environment is set up. From now on, when you start up a new shell and want to work with ComPWA, just activate the venv with the command ``source <PATH_OF_YOUR_VENV>/bin/activate``. It can be deactivated with the command ``deactivate``.

.. _build-compwa-label:

Building ComPWA
^^^^^^^^^^^^^^^

ComPWA uses ``cmake`` as build system. The usual steps to build all libraries and the test executable are the following:


* Activate your python virtual environment (if you want python support) (recommended)
  .. code-block:: shell

       source <PATH_OF_YOUR_VENV>/bin/activate

* Create and enter a build folder (preferably not the ComPWA source folder)
  .. code-block:: shell

       mkdir build
       cd build

* Set your compiler if you do not use the system default compiler
  .. code-block::

       export CC=<path_to_your_compiler> 
       export CXX=<path_to_your_compiler>

* Build the project. You can add ``-DCMAKE_INSTALL_PREFIX=<COMPWA_INSTALL_PATH>`` to specify an install location.
  .. code-block:: shell

       cmake ../<COMPWA_SOURCE_PATH> 
       make
       make install (optional)


.. _finalize-venv-label:

Finalizing the python virtual environment for ComPWA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**IMPORTANT**: It is assumed that you have correctly set up a python virtual environment and **activated** it (see :ref:`Setting up a python virtual environment <setup-venv-label>`).


Install requirements for modules
""""""""""""""""""""""""""""""""
  Each python module of ComPWA contains a requirements.txt file. If you want to use this module simply install the requirements by executing:
  
  .. code-block:: shell

     pip install -r <PATH_TO_COMPWA_PYTHON_MODULE>/requirements.txt
    
  For example: ``pip install -r Physics/ExpertSystem/requirements.txt`` (assuming you are in the `<COMPWA_SOURCE_PATH>` directory)

Modifying the python search paths
"""""""""""""""""""""""""""""""""
  In order to use the ComPWA python modules, some search paths have to be added to python.  If you only called ``make`` and not ``make install``, execute these commands:

  .. code-block:: shell

     source virtualenvwrapper.sh
     add2virtualenv <COMPWA_SOURCE_PATH>/Physics/ExpertSystem
     add2virtualenv <COMPWA_SOURCE_PATH>/Tools
     add2virtualenv <COMPWA_BUILD_DIR>/Tools/PythonInterface

  Here `<COMPWA_SOURCE_PATH>` points to the ComPWA source directory and `<COMPWA_BUILD_DIR>` to the ComPWA build directory (where make was executed).
  
  If you installed ComPWA via make install

  .. code-block:: shell

     source virtualenvwrapper.sh
     add2virtualenv <COMPWA_INSTALL_PATH>/lib/python
     add2virtualenv <COMPWA_INSTALL_PATH>/lib

  Here `<COMPWA_INSTALL_PATH>` points to the install directory

Testing the ComPWA installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can run the test suite via:

.. code-block:: shell
   
   make test

or

.. code-block:: shell
   
   ctest

In case some python tests fail, make sure to install the requirements for these python modules of ComPWA (see :ref:`finalize python venv <finalize-venv-label>`)

Other
^^^^^

* You can also use cmake to create a preconfigured project for an IDE (e.g. `eclipse <https://www.eclipse.org>`_ ):

.. code-block:: shell

       cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_CXX_COMPILER_ARG1=-std=c++14 ../<COMPWA_SOURCE_PATH>

Installation via Docker
^^^^^^^^^^^^^^^^^^^^^^^

A `Dockerfile <https://github.com/ComPWA/ComPWA/blob/master/Dockerfile>`_ for ComPWA is provided. You can use it to build an `docker <https://www.docker.com>`_ image to run ComPWA. Using such an image ComPWA should run on `all systems that are supported by docker <https://docs.docker.com/engine/installation/>`_ including several (commercial) cloud computing services. If you are new to docker you can have a look on `this <https://prakhar.me/docker-curriculum/>`_ tutorial.

System specific notes
^^^^^^^^^^^^^^^^^^^^^

HimsterII / Mogon II
^^^^^^^^^^^^^^^^^^^^

`Mogon2 <https://hpc.uni-mainz.de/>`_ is the supercomputer of the Mainz University. If you work on it you can fulfill the ComPWA `installation requirements <#requirements>`_ by loading a series of modules:

.. code-block:: shell

   module load devel/CMake/3.9.5
   module load toolchain/foss/2017a
   module load devel/Boost/1.65.1-foss-2017a
   module load numlib/GSL/2.4-foss-2017a
   module load ROOT/v6.12-foss-2017a-python3
   export CC=/cluster/easybuild/broadwell/software/compiler/GCCcore/6.3.0/bin/gcc
   export CXX=/cluster/easybuild/broadwell/software/compiler/GCCcore/6.3.0/bin/g++

Now follow :ref:`the build instructions <build-compwa-label>`.

Troubleshooting
---------------

Add content here