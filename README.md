Time Bombs
==========

Inferring time delays from a single light curve, assuming a flaring mixture model. 

### Getting Started

Dependencies include:

* DNest3 (http://github.com/eggplantbren/DNest3)
* RJObject (http://github.com/eggplantbren/RJObject)

Compile with:

    make

Then run with:

    ./main -t 8
The option causes the program to run with 8 threads.

To postprocess results and create the file posterior_sample.txt:

    python2 showresults.py

#### Installation Example

    Building DNest3:

    sudo yum install cmake.x86_64
    sudo yum install gsl-devel.x86_64
    sudo yum install boost-devel.x86_64
    git clone git://github.com/eggplantbren/DNest3.git
    cd DNest3
    <edit the CMakeLists.txt file to set the required cmake version 
     consistent with the yum-installed version (cmake --version)
    >
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX:PATH=/u2/jchiang/local ..
    make
    make install
    cd ..

    export LD_LIBRARY_PATH=/u2/jchiang/local/lib:${LD_LIBRARY_PATH}

    Building RJObject:
    git clone git://github.com/eggplantbren/RJObject.git
    cd RJObject
    <edit Makefile to add path to libdnest3.so and use libboost_thread-mt:
     LIBS = -L/u2/jchiang/local/lib -ldnest3 -lgsl -lgslcblas -lboost_thread-mt -lboost_system 
    >
    <edit the following files to prepend "dnest3" on the relevant includes:
      Distributions/BasicCircular.cpp
      Distributions/ClassicMassInf.cpp
      Distributions/ClassicMassInf1D.cpp
      Distributions/Pareto.cpp
      RJObjectImpl.h
      main.cpp
     e.g., 
      diff --git a/RJObjectImpl.h b/RJObjectImpl.h
      index 0e1492c..69adf1c 100644
      --- a/RJObjectImpl.h
      +++ b/RJObjectImpl.h
      @@ -1,8 +1,8 @@
      #include <cmath>
      #include <iostream>
      #include <iomanip>
     -#include "RandomNumberGenerator.h"
     -#include "Utils.h"
     +#include "dnest3/RandomNumberGenerator.h"
     +#include "dnest3/Utils.h"
 
      template<class Distribution>
      RJObject<Distribution>::RJObject(int num_dimensions, 
    >
    make
    cd ..

    Building TimeBombs:
    git clone git://github.com/jchiang87/TimeBombs.git
    cd TimeBombs
    export RJOBJECT_DIR=../RJObject
    <edit the Makefile to add paths to libdnest3.so and to use libboost_thread-mt:
      CFLAGS = -O2 -Wall -Wextra -ansi -pedantic -I${RJOBJECT_DIR} -I/u2/jchiang/local/include
      LIBS = -L${RJOBJECT_DIR} -lrjobject -L/u2/jchiang/local/lib -ldnest3 \
             -lgsl -lgslcblas \
             -lboost_system -lboost_thread-mt
     >
     export LD_LIBRARY_PATH=${RJOBJECT_DIR}:${LD_LIBRARY_PATH}

### People, Licence, Credits etc

Currently working on this project:

* Brendon Brewer (Auckland)
* Phil Marshall (KIPAC)

Feel free to send us email, tweets, issues, pull requests etc. The code is still under development and not yet written up, so beware. It'd be lovely to hear from you though! 

(c) 2014 Brendon J. Brewer, Philip J. Marshall

LICENCE: GNU GPL v3.

