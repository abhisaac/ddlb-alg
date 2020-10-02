# Install Dependencies
sudo apt install libxml2-dev bison flex build-essential g++

NOTE: make sure g++ version >=9

# (Optional) Install omnetpp libraries
1. Download omnetpp core from https://omnetpp.org/download/
2. Extract the zip file
3. copy omnetpp/configure.user to the unzipped folder and overwrite the file
4. ./configure && make base -j4
5. copy lib/* to omnetp/lib
  
 
# Folder structure
```
.
├── bin
├── data
│   ├── bus_data
│   └── line_data
├── Makefile
├── obj
├── omnetpp	--------------------------------- omnet libraries/headers
├── README.txt
├── src
│   ├── bfs.cpp
│   ├── bfs.h
│   ├── constants.h
│   ├── distAlg.cc
│   ├── distAlg.h
│   ├── gen_fail.cpp
│   ├── gen_inps.cpp
│   ├── loadflow ---------------------------- loadflow code (generated from MATLAB)
│   ├── main.cc	----------------------------- our own main function to invoke the Simulation kernel
│   ├── original_ddlb_alg ------------------- original DDLB algorithm (without main)
│   │   ├── distAlg.cc
│   │   ├── distAlg.h
│   │   └── Makefile
│   └── util_timer.h
├── test ------------------------------------ test gen_inps, gen_fail and distAlg (with main)
│   ├── distAlg.ned
│   └── run.sh
└── test_origalg ---------------------------- test original alg
    ├── distAlg.ned
    └── runalg.sh
```
