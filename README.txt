#Install Dependencies
sudo apt install libxml2-dev
  
#Folder structure
.
├── bin
├── data
│   ├── bus_data
│   └── line_data
├── Makefile
├── obj
├── omnetpp							- omnet libraries/headers
├── README.txt
├── src
│   ├── bfs.cpp
│   ├── bfs.h
│   ├── constants.h
│   ├── distAlg.cc
│   ├── distAlg.h
│   ├── gen_fail.cpp
│   ├── gen_inps.cpp
│   ├── loadflow					- loadflow code (generated from MATLAB)
│   ├── main.cc						- our own main function to invoke the Simulation kernel
│   ├── original_ddlb_alg        	- original DDLB algorithm (without main)
│   │   ├── distAlg.cc
│   │   ├── distAlg.h
│   │   └── Makefile
│   └── util_timer.h
├── test
│   ├── distAlg.ned
│   └── run.sh
└── test_origalg
    ├── distAlg.ned
    └── runalg.sh

