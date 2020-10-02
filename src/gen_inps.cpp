#include "constants.h"

#include <sys/stat.h>
#include <cstdlib>
#include <cstdio>
#include <limits>

#include <cstring>
#include <cmath>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <array>
#include <vector>
#include <tuple>
#include <algorithm>
#include <map>
#include <random>
#include <filesystem>
#include <thread>
#include <mutex>

#include "util_timer.h"


typedef unsigned int uint;
using namespace std;
const uint SETWIDTH   = 14;

mutex cout_lock;
#define trace(x) { scoped_lock<mutex> lock(cout_lock); cout << x << endl; }

extern "C" {
	#include "loadflow/matlabrun.h"
}

const std::string currentDateTime() {
    time_t  now = time(0);

    char       buf[80];
#ifdef __GNUC__
    tm* tstruct = localtime(&now);
    strftime(buf, sizeof(buf), "%Y%m%d", tstruct);
#else
    struct tm  tstruct;
    localtime_s(&tstruct, &now);
    strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);
#endif

    return buf;
}

template<typename T>
bool is_infinite( const T &value )
{
    T max_value = (std::numeric_limits<T>::max)();
    T min_value = - max_value;

    return ! ( min_value <= value && value <= max_value );
}

template<typename T>
bool is_nan( const T &value )
{
    return value != value;
}

template<typename T>
bool is_valid( const T &value )
{
    return (!is_infinite(value) && !is_nan(value) );
}


struct BusVals
{
    int id, type;
    double voltage, activePower, startPower;

    BusVals():  id(0), type(0), voltage(0.0), activePower(0.0), startPower(0.0) {}
};

void printToFile(int BaseInputNo,int TopoNo,BusVals* askedVals);

enum Case {normal, overloaded, underloaded};
//typedef BusVals BusMatrix[BUSNUM];

typedef double LoadFlowArray[BUSNUM*10];
typedef LoadFlowArray LoadFlowInputArray;
typedef LoadFlowArray LoadFlowOutputArray;
typedef double LineArray[LINENUM*7 ];

inline bool exists_test(const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

double fRand(const double fMin, const double fMax){  //does not return 0
    double returnVal,f;
    do{
        f = (double)rand() / RAND_MAX;
        returnVal=fMin + f * (fMax - fMin);
    }while(returnVal == 0.0);
    return returnVal;
}

LineArray linefull;
int getLinevals(const std::string input){
    std::ifstream ifs(input);
    if(!ifs.good()){
        cout << "            ################     error opening line_values_file      ###############" << '\n';
        return 1;
    }

    std::string l1;
    int num=0;
    while(std::getline(ifs, l1)){ // read one line from ifs
        std::istringstream iss(l1); // access line as a stream
        iss >> linefull[num] >> linefull[LINENUM + num] >> linefull[LINENUM*2 + num]
                >> linefull[LINENUM*3 + num] >> linefull[LINENUM*4 + num]
                >> linefull[LINENUM*5 + num] >> linefull[LINENUM*6 + num] ;
        ++num;
    }
    return 0;
}
int getBusvals(const std::string input, BusVals* inputBus){
    std::ifstream ifs(input);
    if(!ifs.good()){
        cout << "            ################     error opening bus_values_file      ###############" << '\n';
        return 1;
    }

    std::string line;
    int num=0; double dummy;
    while(std::getline(ifs, line)){ // read one line from ifs
        std::istringstream iss(line); // access line as a stream
        iss >> inputBus[num].id >> inputBus[num].voltage >> dummy >> dummy >> inputBus[num].activePower
                >> dummy >> dummy >> inputBus[num].type >> inputBus[num].startPower;
        ++num;
    }
    return 0;
}

struct Network {
    LoadFlowInputArray bus_startvals;
    LoadFlowOutputArray bus_sol2;

    double bus[BUSNUM*10];double bus_sol[BUSNUM*10];
    uint    no_Of_Gens      =0,
    no_Of_slack     =0,
    slack_bus_ID    =0,
    Asked_invalid   =0,
    Start_invalid   =0,
    iter            =0,
    iter2           =0,
    _percentage     =20,
    Global_Counter  =0,
    Stats_valid     =0,
    Stats_invalid   =0,
    iter_more       =0;

    double  asked_power[BUSNUM],
    start_power[BUSNUM],
    type[BUSNUM],
    id[BUSNUM],
    voltage[BUSNUM];

    std::stringstream outputStream;
    std::stringstream fileName;
    string directory;

    std::vector<int> myvector;

    double  loadSlack,
    voltageSlack,
    _actpower,
    _reactpower,
    Scale,
    maxvolt=0,
    sum=0,
    vslack;

    bool badValue = false;

    Case checkCase;

    double _normOrOvldFactor    =1;
    int itercheck               =0;
    double absSum               =0;

    uint noOfLoadCorrections,
    noOfGenCorrections;
    std::vector<float> skewArray;
    uint perturbType;
    double skew;
    double endScale;
    uint noOfTopologies,noOfBaseInputs,noOfCasesEach;

    string todayDateTime;

    void findnoofGensAndSlack(int len, int seek,BusVals* inputVals){
        no_Of_slack =0 ;

        for (int i = 0; i < len; ++i){

            if (inputVals[i].type == seek) {
                if(seek == 1){
                    ++no_Of_slack;

                }else if (seek == 2){
                    ++no_Of_Gens;
                }
            }
        }
    }



    bool verify(LoadFlowArray & bus_sol,uint iter,bool startValue){
        double maxvolt =0;
        bool valid=true;
        double volt=0;
        for(uint i=0;i<BUSNUM;++i){

            volt = bus_sol[BUSNUM+i];

            if(fabs(maxvolt) < fabs(volt))
                maxvolt = volt;

            if(!is_valid(bus_sol[BUSNUM*5+i] - bus_sol[BUSNUM*3+i]))
                valid = false;
        }

        if(startValue && ( fabs(maxvolt-1) > 0.05 ) )
            valid=false;

        if(iter > 5){
            ++iter_more;
            return false;
        }

        if(valid==false){
            if(!startValue)
                ++Asked_invalid;
            else
                ++Start_invalid;
            return false;
        }

        return true;
    }



    int convertToMatlabInputAndRun(double* bus2,
                                   BusVals* b, double* bus_sol2){

        double _rpowerPercentage = .20;

        for(uint i=0;i<BUSNUM;i++){
            bus2[i]          = i+1; //convert 0 base of C++ array to 1 base of Matlab**
            bus2[BUSNUM + i]    = b[i].voltage;
            bus2[BUSNUM*9 + i]  = b[i].type;
            bus2[BUSNUM*7 + i] =
                    bus2[BUSNUM*8 + i] =
                    bus2[BUSNUM*2 + i] = 0;

            if (int(b[i].type)  == 3){
                bus2[BUSNUM*5 + i]  =   b[i].activePower;
                bus2[BUSNUM*6 + i]  =   bus2[BUSNUM*5 + i]*_rpowerPercentage;
                bus2[BUSNUM*3 + i]  =   0;
                bus2[BUSNUM*4 + i]  =   0;
            }else{
                bus2[BUSNUM*3 + i]  = - b[i].activePower;
                bus2[BUSNUM*4 + i]  =   bus2[BUSNUM*3 + i]*_rpowerPercentage;
                bus2[BUSNUM*5 + i]  =   0;
                bus2[BUSNUM*6 + i]  =   0;
            }
        }
        return matlabrun(bus2,linefull,no_Of_Gens,
                         slack_bus_ID+1,no_Of_slack,bus_sol2) ;

    }


    void perturb_and_assign_end(const double& Scale,
                                const uint& perturbType,
                                const double& skew,
                                const double& absSum,
                                const uint&noOfLoadCorrections,
                                const uint&noOfGenCorrections,
                                BusVals* startVals,
                                BusVals* askedVals)
    {

        badValue = false;
        double perturbLoad= 0 ,
                perturbGen =0;
        double sumLoads=0,sumGens=0;
        double minValue = 5 ;

        if(perturbType == 3){  // perturb both loads and gens equally
            perturbLoad = (skew*absSum)/(2*noOfLoadCorrections);
            perturbGen 	=  (skew*absSum)/(2*noOfGenCorrections);
        }
        else if(perturbType == 2){  // perturb only Loads
            perturbLoad = (skew*absSum)/(noOfLoadCorrections);
            perturbGen 	=  0;
        }
        else{   // perturb only Gen's
            perturbLoad =  0;
            perturbGen 	=  (skew*absSum)/(noOfGenCorrections);
        }

        for (uint i = 0; i < BUSNUM; ++i) {

            askedVals[i].id        = startVals[i].id;
            askedVals[i].voltage   = startVals[i].voltage;
            askedVals[i].type      = startVals[i].type;
            askedVals[i].startPower= startVals[i].activePower;

            if (int(askedVals[i].type) == 1){
                askedVals[i].activePower  = 0;
                slack_bus_ID = i;
            } else {
                if(startVals[i].activePower >= 0)
                    askedVals[i].activePower =  (startVals[i].activePower + perturbLoad)*Scale;
                else
                    askedVals[i].activePower =  (startVals[i].activePower + perturbGen)*Scale;
            }

            if( askedVals[i].activePower  < 0)
                sumGens  += fabs(askedVals[i].activePower );
            else
                sumLoads += fabs(askedVals[i].activePower );

            if(fabs(bus[BUSNUM*5+i] - bus[BUSNUM*3+i]) < fabs(minValue) && (bus[BUSNUM*9 + i] != slack_bus_ID) )
                minValue = (bus_startvals[BUSNUM*5 + i] - bus_startvals[BUSNUM*3 + i]);

        }

        if(fabs(minValue) < 0.001 || (sumGens > 6) || (sumLoads > 6)){
            //  cout << "badValue" << '\n';
            badValue = true;
        }
    }

    void changeTopologyOf(BusVals* TopoVals){
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle( myvector.begin(), myvector.end(), g); //shuffle topology
        int j=0;
        BusVals* temp = new BusVals[BUSNUM];
        for (std::vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it) {
            int i = int(*it);
            temp[j].id           = TopoVals[i].id;
            temp[j].voltage      = TopoVals[i].voltage;
            temp[j].type         = TopoVals[i].type;
            temp[j].activePower  = TopoVals[i].activePower;
            temp[j].startPower   = TopoVals[i].startPower;
            ++j;
        }
        for (uint i = 0; i < BUSNUM; ++i) {
            TopoVals[i].id           = temp[i].id;
            TopoVals[i].voltage      = temp[i].voltage;
            TopoVals[i].type         = temp[i].type;
            TopoVals[i].activePower  = temp[i].activePower;
            TopoVals[i].startPower   = temp[i].startPower;
        }
        delete[] temp;
    }
    void perturb_and_assign_start(const double Scale, const BusVals* busMatrix ,
                                  BusVals* busMatrixOut,
                                  double& absSum,uint& noOfLoadCorrections, uint& noOfGenCorrections){

        // note : fRand does not return 0.0

        double buffer[BUSNUM]; double correction = 0 ;
        noOfLoadCorrections =0 ; noOfGenCorrections = 0;
        absSum = 0;
        double _randFactor = fRand(0.5,1.0);

        double zero_adjust=0.10;
        for (uint i = 0; i < BUSNUM; ++i) {

            busMatrixOut[i].id          = busMatrix[i].id;
            busMatrixOut[i].voltage     = busMatrix[i].voltage;
            busMatrixOut[i].type        = busMatrix[i].type;

            if((busMatrix[i].activePower) == 0){
                buffer[i] = 0;
                slack_bus_ID = i; //find slack ID
            }
            else{
                buffer[i] =	( busMatrix[i].activePower == 0.0) ? (fRand(-zero_adjust,zero_adjust)*Scale)
                                                               :(busMatrix[i].activePower
                                                                 + fRand(-busMatrix[i].activePower/_randFactor,
                                                                         busMatrix[i].activePower/_randFactor)
                                                                 );
                if(buffer[i] < 0)
                    ++noOfGenCorrections;
                else
                    ++noOfLoadCorrections;
            }
            correction += buffer[i];


        }

        correction /= (noOfLoadCorrections + noOfGenCorrections);
        //correction *= 0.996;

        for (uint j =0 ;j< BUSNUM; ++j) {
            if(int(busMatrixOut[j].type)  == 1)
                busMatrixOut[j].activePower = 0;
            else  // sum of Loads and Gens made same;
                busMatrixOut[j].activePower = buffer[j] - correction ;

            absSum += fabs(busMatrixOut[j].activePower);
        }

    }



    bool makeCase(Case type,BusVals* startVals,BusVals* askedVals
                  ,uint ,uint TopoNo,uint BaseInputNo){
        checkCase=type;
        itercheck =0 ;
        int failCase=0;
        do{

            ++itercheck;
            Scale = fRand(0.8,1.2);
            if(type==normal){
                float skewVals1[] = {0,0};
                skewArray.assign (skewVals1,skewVals1+1);
                _normOrOvldFactor = 1 + 0.2;
            }else if(type==overloaded){
                float skewVals2[] =  {0.1f,0.2f,0.3f};
                skewArray.assign (skewVals2,skewVals2 + 2) ;
                _normOrOvldFactor = 1 + fRand(0.05,0.25);
            }else{
                float skewVals2[] =  {-0.1f,-0.2f,-0.3f};
                skewArray.assign (skewVals2,skewVals2 + 2) ;
                _normOrOvldFactor = 1 + fRand(0.05,0.25);
            }

            skew = skewArray.at(int(rand() % skewArray.size()) );
            perturbType = 1 + rand()%3;
            endScale = Scale*_normOrOvldFactor;

            perturb_and_assign_end( endScale,perturbType,
                                    skew, absSum, noOfLoadCorrections,
                                    noOfGenCorrections,startVals,askedVals);
            if(badValue == true){

                --itercheck;
                if(++failCase > 10)
                    return false;

                continue;
            }

            iter = convertToMatlabInputAndRun(bus,askedVals,bus_sol);
        } while(itercheck > 10 && verify(bus_sol,iter,false) == false );
        if(itercheck > 10)
            return false;


        printToFile(BaseInputNo,TopoNo,askedVals);
        return true;

    }



    void createFileName(int BaseInputNo,int TopoNo){

        fileName << directory << "/";
        fileName << setfill('0') << setw(7) << ++Global_Counter ;
        fileName << std::fixed;
        fileName << std::setprecision(2);
        fileName <<	"_v" << voltageSlack << "_l" << loadSlack
                 << "_sc" << endScale<< "_sk" << skew;

        if 		(perturbType== 3) 	fileName <<	"_B";
        else if (perturbType == 2)	fileName <<	"_L";
        else 						fileName <<	"_G";

        fileName    << "_BN" << BaseInputNo
                    <<  "_TN" << TopoNo;

        if      (checkCase == normal) 		fileName <<	"_N_"	;
        else if (checkCase == overloaded) 	fileName <<	"_O_"	;
        else if (checkCase == underloaded) 	fileName <<	"_U_"	;

        fileName 	<< todayDateTime;

        fileName << ".dat" ;

    }

    void printToFile(int BaseInputNo,int TopoNo,BusVals* askedVals){

        voltageSlack = fRand(0,0.04);
        loadSlack = fRand(0,0.25);
        ++Stats_valid;

        if(1){

            createFileName(BaseInputNo,TopoNo);

            ofstream myfile(fileName.str().c_str(),  ios::out);

            for(uint i=0;i<BUSNUM;++i){

                if(int(bus[BUSNUM*9 + i])==2 || int(bus[BUSNUM*9 + i])==1)
                    vslack = 0;
                else
                    vslack = (askedVals[i].voltage*voltageSlack);

                if(int(bus[BUSNUM*9 + i]) == 1)
                    _actpower=0;
                else _actpower= (askedVals[i].activePower);

                outputStream
                        << fixed << setw(SETWIDTH-6) << right << i+1
                        << setw(SETWIDTH) << right <<  askedVals[i].voltage
                        << setw(SETWIDTH) << right << (askedVals[i].voltage + vslack)
                        << setw(SETWIDTH) << right << (askedVals[i].voltage - vslack)
                        << setw(SETWIDTH) << right << _actpower
                        << setw(SETWIDTH) << right << _actpower + fabs(_actpower*loadSlack)
                        << setw(SETWIDTH) << right << _actpower - fabs(_actpower*loadSlack)
                        << setw(SETWIDTH-6) << right << int(askedVals[i].type )
                        << setw(SETWIDTH) << right << ( askedVals[i].startPower )
                        << "\n";
            }

            myfile << outputStream.rdbuf();
            myfile.close();
        }
        fileName.str("");
        outputStream.str("");
    }

    void resetTopology(){
        std::sort (myvector.begin(), myvector.end());
    }
};

//GLOBALS
uint    g_noOfBaseInputs = 0,
        g_noOfTopologies = 0,
        g_noOfCasesEach  = 0;
string  g_directory;

void generateInputs(uint seed, uint counter, uint eidx, uint nBase, BusVals* inputValsBase) {
    Network nw;

    nw.directory = g_directory;

    nw.noOfBaseInputs = nBase;
    nw.noOfTopologies = g_noOfTopologies;
    nw.noOfCasesEach  = g_noOfCasesEach;

    nw.Global_Counter = counter;

    nw.todayDateTime = currentDateTime();

    BusVals* inputVals = new BusVals[BUSNUM];
    memcpy(inputVals, inputValsBase, sizeof(BusVals)*BUSNUM);

    BusVals* startVals = new BusVals[BUSNUM];
    BusVals* askedVals = new BusVals[BUSNUM];

    for(uint i=0;i<BUSNUM;++i)
    {
        nw.bus[BUSNUM*2 + i]
                = nw.bus[BUSNUM*3 + i]
                = nw.bus[BUSNUM*4 + i]
                = nw.bus[BUSNUM*7 + i]
                = nw.bus[BUSNUM*8 + i]= 0 ;
        nw.bus_startvals[BUSNUM*2 + i]
                = nw.bus_startvals[BUSNUM*3 + i]
                = nw.bus_startvals[BUSNUM*4 + i]
                = nw.bus_startvals[BUSNUM*7 + i]
                = nw.bus_startvals[BUSNUM*8 + i]= 0 ;
    }

    for (uint i=0; i<BUSNUM; ++i) nw.myvector.push_back(i); //fill vector with values

    nw.findnoofGensAndSlack(BUSNUM,2,inputVals);		 //find noofgens
    nw.findnoofGensAndSlack(BUSNUM,1,inputVals);;

    srand(seed);
    rand();

    for(uint j = 0; j < nw.noOfBaseInputs; j++){
        //get one base point
        nw.resetTopology();
        nw.changeTopologyOf(inputVals);
        nw.itercheck=0;
        nw.Scale = (0.7);

        do{
            ++nw.itercheck ;
            nw.perturb_and_assign_start(nw.Scale,inputVals, startVals, nw.absSum,
                                        nw.noOfLoadCorrections,nw.noOfGenCorrections);

            nw.iter2 = nw.convertToMatlabInputAndRun(nw.bus_startvals,startVals,nw.bus_sol2);

        } while(nw.itercheck < 10 && nw.verify(nw.bus_sol2,nw.iter2,true) == false);

        if(nw.itercheck >= 10 ){   --j;  // cout << "start fail" << '\n';
            continue;  }

        for(uint i=0;i<nw.noOfTopologies;i++){
            nw.itercheck=0;
            do{
                ++nw.itercheck ;
                nw.changeTopologyOf(startVals);
                nw.iter2 = nw.convertToMatlabInputAndRun(nw.bus_startvals,startVals,nw.bus_sol2);
            } while(nw.itercheck < 10 && nw.verify(nw.bus_sol2,nw.iter2,true) == false);

            if(nw.itercheck >= 10 ){   --i;  // cout << "start fail" << '\n';
                continue;  }

            for(uint k=0;k < nw.noOfCasesEach ;k++){

                if(nw.makeCase(normal,     startVals,askedVals,k,i,j) &&
                   nw.makeCase(overloaded, startVals,askedVals,k,i,j) &&
                   nw.makeCase(underloaded,startVals,askedVals,k,i,j) )
                {
                    if (nw.Global_Counter == eidx)
                        goto M_EXIT;
                }else{
                    --k;
                    //   cout << "case failed" << '\n';
                    continue;
                }

            }
        }
    }

M_EXIT:
    delete [] inputVals;
    delete [] startVals;
    delete [] askedVals;
    trace( "Valid : " << nw.Stats_valid
         << "	Not Valid : " << (nw.noOfBaseInputs*nw.noOfTopologies*(3*nw.noOfCasesEach)) - nw.Stats_valid
         << "	Iter>5 : " <<  nw.iter_more
         << " 	Asked Invalid : " << nw.Asked_invalid
         << "	Start Invalid : " << nw.Start_invalid);
}

void singleThreaded(BusVals* inputVals) {
    auto beforeStart = GetTimeMs64();
    auto seed = (uint)time(NULL);
    generateInputs(seed, 0, g_noOfBaseInputs * g_noOfTopologies * g_noOfCasesEach * 3, g_noOfBaseInputs, inputVals);

    cout << " -----------------------\n";
    cout << " Elapsed: " << timeit(beforeStart) << endl;
}

void multiThreaded(BusVals* inputVals) {
    vector<thread> threads;

    uint nThreads = thread::hardware_concurrency();
    cout << " Num Threads: " << nThreads << endl;
    cout << " -----------------------\n";
    uint W = g_noOfBaseInputs;
    uint T = g_noOfBaseInputs * g_noOfTopologies * g_noOfCasesEach * 3;
    uint N = 1 + ((W - 1) / nThreads); //ceil

    uint start, end;

    auto seed = (uint)time(NULL);

    auto beforeStart = GetTimeMs64();
    for (uint i = 0; i < nThreads; ++i) {
          start = i*N;
          end = (uint)min(T, start + N);

          seed += i*10;

          uint sidx = start*g_noOfTopologies*g_noOfCasesEach*3;
          uint eidx = (uint) min(sidx + (g_noOfTopologies * g_noOfCasesEach * 3 * N), T);

          cout << " Thread " << i + 1 << ": startIdx:" << sidx << "\t # endIdx:" << eidx << endl;
          threads.push_back(thread(generateInputs, seed, sidx, eidx, N, inputVals));
          if (end == T) break;
    }
    for (auto& t : threads)
        t.join();
    cout << " -----------------------\n";
    cout << " Elapsed: " << timeit(beforeStart) << endl;
}

int main(int argc, char *argv[])
{
    g_noOfBaseInputs=10; //default

    if(argc > 1)
        g_noOfBaseInputs = atoi(argv[1]);

    g_noOfTopologies =10; //default

    if(argc > 2)
        g_noOfTopologies = atoi(argv[2]);

    g_noOfCasesEach =10; //default

    if(argc > 3)
        g_noOfCasesEach = atoi(argv[3]);

    if(argc > 4) {
        g_directory = argv[4];
    }
    if(!exists_test(g_directory)) {
        std::filesystem::create_directories(g_directory.c_str());
    }

    if(getLinevals(argv[6]))             return 1;

    BusVals* inputVals = new BusVals[BUSNUM];
    if(getBusvals(argv[5], inputVals) )  return 1;

    cout << " -----------------------\n";
    cout << "      Gen Inputs\n";
    cout << " -----------------------\n";
    //singleThreaded(inputVals);

    multiThreaded(inputVals);

    return 0;
}

