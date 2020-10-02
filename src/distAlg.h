#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <map>
#include <libgen.h>
#include <memory>
#include <omnetpp.h>
#include "constants.h"

typedef unsigned int uint;

extern "C" {
  #include "loadflow/matlabrun.h"
}

using namespace omnetpp;
using namespace std;

//-----------
//HELPERS
//-----------
template<typename T>
bool is_infinite( const T &value )
{
    T max_value = std::numeric_limits<T>::max();
    T min_value = - max_value;

    return ! ( min_value <= value && value <= max_value );
}

template<typename T>
bool is_nan( const T &value )
{
    // True if NAN
    return value != value;
}

template<typename T>
bool is_valid( const T &value )
{
    return (!is_infinite(value) && !is_nan(value) );
}

struct valmsg : ::omnetpp::cMessage
{
    int src;
    double srcValue = 0;
    double srcVolt = 0;
    double srcLCost = 0;
};

//-----------
//GLOBALS
//-----------
const int RNDS                     = 600;
const std::vector<int> printRounds = {1,150,300,450,600};

const char *bus_filename, *line_filename;

bool  g_doDiffCommunication,
      g_zeroloadslack,
      g_zerovoltslack;


double single_run_slack;
bool exec_once = true;
double bus[BUSNUM*10];
double 	max_volt[BUSNUM],min_volt[BUSNUM],
        max_volt_meas[BUSNUM],min_volt_meas[BUSNUM],
        max_power_meas[BUSNUM],max_power[BUSNUM],
        min_power_meas[BUSNUM],min_power[BUSNUM],
        asked_power[BUSNUM],stable_powers[BUSNUM],
        linefull[LINENUM*7];

double M,L,Slack_error,Slack_error_percentage;
double Va,Vm,V1,V2;
std::vector<int> gens;
std::vector<int>::iterator it;
uint no_Of_Gens=0;  //includes slack bus
uint no_Of_slack=0;
int slack_busm_ID;
uint iter=0;

double maxPerfLoadError_Buffer[RNDS];
double loadAt600Iter =0 ,
       loadAtPresent=0;

double bus_sol[BUSNUM*10];

bool _getBusValsDone=false;
bool _getLineValsDone=false;

uint g_matlabCallCount = 0,
     g_nodesDone = 0,
     g_rounds = 0,
     g_initcalls	= 0;

string st[]={"Voltage |","Voltage Cost |","localloadCost|","loadCost |",
  "React_Powers |","Powers |"," %PowerError"," | %Perf_PowerError"};

double powers[BUSNUM],	rpowers[BUSNUM],
       mat_powers_c[BUSNUM],
       real_volts[BUSNUM],loadCost[BUSNUM];

std::stringstream avgVoltStream;
std::stringstream maxVoltStream;
std::stringstream maxPerfLoadStream;
std::stringstream avgPerfLoadStream;
std::stringstream maxLoadVoltStream;
std::stringstream statsStream;

// ----------------------------------------------
class DistAlg : public cSimpleModule
{
  protected:
	  virtual void initialize() override;
    virtual void handleMessage(cMessage*) override;
	
  public:
	enum BusType {slack=1,gen,load};

  valmsg* msg;
  std::vector<int> nbrToBeDroppedList;

  uint m_ID,
        m_heardFromNeighborsCount,
        m_degree;
  double m_localLoadCost,
          m_vCost,
          m_lCost;

	bool m_processMsgs = 0;
	double m_sumLocalLoadcosts = 0,
         m_sumLocalVoltcosts = 0;

  DistAlg() : cSimpleModule() {}
  ~DistAlg();
	
	virtual valmsg* genMsg(double &Srcval);
	
  void calculateNewLoad(double &load,double &rpower,const double &loadCost,const double &voltCost);
  void calculateLocalLoadcost(double &l_cost,const double &tval,const double &load);
  void calculateVoltagecost(double &v_cost,const double &sumvolts);
  void calculateLoadcost(double &loadCost,const double &sum_local_loadcosts,const double &l_lcost);
  void writetoFileandFinish();
  void resets(double&sum_local_loadcosts,double&sumvolts);
  double degreeBasedAveraging(const double &sumOfOthers, const double &ownValue);
  double Vscale(const double&val);
  double Lscale(const double&val, const double&tval);
  double GenericScale(const double& a,const double&b,const double&c,const double&Delta,
      const double&Delta1,const double&C,double &del1,
      double &del2,double &del3,double &del4,const double&val
      );
  void endlOpenWriteCloseFile(string fileName, std::stringstream &Stream);

  double dynamicDelta(double first, double last);
  void callMatlab();
  int findnoofGensAndSlack(int len, BusType type);
  int getBusvals(const char* input);
  int getLinevals(const std::string &input);
  void initializeBusParameters(double *bus, double *powers, double *rpowers, uint slack_busm_ID, double &genSum, double &loadSum);
  void maxAndAvg(double &max, double value, double &avg);
  void loadFlowSolution();
  void singleRunLoadFlow();
};


DistAlg::~DistAlg()
{}

int DistAlg::findnoofGensAndSlack(int len, BusType type){
  for (int i = 0; i < len; ++i){

    if (bus[BUSNUM*9 + i] == type) {
      if(type == 1){
        ++no_Of_slack;
        gens.push_back(i);
        return i;
      }else{
        ++no_Of_Gens;
        gens.push_back(i);
      }
    }
  }
  return 0;
}

int DistAlg::getBusvals(const char* input){
  if(!_getBusValsDone){
    std::ifstream ifs(input);
    if(!ifs.good()){
      cerr << "			 ################     error opening bus_values_file      ###############" << '\n';
      return 1;
    }
    for(uint i=0;i<BUSNUM;++i){
      bus[BUSNUM*2 + i]=	bus[BUSNUM*3 + i]= bus[BUSNUM*4 + i]=	bus[BUSNUM*7 + i]=	bus[BUSNUM*8 + i]= 0 ;
    }
    std::string line;
    int num=0;
    // double dummy;
    while(std::getline(ifs, line)){ // read one line from ifs
      if (line.length() == 0) continue;

      std::istringstream iavgVoltStream(line); // access line as a stream

      iavgVoltStream >> bus[num] >> bus[BUSNUM + num] >> max_volt[num]
        >> min_volt[num] >> asked_power[num]
        >> max_power[num] >> min_power[num] >> bus[BUSNUM*9 + num] >> stable_powers[num];
		
      max_power_meas[num]=max_power[num];
      min_power_meas[num]=min_power[num];
      max_volt_meas[num]=max_volt[num];
      min_volt_meas[num]=min_volt[num];

      if(!g_zeroloadslack){
        max_power[num] = asked_power[num];
        min_power[num] = asked_power[num];
      }

      if(!g_zerovoltslack){
        max_volt[num] = bus[BUSNUM + num];
        min_volt[num] = bus[BUSNUM + num];
      }

      ++num;
    }
    slack_busm_ID = findnoofGensAndSlack(BUSNUM,slack); //find slack ID
    findnoofGensAndSlack(BUSNUM,gen);				   //find gen ID's
    _getBusValsDone = true;
  }
  return 0;
}

int DistAlg::getLinevals(const std::string& input){
  if(!_getLineValsDone){
    std::ifstream ifs(input);
    if(!ifs.good()){
      cerr << "			 ################     error opening line_values_file      ###############" << '\n';
      return 1;
    }
    std::string l1;
    int num=0;
    while(std::getline(ifs, l1)){ // read one line from ifs
      if (l1.length() == 0) continue; 
      std::istringstream iavgVoltStream(l1); // acceavgVoltStream line as a stream
      iavgVoltStream >> linefull[num] >> linefull[LINENUM + num] >> linefull[LINENUM*2 + num]
        >> linefull[LINENUM*3 + num] >> linefull[LINENUM*4 + num]
        >> linefull[LINENUM*5 + num] >> linefull[LINENUM*6 + num] ;
      ++num;
    }
    _getLineValsDone = true;
  }
  return 0;
}

void DistAlg::initializeBusParameters(double *bus,double *powers,double *rpowers,
    uint slack_busm_ID,double &genSum,double &loadSum){
  genSum = 0;
  loadSum = 0;
  for(uint i=0;i<BUSNUM;i++){
    if (bus[BUSNUM*9 + i] == 3.0){
      bus[BUSNUM*5 + i]	=powers[i];
      bus[BUSNUM*6 + i]	=rpowers[i];
      bus[BUSNUM*3 + i] = bus[BUSNUM*4+i] =0;
    }else{
      bus[BUSNUM*3 + i]	=-powers[i];
      bus[BUSNUM*4 + i]	=-rpowers[i];
      bus[BUSNUM*5 + i] = bus[BUSNUM*6+i] =0;
    }

    if(bus[BUSNUM*5+i]-bus[BUSNUM*3+i] < 0)
      loadSum+=bus[BUSNUM*5+i]-bus[BUSNUM*3+i];
    else{
      if(i != slack_busm_ID)
        genSum +=bus[BUSNUM*5+i]-bus[BUSNUM*3+i];
    }
  }
}

void DistAlg::maxAndAvg(double &max,double value,double& avg){
  if(fabs(max) < fabs(value))
    max = value;
  avg += fabs(value);
}

void DistAlg::loadFlowSolution()
{
  initializeBusParameters(bus,powers,rpowers,slack_busm_ID,M,L);

  matlabrun(bus,linefull,no_Of_Gens,slack_busm_ID+1,
      no_Of_slack,bus_sol);
  ++g_matlabCallCount;

  double	maxVoltError=0,
          avgVoltError=0,
          myVoltError =0 ,
          loadError=0,
          perfLoadError=0,
          maxLoadError=0,
          avgLoadError=0,
          maxPerfLoadError=0,
          avgPerfLoadError=0,
          maxLoadVoltError=0,
          loadVoltError=0,
          avgLoadVoltError=0,
          maxVolt=0,
          avgVolt=0,
          voltCorr=0,
          //reactMax=0,
          reactAvg=0,
          rpower=0;
          //pass;


  for(int i=0;i<BUSNUM;i++){
    mat_powers_c[i] = bus_sol[BUSNUM*5 + i] - bus_sol[BUSNUM*3 + i];
    // rpower = bus_sol[BUSNUM*6 + i] - bus_sol[BUSNUM*4 + i];
    real_volts[i] = bus_sol[BUSNUM + i];
    myVoltError = (real_volts[i]>max_volt_meas[i])?(real_volts[i] - max_volt_meas[i]) : (
        (real_volts[i]<min_volt_meas[i]) ? (real_volts[i] - min_volt_meas[i]) : 0 );
    if(i==slack_busm_ID){
      perfLoadError=0;
      loadError=0;
    }else{
      if(asked_power[i] == 0)
        loadError= (mat_powers_c[i]-asked_power[i]) ;
      else
        loadError= (mat_powers_c[i]-asked_power[i]) / (asked_power[i])*100 ;
      if((max_power_meas[i] * min_power_meas[i]) ==0){
        perfLoadError=(mat_powers_c[i]>max_power_meas[i])
          ?((mat_powers_c[i] - max_power_meas[i] )): (
              (mat_powers_c[i]<min_power_meas[i])
              ? (( mat_powers_c[i] - min_power_meas[i] )) : 0 );
      }else{
        perfLoadError=(mat_powers_c[i]>max_power_meas[i])
          ?((mat_powers_c[i] - max_power_meas[i] )/max_power_meas[i]*100 ): (
              (mat_powers_c[i]<min_power_meas[i])
              ? (( mat_powers_c[i] - min_power_meas[i] )/min_power_meas[i]*100) : 0 );
      }


    }

    loadVoltError = myVoltError*perfLoadError;
    voltCorr = real_volts[i] - 1;
    maxAndAvg(maxVolt,voltCorr,avgVolt);
    maxAndAvg(maxLoadVoltError,loadVoltError,avgLoadVoltError);
    maxAndAvg(maxVoltError,myVoltError,avgVoltError);
    maxAndAvg(maxLoadError,loadError,avgLoadError);
    maxAndAvg(maxPerfLoadError,perfLoadError,avgPerfLoadError);


  }

  reactAvg        /= BUSNUM;
  avgVolt 		/= BUSNUM;
  avgVoltError 	/= BUSNUM;
  avgLoadError 	/= BUSNUM;
  avgPerfLoadError/= BUSNUM;

  Va = avgVoltError;
  Vm = maxVoltError;

  Slack_error= fabs(mat_powers_c[slack_busm_ID]) - fabs(M-L);
  Slack_error_percentage = (Slack_error)*100/max(fabs(M),fabs(L));

  maxPerfLoadError_Buffer[g_matlabCallCount - 1] = maxPerfLoadError;

  maxVoltStream << maxVoltError << ",";  //save to maxvoltStream
  avgVoltStream << avgVoltError << ",";  //save to avgvoltStream

  maxPerfLoadStream << maxPerfLoadError << ",";
  avgPerfLoadStream << avgPerfLoadError << ",";

  maxLoadVoltStream << maxLoadVoltError << ",";    //load*volt error

  if(find (printRounds.begin(), printRounds.end(), g_matlabCallCount) != printRounds.end()){
    statsStream.precision(5);
    statsStream << std::fixed
      //<< std::cout << std::setprecision(5)
      << Slack_error_percentage  				<< ","
      << (mat_powers_c[slack_busm_ID])	<< ","
      << maxVoltError 						      << ","
      << avgVoltError 						      << ","
      << (Vm/V1) 								        << ","
      << (Va/V2)								        << ","
      << maxVolt								        << ","
      << avgVolt 								        << ","
      << maxLoadVoltError						    << ","
      << maxLoadError 						      << ","
      << avgLoadError 						      << ","
      << maxPerfLoadError 					    << ","
      << avgPerfLoadError 					    << ","
      ;
  }

  if(g_matlabCallCount == RNDS)
  {

    double reactAvg = 0;
    for(int i=0;i<BUSNUM;i++)
    {
      mat_powers_c[i] = bus_sol[BUSNUM*5 + i] - bus_sol[BUSNUM*3 + i];
      rpower = bus_sol[BUSNUM*6 + i] - bus_sol[BUSNUM*4 + i];
      if(bus_sol[9+i] != 1)
        reactAvg += fabs(rpower);
    }
    reactAvg /= BUSNUM;

    for(int j=(RNDS-1);j>=0;j--){
      if((fabs(maxPerfLoadError_Buffer[j]))     //5% error in value
          > (0.95 * (fabs(maxPerfLoadError)) ) ){
        statsStream << j << ",";
        statsStream << reactAvg << ",";
        break;
      }
    }
  }
}

void DistAlg::singleRunLoadFlow()
{
  uint _percentage;
  double overallMinimumMaxVolt=20;
  double maxVolt,optimumRpowerPercentage=20;
  double avgVolt=0;
  uint _count=0;
  _percentage = 0;
  double voltCorr;

  double avgLoad=0;
  double load=0;
  double maxVoltError=0, avgVoltError=0, voltError=0;

  do{
    maxVolt = 0 ;
    avgVolt = 0;
    for(uint i=0;i<BUSNUM;i++){
      rpowers[i] = asked_power[i]*_percentage/100;
    }
    initializeBusParameters(bus,asked_power,rpowers,
        slack_busm_ID,M,L);

    matlabrun(bus,linefull,no_Of_Gens,slack_busm_ID+1,
        no_Of_slack,bus_sol);

    avgLoad = 0;
    for(uint i=0;i<BUSNUM;++i){
      load = bus_sol[BUSNUM*5 + i] - bus_sol[BUSNUM*3 + i];
      avgLoad += load;
      voltCorr = bus_sol[BUSNUM + i] - 1;
      maxAndAvg(maxVolt,voltCorr,avgVolt);
    }

    avgLoad/=BUSNUM;

    if(is_valid(maxVolt) && is_valid(avgLoad) && (fabs(overallMinimumMaxVolt) > fabs(maxVolt) ) ){
      overallMinimumMaxVolt = maxVolt;
      optimumRpowerPercentage = _percentage;
    }

    _percentage+=10;

  }while(++_count < 6);

  _percentage = optimumRpowerPercentage ;

  for(uint i=0;i<BUSNUM;i++){
    rpowers[i] = asked_power[i]*_percentage/100;
  }

  initializeBusParameters(bus,asked_power,rpowers,
      slack_busm_ID,M,L);

  iter = matlabrun(bus,linefull,no_Of_Gens,slack_busm_ID+1,
      no_Of_slack,bus_sol);

  single_run_slack = 	bus_sol[BUSNUM*5+slack_busm_ID] - bus_sol[BUSNUM*3+slack_busm_ID];
  Slack_error= fabs(single_run_slack) - fabs(M-L);
  Slack_error_percentage = (Slack_error)*100/max(fabs(M),fabs(L));

  maxVoltError=0, avgVoltError=0, voltError=0;
  avgVolt=0;
  maxVolt=0;

  for(uint i=0;i<BUSNUM;++i){
    voltError = bus_sol[BUSNUM + i] - bus[BUSNUM + i];
    maxAndAvg(maxVoltError,voltError,avgVoltError);
    voltCorr = bus_sol[BUSNUM + i] - 1;
    maxAndAvg(maxVolt,voltCorr,avgVolt);
  }

  avgVolt/=BUSNUM;
  avgVoltError /= BUSNUM;
  V2 = avgVoltError ;
  V1 = maxVoltError ;

  statsStream
    << Slack_error_percentage 	<< ","
    << single_run_slack     	<< ","
    << optimumRpowerPercentage	<< ","
    << V1 						<< ","
    << V2  						<< ","
    << maxVolt 					<< ","
    << avgVolt 					<< ","
    << M 						<< ","
    << L 						<< ","
    ;

  maxVoltStream << V1 << ",";
  avgVoltStream << V2 << ",";
}






