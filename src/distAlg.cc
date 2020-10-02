#include "distAlg.h"

Define_Module(DistAlg)

void DistAlg::initialize() {
  if (exec_once){ //run only once for entire network
    
    if(getLinevals(line_filename)){
      callFinish();
      endSimulation();
    }

    if(getBusvals(bus_filename)){
      callFinish();
      endSimulation();
    }

    for(int i =0;i < BUSNUM;++i){
      loadCost[i]=0;
    }

    statsStream << basename((char*)bus_filename) << ",";

    singleRunLoadFlow(); // run once
    exec_once = false;
  }
  m_ID = getIndex();
  m_heardFromNeighborsCount=0;
  m_localLoadCost=0;
  m_vCost=0;
  m_lCost=0;
  m_degree = gateSize("gate$o");
  powers[m_ID]= stable_powers[m_ID];
  
  ++g_initcalls;

	rpowers[m_ID]=powers[m_ID]*0.2;  // 20% of active power set as reactive power initially

	calculateLocalLoadcost(m_lCost,asked_power[m_ID],powers[m_ID]);

	if(g_initcalls==BUSNUM){
	  callMatlab();
	  ++g_rounds;
	}

	for(uint k=m_degree;k--;)
	  send(DistAlg::genMsg(powers[m_ID]), "gate$o", k);
	
	m_processMsgs = 1;
}

void DistAlg::handleMessage(cMessage* in_msg)
{

  valmsg* msg = (valmsg*)in_msg;

    powers[m_ID] = mat_powers_c[m_ID];
    calculateLocalLoadcost(m_lCost,asked_power[m_ID],powers[m_ID]);

    if(m_processMsgs==0){
      m_processMsgs=1;
      for(uint k=m_degree;k--;)
        send(genMsg(powers[m_ID]), "gate$o", k);
    }

    if( m_processMsgs ) { 	//hear from neighbours

      if(g_doDiffCommunication){
        if(g_rounds==1){ //coin toss for diff comm topology edge drops
          if(rand()%9 < 5){             //33% drop rate
            m_sumLocalVoltcosts += Vscale(msg->srcVolt);
            m_sumLocalLoadcosts += (msg->src==slack_busm_ID)
              ? 0 : msg->srcLCost;
          }else{
            nbrToBeDroppedList.push_back(msg->src);
          }
        }else{
          if( find (nbrToBeDroppedList.begin(), nbrToBeDroppedList.end(),
                msg->src)
              != nbrToBeDroppedList.end())
          {}else{
                  m_sumLocalVoltcosts += Vscale(msg->srcVolt);
                  m_sumLocalLoadcosts += (msg->src==slack_busm_ID)
                    ? 0 : msg->srcLCost;
                }
        }
      }else{
        m_sumLocalVoltcosts += Vscale(msg->srcVolt);
        m_sumLocalLoadcosts += (msg->src==slack_busm_ID)
          ? 0 : msg->srcLCost;
      }
      ++m_heardFromNeighborsCount;

      if(m_heardFromNeighborsCount == m_degree) {// heard from all nbrs
        m_heardFromNeighborsCount=0;

        calculateLoadcost(loadCost[m_ID],m_sumLocalLoadcosts,m_lCost);

        calculateVoltagecost(m_vCost,m_sumLocalVoltcosts);

        calculateNewLoad(powers[m_ID],rpowers[m_ID],loadCost[m_ID],m_vCost);

		    m_processMsgs = 0;
        resets(m_sumLocalLoadcosts,m_sumLocalVoltcosts);

        ++g_nodesDone;
      }

      if(g_nodesDone == BUSNUM) {  // all nodes compute done
        callMatlab();
        ++ g_rounds ;
		
		    //cout << g_rounds <<  endl;
        
		if(g_rounds == RNDS) {
		  delete msg;
		  writetoFileandFinish();
			
			callFinish();
			endSimulation();
        }
        m_processMsgs=1;
        calculateLocalLoadcost(m_lCost,asked_power[m_ID],
            powers[m_ID]);
		
		for(uint k=m_degree;k--;)
			send(genMsg(mat_powers_c[m_ID]), "gate$o", k);
		
        g_nodesDone = 0;
      }
    }
    delete msg;
}

valmsg* DistAlg::genMsg(double &Srcval){

  valmsg* msg = new valmsg;
  msg->src = m_ID;
  msg->srcValue = Srcval;
  msg->srcVolt = real_volts[m_ID];
  msg->srcLCost = m_lCost;
  return msg;
}

void DistAlg::calculateNewLoad(double &load, double &rpower,
    const double &loadCost, const double &voltCost){


  double _loadScale,_voltScale,_closeToZero,
  _maxLoadChangeFactor;
  _loadScale              = 50;
  _voltScale              = 30;
  _maxLoadChangeFactor    = 0.1;
  _closeToZero            = 0.0000001 ;

  if(loadCost < 0) {
    if(load>= -_closeToZero && load <=_closeToZero)
      load = _closeToZero*10;
    load += ( fabs(load) * min(_maxLoadChangeFactor,(fabs(loadCost)/_loadScale))) ;
  }

  if(loadCost >= 0) {
    if(load>= -_closeToZero && load <=_closeToZero)
      load = -_closeToZero*10;
    load -= (fabs(load) *  min(_maxLoadChangeFactor,(fabs(loadCost)/_loadScale)) );
  }

  if(voltCost > 0){
    if(rpower>= -_closeToZero && rpower <=_closeToZero)
      rpower = _closeToZero*10;
    rpower += fabs(rpower)*(fabs(voltCost)/_voltScale) ;
  }

  if(voltCost <= 0){
    if(rpower>= -_closeToZero && rpower <=_closeToZero)
      rpower = -_closeToZero*10;
    rpower -= fabs(rpower)*(fabs(voltCost)/_voltScale) ;
  }

}

void DistAlg::calculateLocalLoadcost(double &l_cost, const double &tval, const double &load){

  l_cost = Lscale(load,tval) ;
}

void DistAlg::calculateVoltagecost(double &v_cost, const double &sum_local_voltcosts){

  if(std::find(gens.begin(), gens.end(), m_ID ) != gens.end())
    v_cost = 0;
  else
    v_cost = degreeBasedAveraging(sum_local_voltcosts,
        Vscale(real_volts[m_ID]));
}

void DistAlg::calculateLoadcost(double &loadCost, const double &sum_local_loadcosts, const double &l_lcost){

  loadCost = degreeBasedAveraging(sum_local_loadcosts,l_lcost);
}

void DistAlg::endlOpenWriteCloseFile(string fileName, std::stringstream& Stream){
  Stream << '\n';
  ofstream myfile(fileName,  ios::out | ios::app);
  myfile << Stream.rdbuf();
  myfile.close();
}

void DistAlg::writetoFileandFinish(){

  //write to file
  endlOpenWriteCloseFile("Stats_collection.csv",statsStream);
  endlOpenWriteCloseFile("maxVoltError.csv",maxVoltStream);
  endlOpenWriteCloseFile("avgVoltError.csv",avgVoltStream);
  endlOpenWriteCloseFile("maxPerfLoadError.csv",maxPerfLoadStream);
  endlOpenWriteCloseFile("avgPerfLoadError.csv",avgPerfLoadStream);
  endlOpenWriteCloseFile("maxLoadVoltError.csv",maxLoadVoltStream);

  callFinish();
  endSimulation();
}

double DistAlg::degreeBasedAveraging(const double &sumOfOthers, const double &ownValue){
  double a = 27/2,
         b = 3/2,
         fm_degree = 1 + a/pow(b,m_degree);

  return (sumOfOthers + ( fm_degree*ownValue ))/(m_degree + fm_degree);
}

double DistAlg::Vscale(const double &val){ // local voltage cost actually

  double a,b,c,Delta,Delta1,C,del,del1,del2,del3,del4;
  C=2; //scaling from slope 2 to 3
  del=0.01;
  del1 = del2 = del3 = del4 =del;
  c=max_volt[m_ID];
  a=min_volt[m_ID];
  b=bus[BUSNUM + m_ID]; //ideal voltages

  Delta = dynamicDelta(0.5,1.2);
  Delta1 = dynamicDelta(0.25,0.6);

  return GenericScale(a,b,c,Delta,Delta1,C,del1,del2,del3,del4,val);
}

double DistAlg::Lscale(const double&val,const double&tval){ // local LOAD cost
  double a,b,c,Delta,Delta1,C,alpha,del1,del2,del3,del4;
  C=2;
  alpha = 0.10;
  c=max_power[m_ID];
  a=min_power[m_ID];

  if(a==0)
    del1 = del2 = fabs(c-a)/10;
  else
    del1 = del2 = fabs(alpha*a);
  if(c==0)
    del3 = del4 = fabs(c-a)/10;
  else
    del3 = del4 = fabs(alpha*c);

  b=tval; 	   //ideal voltages

  Delta = dynamicDelta(0.05,0.2);
  Delta1 = dynamicDelta(0.2,0.6);

  return GenericScale(a,b,c,Delta,Delta1,C,del1,del2,del3,del4,val);
}

double DistAlg::GenericScale(const double& a,const double&b,const double&c,
    const double&Delta,const double&Delta1,
    const double&C,
    double &del1,double &del2,double &del3,double &del4,
    const double&val
    ) {

  double out,a1,c1,s1,ya,yc,s2,gammab,gammac,gammad,gammae,s3,s4,s5;

  a1=a - del2;
  c1=c + del3;

  if(c==b && b==a)
    s1 =0;
  else
    s1=fabs(Delta/max(fabs(c-b),fabs(b-a))); //line b
  ya=s1*(b-a);
  yc=s1*(c-b);

  if(del1 == 0)
    del1 = 0.1;
  s2 = (Delta1 / del1); //line a
  gammab = -ya;

  s3 = (C*s2) ; //line c
  gammac = -ya - Delta1;

  if(del3 == 0)
    del3 = 0.1;
  s4 = (Delta1 / del3);
  gammad = yc ;

  s5 = (C*s4);
  gammae = yc + Delta1 ;

  if(val < a1)
    out=s3*((val)-a1) + (gammac);
  else if(val >= a1 && val < a)
    out=s2*((val)-a) + (gammab);
  else if(val >= a && val < c)
    out =s1*((val)-b);
  else if(val >= c && val < c1)
    out=s4*((val)-c) + (gammad);
  else // > c1
    out=s5*((val)-c1) + (gammae);

  return out;
}

void DistAlg::resets(double&sum_local_loadcosts, double&sumvoltcosts){
  //m_processMsgs=0;
  sum_local_loadcosts = 0;
  sumvoltcosts = 0 ;
}

double DistAlg::dynamicDelta(double first, double last){
  return (last-first)*(g_matlabCallCount+1)/RNDS + first;
}

void DistAlg::callMatlab() {
  loadFlowSolution();
}
