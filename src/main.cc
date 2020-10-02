#include <omnetpp.h>

using namespace omnetpp;

/**
 * Emulates an empty omnetpp.ini, causing config options' default values
 * to be used.
 */
class EmptyConfig : public cConfiguration
{
  protected:
    class NullKeyValue : public KeyValue {
      public:
        virtual const char *getKey() const override {return nullptr;}
        virtual const char *getValue() const override {return nullptr;}
        virtual const char *getBaseDirectory() const override {return nullptr;}
    };
    NullKeyValue nullKeyValue;

  protected:
    virtual const char *substituteVariables(const char *value) const override {return value;}

  public:
    virtual const char *getConfigValue(const char *key) const override {return nullptr;}
    virtual const KeyValue& getConfigEntry(const char *key) const override {return nullKeyValue;}
    virtual const char *getPerObjectConfigValue(const char *objectFullPath, const char *keySuffix) const override {return nullptr;}
    virtual const KeyValue& getPerObjectConfigEntry(const char *objectFullPath, const char *keySuffix) const override {return nullKeyValue;}
};

/**
 * Defines a minimal environment for the simulation. Module parameters get
 * initialized to their default values (causing an error if there's no
 * default); module log messages (EV<<) get printed on the standard output;
 * calls to record results will be ignored.
 */
class MinimalEnv : public cNullEnvir
{
  public:
    // constructor
    MinimalEnv(int ac, char **av, cConfiguration *c) : cNullEnvir(ac, av, c) {}

    // model parameters
    virtual void readParameter(cPar *par) override {
        if (par->containsValue())
            par->acceptDefault();
        else
            throw cRuntimeError("no value for parameter %s", par->getFullPath().c_str());
    }

    // see module output
    virtual void sputn(const char *s, int n) {
        (void) ::fwrite(s,1,n,stdout);
    }

};

void simulate(const char *networkName, simtime_t limit)
{
    // set up the network
    cModuleType *networkType = cModuleType::find(networkName);
    if (networkType == nullptr) {
        printf("No such network: %s\n", networkName);
        return;
    }
    getSimulation()->setupNetwork(networkType); //XXX may throw exception
    getSimulation()->setSimulationTimeLimit(limit);

    // prepare for running it
    getSimulation()->callInitialize();

    // run the simulation
    bool ok = true;
    try {
        while (true) {
            cEvent *event = getSimulation()->takeNextEvent();
            if (!event)
                break;  //XXX
            getSimulation()->executeEvent(event);
        }
    }
    catch (cTerminationException& e) {
        //printf("Finished: %s\n", e.what());
    }
    catch (std::exception& e) {
        ok = false;
        printf("ERROR: %s\n", e.what());
    }

    if (ok)
        getSimulation()->callFinish();  //XXX may throw exception

    // finish the simulation and clean up the network
    getSimulation()->deleteNetwork();
}


extern const char *bus_filename, *line_filename;
extern bool g_doDiffCommunication, g_zeroloadslack, g_zerovoltslack;
		 
int main(int argc, char *argv[])
{
    // the following line MUST be at the top of main()
    cStaticFlag dummy;
	
	bus_filename = argv[1];
	line_filename = argv[2];
	g_doDiffCommunication      = atoi(argv[3]);
    g_zeroloadslack            = atoi(argv[4]);
    g_zerovoltslack            = atoi(argv[5]);

    // initializations
    CodeFragments::executeAll(CodeFragments::STARTUP);
    SimTime::setScaleExp(-12);

    // set up an environment for the simulation
    cEnvir *env = new MinimalEnv(argc, argv, new EmptyConfig());
    cSimulation *sim = new cSimulation("simulation", env);
    cSimulation::setActiveSimulation(sim);

    // load NED files
    cSimulation::loadNedSourceFolder(".");
    cSimulation::doneLoadingNedFiles();

    // set up and run a simulation model
    simulate("Custom1", 1000);

    // exit
    cSimulation::setActiveSimulation(nullptr);
    delete sim;

    // deallocate registration lists, loaded NED files, etc.
    CodeFragments::executeAll(CodeFragments::SHUTDOWN);
    return 0;
}

