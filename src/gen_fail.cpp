#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <thread>

#include <sys/stat.h>
#include <filesystem>
#include "constants.h"
#include "util_timer.h"
#include "bfs.h"
using namespace std;
typedef unsigned int uint;

inline bool exists_test(const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

struct Bus
{
    int	nodeID;
    double
    voltage,
    voltageMax,
    voltageMin,
    activePower,
    activePowerMax,
    activePowerMin,
    busType,
    startActivePower ;
};

struct Line
{
    int node1,
    node2;
    double
    dummy1,
    dummy2,
    dummy3,
    dummy4,
    dummy5;
};

class FileOperations
{

    stringstream outFileName;
public:
    //  FileOperations(int BUSNUM) { this->BUSNUM = BUSNUM;}
    void appendToFileName(string, int, int, string, int iterNum);
    void printToFile(int, int, string, Bus*, string, int iterNum);
    void removesubstr(string&,const string&);
    int readInputFiles(Line*, Bus*, string, string lineFileName);
};

void FileOperations::removesubstr(string& s,const string& remove){
    s.erase(s.find(remove),remove.length());
}

void FileOperations::appendToFileName(string fileName,int radius,int noOfGensRemoved,string directory,int iterNum){

    removesubstr(fileName,".dat");
    outFileName << directory << "/" << fileName << "_R" << radius << "_i" << iterNum
                << "_nG" << noOfGensRemoved << ".dat";
}

void FileOperations::printToFile(int radius,int noOfGensRemoved,string fileName,
                                 Bus* bus,string directory,int iterNum){

    outFileName.str("");
    const int SETWIDTH = 14;
    appendToFileName(fileName,radius,noOfGensRemoved,directory,iterNum);
    ofstream myfile(outFileName.str().c_str(), ios::out);
    //std::cout << outFileName.str().c_str() << std::'\n';
    stringstream outputStream;
    for(int i=0;i<BUSNUM;++i){

        outputStream
                << fixed << setw(SETWIDTH-6) << right << bus[i].nodeID
                << setw(SETWIDTH) << right <<  bus[i].voltage
                << setw(SETWIDTH) << right <<  bus[i].voltageMax
                << setw(SETWIDTH) << right <<  bus[i].voltageMin
                << setw(SETWIDTH) << right <<  bus[i].activePower
                << setw(SETWIDTH) << right <<  bus[i].activePowerMax
                << setw(SETWIDTH) << right <<  bus[i].activePowerMin
                << setw(SETWIDTH-6) << right << int(bus[i].busType)
                << setw(SETWIDTH) << right << bus[i].startActivePower
                << "\n";
    }

    myfile << outputStream.rdbuf();
    myfile.close();

}

Line line[LINENUM];
int getLineData(const string& lineFileName) {

    ifstream lineFile(lineFileName);
    if(!lineFile.good()){
        cout << "			 ################     error opening line_values_file      ###############" << '\n';
        return 1;
    }

    std::string lineIn;
    int num=0;

    while(std::getline(lineFile, lineIn)){ // read one line from ifs
        std::istringstream iStream(lineIn); // access line as a stream
        iStream >> line[num].node1 >> line[num].node2 >> line[num].dummy1 >> line[num].dummy2 >> line[num].dummy3 >>
                line[num].dummy4 >> line[num].dummy5;
        ++num;
    }
    return 0;
}

int getBusData(Bus* bus, string busFileName){
    ifstream busFile(busFileName);
    if(!busFile.good()){
        cout << "			 ################     error opening bus_values_file      ###############" << '\n';
        return 1;
    }

    std::string lineIn;
    int num=0;
    while(std::getline(busFile, lineIn)){ // read one line from ifs
        std::istringstream iStream(lineIn); // access line as a stream
        iStream >> bus[num].nodeID >> bus[num].voltage >> bus[num].voltageMax >> bus[num].voltageMin
                >> bus[num].activePower >> bus[num].activePowerMax >> bus[num].activePowerMin
                >> bus[num].busType >> bus[num].startActivePower ;
        ++num;
    }

    return 0;
}
string inputFolderName;
string outputDirectory;
string lineFileName;
uint numInstances;

void genFailCases(uint seed, string& fileName) {
    Bus bus[BUSNUM] ;
    Bus outBus[BUSNUM];
    int noOfGenRemoved;int randomNode;

    FileOperations f;
    Graph g(BUSNUM);

    if (getBusData(bus, fileName)) return ;

    fileName = filesystem::path(&fileName[0u]).filename().string();
    for (int i = 0; i < LINENUM; i++) {
        g.addEdge(line[i].node1-1, line[i].node2-1);
    }

    vector<int> neighbors;
    int radiusiter = 1;
    int radius;

    srand(seed);
    rand();
    uint iterNum = 0;

    while(iterNum < numInstances){

        radiusiter = !radiusiter;
        radius = radiusiter + 1; //toggle between 1 and 2
        memcpy(&outBus, &bus, sizeof(bus)); //reinitilize to original

        neighbors.clear();

        randomNode = rand() % BUSNUM;
        g.BFS(randomNode,radius,neighbors);

        noOfGenRemoved = 0;

        for(auto i: neighbors)
        {
            if(int(bus[i].busType) == 2){
                ++noOfGenRemoved;
                outBus[i].activePower   = 0;
                outBus[i].activePowerMax= 0;
                outBus[i].activePowerMin= 0;
            }
        }

        f.printToFile(radius,noOfGenRemoved,fileName,outBus,outputDirectory,iterNum);
        ++iterNum;
    }
}
void processFiles(vector<string>& files, uint seed, uint start, uint end) {
    for (uint i = start; i < end; ++i) {
        genFailCases(seed, files[i]);
    }
}
void singleThreaded() {
    auto now = GetTimeMs64();
    auto seed = (uint)time(NULL);
    for(auto& p: std::filesystem::directory_iterator(inputFolderName)) {
        string busFileName = p.path().string();
        genFailCases(seed, busFileName);
    }
    cout << "Elapsed: " << timeit(now) << endl;
}
void multiThreaded() {
    vector<thread> threads;

    vector<string> files;
    for(auto& p: std::filesystem::directory_iterator(inputFolderName)) {
        files.push_back(p.path().string());
    }

    uint nThreads = thread::hardware_concurrency();
    cout << " Num Threads: " << nThreads << endl;
    cout << " -----------------------\n";
    uint W = (uint)files.size();

    uint N = 1 + ((W - 1) / nThreads); //ceil
    auto seed = (uint)time(NULL);
    uint start, end;
    auto now2 = GetTimeMs64();

    for (uint i = 0; i < nThreads; ++i) {
          start = i*N;
          end = (uint)min(W, start + N);
          seed += i*10;
          cout << " Thread " << i + 1 << ": " << start << " to " << end << " files" << endl;
          threads.push_back(thread(processFiles, ref(files), seed, start, end));
          if (end == W)
              break;
    }
    for (auto& t : threads)
        t.join();
    cout << " -----------------------\n";
    cout << " Elapsed: " << timeit(now2) << endl;
}
int main(int , char *argv[])
{
    inputFolderName  = argv[1];
    outputDirectory  = argv[2];
    lineFileName     = argv[3];
    numInstances     = 2*atoi(argv[4]); // radius is toggled between 1&2 thus times 2

    if (!exists_test(outputDirectory))
        filesystem::create_directories(outputDirectory.c_str());

    if (getLineData(lineFileName)) return 1;

    cout << " -----------------------\n";
    cout << "      Gen Fails\n";
    cout << " -----------------------\n";
    //    singleThreaded();

    multiThreaded();

    cout << " Done!\n";

    return 0;
}
