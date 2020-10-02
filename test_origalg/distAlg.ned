
// this file is auto-generated using bin/updateconstants.sh
// do not modify directly. Please update bin/updateconstants.sh.

simple DistAlg
{
parameters:

  string bus_filename=default("");
  string line_filename=default("");
  int dodiffCommunication=default(0);
  int zeroloadslack=default(0);
  int zerovoltslack=default(0);

  @display("i=block/routing");

gates:
  inout gate[];
}

network Custom1
{
types:

  channel delay1 extends ned.DelayChannel  {
    delay = 100	ms;// uniform(0ms,100ms) ;
  }

submodules:

  node[39]:DistAlg;

connections:

node[0].gate++ <--> delay1 <--> node[1].gate++ ;
node[0].gate++ <--> delay1 <--> node[38].gate++ ;
node[1].gate++ <--> delay1 <--> node[2].gate++ ;
node[1].gate++ <--> delay1 <--> node[24].gate++ ;
node[2].gate++ <--> delay1 <--> node[3].gate++ ;
node[2].gate++ <--> delay1 <--> node[17].gate++ ;
node[3].gate++ <--> delay1 <--> node[4].gate++ ;
node[3].gate++ <--> delay1 <--> node[13].gate++ ;
node[4].gate++ <--> delay1 <--> node[5].gate++ ;
node[4].gate++ <--> delay1 <--> node[7].gate++ ;
node[5].gate++ <--> delay1 <--> node[6].gate++ ;
node[5].gate++ <--> delay1 <--> node[10].gate++ ;
node[6].gate++ <--> delay1 <--> node[7].gate++ ;
node[7].gate++ <--> delay1 <--> node[8].gate++ ;
node[8].gate++ <--> delay1 <--> node[38].gate++ ;
node[9].gate++ <--> delay1 <--> node[10].gate++ ;
node[9].gate++ <--> delay1 <--> node[12].gate++ ;
node[12].gate++ <--> delay1 <--> node[13].gate++ ;
node[13].gate++ <--> delay1 <--> node[14].gate++ ;
node[14].gate++ <--> delay1 <--> node[15].gate++ ;
node[15].gate++ <--> delay1 <--> node[16].gate++ ;
node[15].gate++ <--> delay1 <--> node[18].gate++ ;
node[15].gate++ <--> delay1 <--> node[20].gate++ ;
node[15].gate++ <--> delay1 <--> node[23].gate++ ;
node[16].gate++ <--> delay1 <--> node[17].gate++ ;
node[16].gate++ <--> delay1 <--> node[26].gate++ ;
node[20].gate++ <--> delay1 <--> node[21].gate++ ;
node[21].gate++ <--> delay1 <--> node[22].gate++ ;
node[22].gate++ <--> delay1 <--> node[23].gate++ ;
node[24].gate++ <--> delay1 <--> node[25].gate++ ;
node[25].gate++ <--> delay1 <--> node[26].gate++ ;
node[25].gate++ <--> delay1 <--> node[27].gate++ ;
node[25].gate++ <--> delay1 <--> node[28].gate++ ;
node[27].gate++ <--> delay1 <--> node[28].gate++ ;
node[11].gate++ <--> delay1 <--> node[10].gate++ ;
node[11].gate++ <--> delay1 <--> node[12].gate++ ;
node[5].gate++ <--> delay1 <--> node[30].gate++ ;
node[9].gate++ <--> delay1 <--> node[31].gate++ ;
node[18].gate++ <--> delay1 <--> node[32].gate++ ;
node[19].gate++ <--> delay1 <--> node[33].gate++ ;
node[21].gate++ <--> delay1 <--> node[34].gate++ ;
node[22].gate++ <--> delay1 <--> node[35].gate++ ;
node[24].gate++ <--> delay1 <--> node[36].gate++ ;
node[1].gate++ <--> delay1 <--> node[29].gate++ ;
node[28].gate++ <--> delay1 <--> node[37].gate++ ;
node[18].gate++ <--> delay1 <--> node[19].gate++ ;

}
