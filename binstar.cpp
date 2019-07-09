/* Integrate the 2-body problem under the influence
of many stars and a central BH.

Units of AU, yr, Msun. 
All stars 1 M_sun.

Usage:
./binstar RUN_TYPE(1,2,3) M_SMBH(double) N_FILES(int) ERR_TOL(double)
RUN_TYPE 1 = normal start: no further input needed. Reads from fulllist.dat
RUN_TYPE 2 = restart at last output: Reads from starout_old.dat. Additionally requires...
     NB6_STARTFILE_NUM(int) NB6_OUTNUM(int)
RUN_TYPE 3 = custom: additionally requires...
     NB6_STARTFILE_NUM(int) NB6_OUTNUM(int) INPUT_FILENAME(string) INPUT_LINENUM(int)

Requires Nbody6 output files "OUT3_1,OUT3_2..."

Binary initial conditions format:
    Mass1 x1 y1 z1 vx1 vy1 vz1
    Mass2 x2 y2 z2 vx2 vy2 vz2
For rerun (assumes masses are 10 M_sun):
    Time x1 y1 z1 vx1 vy1 vz1
    Time x2 y2 z2 vx2 vy2 vz2

Each NBODY6 file: 1st and last output skipped (restart goes back 2 steps).
1st Nbody6 output is output 0 so N_out-1 outputs read.
NB6_OUTNUM is the NEXT nb6 output!

INPUT_LINENUM = 2[NB6_OUTNUM + (N_previous-1) + (N_previous-1) + ...] - 1

Numerical recipes output structure does not work well with current setup for high precision. Fixed (possibly).
*/

#define DIST_CONTROL 0.000001 //Exit integration if binary separation less than this in AU or greater than 100 AU
#define PN_1 //Include 1 PN term in gravity calculation 
#define PN_25 //Include 2.5 PN term in gravity calculation
//#define FIN 8 //Alternate number of steps to run for
#define HOME //Alternate location of OUT3 files (for personal tests)
//#define CLOSE 0.2 //Print close approaches within this distance to output

#define SQ(x) ((x)*(x))
#define G 39.47841760435743 /* Newton's constant for AU,yr,Msun*/
#define c_speed 6.3241e4 /*Speed of light for Au,yr,Msun */
#define PI 3.141592653589793

#include "nr3_precise.h" //Numerical recipes c++ header
#include "stepper.h" //NR stepper header
#include "odeint.h" //NR wrapper to the integrator
#include "stepperdopr853.h" //NR Runge-Kutta integrator
#include "nb6_new.h" //Reads Nbody6 output
#include "gravity_GR.h" //Calculates gravity
using namespace std;

int main(int argc,char *argv[ ]) {
  if (argc <= 4) {
    cerr << "Usage: " << argv[0] << " RUN_TYPE(1,2,3) M_SMBH(double) N_FILES(int) ERR_TOL(double)" << endl;
    cerr << "RUN_TYPE 1 = normal start: no further input needed." << endl;
    cerr << "RUN_TYPE 2 = restart at last output: additionally requires..." << endl;
    cerr << "     NB6_STARTFILE_NUM(int) NB6_OUTNUM(int)" << endl;
    cerr << "RUN_TYPE 3 = custom: additionally requires..." << endl;
    cerr << "     NB6_STARTFILE_NUM(int) NB6_OUTNUM(int) INPUT_FILENAME(string) INPUT_LINENUM(int)" << endl;
    return 1;
  }

  cout << "Running..." << endl;
  int RUN_TYPE = atoi(argv[1]); //Normal run or rerun
  double M_SMBH = atof(argv[2]); //Mass of SMBH
  int N_FILES = atoi(argv[3]); //Number of NBODY6 output files
  double ERR = atof(argv[4]); //Integrator error tolerance
  cout << "RUN_TYPE = " << RUN_TYPE << ", M_SMBH = " << M_SMBH << ", N_FILES = " << N_FILES << ", ERR = " << ERR << endl;

  int i,j,k,beg,fin,BHB_NAME;
  ostringstream file;
  VecDoub mass(2); //Masses of binary members
  VecDoub var(12); //Contains [x1,x2,y1,y2,z1,z2,vx1,vx2,vy1,vy2,vz1,vz2]
  VecDoub tempcm(6);
  double time;
  Doub atol = 0.0; //Integrator absolute error tolerance
  Doub rtol = ERR; //Integrator relative error tolerance
  Doub step = 0.1; //Initial integrator step size
  Int count = 1; //Use Output out(count) for more output
  Output out(count); //Numerical recipes output structure, currently saves first and last output
  string BIN_INIT; //Binary initial conditions file
  int line_num; //Read BIN_INIT from this line #
  int NB6_INITFILE; //NBODY6 file # to start from
  int start; //NBODY6 output to start from
  ifstream init;

  clock_t nb6_t, gravity_t, integrate_t, total_t;
  double nb6_total=0, gravity_total=0, integrate_total=0, total_total=0;

//Normal start
  if (RUN_TYPE == 1) {
    if (argc != 5) {
      cerr << "Incorrect number of inputs!" << endl;
      return 1;
    }
    BIN_INIT = "fulllist.dat";
    line_num = 1;
    NB6_INITFILE = start = 1;
    //Read in binary initial conditions
    init.open(BIN_INIT.c_str());
    if(init.fail()) {
      cerr << "Unable to open" << BIN_INIT << "!" << endl;
      return 1;
    }
    i = 0;
    time = 0.0;
    while(i<2) {
      init >> mass[i] >> var[i] >> var[2+i] >> var[4+i] >> var[6+i] >> var[8+i] >> var[10+i];
      ++i;
    }
  }

  //Restart at last output
  else if (RUN_TYPE == 2) {
    if (argc != 7) {
      cerr << "Incorrect number of inputs!" << endl;
      return 1;
    }
    BIN_INIT = "starout_old.dat";
    line_num = 0;
    NB6_INITFILE = atoi(argv[5]);
    start = atoi(argv[6]);
    //Read in binary initial conditions
    init.open(BIN_INIT.c_str());
    if(init.fail()) {
      cerr << "Unable to open" << BIN_INIT << "!" << endl;
      return 1;
    }
    //Read last 2 lines of file
    string temp;
    while (getline(init,temp)) {
        ++line_num;
    }
    line_num -= 1;
    init.clear();
    init.seekg(0,ios::beg);
    for (i=1;i<line_num;i++) {
      getline(init,temp);
    }
    i = 0;
    while(i<2) {
      init >> time >> var[i] >> var[2+i] >> var[4+i] >> var[6+i] >> var[8+i] >> var[10+i];
      mass[i] = 10;
      ++i;
    }
  }

  //Custom restart
  else if (RUN_TYPE == 3) {
    if (argc != 9) {
      cerr << "Incorrect number of inputs!" << endl;
      return 1;
    }
    BIN_INIT = argv[7];
    line_num = atoi(argv[8]);
    NB6_INITFILE = atoi(argv[5]);
    start = atoi(argv[6]);
    //Read in binary initial conditions
    init.open(BIN_INIT.c_str());
    if(init.fail()) {
      cerr << "Unable to open" << BIN_INIT << "!" << endl;
      return 1;
    }
    string temp;
    for (i=1;i<line_num;i++) {
      getline(init,temp);
    }
    i = 0;
    while(i<2) {
      init >> time >> var[i] >> var[2+i] >> var[4+i] >> var[6+i] >> var[8+i] >> var[10+i];
      mass[i] = 10;
      ++i;
    }
  }
  cout << "Read " << BIN_INIT << " from output " << line_num << endl;
  init.close();

  //Write initial conditions to output
  ofstream output("starout.dat");
  output.precision(16); //Set precision of output
  output << time << " " << var[0] << " " << var[2] << " " << var[4] << " " 
	 << var[6] << " " << var[8] << " " << var[10] << endl;
  output << time << " " << var[1] << " " << var[3] << " " << var[5] << " " 
	 << var[7] << " " << var[9] << " " << var[11] << endl;

  //Cycle through NBODY6 output files
  for (k=NB6_INITFILE;k<=N_FILES;k++) {
    file.str("");
#ifdef HOME
    file << "OUT3_" << k;
#else
    file << "../OUT3_" << k;
#endif
    cout << "Opening " << file.str() << endl;
    nbody6 nb6(file.str());
    //Determine number of outputs in Nbody6 file
    nb6(-100,0);
    cout << nb6.N_out << " outputs" << endl;
    //Set BHB_NAME to N_obj unless otherwise required
    BHB_NAME = nb6.N_obj;
    //For each Nbody6 output
    if (k == NB6_INITFILE) beg = start;
    else beg = 1;
#ifdef FIN
    fin = FIN;
#else
    fin = nb6.N_out;  
#endif
    cout << "Starting at output " << beg << endl;
    for (i=beg;i<fin;i++) {
      total_t = clock();
      nb6_t = clock();
      //Find positions of stars at beginning and end of step
      //Also find current NBODY6 COM position of BHB
      nb6(i,BHB_NAME);
      nb6_t = clock()-nb6_t;
      //Instantiate gravity structure with star data
      gravity_t = clock();
      gravity g(nb6.N_cur,M_SMBH,nb6.T_init,nb6.T_fin,mass,nb6.part_1,nb6.part_2);
      //After each step replace the BHB at the NBODY6 COM position
      for (j=0;j<6;j++) {
	//COM pos and vel of BHB
	tempcm[j]=(mass[0]*var[2*j]+mass[1]*var[2*j+1])/(mass[0]+mass[1]);
	//Reposition
	var[2*j]+=nb6.BHB[j+1]-tempcm[j];
	var[2*j+1]+=nb6.BHB[j+1]-tempcm[j];
      }
      //Instantiate integrator
      Odeint<StepperDopr853<gravity> > ode(var,nb6.T_init,nb6.T_fin,atol,rtol,step,0.0,out,g);
      gravity_t = clock()-gravity_t;
      //Integrate!
      integrate_t = clock();
      ode.integrate();
      integrate_t = clock()-integrate_t;
      //New trial step is last integration step
      step = ode.h;
      //Print to output
      output << ode.x << " " << var[0] << " " << var[2] << " " << var[4] << " " 
	     << var[6] << " " << var[8] << " " << var[10] << endl;
      output << ode.x << " " << var[1] << " " << var[3] << " " << var[5] << " " 
	     << var[7] << " " << var[9] << " " << var[11] << endl;
      cout << "NB6_OUTFILE_NUM = " << k << " NB6_OUTNUM = " << i+1 << endl;
      //Exit loop if integrator returns before normal end time
      if(ode.x < nb6.T_fin) {
	cout << "Integrator exited early!" << endl;
	return 0;
      }
      total_total += (double)total_t;
      nb6_total += (double)nb6_t;
      gravity_total += (double)gravity_t;
      integrate_total += (double)integrate_t;
      cout << "Total: " << total_total/CLOCKS_PER_SEC << " NB6: " << nb6_total/CLOCKS_PER_SEC << " Gravity: " << gravity_total/CLOCKS_PER_SEC << " Integrate: " << integrate_total/CLOCKS_PER_SEC << endl;
    }
  }
  output.close();
  cout << "Exited normally" << endl;

  return 0;
}
