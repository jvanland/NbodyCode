/* 
To instantiate give an nbody6 output filename to read.
To run give:
   An integer for the output step to return values for.
   The integer NBODY6 name corresponding to the BHB particle being tracked.
If out < -2: Fills N_out, N_obj
Otherwise returns:
   Positions of stars at the output step and the next.
   Position and velocity of the BHB at the output step.
   Time in years at output step and output step+1.

When doing a run with multiple parts, will the different number of starting stars matter?
 */
struct nbody6 {
  string input; //Input file name
  int N_obj; //Initial number of objects (including BH's)
  int N_cur; //Current number of stars
  Doub T_init,T_fin; //Time at out and out+1 in yrs
  int N_out; //Total number of outputs in file
  VecDoub part_1; //position of stars at out in AU relative to SMBH
  VecDoub part_2; //position of stars at out+1 in AU relative to SMBH
  VecDoub BHB; //NEW position and velocity of BHB at out relative to SMBH
  nbody6(string inputt) : input(inputt) {}
  void operator() (int out, int BHB_NAME) {
    //Mtot = rec[3], Rconv = 1/(2.063e5*rec[2]), Vconv = 1/2*PI*sqrt(Rconv)*sqrt(Mtot), Tconv = Rconv/Vconv
    vector<int> first;
    vector<int> name;
    vector<float> rec;
    vector<float> m;
    vector<float> pos;
    NRvector<float> vel;
    VecDoub SMBH;
    BHB.assign(7,0); //Assign size of BHB vector
    SMBH.assign(6,0);
    first.assign(7,0);
    name.assign(1,0);
    rec.assign(20,0);
    m.assign(1,0);
    pos.assign(1,0);
    vel.assign(1,0);
    //vector<int> first(7);
    //vector<float> rec(20);
    //vector<float> m(1);
    //vector<float> pos(1);
    //vector<float> vel(1);
    //vector<int> name(1);
    char buff[60000];
    streampos pos1,pos2;
    int N_tot,N_pairs,i,j;
    Doub Mtot,Rconv,Tconv,Vconv;
    //VecDoub SMBH(6);
    VecDoub pos1_temp,pos2_temp;

    //Open input file
    ifstream initfile;
    initfile.open(input.c_str(), ios::binary);
    //Check to see if file exists
    if(initfile.fail()) {
      cerr << "Unable to open NBODY6 file " << input << "!" << endl; 
    }
    N_out = 0;
    initfile.seekg(0,ios::beg);
    while(initfile) {
      //Read first 7 values into first
      initfile.read(buff, 7*sizeof(float));
      for (i=0;i<7;i++) {
	first[i] = *reinterpret_cast<int*>(buff+i*sizeof(float));
      }
      N_tot = first[1];
      //Set N_obj = N_tot at the first output
      if (N_out == 0) {
	N_obj = N_tot; 
	pos1_temp.assign(3*N_obj,0);
	pos2_temp.assign(3*N_obj,0);
      }
      //resize vectors to current number of objects
      m.resize(N_tot);
      pos.resize(3*N_tot);
      vel.resize(3*N_tot);
      name.resize(N_tot);
      //Fill values at selected output step
      if (N_out == out) {
	//Read next 20 values into rec
	initfile.read(buff,20*sizeof(float));
	for (i=0;i<20;i++) {
	  rec[i] = *reinterpret_cast<float*>(buff+i*sizeof(float));
	}
	N_pairs = (int) rec[1];
	//Current number of stars not including BH's or pairs
	N_cur = N_tot-N_pairs-2; 
	Mtot = double(round(rec[3]));
	Rconv = 1/(2.063e5*double(rec[2]));
	Vconv = 1.0/(2.0*PI*sqrt(Rconv*Mtot));
	Tconv = Rconv/Vconv;
	//Time at out in yrs
	T_init = double(rec[0])/Tconv;
	//Read masses of objects
      	initfile.read(buff,N_tot*sizeof(float));
	for (i=0;i<N_tot;i++) {
	  m[i] = *reinterpret_cast<float*>(buff+i*sizeof(float));
	}
	//Read position of objects
	initfile.read(buff,3*N_tot*sizeof(float));
	for (i=0;i<3*N_tot;i++) {
	  pos[i] = *reinterpret_cast<float*>(buff+i*sizeof(float));
	}
	//Read velocities of objects
	initfile.read(buff,3*N_tot*sizeof(float));
	for (i=0;i<3*N_tot;i++) {
	  vel[i] = *reinterpret_cast<float*>(buff+i*sizeof(float));
	}
	//Read names of objects
	initfile.read(buff,N_tot*sizeof(float));
	for (i=0;i<N_tot;i++) {
	  name[i] = *reinterpret_cast<int*>(buff+i*sizeof(float));
	}
	initfile.seekg(4,ios::cur);
	for (i=0;i<N_tot;i++) {
	  //SMBH position and velocity
	  if(name[i] == 1) {
	    SMBH[0] = double(pos[3*i]);
	    SMBH[1] = double(pos[3*i+1]);
	    SMBH[2] = double(pos[3*i+2]);
	    SMBH[3] = double(vel[3*i]);
	    SMBH[4] = double(vel[3*i+1]);
	    SMBH[5] = double(vel[3*i+2]);
	  }
	  //BHB position and velocity
	  if(name[i] == BHB_NAME) {
	    BHB[0] = double(m[i]/m[i-1]);
	    BHB[1] = double(pos[3*i]-SMBH[0])/Rconv;
	    BHB[2] = double(pos[3*i+1]-SMBH[1])/Rconv;
	    BHB[3] = double(pos[3*i+2]-SMBH[2])/Rconv;
	    BHB[4] = double(vel[3*i]-SMBH[3])/Vconv;
	    BHB[5] = double(vel[3*i+1]-SMBH[4])/Vconv;
	    BHB[6] = double(vel[3*i+2]-SMBH[5])/Vconv;
	  }
	}
	for (i=0;i<N_tot;i++) {
	  if (name[i] < N_obj && name[i] != 1) { //Pairs and BH's not allowed
	    //Converts to (x1,x2...xn,y1,y2...) from (x1,y1,z1,x2,y2...)
	    //Converts to name order from nbody6 order
	    //Converts to AU relative to SMBH
	    pos1_temp[name[i]] = (double(pos[3*i])-SMBH[0])/Rconv;
	    pos1_temp[N_obj+name[i]] = (double(pos[3*i+1])-SMBH[1])/Rconv;
	    pos1_temp[2*N_obj+name[i]] = (double(pos[3*i+2])-SMBH[2])/Rconv;
	  }
	}
      }
      //Fill values at step after selected output step
      else if (N_out == out+1) {
	//Read next 20 values into rec 
	initfile.read(buff,20*sizeof(float));
	for (i=0;i<20;i++) {
	  rec[i] = *reinterpret_cast<float*>(buff+i*sizeof(float));
	}
	N_pairs = (int) rec[1];
	//Current number of stars not including BH's or pairs
	N_cur = N_tot-N_pairs-2;
	Mtot = double(round(rec[3]));
	Rconv = 1/(2.063e5*double(rec[2]));
	Tconv = 2.0*PI*Rconv*sqrt(Rconv)*sqrt(Mtot);
	//Time at out+1 in yrs
	T_fin = double(rec[0])/Tconv;
	initfile.seekg(N_tot*sizeof(float),ios::cur);
	//Read position of objects
	initfile.read(buff,3*N_tot*sizeof(float));
	for (i=0;i<3*N_tot;i++) {
	  pos[i] = *reinterpret_cast<float*>(buff+i*sizeof(float));
	}
	initfile.seekg(3*N_tot*sizeof(float),ios::cur);
	//Read names of objects
	initfile.read(buff,N_tot*sizeof(float));
	for (i=0;i<N_tot;i++) {
	  name[i] = *reinterpret_cast<int*>(buff+i*sizeof(float));
	}
	initfile.seekg(4,ios::cur);
	//SMBH position
	for (i=0;i<N_tot;i++) {
	  if(name[i] == 1) {
	    SMBH[0] = double(pos[3*i]);
	    SMBH[1] = double(pos[3*i+1]);
	    SMBH[2] = double(pos[3*i+2]);
	  }
	}
	for (i=0;i<N_tot;i++) {
	  if (name[i] < N_obj && name[i] != 1) { //Pairs and BH's not allowed
	    //Converts to (x1,x2...xn,y1,y2...) from (x1,y1,z1,x2,y2...)
	    //Converts to name order from nbody6 order
	    //Converts to AU relative to SMBH
	    pos2_temp[name[i]] = (double(pos[3*i])-SMBH[0])/Rconv;
	    pos2_temp[N_obj+name[i]] = (double(pos[3*i+1])-SMBH[1])/Rconv;
	    pos2_temp[2*N_obj+name[i]] = (double(pos[3*i+2])-SMBH[2])/Rconv;
	  }
	}
	part_1.assign(3*N_cur,0);
	part_2.assign(3*N_cur,0);
	j=0;
	//Weeds out particles that have been ejected
	for (i=0;i<N_obj;i++) {
	  if (pos2_temp[i] != 0) {
	    part_1[j] = pos1_temp[i];
	    part_1[N_cur+j] = pos1_temp[N_obj+i];
	    part_1[2*N_cur+j] = pos1_temp[2*N_obj+i];
	    part_2[j] = pos2_temp[i];
	    part_2[N_cur+j] = pos2_temp[N_obj+i];
	    part_2[2*N_cur+j] = pos2_temp[2*N_obj+i];
	    j++;
	  }
	}
      }
      //Scan through outputs to get to out or to find N_out
      else if (N_out < out || out < 0) {
	//	initfile.seekg(20*sizeof(float),ios::cur);
	initfile.read(buff,20*sizeof(float));
	initfile.seekg(N_tot*sizeof(float),ios::cur);
	initfile.seekg(3*N_tot*sizeof(float),ios::cur);
	initfile.seekg(3*N_tot*sizeof(float),ios::cur);
	initfile.seekg(N_tot*sizeof(float),ios::cur);
	//Check to see how close to the end of the file we are
	pos1 = initfile.tellg();
	initfile.seekg(0,ios::end);
	pos2 = initfile.tellg();
	//If we are at the end of the file break out of the loop
	if (pos2-pos1 <= 4) break;
	//Otherwise keep going
	initfile.seekg(pos1,ios::beg);
	initfile.seekg(4,ios::cur);
      }
      //After reading out and out+1 break
      else break;
      ++N_out;
    }
    //free (first);
    //free (rec);
    //free (m);
    //free (pos);
    //free (vel);
    //free (name);
    //free (SMBH);
    initfile.close();
  }
};

