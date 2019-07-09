/* 
Returns the velocities and accelerations of members of a binary due to the gravity from many stars
    and a SMBH.
To instantiate give:
    The current number of stars, the mass of the SMBH, the time at beginning and end of the current run,
    the masses of the members of the binary, and the position of stars at beginning and end of the current run.
To run give:
    The current time, the positions and velocities of the members of the binary,
    and the structure to place the derivatives of the positions and velocities into.
 */
struct gravity {
  int N_cur; //Current number of stars
  Doub M_SMBH; //Mass of SMBH
  Doub T_init; //Time at beginning of current run
  Doub T_fin; //Time at end of current run
  VecDoub mass; //Masses of members of binary
  VecDoub part_1; //Positions of stars at beginning of current run
  VecDoub part_2; //Positions of stars at end of current run
  gravity(int N_curr, Doub M_SMBHH, Doub T_initt, Doub T_finn, VecDoub mmass, VecDoub ppart_1, VecDoub ppart_2) : 
    N_cur(N_curr), M_SMBH(M_SMBHH), T_init(T_initt), T_fin(T_finn), mass(mmass), part_1(ppart_1), part_2(ppart_2) {}
  void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
    //x = Current time, y = positions and velocities of binary members, dydx = derivatives
    int i,j;
    VecDoub posx(N_cur),posy(N_cur),posz(N_cur);
    VecDoub xforce(2,0.0),yforce(2,0.0),zforce(2,0.0);
    Doub M_star = 1.0; //Stars are 1 M_sun
    Doub Rsq,Rcube,d_frac;
    Doub Lx,Ly,Lz,Lsq,Rfive;

    //GR precession
    //From potential V = - G(M+m)L^2/c^2*mu*r^3
    //Force/mass = -G(M+m)((r x v)^2/c^2)(3*rvec/r^5)
    //c.o.m.?

    //Interpolate positions of stars
    d_frac = (x-T_init)/(T_fin-T_init);
    for (i=0;i<N_cur;i++) {
      posx[i] = part_1[i]+(part_2[i]-part_1[i])*d_frac;
      posy[i] = part_1[N_cur+i]+(part_2[N_cur+i]-part_1[N_cur+i])*d_frac;
      posz[i] = part_1[2*N_cur+i]+(part_2[2*N_cur+i]-part_1[2*N_cur+i])*d_frac;
    }

    for (i=0;i<6;i++) {
      //set x, y, and z derivatives equal to their velocities for both particles
      dydx[i] = y[6+i];
      //set vx, vy, and vz derivatives equal to 0 for both particles
      dydx[6+i] = 0;
    }
    //Gravitational force from other member of binary
    Rsq = SQ(y[0]-y[1])+SQ(y[2]-y[3])+SQ(y[4]-y[5]);
    Rcube = Rsq*sqrt(Rsq);
    xforce[0] += -G*mass[0]*mass[1]*(y[0]-y[1])/Rcube;
    yforce[0] += -G*mass[0]*mass[1]*(y[2]-y[3])/Rcube;
    zforce[0] += -G*mass[0]*mass[1]*(y[4]-y[5])/Rcube;
    //GR precessional force from other member of binary
    /*Lx = (y[2]-y[3])*(y[10]-y[11])-(y[4]-y[5])*(y[8]-y[9]);
    Ly = (y[4]-y[5])*(y[6]-y[7])-(y[0]-y[1])*(y[10]-y[11]);
    Lz = (y[0]-y[1])*(y[8]-y[9])-(y[2]-y[3])*(y[6]-y[7]);
    Lsq=Lx*Lx+Ly*Ly+Lz*Lz;
    Rfive = Rcube*Rsq;
    xforce[0] += -3*G*mass[0]*(mass[0]+mass[1])*Lsq*(y[0]-y[1])/(c_speed*c_speed*Rfive);
    yforce[0] += -3*G*mass[0]*(mass[0]+mass[1])*Lsq*(y[2]-y[3])/(c_speed*c_speed*Rfive);
    zforce[0] += -3*G*mass[0]*(mass[0]+mass[1])*Lsq*(y[4]-y[5])/(c_speed*c_speed*Rfive);*/
    //These 2 forces are equal and opposite
    xforce[1] += -xforce[0];
    yforce[1] += -yforce[0];
    zforce[1] += -zforce[0];
    //For each member of the binary
    for (i=0;i<2;i++) {
      //Force from SMBH (at [0,0])
      Rsq = SQ(y[i])+SQ(y[2+i])+SQ(y[4+i]);
      Rcube = Rsq*sqrt(Rsq);
      xforce[i] += -G*mass[i]*M_SMBH*y[i]/Rcube;
      yforce[i] += -G*mass[i]*M_SMBH*y[2+i]/Rcube;
      zforce[i] += -G*mass[i]*M_SMBH*y[4+i]/Rcube;
      //Force from stars
      for (j=0;j<N_cur;j++){
	Rsq = SQ(y[i]-posx[j])+SQ(y[2+i]-posy[j])+SQ(y[4+i]-posz[j]);
	Rcube = Rsq*sqrt(Rsq);
	xforce[i] += -G*mass[i]*M_star*(y[i]-posx[j])/Rcube;
	yforce[i] += -G*mass[i]*M_star*(y[2+i]-posy[j])/Rcube;
	zforce[i] += -G*mass[i]*M_star*(y[4+i]-posz[j])/Rcube;
      }
      //Set the derivatives of the velocities equal to the sum of the accelerations
      dydx[6+i] = (1/mass[i])*xforce[i];
      dydx[8+i] = (1/mass[i])*yforce[i];
      dydx[10+i] = (1/mass[i])*zforce[i];
    }
  }
};
