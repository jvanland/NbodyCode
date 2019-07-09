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
    Doub x_rel,y_rel,z_rel,vx_rel,vy_rel,vz_rel;
    Doub R,Rsq,Rcube,d_frac;
    Doub pot_0,pot_1,fac_0,fac_star,fac_SMBH;

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
    //From Itoh + Futumase + Asada 2001 Eq. 7.6
    x_rel = y[0]-y[1];
    y_rel = y[2]-y[3];
    z_rel = y[4]-y[5];
    vx_rel = y[6]-y[7];
    vy_rel = y[8]-y[9];
    vz_rel = y[10]-y[11];
    Rsq = SQ(x_rel)+SQ(y_rel)+SQ(z_rel);
    R = sqrt(Rsq);
    Rcube = Rsq*R;
    pot_0 = G*mass[0]/R;
    pot_1 = G*mass[1]/R;
    //0 PN term
    fac_0 = -G*mass[0]*mass[1]/Rcube;
    xforce[0] += fac_0*x_rel;
    yforce[0] += fac_0*y_rel;
    zforce[0] += fac_0*z_rel;
    //These 2 forces are equal and opposite
    xforce[1] += -xforce[0];
    yforce[1] += -yforce[0];
    zforce[1] += -zforce[0];
#ifdef PN_1
    //1 PN term
    Doub vsq_0,vsq_1,v0dotv1,rdotv0,rdotv1,fac_1;
    Doub rad1_term_0,rad1_term_1,vel1_term_0,vel1_term_1;
    vsq_0 = y[6]*y[6]+y[8]*y[8]+y[10]*y[10];
    vsq_1 = y[7]*y[7]+y[9]*y[9]+y[11]*y[11];
    v0dotv1 = y[6]*y[7]+y[8]*y[9]+y[10]*y[11];
    rdotv0 = (x_rel*y[6]+y_rel*y[8]+z_rel*y[10]);
    rdotv1 = (x_rel*y[7]+y_rel*y[9]+z_rel*y[11]);
    fac_1 = G*mass[0]*mass[1]/(c_speed*c_speed*Rcube);
    rad1_term_0 = -vsq_0-2.0*vsq_1+4.0*v0dotv1+1.5*rdotv1*rdotv1/Rsq+5.0*pot_0+4.0*pot_1;
    rad1_term_1 = -vsq_1-2.0*vsq_0+4.0*v0dotv1+1.5*rdotv0*rdotv0/Rsq+5.0*pot_1+4.0*pot_0;
    vel1_term_0 = 4.0*rdotv0-3.0*rdotv1;
    vel1_term_1 = 4.0*rdotv1-3.0*rdotv0;
    xforce[0] += fac_1*(rad1_term_0*x_rel+vel1_term_0*vx_rel);
    yforce[0] += fac_1*(rad1_term_0*y_rel+vel1_term_0*vy_rel);
    zforce[0] += fac_1*(rad1_term_0*z_rel+vel1_term_0*vz_rel);
    xforce[1] -= fac_1*(rad1_term_1*x_rel+vel1_term_1*vx_rel);
    yforce[1] -= fac_1*(rad1_term_1*y_rel+vel1_term_1*vy_rel);
    zforce[1] -= fac_1*(rad1_term_1*z_rel+vel1_term_1*vz_rel);
#endif
#ifdef PN_25
    //2.5 PN term
    Doub vdotr,Vsq,c_five,fac_25;
    Doub rad25_term_0,rad25_term_1,vel25_term_0,vel25_term_1;
    vdotr = x_rel*vx_rel+y_rel*vy_rel+z_rel*vz_rel;
    Vsq = SQ(vx_rel)+SQ(vy_rel)+SQ(vz_rel);
    c_five = c_speed*c_speed*c_speed*c_speed*c_speed;
    fac_25 = 4.0*G*G*mass[0]*mass[1]/(5.0*c_five*Rcube);
    rad25_term_0 = (-6.0*pot_0+(52.0/3.0)*pot_1+3.0*Vsq)*vdotr/Rsq;
    rad25_term_1 = (-6.0*pot_1+(52.0/3.0)*pot_0+3.0*Vsq)*vdotr/Rsq;
    vel25_term_0 = 2.0*pot_0-8.0*pot_1-Vsq;
    vel25_term_1 = 2.0*pot_1-8.0*pot_0-Vsq;
    xforce[0] += mass[0]*fac_25*(rad25_term_0*x_rel+vel25_term_0*vx_rel);
    yforce[0] += mass[0]*fac_25*(rad25_term_0*y_rel+vel25_term_0*vy_rel);
    zforce[0] += mass[0]*fac_25*(rad25_term_0*z_rel+vel25_term_0*vz_rel);
    xforce[1] -= mass[1]*fac_25*(rad25_term_1*x_rel+vel25_term_1*vx_rel);
    yforce[1] -= mass[1]*fac_25*(rad25_term_1*y_rel+vel25_term_1*vy_rel);
    zforce[1] -= mass[1]*fac_25*(rad25_term_1*z_rel+vel25_term_1*vz_rel);
#endif
    //For each member of the binary
    for (i=0;i<2;i++) {
      //Force from SMBH (at [0,0])
      Rsq = SQ(y[i])+SQ(y[2+i])+SQ(y[4+i]);
      Rcube = Rsq*sqrt(Rsq);
      fac_SMBH = -G*mass[i]*M_SMBH/Rcube;
      xforce[i] += fac_SMBH*y[i];
      yforce[i] += fac_SMBH*y[2+i];
      zforce[i] += fac_SMBH*y[4+i];
      //Force from stars
      for (j=0;j<N_cur;j++){
	Rsq = SQ(y[i]-posx[j])+SQ(y[2+i]-posy[j])+SQ(y[4+i]-posz[j]);
#ifdef CLOSE
	if (sqrt(Rsq) < CLOSE) {
	  cout << "Star " << j << " within " << sqrt(Rsq) << " AU of binary member " << i << "!" << endl;
	}
#endif
	Rcube = Rsq*sqrt(Rsq);
	fac_star = -G*mass[i]*M_star/Rcube;
	xforce[i] += fac_star*(y[i]-posx[j]);
	yforce[i] += fac_star*(y[2+i]-posy[j]);
	zforce[i] += fac_star*(y[4+i]-posz[j]);
      }
      //Set the derivatives of the velocities equal to the sum of the accelerations
      dydx[6+i] = (1/mass[i])*xforce[i];
      dydx[8+i] = (1/mass[i])*yforce[i];
      dydx[10+i] = (1/mass[i])*zforce[i];
    }
  }
};
