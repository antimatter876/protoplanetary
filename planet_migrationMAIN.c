/* MODELS PLANET MIGRATION IN alpha-DISC, WITH PHOTOEVAPORATION   */
/*                                                                */
/* (c) Richard Alexander - 2004-2019                              */
/*                                                                */
/* -------------------------------------------------------------  */
/*                                                                */
/* SOLVES FOLLOWING EQUATIONS:                                    */
/*                                                                */
/* (1)  dSigma/dt = (3/R)d/dR[ sqrt(R)d/dR(nu Sigma sqrt(R))      */
/*          - (2 Lambda Sigma R^1.5)/(sqrt(G*M)) ] - Sigmadot(R)  */
/*   nu = alpha c_s H                                             */
/*                                                                */
/* CODE UNITS ARE Msun, AU & yr   (G = 4pi^2)                     */
/*                                                                */
/* GRID IS EQUISPACED IN R^(1/2)                                  */
/* Eq (1) SOLVED USING FIRST-ORDER EXPLCIT METHOD                 */
/*                                                                */
/* VARIABLE CHANGES:                                              */
/*        X = 2*sqrt(R)                                           */
/*        Y = 3*nu*Sigma*sqrt(R)                                  */
/*        S = Sigma*R^(3/2)                                       */
/*                                                                */
/*  => Eq (1) becomes   dS/dt = d2Y/dX2 + (torque) + (wind)       */
/*                                                                */
/* Nov 2019: THIS IS A STRIPPED-DOWN VERSION DESIGNED TO HANDLE   */
/*           TYPE I MIGRATION, NOT TYPE II. NO BACK-REACTION OF   */
/*           PLANET ON DISC (i.e., Lambda = 0 in Eq (1) ).        */
/*                                                                */
/* COMPILE WITH gcc, SYNTAX:                                      */
/*           gcc planet_migration.c -o planet_migration.e -lm     */
/* THEN RUN VIA: ./planet_migration.e                             */
/*                                                                */
/* OUTPUT FILES ARE:                                              */
/*     sigma_???_x.dat - R, Sigma(R) AND Lambda(R) FOR DUMP x     */
/*     mdot_???.dat    - GLOBAL VARIABLES, ONE LINE PER DUMP      */
/*     planet_ICs.dat  - PLANET INTIAL CONDITIONS                 */
/*     planet_data.dat - PLANET FINAL CONDITIONS                  */ 
/*                                                                */
/*                                                                */
/* -------------------------------------------------------------  */


/* STANDARD LIBRARY FILES */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

/* MAX AND MIN FUNCTIONS */
double max(double a,double b)
{
  double result;
  if (a < b) 
    result = b;
  else
    result = a;

  return result;
}

double min(double a,double b)
{
  double result;
  if (a > b) 
    result = b;
  else
    result = a;

  return result;
}

/* RANDOM NUMBER GENERATOR */
double rand_num(double min, double max)
{
    double random_num;
       
    random_num = ( (double)rand() / ((double)(RAND_MAX)) );   
   
    random_num = (random_num * (max-min)) + min;   // output # from [r_min, r_max]
   
    return (random_num);
}


/* MAIN PROGRAM */
int main()
{
  int i,j,k,t,in,out,mark,stop_flag,seed,gauss_flag;
  int in_low,out_low,in_high,out_high,res_switch,counter,factor;
  int edge,direct_flag,oneAU,pointoneAU;
  int dump_counter;
  /* CONSTANTS */
  static double M_sun=1.989E33, AU=1.496E13, year=31536000.0;
  /* MODEL PARAMETERS */ 
  static double alpha=0.01;
  static double H_R=0.0333333333333333;    /* H/R at 1AU */
  /* DISC ICs */
  static double gamma = 1.0;     /* Slope of initial Sigma(r)  */
  static double M_star = 1.0;    /* Msun - central object mass */
  static double R_s = 10.0;      /* AU - initial disc radius   */
  static double q=1.25;          /* Flaring parameter          */
  static double M_jup = 0.0009548;       /* SOLAR MASSES */
  static double log_M_d_mean = -1.5;     /* SOLAR MASSES, LOG UNITS */
  double M_d_0,M_d_0_log,x,f;
  /* PLANET PARAMETERS */
  double R_planet,M_planet,t_planet;
  double t_s,M_d_min,t_p_max,M_dot_min,t_p_max_alt;
  int planet_number,planet_flag;  
  /* WIND PARAMETERS */
  static double sigma_crit_cgs=1E-5;
  static double R_g=1.4;             /* AU */
  static double R_p_min=0.15;        /* AU */
  static double a = 0.3423, b = -0.3612, d = 0.2457;
  static double R_G = 8.9;         /* AU */
  double rho_g, r, v, rho, rho_in;
  static double c_s = 2.10802;
  static double mu=1.35;
  static double CD = 0.235, aa = -2.42;
  static double Phi = 10.0;        /* 10^41 photon/s */
  /* OUTER BOUNDARY CONDITION:            */
  /* 0 = ZERO-TORQUE ; 1 = REFLECTING     */
  static int BC_flag=0;
  /* COMPUTATIONAL PARAMS */
  static double epsilon=1E-10;
  static double Courant=0.05;
  static double T_min = 10.0;        /* K */
  double Courant_fix;
  /* GRID & VISC */
  double Rz[10001],X[10001],nu[10001],Rf[10001],XX[10001];
  /* GAS */
  double Y[10001][2],sigma[10001][2],S[10001][3];
  double H[10001],Omega[10001],sigmadot[10001];
  /* PLANET */
  double R_p;
  double alpha0; /*New */
  double beta0; /*New*/
  double Il;   /*New*/
  double c0;   /*New*/
  static double gammaI=1.4; /*New*/
  double Lambda[10001],sigma_visc[10001],Delta_p,M_p,Q,G;
  double da_dt,dR,I,dS_dt_p,dS_dt_nu,v_r_p,v_r,v_r_nu,Lambda_max,v_R_max;
  /* OTHERS */
  double tau,ttt,M_d_t,R_0;
  double M_dot[3],nu_s,M_acc,M_out,M_wind;
  double two_gamma,rad;
  double dX,dt,minstep,gstep,age,tstep;
  double disc_mass,wind_rate,sigma_crit,M_d;
  double Time,time_to_dump;
  double dX_low,dX_high,Rz_low,X_low,slope;
  double t_p_acc,mass_cell,Rz_high,Rf_high,Rf_high_p1;
  int transition_flag;
  double sigma_trans;
  static double sigma_trans_cgs=-1.5;   /* log UNITS */
  /* DATA DUMP FREQ IN yr */
  static double step=1.0E4;
  /* STOP TIME IN yr */
  static double t_max=3.0E7;
  char sigfile[100],mfile[100],root[100],resfile[100];
  FILE *sigmaout,*mdot,*restart,*data,*planet_data;
  /* NEW VARIABLE */
  int i_p;

  /* SET TIMESTEP CORRECTION FACTOR */
  Courant_fix = alpha/0.01;

  /* PLANET MASS AND INITIAL RADIUS */
  M_p = 0.0;
  R_p = 1E-15;
  Q = M_p/M_star;
  planet_number = 0;
  planet_flag=0;

  t_p_acc = 1.0E30;
 
  /* RANDOMLY SAMPLE PLANET PARAMETERS */
  /* INITIALISE RANDOM NUMBER GENERATOR WITH TIME OF DAY */
  seed = time(0);
  srand(seed);

  /* SET INITIAL DISC MASS */
  M_d_0 = pow(10.0,log_M_d_mean);
  
  /* PLANET FORMATION RADIUS (AU) */
  R_planet = 5.0;
  
  /* PLANET MASS (M_Jup) */
  M_planet = 3E-2;   /* ~10M_e */
      
  /* FORMATION TIME (yr) */
  t_planet = 4.6E6;


  /* WRITE INITIAL CONDITIONS TO DATA FILE */
  planet_data = fopen("planet_ICs.dat","a");
  fprintf(planet_data,"%g %g %g %g\n",M_planet,R_planet,t_planet,M_d_0);
  fclose(planet_data);


  /* GRAV CONSTANT */
  G = 4.0*(double)M_PI*(double)M_PI;

  /* "BASE" FOR MAKING OUTPUT FILENAMES */
  sprintf(root,"basic_model");

  /* SET FILENAMES */
  sprintf(mfile,"Files/mdot_%s.dat",root);  
  
  /* INITIALISE PROBLEM */
  /* OPEN OUTPUT FILES */
  mdot=fopen(mfile,"w");
  
  /* GRID SCALE  - dX=2dr^1/2  */
  /* INNER FACE IS AT 0.04AU (X=0.4); outer face at 10000AU (X=200) */
  dX = 0.2;
  in = 1;
  out = 1000;


  /* MAKE GRID */
  /* ZONE CENTRED */
  for (i=0; i<out+1; i++)
    {
      X[i]=(dX*((double)i+0.5));   /* UNITS ARE AU^0.5 */
      Rz[i]=(X[i]*X[i])/4.0;
    }
  /* FACE CENTRED */
  for (i=0; i<out+1; i++)
    {
      XX[i]=(dX*((double)i));  /* AU^0.5 */
      Rf[i]=(XX[i]*XX[i])/4.0;
    }
  
  
  /* DEFINE ORBITAL FREQUENCY */
  for (i=in; i<out+1; i++)
    {
      Omega[i] = (2.0*(double)M_PI)/(sqrt((Rz[i]*Rz[i]*Rz[i])/M_star)); /* yr^-1 */
    }
  
  
  /* PRINT RANGES */
  printf("Inner boundary (first non-zero cell) at %gAU\n",Rf[in+1]);
  if (BC_flag==0)
    printf("Zero-torque outer boundary (first non-zero cell) at %gAU\n",Rf[out]);
  if (BC_flag==1)
    printf("Reflecting outer boundary (first non-zero cell) at %gAU\n",Rf[out]);
  
  
  /* SLOPE OF INITIAL DISC */
  two_gamma = 2.0-gamma;
  
  /* SURFACE DENSITY */
  for (i=in+1; i<out+1; i++)
    {
      rad = Rz[i]/R_s;
      
      sigma[i][0]=( (M_d_0*two_gamma)/(2.0*(double)M_PI*R_s*R_s*pow(rad,gamma))) * exp(-1.0*pow(rad,two_gamma));  /* SOLAR MASSES/AU^2 */
      
    }
  
  /* MAKE DISC MASS EXACT */
  disc_mass=0.0;
  for (i=1; i<out; i++)
    {
      disc_mass = disc_mass + (2.0*(double)M_PI*Rz[i]*(Rf[i+1]-Rf[i])*sigma[i][0]);
    }
  for (i=in+1; i<out+1; i++)
    sigma[i][0] = sigma[i][0] * (M_d_0/disc_mass);
  
  
  /* SET VISCOSITY LAW */
  for (i=in; i<out+1; i++)
    {
      H[i] = H_R * pow(Rz[i],q);
      nu[i]= alpha * H[i] * H[i] * Omega[i];  /* AU^2/yr */
    }
  
  /* CREATE WIND PROFILE */
  wind_rate = 0.0;
  rho_g = 2.75E4 * sqrt(Phi) * mu * 2.81E-18 * pow((R_G/8.9),-1.5);   /* Msun/AU^3 */ 
  for (i=in; i<out+1; i++) 
    { 
      r = Rz[i]/R_G;
      if (r >= 0.1)
	{
	  rho = rho_g * pow( 2.0/(pow(r,7.5) + pow(r,12.5)) , 0.2);  /* Msun/AU^3 */
	  v = a * exp(b*(r-0.1)) * pow( (r-0.1) , d) * c_s;  /* AU/yr */   
	  sigmadot[i] = 2.0*rho*v;
	}
      else
	sigmadot[i] = 0.0;
      wind_rate = wind_rate + (2.0*(double)M_PI*Rz[i]*(Rf[i+1]-Rf[i])*sigmadot[i]);
    }
  printf("Wind rate = %gMun/yr\n",wind_rate);
  /* SET FLAG */
  direct_flag = 0;
  sigma_crit = sigma_crit_cgs * (AU*AU/M_sun);
  pointoneAU = 0;
  while (Rf[pointoneAU+1] < 0.05*R_g)
    pointoneAU++;
  oneAU = 0;
  while (Rf[oneAU+1] < 0.95*R_g)
    oneAU++;


  /* Y & S */
  for (i=in+1; i<out; i++)
    {
      Y[i][0]=3.0*nu[i]*sqrt(Rz[i])*sigma[i][0];
      S[i][0]=pow(Rz[i],1.5)*sigma[i][0];
    }
  
  
  /* INNER BC */
  sigma[in][0]=0.0;
  S[in][0]=0.0;
  Y[in][0]=0.0;
  
  /* OUTER BC */
  /* ZERO-TORQUE - Sigma=0 */
  if (BC_flag==0)
    {
      sigma[out][0]=0.0;
      S[out][0]=0.0;
      Y[out][0]=0.0;
    }
  /* REFLECTIVE - v_r=0 */
  if (BC_flag==1)
    {
      Y[out][0]=Y[out-1][0];
    }
    
  
  /* INITIAL ACCRETION RATE */
  M_dot[0]=2.0*(double)M_PI*((Y[in+1][0]-Y[in][0])/dX);
  /* INITIAL MASSES */
  disc_mass=0.0;
  /* ADD UP DISC MASS */
  for (i=1; i<out; i++)
    {
      disc_mass = disc_mass + (2.0*(double)M_PI*Rz[i]*(Rf[i+1]-Rf[i])*sigma[i][0]);
    }
  M_d = disc_mass;
  
  M_acc = 0.0;
  M_out = 0.0;
  M_wind = 0.0;

  /* OUTPUT ARRAYS */
  dump_counter=0;
  /* CREATE FILENAME */
  sprintf(sigfile,"Files/sigma_%s_%.5i.dat",root,dump_counter);
  /* OPEN FILE */
  sigmaout=fopen(sigfile,"w");
  /* DUMP OUTPUT */
  fprintf(sigmaout,"0 %g\n",R_p);
  fprintf(mdot,"0.0 %g %g %g 0.0 1.0\n",M_dot[0],disc_mass,R_p);
  for (i=0; i<out+1; i++)
    {
      fprintf(sigmaout,"%g  %g  0\n",Rz[i],sigma[i][0]);
    }      
  /* CLOSE FILE */
  fclose(sigmaout);
  /* INCREMENT COUNTER */
  dump_counter++;
 
  
  /* TIME EVOLUTION STARTS HERE */
  /* SET TIME VARIABLES */
  t=1;
  Time=0.0;
  time_to_dump=step;
  
/*   printf("Starting time evolution...\n"); */


  /* stop_flag ALLOWS CONDITIONAL ABORT */
  stop_flag=0;
  while((Time < t_max) && (stop_flag == 0))
    {

      /* TIME INTEGRATION LOOP */
      /* USES ARRAY STEP 0 TO FIND STEP 1 */
      /* THEN COPIES 1 TO 0 AT THE END AND ITERATES */
            
      /* ------------ STEP 1 - FIND NEW TIMESTEP ------------ */	  
      minstep=1E10;
      for (i=in+1; i<out; i++)
	{
	  /* FROM VISCOSITY */
	  gstep=Courant*dX*dX*(S[i][0]/Y[i][0]);
	  if (gstep < minstep)
	    minstep=gstep;	   
	}
       

      /* SET TIMESTEP */
      dt = minstep;

      /* UPDATE TIME COUNTER */
      Time=Time+dt;
      time_to_dump = time_to_dump - dt;

      
      /* ------------ STEP 2 - UPDATE WIND PROFILE -------- */
      if ((sigma[oneAU][0] < sigma_crit) && (sigma[pointoneAU][0] < sigma_crit))
	direct_flag = 1;
      else
	direct_flag = 0;

      /* DIFFUSE WIND PROFILE */
      if (direct_flag == 0)
	{
	  wind_rate = 0.0;
	  for (i=in; i<out+1; i++) 
	    { 
	      r = Rz[i]/R_G;
	      if (r >= 0.1)
		{
		  rho = rho_g * pow( 2.0/(pow(r,7.5) + pow(r,12.5)) , 0.2);  /* Msun/AU^3 */
		  v = a * exp(b*(r-0.1)) * pow( (r-0.1) , d) * c_s;  /* AU/yr */   
		  sigmadot[i] = 2.0*rho*v;
		}
	      else
		sigmadot[i] = 0.0;
	      wind_rate = wind_rate + (2.0*(double)M_PI*Rz[i]*(Rf[i+1]-Rf[i])*sigmadot[i]);
	    }
	}
      
      /* DIRECT WIND PROFILE */
      if (direct_flag == 1)
	{
	  
	  /* FIND INNER EDGE */
	  edge=in;
	  while (sigma[edge][0] < sigma_crit)
	    edge++;
	  /* STOP AT OUTER EDGE */
	  if (edge >= out-2)
	    stop_flag=1;

	  /* CREATE WIND PROFILE */
	  wind_rate = 0.0;
	  for (i=in; i<out+1; i++)
	    {
	      if (i >= edge)
		{
		  rho_in = 2.60E6 * CD * mu * 2.81E-18 * pow((3.0/Rz[edge]),1.5);   /* Msun/AU^3 */ 
		  sigmadot[i] = 2.0 * rho_in * c_s * pow((Rz[i]/Rz[edge]),aa);  /* AU/yr */
		  wind_rate = wind_rate + (2.0*(double)M_PI*Rz[i]*(Rf[i+1]-Rf[i])*sigmadot[i]);
		}
	    }
	  
	}
      
      
      
      
      
      /* ------------ STEP 3 - SURFACE DENSITY ------------ */
      /*  ALSO DOES PLANET TORQUE USING OPERATOR SPLITTING  */
      /* LOOP THROUGH GRID */
      for (i=in+1; i<out; i++)
	{

	  /* 3a - EVALUATE PLANET TORQUE ON DISC */
	  
	  
	  /* CURRENTLY BLANK (WAS ORIGINALLY WRITTEN FOR TYPE II) */
	  /****************************************/
	  /****************************************/
	  /* Lambda[i] = 0.0; */
	  /****************************************/
	  /****************************************/
	  


	  /* 3b - UPDATE Sigma & S DUE TO DIFFUSION */
	  /* USES OPERATOR SPLTTING METHOD */
	  dS_dt_nu = (1.0/(dX*dX)) * (Y[i-1][0]+Y[i+1][0] - 2.0*Y[i][0]);
	  S[i][1] = S[i][0] + (dS_dt_nu*dt);

	  /* UPDATE Sigma */
	  sigma_visc[i] = S[i][1]/(pow(Rz[i],1.5));

	  /* 3c - ADD MASS-LOSS DUE TO WIND */
	  S[i][1] = S[i][1] - (sigmadot[i] * pow(Rz[i],1.5) * dt);
	  /* FLOOR FUNCTION (to prevent negative densities) */
	  if (S[i][1] < 0.0)
	    S[i][1] = epsilon*S[i][0];
	
	}


      /* 3d - COMPUTE (INTEGRATED) TORQUE ON PLANET */

      /****************************************/
      /****************************************/
      /****************************************/
      /****************************************/
      I=0.0;
      /* /\* i_p = while (i <= R_planet); *\/ */
      /* /\* print(i_p) *\/ */

      if (planet_flag == 1)
	{

      /* FIND PLANET - START AT INNER EDGE */
      i_p = in;
      /* LOOP UNTIL PLANET IS FOUND */
      while (Rf[i_p] <= R_p)
	{
	  i_p++;
	}
      /* ONE STEP BACK */
      i_p--;
      /* printf("%i  %g  %g  %g\n",i_p,Rf[i_p],R_planet,Rf[i_p+1]); */


	
      I = ((Q*Q)/((H[i_p]/Rz[i_p])*(H[i_p]/Rz[i_p])))*sigma[i_p][0]*(Rz[i_p]*Rz[i_p]*Rz[i_p]*Rz[i_p])*(Omega[i_p]*Omega[i_p]);
      
      alpha0 = -(log(sigma[i_p+1][0])-log(sigma[i_p][0]))/(log(Rz[i_p+1])-log(Rz[i_p]));
      beta0 = -(log((H[i_p+1]*H[i_p+1])*1/(Rz[i_p+1]*Rz[i_p+1]*Rz[i_p+1]))-log((H[i_p]*H[i_p])*1/(Rz[i_p]*Rz[i_p]*Rz[i_p])))/(log(Rz[i_p+1])-log(Rz[i_p]));
      
      if (alpha0 > 10.0)
	alpha0=10.0;
      if (alpha0 < -10.0)
	alpha0=-10.0;
      
          if ((S[i_p][0] < 10.0*epsilon) || (S[i_p+1][0] < 10.0*epsilon))
            c0=0;
          else
            c0 = (1/gammaI)*(-2.5-(1.7*beta0)+(0.1*alpha0));

	  /*c0 = (1/gammaI)*(-2.5-(1.7*beta0)+(0.1*alpha0));*/
      /*printf("%g",Il);*/
      /* printf("%g",alpha0);*/
      /*  printf("%g",beta0);*/

	}
	/*Rz[i_p]*Rz[i_p]*Rz[i_p]*Rz[i_p])*/
      /****************************************/
      /****************************************/
      /****************************************/
      /****************************************/
      /*Gonna delete -1.0 as it was used to correct the torque direction*/
      
      /* 3c - MOVE PLANET */
      /* CALCULATE MIGRATION RATE */
      da_dt = sqrt(R_p/(G*M_star)) * ((4.0*(double)M_PI)/M_p) * I*c0;
      
      /* CHANGE RADIUS */
      if (planet_flag == 1)
	R_p = R_p + (da_dt * dt);

      /* STOP PLANET AT INNER BOUNDARY */
      if ((R_p <= R_p_min) && (planet_flag == 1))
	{
	  R_p = 1.0E-15;
	  Q = 0.0;
	  da_dt=0.0;
	  planet_flag = 0;
	  t_p_acc = Time;
	}
      
      
      
      /* LOOP GRID TO SET OTHER VARIABLES */
      for (i=in+1; i<out; i++)
	{
      	  /* MAKE SURFACE DENSITY */
	  sigma[i][1]=S[i][1]/(pow(Rz[i],1.5));
	  
	  /* MAKE Y */
	  Y[i][1]=3.0*nu[i]*sqrt(Rz[i])*sigma[i][1];

	}
      

      
      /* ------------ STEP 4 - BOUNDARY CONDITIONS ------------ */	  
      /* SET INNER ZONE TO ZERO: ZERO-TORQUE BC */
      sigma[in][1]=0.0;
      Y[in][1]=0.0;
      S[in][1]=0.0;
      
      /* OUTER BC */
      /* ZERO-TORQUE - Sigma=0 */
      if (BC_flag==0)
	{
	  sigma[out][1]=0.0;
	  S[out][1]=0.0;
	  Y[out][1]=0.0;
	}
      /* REFLECTIVE - v_r=0 */
      if (BC_flag==1)
	{
	  Y[out][1]=Y[out-1][1];
	}
      
      /* ACCRETION RATE AT ORIGIN */
      /* FOR CONSERVATION TRACKING */
      M_dot[1]=2.0*(double)M_PI*((Y[in+1][0]-Y[in][0])/dX);
      M_acc = M_acc + M_dot[1]*dt;

 
      M_out = M_out + (dt*2.0*(double)M_PI*((Y[out-1][0]-Y[out][0])/dX));

      M_wind = M_wind + (wind_rate*dt);

      /* MAKE "FAKE" Mdot IN HIGH-RES CASE (FOR CONSISTENT OUTPUT) */
      if (res_switch == 1)
	M_dot[2]=2.0*(double)M_PI*((Y[in+factor+1][1]+Y[in+factor+2][1])/(2.0*dX_low));
      

      /* ADD UP DISC MASS */
      disc_mass=0.0;
      for (i=in; i<out; i++)
	{
	  disc_mass = disc_mass + (2.0*(double)M_PI*Rz[i]*(Rf[i+1]-Rf[i])*sigma[i][1]);
	}
      M_d=disc_mass;

      /* ------------ STEP 4 - DATA I/O & BOOK-KEEPING ------------ */	        
      /* OUTPUT DUMP */
      if ((time_to_dump <= 0.0) || (stop_flag==1))
	{
	  /* ECHO TO SCREEN */
	  printf("%i %gyr  %gMsun/yr  %gAU  %gMsun  %g  %gMsun/yr  dt=%gyr   %gM_J  %gAU/yr\n",t,Time,M_dot[1],R_p,disc_mass,(disc_mass+M_acc+M_out+M_wind)/M_d_0,wind_rate,dt,M_p/M_jup,da_dt);

	  /* DUMP TO FILES */ 
	  fprintf(mdot,"%g %g %g %g %g %g\n",Time,M_dot[1],disc_mass,R_p,M_p/M_jup,(disc_mass+M_acc+M_out+M_wind)/M_d_0);


	  /* OUTPUT ARRAYS */
	  /* CREATE FILENAME */
	  sprintf(sigfile,"Files/sigma_%s_%.5i.dat",root,dump_counter);
	  /* OPEN FILE */
	  sigmaout=fopen(sigfile,"w");
	  /* OUTPUT */
	  fprintf(sigmaout,"%g %g\n",Time,R_p);
	  /* DUMP SURFACE DENSITY & TORQUE */
	  for (i=0; i<out+1; i++)
	    {
	      fprintf(sigmaout,"%g  %g  %g\n",Rz[i],sigma[i][1],Lambda[i]);
	    }      
	  /* CLOSE FILE */
	  fclose(sigmaout);
	  /* INCREMENT COUNTER */
	  dump_counter++;

	  /* RESET COUNTER */
	  time_to_dump = step - ( Time - (round(Time/step)*step) );
	  fflush(sigmaout);
	  fflush(mdot);
	  fflush(stdout);
	}


      /* COPY ARRAYS TO 0 FOR ITERATION */
      for (i=0; i<out+1; i++)
	{
	  S[i][0]=S[i][1];
	  Y[i][0]=Y[i][1];
	  sigma[i][0]=sigma[i][1];
	}


      /* INCREMENT COUNTER */
      t++;

      /* CHECK AND ADD NEW PLANET */
      if ( (Time > t_planet) && (planet_flag == 0) && (planet_number==0) )
	{
	  planet_number++;
	  planet_flag=1;
	  M_p = M_planet * M_jup;
	  R_p = R_planet;
	  Q = M_p/M_star;
	  printf("NEW PLANET ADDED\n");
	  
	}
      

    }
printf("Made it to HERE,THANK YOU\n");
  
  
  /* CLOSE DATA FILES */
  fclose(mdot);


  /* WRITE TO PLANET DATA FILE */
  planet_data = fopen("planet_data.dat","a");
  fprintf(planet_data,"%i %g %g %g %g %g %g\n",planet_number,M_planet,M_p/M_jup,R_planet,R_p,t_planet,M_d_0);
  fclose(planet_data);
  
  exit(EXIT_SUCCESS);

}


