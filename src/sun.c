/* sun.c
 * 2003/02/15 Kenlo Nishida
 * 2017/06/21 Kenlo Nishida Nasahara
 * compile: $ gcc sun.c -o sun -lm
 */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# define NV 4                // number of (dependent) variables for Runge-Kutta solver
# define M 1.989e30          // mass of sun (in kg)
# define G 6.67259e-11       // gravity constant in SI unit
# define START_DOY 185       // July 4th (earth passes the far point (2003))
# define PI 3.14159265
# define AXIS_INCL 0.408407  // inclination of axis of earth (in radian) = 23.4 degree

void func_gravity(double t, double *x, double *f)
	// t: time
	// (x[0], x[1]): location of earth; (x[2], x[3]): velocity of earth
	// f is output
{double r;
 r=sqrt(x[0]*x[0] + x[1]*x[1]);   // distance between sun & earth
 f[0]=x[2];                       // dx/dt
 f[1]=x[3];                       // dy/dt
 f[2]=x[0]*(-1)*G*M/(r*r*r);      // du/dt  (Newton equation)
 f[3]=x[1]*(-1)*G*M/(r*r*r);      // dv/dt  (Newton equation)
}


// calculatitudee location of earth with Runge-Kutta
void revolution(double doy, double *distance, double *phase)
{double x[NV];                           // (x[0], x[1]): earth location; (x[2], x[3]): earth velocity
 double x1[NV], x2[NV], x3[NV];          // intermediate solutionis in Runge-Kutta solver
 double f[NV];                           // d/dt x[]
 double f0[NV], f1[NV], f2[NV], f3[NV];  // intermediate solutionis in Runge-Kutta solver
 double t, t_start, t_end;               // time
 double dt=0;                            // time step
 int i;

// setting time frame
 t_start=(double)START_DOY*24*60*60;
 t_end=doy*24*60*60;
 if (t_end<t_start) t_end=t_end+365.25*24*60*60;
 dt=10*60;

// setting initial condition
 x[0]=1.5210e11;                              // x ... the longest distance
 x[1]=0;                                      // y
 x[2]=0;                                      // u
 x[3]=1.47122e11*2*3.1415/(365.25*24*60*60);  // v ... manually adjusted

// Runge-Kutta solver
 for (t=t_start; t<=t_end; t=t+dt)
{
 func_gravity(t, x,  f0);
 for (i=0; i<NV; i++) x1[i]= x[i]+f0[i]*dt/2;
 func_gravity(t+dt/2, x1, f1);
 for (i=0; i<NV; i++) x2[i]=x[i]+f1[i]*dt/2;
 func_gravity(t+dt/2, x2, f2);
 for (i=0; i<NV; i++) x3[i]=x[i]+f2[i]*dt;
 func_gravity(t+dt,   x3, f3);

 for (i=0; i<NV; i++) f[i]=(f0[i] + 2*f1[i] + 2*f2[i] + f3[i])/6;
 for (i=0; i<NV; i++) x[i]=x[i]+f[i]*dt;
}
*distance=sqrt(x[0]*x[0] + x[1]*x[1]);
*phase=acos(x[0] / *distance);
if (x[1]<0) *phase = 2 * PI - *phase;
}


void sun_position(double lat, double doy, double *answer)
{double doy_;              // day of year with respect to the base DOY
 double doy_base;
 double soy;               // second of year (doy*24*60*60)
 double distance;          // between sun & earth
 double phase_revol;
 double phase_revol_base;  // standard of revolution angle (the longest point)
 double phase_rotate;
 double latitude;
 double px0, py0, pz0;     // ground location (static geographic system)
 double px1, py1, pz1;     // ground location (after diurnal rotation)
 double px2, py2, pz2;     // ground location (after axis-inclination)
 double ax,  ay,  az;      // rotation axis 
 double sx,  sy,  sz;      // earth-center to sun
 double ex,  ey,  ez;      // direction of east on the ground
 double nx,  ny,  nz;      // direction of north on the ground
 double se,  sn;           // projection of s on ground east-north coordinate
 double sol_zenith, sol_azimuth;
 double r;

 latitude=lat*PI/180;
 doy_base=-9; 
 revolution(doy_base, &distance, &phase_revol_base);  // winter solstice
 revolution(doy, &distance, &phase_revol);
 doy_=doy-doy_base;
 soy=doy_*24*60*60;
 phase_revol = (phase_revol-phase_revol_base + PI);  // direction of sun from earth
 phase_rotate = (2 * PI * 366.25 * doy_ / (365.25));

 while (phase_revol>2*PI)  phase_revol=phase_revol-2*PI;
 while (phase_revol<0)     phase_revol=phase_revol+2*PI;
 while (phase_rotate>2*PI) phase_rotate=phase_rotate-2*PI;
 while (phase_rotate<0)    phase_rotate=phase_rotate+2*PI;
  
// initial condition (winter solstice)
 px0=cos(latitude);
 py0=0;
 pz0=sin(latitude);

// rotate aroud the axis
 px1 = px0 * cos(phase_rotate) - py0 * sin(phase_rotate);
 py1 = px0 * sin(phase_rotate) + py0 * cos(phase_rotate);
 pz1 = pz0;

// inclination of axis
 px2 = px1 * cos(AXIS_INCL) +    pz1 * sin(AXIS_INCL);
 py2 = py1;
 pz2 = px1 * sin(-1*AXIS_INCL) + pz1 * cos(AXIS_INCL);

// unit vector of rotation axis
 ax = sin(AXIS_INCL);
 ay = 0;
 az = cos(AXIS_INCL);

// direction to the sun
 sx = cos(phase_revol);
 sy = sin(phase_revol); 
 sz = 0;

// direction to east
 ex = ay*pz2 - az*py2;
 ey = az*px2 - ax*pz2;
 ez = ax*py2 - ay*px2;
 r = sqrt(ex*ex + ey*ey + ez*ez);
 ex=ex/r; ey=ey/r; ez=ez/r;

// direction to north
 nx = py2*ez - pz2*ey;
 ny = pz2*ex - px2*ez;
 nz = px2*ey - py2*ex;
 
// solar zenith
 sol_zenith=acos(sx*px2 + sy*py2 + sz*pz2);

// solar azimuth
 se = sx*ex + sy*ey + sz*ez;
 sn = sx*nx + sy*ny + sz*nz;
 r = sqrt(se*se + sn*sn);
 se = se/r; sn = sn/r; 
 sol_azimuth=acos(se);
 if (sn>0) sol_azimuth=2*PI-sol_azimuth;   // 0 to PI in southward, PI to 2PI in northward

 answer[0]=180*sol_azimuth/PI;
 answer[1]=180*sol_zenith/PI;
 answer[2]=distance;
 answer[3]=180*phase_revol/PI;
}


int main(int argc, char **argv)
{double doy; 
 double lat, lon;           // in degree
 double answer[4];

 if (argc<2) {
     printf("usage: $ sun latitude longitude DOY\n");
     printf("note: If DOY is in UTC, you should set longitude as the real value.\n");
     printf("note: If DOY is in local time, you should set longitude as 0.\n");
     printf("output: solar_azimuth solar_zenith. Both in degree, not radian.\n"); 
     printf("note: solar_azimuth is 0 at east, 90 deg at south, 180 deg ast west.\n"); 
     exit(0);
 } 
 lat=atof(argv[1]);
 lon=atof(argv[2]);
 doy=atof(argv[3]);

 doy=doy+(lon/15.0)/24.0;
 sun_position(lat, doy, answer);
 printf("%lf %lf\n", answer[0], answer[1]);
}
