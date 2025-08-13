#include "read_hdf5.h"

#define M_PI 3.14159265358979323846
#define Mpc_to_cm 3.0856e24 /*Conversion factor from Mpc to cm */
#define Msun_to_g 1.989e33 /*Conversion factor from Msun to grams*/
#define G 6.672e-8 /*Universal Gravitational Constant in cgs units*/
#define Hubble 3.2407789e-18 /*Hubble's constant h/sec*/
#define rho_crit (3.0*pow(Hubble, 2.0) / (8.0 * M_PI * G)) * (pow(Mpc_to_cm, 3.0) / Msun_to_g) /*Cosmological Critical Density in Msun h^2 / Mpc^3 */
#define INDEX10(i,j) (i*10 + j)

/* These functions populate halos with HOD galaxies from an HDF5 halo catalog */
/* The outputs specify galaxy positions, and velocities (host mass and satellite/central identification have been dropped for speed and space considerations) */

hostDMH * find_galaxy_hosts(struct halo halos[], halo_metadata * env, double siglogM, double logMmin, double q_cen, double alpha_cen, double redshift, unsigned long int N_h, unsigned long int *Ncen, gsl_rng *r)
{
  /*This function uses the Zehavi 2011 prescription to find the halos that host central galaxies*/

  float * hosts = malloc(N_h * 10 * sizeof(float));
  int * IDs = malloc(N_h * sizeof(int));
  int i;
  unsigned long j = 0;

  float f_logMmin_0 = (float)logMmin;
  float f_siglogM = (float)siglogM;
  float f_q_cen = (float)q_cen;
	
  /*#pragma omp parallel for*/
  for(i = 0; i < N_h; i++)
    {
      float logM = (float)log10(halos[i].mass);
      float env_rank = env[i].percentile;
      double vrms = halos[i].vrms;
	      
      /* compute modified logMmin for this halo's environment */
      float f_logMmin = f_logMmin_0 + f_q_cen * (env_rank - 0.5f);

      /*Mean central occupation or the probability of hosting a central*/
      float prob = 0.5 * (1.0 + erf( (logM - f_logMmin) / f_siglogM) );

      if(prob > gsl_rng_uniform(r))
	{
	  hosts[INDEX10(j,0)] = halos[i].X;
	  hosts[INDEX10(j,1)] = halos[i].Y;
	  hosts[INDEX10(j,2)] = halos[i].Z;
	  hosts[INDEX10(j,3)] = halos[i].vx + gsl_ran_laplace(r, alpha_cen*vrms*pow(2, -1.0/2.0)); 
	  hosts[INDEX10(j,4)] = halos[i].vy + gsl_ran_laplace(r, alpha_cen*vrms*pow(2, -1.0/2.0));
	  hosts[INDEX10(j,5)] = halos[i].vz + gsl_ran_laplace(r, alpha_cen*vrms*pow(2, -1.0/2.0));
	  hosts[INDEX10(j,6)] = halos[i].mass;
	  hosts[INDEX10(j,7)] = vrms;
	  hosts[INDEX10(i,8)] = redshift;
	  hosts[INDEX10(i,9)] = env[i].percentile; 
	  IDs[j] = halos[i].HID;
	  j ++;
	}
    }
  
  *Ncen = j;

  hostDMH *host_coords = malloc(j*sizeof(hostDMH));

  for(i=0;i<j;i++)
    {
      host_coords[i].X              = hosts[INDEX10(i,0)];
      host_coords[i].Y              = hosts[INDEX10(i,1)];
      host_coords[i].Z              = hosts[INDEX10(i,2)];
      host_coords[i].vx             = hosts[INDEX10(i,3)];
      host_coords[i].vy             = hosts[INDEX10(i,4)];
      host_coords[i].vz             = hosts[INDEX10(i,5)];
      host_coords[i].mass           = hosts[INDEX10(i,6)];
      host_coords[i].vrms           = hosts[INDEX10(i,7)];
      host_coords[i].redshift       = hosts[INDEX10(i,8)];
      host_coords[i].env_percentile = hosts[INDEX10(i,9)];
      host_coords[i].HID = IDs[i];
    }

  return host_coords; 
}

int * find_satellites(struct halo halos[], double siglogM, double logMmin, double logM0, double logM1, double alpha, double q_sat, unsigned long int N_h, unsigned long int *Nsat, gsl_rng *r)
{
  /*This function determines how many satellite galaxies each halo has using the same Zehavi 2011 prescription*/

  int i;
  unsigned long j =0;
  int * satellites = malloc(N_h*sizeof(int));

  /*#pragma omp parallel for*/
  for(i=0; i < N_h; i++) 
    {
      double env_rank = (double)halos[i].env_percentile;	  
      double M0 = pow(10.0, logM0);
      double logM1_init = logM1;
      float f_q_sat = (float)q_sat;
      double f_logM1 = logM1_init + f_q_sat * (env_rank - 0.5f);
      double M1 = pow(10.0, f_logM1);
      
      float logM = log10(halos[i].mass); //Come back to this once hdf5 is figured out
      /* double mean_cen = 0.5 * (1.0 + erf( (logM - logMmin) / siglogM) ); */
      double mean_cen = 1.0;
      double mean_sat;
      if(logM < logM0)
	{
	  mean_sat = 0.0; /* Enforcing the satellite cutoff */
	}
      else
	{
	  mean_sat = mean_cen * pow( ( ( halos[i].mass - M0 ) / M1 ), alpha );
	}
      unsigned int N_sat_this_halo = gsl_ran_poisson(r, mean_sat); /* Drawing from Poisson distribution centered at mean_sat to determine Nsat */
      satellites[i] =  N_sat_this_halo;
      j = j + N_sat_this_halo;
    }
  
  *Nsat = j;

  return satellites;
}

galaxy * pick_NFW_satellites(struct halo host, const int N_sat, double alpha_sat, double O_m0, double del_gamma, double A_con, gsl_rng *r)
{
  /* This function determines the spatial distribution of the satellite galaxies */
  /* Galaxies are NFW-distributed using results from Correa et al. 2015 */
  /* We allow for a deviation from NFW in the galaxy distribution by a factor of r^delta_gamma */

  galaxy * coords = malloc(N_sat * sizeof(galaxy));

  double z = host.redshift;
  double alpha = 1.62774 - 0.2458*(1.0 + z) + 0.01716*pow(1.0 + z, 2.0);
  double beta = 1.66079 + 0.00359*(1.0 + z) - 1.6901*pow(1.0 + z, 0.00417);
  double gamma = -0.02049 + 0.0253*pow(1.0 + z ,-0.1044);
  
  double D_vir = 200.0, rho_u = O_m0 * rho_crit * pow(1+z, 3.0);

  double logM = log10(host.mass), x0 = host.X, y0 = host.Y, z0 = host.Z, vx0 = host.vx, vy0 = host.vy, vz0 = host.vz, vrms = host.vrms;
  int HID = host.HID;
  double exponent = alpha + beta*logM*(1.0 + gamma*pow(logM, 2.0)); /* Fit from Correa et al. 2015 */
  double cvir = A_con * sqrt(2.0)*pow(10.0, exponent); /* Approximate factor to rescale Rvir between crit, matter */
  double R_vir = pow((3.0/4.0)*(1.0/M_PI)*(host.mass/(D_vir*rho_u)), 1.0/3.0);
  
  int j;

  /* pre-compute NFW profile for this satellite */
  float CDF[1000];
  size_t i;
  double prefac = (1.0 + cvir) / ( 2.0 + del_gamma - (1.0 + del_gamma)*gsl_sf_hyperg_2F1(1.0, 1.0, 3.0 + del_gamma, cvir/(1.0 + cvir)) );
  /* double prefac = 1.0 / ( log( 1.0 + cvir ) - (cvir / ( 1.0 + cvir )) ); */ /* Prefactor 1/A(c_vir) */
  float f_c_vir = (float)cvir;
	
 #pragma omp simd
 for(i=0; i<1000; i++)
    {
      float x = (float)i / 1000.0;
      /* CDF[i] = prefac * ( log( 1.0 + x * f_c_vir ) - (x * f_c_vir / ( 1.0 + x*f_c_vir )) ); */ 
      CDF[i] = prefac * pow( x, (2.0 + del_gamma) ) * ( 2.0 + del_gamma - (1.0 + del_gamma)*gsl_sf_hyperg_2F1(1.0, 1.0, 3.0 + del_gamma, x*f_c_vir/(x*f_c_vir + 1.0)) ) / (1 + f_c_vir*x);
    }
  
  for(j=0; j<N_sat; j++)
    {
      double frac = NFW_CDF_sampler(&CDF[0], r);
      double R = R_vir * frac;
      double phi = 2.0*M_PI*gsl_rng_uniform(r), costheta = 2.0*gsl_rng_uniform(r) - 1.0; /* Sphere point picking */
      double sintheta = sqrt(1.0 - costheta*costheta);
      double x = R*sintheta*cos(phi)+x0 , y = R*sintheta*sin(phi)+y0 , z = R*costheta+z0; /* Satellite Galaxy Coordinates */
      coords[j].X = x;
      coords[j].Y = y;
      coords[j].Z = z;
      coords[j].vx = vx0 + alpha_sat*vrms*gsl_ran_gaussian(r, 1.0);
      coords[j].vy = vy0 + alpha_sat*vrms*gsl_ran_gaussian(r, 1.0);
      coords[j].vz = vz0 + alpha_sat*vrms*gsl_ran_gaussian(r, 1.0);
      coords[j].cen_flag = 1; /*cen/sat flag, 1 for satellite*/
      coords[j].host_id = HID;
      coords[j].env_percentile = host.env_percentile;
      coords[j].host_mass = host.mass;
    }

  return coords;
}

inline double wrap_periodic(double x, double Lbox)
{
  if((x < Lbox) && (x >= 0.)) {
    return x;
  } else if (x >= Lbox) {
    return (x-Lbox);
  } else if (x < 0.) {
    return (x+Lbox);
  }
}

void populate_hod(double siglogM, double logMmin, double logM0, double logM1, double alpha, double q_cen, double q_sat, double A_con, double del_gamma, double f_cen, double alpha_cen, double alpha_sat_l, double alpha_sat_h, double M_pivot_as, unsigned long int seed, double Omega_m0, double redshift, double Lbox, char *input_fname, char *output_fname, char *env_fname)
{
  herr_t status;
  size_t NumData,i;
  hostDMH *data;
  
  unsigned long Ncen;
  unsigned long Nsat;
  unsigned long fNcen = 0;

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed); /* Seeding random distribution */
	
  data = read_halo_hdf5(input_fname,"halos",&NumData);

  halo_metadata *env;
  size_t NumEnv;
  env = read_env_hdf5(env_fname,"halos",&NumEnv);
  	
  /* check that NumData == NumEnv */
  if(!(NumEnv == NumData)) exit(1);

  /* compute HOD parameters from number density, mass function, environment density */  
	
  hostDMH *cenhalos; //Central Coordinates
  cenhalos = find_galaxy_hosts(data, env, siglogM, logMmin, q_cen, alpha_cen, redshift, NumData, &Ncen, r);
  galaxy * cens = malloc(Ncen*sizeof(galaxy));
	
  for(i=0; i<Ncen; i++){
    if(gsl_rng_uniform(r) < f_cen)
    {
      cens[fNcen].X = cenhalos[i].X;
      cens[fNcen].Y = cenhalos[i].Y;
      cens[fNcen].Z = cenhalos[i].Z;
      cens[fNcen].vx = cenhalos[i].vx;
      cens[fNcen].vy = cenhalos[i].vy;
      cens[fNcen].vz = cenhalos[i].vz;
      cens[fNcen].cen_flag = 0; /*zero identifies as central */
      cens[fNcen].host_id = cenhalos[i].HID;
      cens[fNcen].host_mass = cenhalos[i].mass;
      cens[fNcen].env_percentile = cenhalos[i].env_percentile;
      fNcen ++;
    }
  }

  int *sats; //Number of Satellites for each halo
  sats = find_satellites(cenhalos, siglogM, logMmin, logM0, logM1, alpha, q_sat, Ncen, &Nsat, r);
  galaxy * coords  = malloc(Nsat*sizeof(galaxy)); //Satellite Coordinates
  int j,k,l=0;

  /*#pragma omp parallel for*/
  for(j=0;j<Ncen;j++)
  {
    if(sats[j]>0){
      galaxy * halosats = malloc(sats[j] * sizeof(galaxy));
	  if(cenhalos[j].mass > M_pivot_as){
      halosats = pick_NFW_satellites(cenhalos[j], sats[j], alpha_sat_h, Omega_m0, del_gamma, A_con, r);
	  }
	  else {
      halosats = pick_NFW_satellites(cenhalos[j], sats[j], alpha_sat_l, Omega_m0, del_gamma, A_con, r);	  
	  }
      for(k=0; k<sats[j]; k++)
	{
	  coords[l].X = wrap_periodic(halosats[k].X,Lbox);
	  coords[l].Y = wrap_periodic(halosats[k].Y,Lbox);
	  coords[l].Z = wrap_periodic(halosats[k].Z,Lbox);
	  coords[l].vx = halosats[k].vx;
	  coords[l].vy = halosats[k].vy;
	  coords[l].vz = halosats[k].vz;
	  coords[l].cen_flag = 1; /*one identifies as sat*/
 	  coords[l].host_id = halosats[k].host_id;
	  coords[l].host_mass = halosats[k].host_mass;
	  coords[l].env_percentile = halosats[k].env_percentile;
	  l++;
	}
      free(halosats);
    }
  }
  
  free(cenhalos);
  free(sats);

  unsigned long int len = Nsat + fNcen;

  galaxy *HODgals = malloc(len*sizeof(galaxy));

  for(i=0; i<fNcen;i++)
    {
      HODgals[i].X = cens[i].X;
      HODgals[i].Y = cens[i].Y;
      HODgals[i].Z = cens[i].Z;
      HODgals[i].vx = cens[i].vx;	 
      HODgals[i].vy = cens[i].vy;	 
      HODgals[i].vz = cens[i].vz;
      HODgals[i].cen_flag = (int) cens[i].cen_flag;
      HODgals[i].host_id = (int) cens[i].host_id;	
      HODgals[i].env_percentile = cens[i].env_percentile;
      HODgals[i].host_mass = cens[i].host_mass;
    }
  free(cens);
  
  for(i=0; i<Nsat; i++)
  {
    HODgals[i+fNcen].X = coords[i].X;
    HODgals[i+fNcen].Y = coords[i].Y;
    HODgals[i+fNcen].Z = coords[i].Z;
    HODgals[i+fNcen].vx = coords[i].vx;	 
    HODgals[i+fNcen].vy = coords[i].vy;	 
    HODgals[i+fNcen].vz = coords[i].vz;
    HODgals[i+fNcen].cen_flag = (int) coords[i].cen_flag;
    HODgals[i+fNcen].host_id = (int) coords[i].host_id;
    HODgals[i+fNcen].env_percentile = coords[i].env_percentile;
    HODgals[i+fNcen].host_mass = coords[i].host_mass;	  
  }
  free(coords);

  printf("Satellites Found. Writing to HDF5 file: %s\n", output_fname);
  status = write_gal_hdf5(output_fname, "particles", (size_t)len, HODgals);
  
  free(HODgals);
 
  gsl_rng_free(r);
}
