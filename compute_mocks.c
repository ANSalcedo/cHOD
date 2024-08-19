#include "read_hdf5.h"

int main(int argc, char *argv[])
{
  if(argc != 20) {
    printf("%d usage: ./compute_mocks Omega_m0 redshift siglogM logMmin logM0 logM1 alpha q_cen q_sat del_gamma A_con f_cen alpha_cen alpha_sat boxsize seed [halo catalog file] [galaxy mock file] [halo environment file]\n",argc);
    return -1;
  }

  double Omega_m0 = strtod(argv[1], NULL);
  double redshift = strtod(argv[2], NULL);

  double siglogM = strtod(argv[3], NULL);
  double logMmin = strtod(argv[4], NULL);
  double logM0 = strtod(argv[5], NULL);
  double logM1 = strtod(argv[6], NULL);
  double alpha = strtod(argv[7], NULL);
  double q_cen = strtod(argv[8], NULL);
  double q_sat = strtod(argv[9], NULL);
  double del_gamma = strtod(argv[10], NULL);
  double A_con = strtod(argv[11], NULL);
  double f_cen = strtod(argv[12], NULL);
  double alpha_cen = strtod(argv[13], NULL);
  double alpha_sat = strtod(argv[14], NULL);

  double boxsize = strtod(argv[15], NULL);

  fprintf(stderr,"\tboxsize: %lf\n", boxsize);

  int seed = atoi(argv[16]);

  char *halo_file, *output_file, *env_file;
  size_t halo_file_ssize, output_file_ssize, env_file_ssize;

  halo_file_ssize = sizeof(char)*(strlen(argv[17]) + 1);
  output_file_ssize = sizeof(char)*(strlen(argv[18]) +1 );
  env_file_ssize = sizeof(char)*(strlen(argv[19]) +1 );

  halo_file = malloc(halo_file_ssize);
  output_file = malloc(output_file_ssize);
  env_file = malloc(env_file_ssize);
  snprintf(halo_file, halo_file_ssize, "%s", argv[17]);
  snprintf(output_file, output_file_ssize, "%s", argv[18]);
  snprintf(env_file, env_file_ssize, "%s", argv[19]);

  fprintf(stderr,"Computing HOD from %s\n", halo_file);
  fprintf(stderr,"Reading environment density from %s\n", env_file);
  fprintf(stderr,"Saving to output file %s\n", output_file);

  populate_hod(siglogM, logMmin, logM0, logM1, alpha, q_cen, q_sat, del_gamma, A_con, f_cen, alpha_cen, alpha_sat, \
	       seed, Omega_m0, redshift, boxsize, halo_file, output_file, env_file);

  return 0;
}
