#include "read_hdf5.h"

int main(int argc, char *argv[])
{
  if(argc != 17) {
    printf("%d usage: ./compute_mocks Omega_m0 siglogM logMmin logM0 logM1 alpha q_env del_gamma f_cen alpha_cen alpha_sat boxsize seed [halo catalog file] [galaxy mock file] [halo environment file]\n",argc);
    return -1;
  }

  double Omega_m0 = strtod(argv[1], NULL);

  double siglogM = strtod(argv[2], NULL);
  double logMmin = strtod(argv[3], NULL);
  double logM0 = strtod(argv[4], NULL);
  double logM1 = strtod(argv[5], NULL);
  double alpha = strtod(argv[6], NULL);
  double q_env = strtod(argv[7], NULL);
  double del_gamma = strtod(argv[8], NULL);
  double f_cen = strtod(argv[9], NULL);
  double alpha_cen = strtod(argv[10], NULL);
  double alpha_sat = strtod(argv[11], NULL);

  double boxsize = strtod(argv[12], NULL);

  fprintf(stderr,"\tboxsize: %lf\n", boxsize);

  int seed = atoi(argv[13]);

  char *halo_file, *output_file, *env_file;
  size_t halo_file_ssize, output_file_ssize, env_file_ssize;

  halo_file_ssize = sizeof(char)*(strlen(argv[14]) + 1);
  output_file_ssize = sizeof(char)*(strlen(argv[15]) +1 );
  env_file_ssize = sizeof(char)*(strlen(argv[16]) +1 );

  halo_file = malloc(halo_file_ssize);
  output_file = malloc(output_file_ssize);
  env_file = malloc(env_file_ssize);
  snprintf(halo_file, halo_file_ssize, "%s", argv[14]);
  snprintf(output_file, output_file_ssize, "%s", argv[15]);
  snprintf(env_file, env_file_ssize, "%s", argv[16]);

  fprintf(stderr,"Computing HOD from %s\n", halo_file);
  fprintf(stderr,"Reading environment density from %s\n", env_file);
  fprintf(stderr,"Saving to output file %s\n", output_file);

  populate_hod(siglogM, logMmin, logM0, logM1, alpha, q_env, del_gamma, f_cen, alpha_cen, alpha_sat, \
	       seed, Omega_m0, boxsize, halo_file, output_file, env_file);

  return 0;
}
