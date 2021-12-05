#include "bbox.h"
#include "bvh.h"

#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define MAXTHREADS 128

static int _argc;
static const char **_argv;


static void show_help(const char *program_path) {
  printf("Usage: %s OPTIONS\n", program_path);
  printf("\n");
  printf("OPTIONS:\n");
  printf("\t-f <input_filename> (required)\n");
  printf("\t-n <num_of_threads> (required)\n");
}

const char *get_option_string(const char *option_name,
                              const char *default_value) {
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return _argv[i + 1];
  return default_value;
}

int get_option_int(const char *option_name, int default_value) {
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return atoi(_argv[i + 1]);
  return default_value;
}

float get_option_float(const char *option_name, float default_value) {
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return (float)atof(_argv[i + 1]);
  return default_value;
}

int main(int argc, const char *argv[]) {

    _argc = argc - 1;
    _argv = argv + 1;

    const char *input_filename = get_option_string("-f", NULL);
    int num_of_threads = get_option_int("-n", 1);

    int error = 0;

    if (input_filename == NULL) {
        printf("Error: You need to specify -f.\n");
        error = 1;
    }

    if (error) {
        show_help(argv[0]);
        return 1;
    }


    printf("Number of threads: %d\n", num_of_threads);
    printf("Input file: %s\n", input_filename);

    omp_set_nested(1);
    omp_set_num_threads(num_of_threads);

    //Setup finished, actual code goes here!
}

