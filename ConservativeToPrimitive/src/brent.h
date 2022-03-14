#include <math.h>
#include <stdbool.h>

#define par_q 0
#define par_r 1
#define par_s 2
#define par_t 3
#define conD  4
#define conYE  5




struct brent_report {
  bool failed;
  char err_msg[200];
  int count;
};

