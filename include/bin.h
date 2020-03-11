#include "localimports.h"
#include "run.h"
#include "region.h"
#ifndef BIN
#define BIN
namespace jemalloc{
class bin{
  region regions;
  vector<run> runs;
  run curr_run;
};
}
#endif
