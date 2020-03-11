#include "localimports.h"
#include "bins.h"
#include "regions.h"
#include "page.h"

#ifndef RUN
#define RUN
namespace jemalloc {
class run{
  bin bins;
  region curr_region;
  vector<page> pages;
  int numFreeRegions;
  int totalFreeRegions;
  bool isRegionFree[MAX_REGIONS];  // randomly assigning 200 pages for each run
};
}
#endif
