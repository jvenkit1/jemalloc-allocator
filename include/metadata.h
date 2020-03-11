#include<vector>
#include "thread.h"
#include "chunk.h"
using namespace std;

#define MAXREGIONS 200

class Thread{
  int tid;  // defines the thread id
};

class region{
  int allocationSize;
  // add enum to include the current type
};

class run{
  bin bins;
  region curr_region;
  vector<page> pages;
  int numFreeRegions;
  int totalFreeRegions;
  bool isRegionFree[MAXREGIONS];  // randomly assigning 200 pages for each run
};

class chunk{
  arena curr_arena; // arena for the current chunk;
  bool isDirty;  // indicates if the current chunk is dirty.
  vector<run> runs;

  map<region, page> freePages;  // contains a mapping for all free pages for a particular region
};


class bin{
  region regions;
  vector<run> runs;
  run curr_run;
};

class arena{
  vector<thread> threads;
  vector<chunk> chunks;

  set<run> cleanRunsAvailable;
  set<run> dirtyRunsAvailable;

  vector<bin> bins;
};
