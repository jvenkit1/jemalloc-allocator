#include "localimports.h"
#include "thread.h"
#include "chunk.h"
#include "run.h"
#include "bin.h"

#ifndef ARENA
#define ARENA
namespace jemalloc{

class Arena{
  vector<thread> threads;
  vector<chunk> chunks;

  set<run> cleanRunsAvailable;
  set<run> dirtyRunsAvailable;

  vector<bin> bins;
};

}

#endif
