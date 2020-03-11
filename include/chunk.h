#include "localimports.h"
#include "run.h"
#include "region.h"
#include "page.h"

#ifndef CHUNK
#define CHUNK
namespace jemalloc{

class chunk{
  arena curr_arena; // arena for the current chunk;
  bool isDirty;  // indicates if the current chunk is dirty.
  vector<run> runs;

  map<region, page> freePages;  // contains a mapping for all free pages for a particular region
};

}
#endif
