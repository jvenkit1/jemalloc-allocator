#include "localimports.h"

#ifndef THREAD
#define THREAD
namespace jemalloc {

class Thread{
  int tid;  // thread_id for the current thread
};

}
#endif
