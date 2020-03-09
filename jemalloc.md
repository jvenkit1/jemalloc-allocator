# Title:

A Scalable Concurrent malloc implementation for FreeBSD -- Jemalloc

# Introduction:

SMP systems have advanced to a state where malloc has become a scalability bottleneck for multi-threaded applications

* ***Allocator Performance***: <br/>

Typically, allocator performance is measured by a combination of the application execution time and average/peak application memory usage.<br/>

Memory layout can have a huge impact on this performance of the application, owing to the latency associated with CPU caches, RAM and virtual memory paging. <br/>

This factor can have a significant impact on **measuring** the performance, since an allocator might perform poorly for certain memory allocation patterns, and if it isn't tested on any such pattern, then the allocator may seem to perform well. <br/>

This demands an allocator design that minimizes the number and severity of degenerate edge cases.<br/>

* ***Fragmentation***: <br/>

Two types of fragmentation exist, Internal and External fragmentations. <br/>

Internal fragmentation is a measure of wasted space associated with individual allocations, due to unusable leading/trailing spaces. <br/>

External fragmentation is a measure of space physically backed by the virtual memory system. <br/>

It is important for a memory allocator to make tradeoffs that impact the amount of occurence of each type of fragmentation. <br/>

* ***Cache Locality and memory usage***: <br/>

An allocator that uses lesser memory does not exhibit better cache locality, since if data pertinent to an application does not fit in the cache, performance will still improve if the data is tightly packed in memory. <br/>

Object allocated together also tend to be used together. Hence, if an allocator allocates data contiguously, there is potential for locality improvement. <br/>

Thus, **jemalloc** tries to minimize memory usage and allocates contiguously, only when it does not conflict with the first goal. <br/>

* ***False Sharing***: <br/>

This phenomenon, as observed in Hoard allocator, occurs when 2 threads, running on separate processors use different data that occurs on the same cache line. In such a case, the processor must arbitrate the ownership of cache lines (cache coherency). <br/>

To fix this, jemalloc relies on multiple allocation arenas to reduce the problem, and allows the programmer to pad allocations to avoid false cache line sharing in performance critical code. <br/>

* ***Lock contention*** <br/>

An important goal for this allocator was to reduce the lock contention for multi threaded applications running on multi processor systems. <br/>

An interesting method on how this can be achieved is to have multiple allocator locks. Each free list can be assigned its own lock. However this approach does not scale well despite minimal lock contention. <br/>

Cache sloshing: quick migration of cached data among processors during the manipulation of allocator data structures. <br/>

Jemalloc uses multiple arenas but combined with a mechanism different to hashing for assignment of threads to these arenas. <br/>

# Algorithm:

Each application is configured at run time to have a fixed number of arenas. Thus threads are assigned to an arena based on a round robin scheme. This technique reduces contention among the different threads while requesting memory. The number of arenas depends on the number of processors. <br/>

The general rule of thumb is to have at least 4 times as many arenas as there are processors in a multi processor environment, while there is only 1 arena in a single processor environment. <br/>
Thread local storage is essential in this round robin implementation, since each thread's arena assignment needs to be stored somewhere.

* Memory Management:

Memory is requested from kernel via sbrk or mmap. All this memory is managed in multiples of chunk size. This means the base address is a multiple of the chunk size, which facilitates constant time calculation of a chunk associated with an allocation.<br/>
Chunks are managed by particular arenas and the chunk size is ***2 Mb*** by default.<br/>

Allocation size classes are of 3 categories - small, large and huge. Each allocation request is rounded up to the nearest size class boundary. Huge allocations are larger than half of a chunk and are directly backed by dedicated chunks(i.e entire chunk for an allocation).<br/>
Metadata about huge allocations are stored in a single red-black tree.<br/>

For small and large allocations, chunks are carved into page runs using the **binary buddy algorithm**. Blocks are repeatedly split into smaller sized blocks to cater to a allocation request, but can only be coalesced in ways that reverse the splitting process.

Small allocations fall into 3 subcategories - tiny, quantum-spaced and sub-page. This particular distribution of the splits reduces average internal fragmentation.
