#include "strassen.h"

boost::mutex counter_mutex;

size_t Strassen :: _maxMemoryUse = 0;
size_t Strassen :: _nowMemoryUse = 0;