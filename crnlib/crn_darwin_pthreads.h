#pragma once

#include "crn_core.h"

#include <semaphore.h>

int sem_timedwait(sem_t* sem, const struct timespec* abs_timeout);
