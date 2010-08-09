#ifndef __LOG_H_

#define __LOG_H_

#include "globals.h"

#include <cstdio>

void SetLogFile(std::FILE *newLogFile);
void LogMessage(const char *fmt, ...);
void LogDebug(const char *fmt, ...);
void LogError(const char *fmt, ...);

#endif