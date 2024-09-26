#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stddef.h>
#include <stdio.h> // IWYU pragma: export

typedef double Scalar;
typedef ptrdiff_t Index;

#define LOG_INFO(...) fprintf(stderr, __VA_ARGS__)

#endif // __CONFIG_H__
