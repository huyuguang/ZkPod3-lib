#ifndef __MSVC_HACK_H__
#define __MSVC_HACK_H__

#include <stdlib.h>

#if defined(_MSC_VER)
#include <intrin.h>

#if defined(_WIN64)
typedef __int64 ssize_t;
#define builtin_clz(x) __lzcnt(x)
#define builtin_clzl(x) __lzcnt64(x)
#else
typedef _W64 __int32 ssize_t;
#define buildin_clz(x) __lzcnt(x)
#define buildin_clzl(x) __lzcnt(x)
#endif  // _WIN64

#define __noinline__ __declspec(noinline)

inline int setenv(const char *name, const char *value, int overwrite) {
  int errcode = 0;
  if (!overwrite) {
    size_t envsize = 0;
    errcode = getenv_s(&envsize, NULL, 0, name);
    if (errcode || envsize) return errcode;
  }
  return _putenv_s(name, value);
}

#endif  // _MSC_VER

#ifdef __GNUC__
#define builtin_clz(x) __builtin_clz(x)
#define builtin_clzl(x) __builtin_clzl(x)
#define __noinline__ __attribute__((noinline))
#endif  // __GNUC__

#ifdef _DEBUG
#define DEBUG
#endif

#endif  // __MSVC_HACK_H__