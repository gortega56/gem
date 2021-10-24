#pragma once

#if defined(_MSC_VER)
#define GEM_INLINE __forceinline
#else
#define GEM_INLINE inline
#endif

#ifdef GEM_VECTORCALL_ENABLE
#if defined(_MSC_VER)
#define GEM_VECTORCALL __vectorcall
#else
#define GEM_VECTORCALL
#endif
#else
#define GEM_VECTORCALL
#endif