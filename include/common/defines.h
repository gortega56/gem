#pragma once
#if defined(_MSC_VER)
#define GEM_INLINE __forceinline
#define GEM_VECTORCALL __vectorcall
#else
#define GEM_INLINE inline
#define GEM_VECTORCALL
#endif