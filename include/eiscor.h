#ifndef STDERR 
#define STDERR (0)
#endif 

#ifdef DEBUG 
#undef DEBUG
#define DEBUG (.TRUE.)
#else 
#define DEBUG (.FALSE.)
#endif   

#ifdef VERBOSE 
#undef VERBOSE
#define VERBOSE (.TRUE.)
#else 
#define VERBOSE (.FALSE.)
#endif 

#ifdef EISCOR_DBL_EPS 
#undef EISCOR_DBL_EPS
#define EISCOR_DBL_EPS (epsilon(1d0))
#else 
#define EISCOR_DBL_EPS (epsilon(1d0))
#endif

#ifdef EISCOR_DBL_INF 
#undef EISCOR_DBL_INF
#define EISCOR_DBL_INF (huge(1d0))
#else 
#define EISCOR_DBL_INF (huge(1d0))
#endif

#ifdef EISCOR_DBL_PI 
#undef EISCOR_DBL_PI
#define EISCOR_DBL_PI (3.141592653589793239d0)
#else 
#define EISCOR_DBL_PI (3.141592653589793239d0)
#endif

