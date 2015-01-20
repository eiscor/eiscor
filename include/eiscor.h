#ifdef DEBUG 
#undef DEBUG
#define DEBUG (.TRUE.)
#else 
#define DEBUG (.FALSE.)
#endif   

#ifdef STDERR 
#undef STDERR
#define STDERR (0)
#else 
#define STDERR (0)
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

