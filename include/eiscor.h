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


