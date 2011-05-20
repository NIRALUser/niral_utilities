#ifndef _MINMAX_H 
#define _MINMAX_H 
//#ifdef WIN32 
// needed to cope with bug in MS library: 
// it fails to define min/max 
namespace std{
template <class T> inline T Max(const T& a, const T& b) 
{ 
    return (a > b) ? a : b;
} 

template <class T> inline T Min(const T& a, const T& b) 
{ 
    return (a < b) ? a : b;
} 
}

//#endif 

#endif 
