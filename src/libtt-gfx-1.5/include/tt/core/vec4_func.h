#ifndef vec4_func_h
#define vec4_func_h

template<typename T1, typename T2> inline
void vec4_copy(T1 to[4], const T2 from[4])
{
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
	to[3] = from[3];
}

template<typename T1, typename T2> inline
void vec4_set(T1 v[3], T2 x, T2 y, T2 z, T2 w)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
	v[3] = w;
}

#endif  // vec4_func_h

