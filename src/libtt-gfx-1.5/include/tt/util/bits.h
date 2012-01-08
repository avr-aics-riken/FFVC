#ifndef bits_h
#define bits_h

#include <string>

template<typename T>
inline void tt_bits_on(T& bits, int i)
{
	bits |= (1 << i);
}

template<typename T>
inline void tt_bits_off(T& bits, int i)
{
	bits &= ~(1 << i);
}

template<typename T>
inline void tt_bits_toggle(T& bits, int i)
{
	T mask = (1 << i);
	bits = (bits ^ mask) & mask | bits & ~mask;
}

template<typename T>
inline void tt_bits_set(T& bits, int i, int v)
{
	(v == 0) ? tt_bits_off(bits, i) : tt_bits_on(bits, i);
}

template<typename T>
inline int tt_bits_get(T bits, int i)
{
	return (bits & (1 << i)) >> i;
}

template<typename T>
inline void tt_bits_setn(T& bits, int i, int nbits, T v)
{
	T mask = 0;
	for (int j=0; j<nbits; j++) mask |= 0x1 << j;
	mask <<= i;
	bits = bits & ~mask | v << i;
}

template<typename T>
inline T tt_bits_getn(T bits, int i, int nbits)
{
	T mask = 0;
	for (int j=0; j<nbits; j++) mask |= 0x1 << j;
	mask <<= i;
	return (bits & mask) >> i;
}

template<typename T>
inline std::string tt_bits_string(T bits)
{
	int nbits = sizeof(T) * 8;
	std::string s;
	for (int i=nbits-1; i>=0; --i)
	{
		s += (bits >> i & 1) ? '1' : '0';
	}
	return s;
}

#endif  // tt_bits_h

