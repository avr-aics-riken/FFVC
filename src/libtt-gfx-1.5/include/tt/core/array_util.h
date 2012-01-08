#ifndef array_util_h
#define array_util_h

#include <vector>
#include <list>

template<class T> inline
void copy_data(TtArray<T>& ar, const std::vector<T>& vec)
{
	int sz = vec.size();
	ar.resize(sz);
	for (int i=0; i<sz; i++)
	{
		ar[i] = vec[i];
	}
}

template<class T> inline
void copy_data(TtArray<T>& ar, const std::list<T>& lst)
{
	typename std::list<T>::const_iterator itr = lst.begin();
	int sz = lst.size();
	ar.resize(sz);
	for (int i=0; i<sz; i++, itr++)
	{
		ar[i] = *itr;
	}
}

template<class T> inline
void copy_data(std::vector<T>& vec, const std::list<T>& lst)
{
	typename std::list<T>::const_iterator itr = lst.begin();
	int sz = lst.size();
	vec.resize(sz);
	for (int i=0; i<sz; i++, itr++)
	{
		vec[i] = *itr;
	}
}

#endif  // array_util_h

