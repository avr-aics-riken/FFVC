#ifndef ptr_set_h
#define ptr_set_h

#include <string>
#include <vector>

//=========================================================================
// PointerSetEntry
//=========================================================================

class PointerSetEntry
{
public:
	void* m_data;
};

//=========================================================================
// PointerSet
//=========================================================================

template<class T>
class PointerSet
{
public:
	PointerSet() {}
	virtual ~PointerSet() { clearData(); }

	int size() const { return m_list.size(); }

	T* getData(int i) const;
	void addData(T* data);

	void deleteData(int i);
	void clearData();

private:
	std::vector<PointerSetEntry> m_list;
};

//=========================================================================

template<class T>
inline T*
PointerSet<T>::getData(int i) const
{
	int sz = m_list.size();
	if (sz > 0 && i >= 0 && i < sz)
		return static_cast<T*>(m_list[i].m_data);
	else
		return 0;
}

template<class T>
inline void
PointerSet<T>::addData(T* data)
{
	PointerSetEntry entry;
	entry.m_data = data;
	m_list.push_back(entry);
}

template<class T>
inline void
PointerSet<T>::deleteData(int i)
{
	int sz = m_list.size();
	if (sz > 0 && i >= 0 && i < sz)
	{
		delete getData(i);
		m_list.erase(m_list.begin() + i);
	}
}

template<class T>
inline void
PointerSet<T>::clearData()
{
	int sz = m_list.size();
	for (int i=0; i<sz; i++) delete getData(i);
	m_list.clear();
}

#endif  // ptr_set_h

