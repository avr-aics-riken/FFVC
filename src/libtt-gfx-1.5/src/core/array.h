#ifndef array_h
#define array_h

#include <assert.h>

//=========================================================================
// class TtArray<T>
//=========================================================================

template<class T>
class TtArray
{
public:
	TtArray() : m_data(0), m_size(0), m_capacity(0) {}
	TtArray(int size);
	TtArray(int size, int capacity);
	~TtArray();

	// destroy the original data
	void reallocate(int size, int capacity);
	void resize(int size);

	// keep the original data
	void resizeConservatively(int size);

	void pack();

	void push_back(const T& v);

	void clear();
	void clear(const T& v);

	void setSize(int size) {
		assert(size >= 0 && size <= m_capacity);
		m_size = size;
	}
	int  size() const { return m_size; }
	int  capacity() const { return m_capacity; }

	T& operator [](int i) {
		assert(i >= 0 && i < m_size);
		return m_data[i];
	}
	const T& operator [](int i) const {
		assert(i >= 0 && i < m_size);
		return m_data[i];
	}

	operator       T*()       { return m_data; }
	operator const T*() const { return m_data; }
	      T* ptr()         { return m_data; }
	const T* ptr()   const { return m_data; }

	      T* begin()       { return m_data; }
	const T* begin() const { return m_data; }
	      T* end()         { return m_data + m_size; }
	const T* end()   const { return m_data + m_size; }

	      T& front()       { return m_data[0]; }
	const T& front() const { return m_data[0]; }
	      T& back()        { return m_data[m_size - 1]; }
	const T& back()  const { return m_data[m_size - 1]; }

private:
	T*   m_data;
	int  m_size;
	int  m_capacity;

	// prohibit the copy constructor
	TtArray(const TtArray& o) {}
	// prohibit the copy operator
	TtArray& operator =(const TtArray& o) { return *this; }

	void __destroy();
	void __allocate(int capacity);
};

//=========================================================================

template<class T> inline
TtArray<T>::TtArray(int size)
: m_data(0), m_size(0), m_capacity(0) {
	__allocate(size);
	m_size = size;
}

template<class T> inline
TtArray<T>::TtArray(int size, int capacity)
: m_data(0), m_size(0), m_capacity(0) {
	__allocate(capacity);
	m_size = size;
}

template<class T> inline
TtArray<T>::~TtArray() {
	__destroy();
}

template<class T> inline
void TtArray<T>::__allocate(int capacity) {
	if (capacity > 0) {
		m_data = new T[capacity];
		m_capacity = capacity;
	}
}

template<class T> inline
void TtArray<T>::__destroy() {
	if (m_data) {
		delete [] m_data;
		m_data = 0;
		m_size = 0;
		m_capacity = 0;
	}
}

template<class T> inline
void TtArray<T>::reallocate(int size, int capacity) {
	if (capacity > 0) {
		__destroy();
		__allocate(capacity);
		m_size = size;
	}
}

template<class T> inline
void TtArray<T>::resize(int size) {
	if (size > m_capacity) {
		reallocate(size, size);
	} else {
		m_size = size;
	}
}

template<class T> inline
void TtArray<T>::resizeConservatively(int size) {
	if (size > m_capacity) {
		T* old_data = m_data;
		int old_size = m_size;

		__allocate(size);
		m_size = size;
		for (int i=0; i<old_size; i++) {
			m_data[i] = old_data[i];
		}
		delete [] old_data;
	} else {
		m_size = size;
	}
}

template<class T> inline
void TtArray<T>::push_back(const T& v) {
	if (m_size == m_capacity) {
		if (m_size == 0)
		{
			__allocate(1);
		}
		else
		{
			T* old_data = m_data;
			int old_size = m_size;

			__allocate(m_size * 2);
			for (int i=0; i<old_size; i++) {
				m_data[i] = old_data[i];
			}
			delete [] old_data;
		}
	}
	m_data[m_size] = v;
	m_size++;
}

template<class T> inline
void TtArray<T>::clear() {
	for (int i=0; i<m_size; i++) m_data[i] = T();
}

template<class T> inline
void TtArray<T>::clear(const T& v) {
	for (int i=0; i<m_size; i++) m_data[i] = v;
}

template<class T> inline
void TtArray<T>::pack() {
	if (m_capacity > m_size)
	{
		T* old_data = m_data;
		int old_size = m_size;

		__allocate(m_size);
		for (int i=0; i<old_size; i++) {
			m_data[i] = old_data[i];
		}
		delete [] old_data;
	}
}

#endif  // array_h

