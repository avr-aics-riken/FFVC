#include <stdarg.h>
#include "file.h"

#define BUFFER_SIZE   1024

using namespace std;


TtEndianType
tt_check_machine_endian()
{
	int v = 1;
	char* p = (char*)&v;

	if (p[0])
	{
		return TT_LITTLE_ENDIAN;
	}
	else if (p[sizeof(int)-1])
	{
		return TT_BIG_ENDIAN;
	}
	else
	{
		return TT_OTHER_ENDIAN;
	}
}

void 
tt_invert_byte_order(void* _mem, int size, int n)
{
	if (size == 1)
	{
		return;
	}

	char* mem = (char*)_mem;
	char c;
	int i;

	while (n)
	{
		for(i=0; i<size/2; i++)
		{
			c = mem[i];
			mem[i] = mem[size-1-i];
			mem[size-1-i] = c;
		}
		mem += size;
		n--;
	}
}

size_t
tt_read(void* data, size_t size, size_t n, FILE* fp, int inv)
{
    size_t r = fread(data, size, n, fp);

	if (inv)
	{
		tt_invert_byte_order(data, size, n);
	}

	return r;
}

size_t
tt_write(const void* _data, size_t size, size_t n, FILE* fp, int inv)
{
	const char* data = (const char*)_data;
	char* tmp = 0;

	if (inv)
	{
		int sz = size * n;
		tmp = new char[sz];
		for (int i=0; i<sz; i++) tmp[i] = data[i];
		tt_invert_byte_order(tmp, size, n);
		data = tmp;
	}

    size_t r = fwrite(data, size, n, fp);

	if (inv)
	{
		delete [] tmp;
	}

	return r;
}

void
tt_read(std::istream& is, void* _data, int size, int n, int inv)
{
	char* data = (char*)_data;

    is.read(data, size * n);

	if (inv)
	{
		tt_invert_byte_order(data, size, n);
	}
}

void
tt_write(std::ostream& os, const void* _data, int size, int n, int inv)
{
	const char* data = (const char*)_data;
	char* tmp = 0;

	if (inv)
	{
		int sz = size * n;
		tmp = new char[sz];
		for (int i=0; i<sz; i++) tmp[i] = data[i];
		tt_invert_byte_order(tmp, size, n);
		data = tmp;
	}

    os.write(data, size * n);

	if (inv)
	{
		delete [] tmp;
	}
}

char
tt_read_until_chars(std::istream& is, std::string& str, const char* s)
{
	char c;
	while ((c = is.get()))
	{
		int len = strlen(s);
		for (int i=0; i<len; i++)
		{
			if (c == s[i]) goto END;
		}
		str.push_back(c);
	}
END:
	return c;
}

char
tt_skip_until_chars(std::istream& is, const char* s)
{
	char c;
	while ((c = is.get()))
	{
		int len = strlen(s);
		for (int i=0; i<len; i++)
		{
			if (c == s[i]) goto END;
		}
	}
END:
	return c;
}

std::istream&
tt_skip_until(std::istream& is, std::string s)
{
	std::string token;
	while (is >> token && token != s) {}
	return is;
}

std::istream&
tt_skip_line(std::istream& is)
{
	return is.ignore(BUFFER_SIZE, '\n');
}

std::istream&
tt_skip_spaces(std::istream& is)
{
	char c;
	while (is.get(c))
	{
		if (!isspace(c))
		{
			is.putback(c);
			break;
		}
	}
	return is;
}

//=========================================================================
// class TtFileBase
//=========================================================================

TtFileBase::TtFileBase()
: m_fp(0), m_inv(0), m_flag(0)
{
}

TtFileBase::~TtFileBase()
{
	close();
}

int
TtFileBase::close()
{
	if (!m_fp) return 0;

	int r;

#ifndef WIN32
	if (isPipe())
	{
		r = pclose(m_fp);
	}
	else
#endif // !WIN32
	{
		r = fclose(m_fp);
	}

	m_fp = 0;

	return r;
}

int
TtFileBase::seek(long offset, int whence)
{
	int r = fseek(m_fp, offset, whence);

	if (r == 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

long
TtFileBase::seek()
{
	return ftell(m_fp);
}

//=========================================================================
// class TtInFile
//=========================================================================

TtInFile::TtInFile(const TtFileName& fname)
{
	open(fname);
}

FILE*
TtInFile::open(const TtFileName& fname)
{
	string ext = fname.suffix();

#ifndef WIN32
	if (ext == ".gz")
	{
		string cmd = string("gzip -dc ") + fname.path();
		m_fp = popen(cmd.c_str(), "r");
		isPipe(1);
	}
	else if (ext == ".bz2")
	{
		string cmd = string("bzip2 -dc ") + fname.path();
		m_fp = popen(cmd.c_str(), "r");
		isPipe(1);
	}
	else
#endif // !WIN32
	{
		m_fp = fopen(fname.c_str(), "rb");
	}

	return m_fp;
}

int
TtInFile::getline(char* data, int size)
{
	char* s = fgets(data, size, m_fp);

	// EOF or error
	if (s == NULL) return 0;

	return 1;
}

int
TtInFile::read(void* data, size_t size, size_t n)
{
	size_t r = tt_read(data, size, n, m_fp, m_inv);

	// error
	if (r != n) return 0;

	return 1;
}

//=========================================================================
// class TtOutFile
//=========================================================================

TtOutFile::TtOutFile(const TtFileName& fname)
{
	open(fname);
}

FILE*
TtOutFile::open(const TtFileName& fname)
{
	string ext = fname.suffix();

#ifndef WIN32
	if (ext == ".gz")
	{
		string cmd = string("gzip > ") + fname.path();
		m_fp = popen(cmd.c_str(), "w");
		isPipe(1);
	}
	else if (ext == ".bz2")
	{
		string cmd = string("bzip2 > ") + fname.path();
		m_fp = popen(cmd.c_str(), "w");
		isPipe(1);
	}
	else
#endif // !WIN32
	{
		m_fp = fopen(fname.c_str(), "wb");
	}

	return m_fp;
}

int
TtOutFile::putline(const char* data)
{
	int r = fputs(data, m_fp);

	// error
	if (r <= 0) return 0;

	return 1;
}

int
TtOutFile::putlinef(const char* fmt, ...)
{
	char buffer[BUFFER_SIZE];
	va_list ap;

	va_start(ap, fmt);
	vsprintf(buffer, fmt, ap);
	va_end(ap);

	return putline(buffer);
}

int
TtOutFile::write(const void* data, size_t size, size_t n)
{
	size_t r = tt_write(data, size, n, m_fp, m_inv);

	// error
	if (r != n) return 0;

	return 1;
}

//=========================================================================
// class TtIOFile
//=========================================================================

TtIOFile::TtIOFile(const TtFileName& fname)
{
	open(fname);
}

FILE*
TtIOFile::open(const TtFileName& fname)
{
	string ext = fname.suffix();

	m_fp = fopen(fname.c_str(), "r+b");
	if (m_fp == NULL)
	{
		m_fp = fopen(fname.c_str(), "w+b");
	}

	return m_fp;
}

