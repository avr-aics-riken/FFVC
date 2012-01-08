#ifndef tt_file_h
#define tt_file_h

#include <iostream>
#include <stdio.h>
#include <string>
#include "filename.h"

enum TtEndianType {TT_OTHER_ENDIAN, TT_LITTLE_ENDIAN, TT_BIG_ENDIAN};

TtEndianType tt_check_machine_endian();
void tt_invert_byte_order(void* mem, int size, int n);

size_t tt_read(void* data, size_t size, size_t n, FILE* fp, int inv = 0);
size_t tt_write(const void* data, size_t size, size_t n, FILE* fp, int inv = 0);

void tt_read(std::istream& is, void* data, int size, int n, int inv = 0);
void tt_write(std::ostream& os, const void* data, int size, int n, int inv = 0);

char tt_read_until_chars(std::istream& is, std::string& str, const char* s);
char tt_skip_until_chars(std::istream& is, const char* s);
std::istream& tt_skip_until(std::istream& is, std::string s);
std::istream& tt_skip_line(std::istream& is);
std::istream& tt_skip_spaces(std::istream& is);

//=========================================================================
// class TtFileBase
//=========================================================================

class TtFileBase
{
public:
	TtFileBase();
	virtual ~TtFileBase();

	virtual FILE* open(const TtFileName& fname) = 0;
	virtual int   close();

	int  seek(long offset, int whence = SEEK_SET);
	long seek();

	FILE* fp() { return m_fp; }
	int  getInverseFlag() const { return m_inv; }
	void setInverseFlag(int b) { m_inv = b; }

protected:
	int  isPipe() const { return m_flag; }
	void isPipe(int b)  { m_flag = b; }

	FILE* m_fp;
	int m_inv;

private:
	int m_flag;
};

//=========================================================================
// class TtInFile
//=========================================================================

class TtInFile : virtual public TtFileBase
{
public:
	TtInFile() {}
	TtInFile(const TtFileName& fname);

	virtual FILE* open(const TtFileName& fname);

	// read a line
	int getline(char* data, int size);

	// read binary data
	int read(void* data, size_t size, size_t n);
};

//=========================================================================
// class TtOutFile
//=========================================================================

class TtOutFile : virtual public TtFileBase
{
public:
	TtOutFile() {}
	TtOutFile(const TtFileName& fname);

	virtual FILE* open(const TtFileName& fname);

	// write a line
	int putline(const char* data);
	int putlinef(const char* fmt, ...);

	// write binary data
	int write(const void* data, size_t size, size_t n);
};

//=========================================================================
// class TtIOFile
//=========================================================================

class TtIOFile : public TtInFile, public TtOutFile
{
public:
	TtIOFile() {}
	TtIOFile(const TtFileName& fname);

	FILE* open(const TtFileName& fname);
};

#endif  // tt_file_h

