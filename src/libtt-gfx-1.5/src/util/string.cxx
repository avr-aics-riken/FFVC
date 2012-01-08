#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <algorithm>   // find()
#include "string.h"

#define BUFFER_SIZE 1024

#ifndef whitespace
#define whitespace(c) (((c) == ' ') || ((c) == '\t'))
#endif

//=========================================================================

std::string
va2str(const char* fmt, ...)
{
	char buffer[BUFFER_SIZE];
	va_list ap;

	va_start(ap, fmt);
	vsprintf(buffer, fmt, ap);
	va_end(ap);

	return std::string(buffer);
}

//=========================================================================

char*
tt_remove_newline(char* c_str)
{
	int n = strlen(c_str) - 1;
	if (c_str[n] == '\n')
	{
		c_str[n] = '\0';
	}
	return c_str;
}

std::string&
tt_remove_newline(std::string& str)
{
	int n = str.size() - 1;
	if (str[n] == '\n')
	{
		str.erase(n, 1);
	}
	return str;
}

// Strip whitespace from the start and end of STRING.
// Return a pointer into STRING.
char*
tt_strip_ws(char* c_str)
{
  char *s, *t;

  for (s = c_str; whitespace(*s); s++) ;

  if (*s == 0)
    return s;

  t = s + strlen(s) - 1;
  while (t > s && whitespace(*t))
    t--;
  *++t = '\0';

  return s;
}

void
tt_split(std::vector<std::string>& s, const char* cs, char separator)
{
	const char* p   = cs;
	const char* end = cs + strlen(cs);

	const char* q;
	while ((q = std::find(p, end, separator)) < end)
	{
		s.push_back(std::string(p, q));
		p = q + 1;
	}
	s.push_back(std::string(p, end));
}

std::string
tt_zero2(int n)
{
	char buffer[20];
	sprintf(buffer, "%02d", n);
	return std::string(buffer);
}

std::string
tt_zero4(int n)
{
	char buffer[20];
	sprintf(buffer, "%04d", n);
	return std::string(buffer);
}

std::string
tt_to_string(int n)
{
	char buffer[20];
	sprintf(buffer, "%d", n);
	return std::string(buffer);
}

std::string
tt_to_string(float n)
{
	char buffer[20];
	sprintf(buffer, "%f", n);
	return std::string(buffer);
}

std::string
tt_to_string(double n)
{
	char buffer[20];
	sprintf(buffer, "%f", n);
	return std::string(buffer);
}

std::string
tt_to_string(const std::type_info& info, const std::string& format)
{
	int f = (format == "short") ? 1 : 0;
	std::string s = "other";

	if      (info == typeid(char))            s = f ? "c"  : "char";
	else if (info == typeid(short))           s = f ? "s"  : "short";
	else if (info == typeid(int))             s = f ? "i"  : "int";
	else if (info == typeid(long))            s = f ? "l"  : "long";
	else if (info == typeid(float))           s = f ? "f"  : "float";
	else if (info == typeid(double))          s = f ? "d"  : "double";
	else if (info == typeid(unsigned char))   s = f ? "uc" : "unsigned char";
	else if (info == typeid(unsigned short))  s = f ? "us" : "unsigned short";
	else if (info == typeid(unsigned int))    s = f ? "ui" : "unsigned int";
	else if (info == typeid(unsigned long))   s = f ? "ul" : "unsigned long";
	else if (info == typeid(long double))     s = f ? "ld" : "long double";

	return s;
}

std::string
tt_bytes_to_string(unsigned long long bytes)
{
	static unsigned long long KB = 1024;
	static unsigned long long MB = 1024*1024;
	static unsigned long long GB = 1024*1024*1024;
	static unsigned long long TB = 1024LL*1024LL*1024LL*1024LL;
	static unsigned long long PB = 1024LL*1024LL*1024LL*1024LL*1024LL;

	char buffer[128];

	if (bytes < KB)
	{
		sprintf(buffer, "%lld bytes", bytes);
	}
	else if (bytes < MB)
	{
		sprintf(buffer, "%.1f KB", double(bytes)/KB);
	}
	else if (bytes < GB)
	{
		sprintf(buffer, "%.1f MB", double(bytes)/MB);
	}
	else if (bytes < TB)
	{
		sprintf(buffer, "%.1f GB", double(bytes)/GB);
	}
	else if (bytes < PB)
	{
		sprintf(buffer, "%.1f TB", double(bytes)/TB);
	}
	else
	{
		sprintf(buffer, "%.1f PB", double(bytes)/PB);
	}

	return std::string(buffer);
}

