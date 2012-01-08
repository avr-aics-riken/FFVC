#ifndef tt_string_h
#define tt_string_h

#include <sys/types.h>
#include <typeinfo>
#include <vector>
#include <string>

std::string va2str(const char* fmt, ...);

char* tt_remove_newline(char* c_str);
std::string& tt_remove_newline(std::string& str);
char* tt_strip_ws(char* c_str);
void tt_split(std::vector<std::string>& s, const char* cs, char separator = ' ');
std::string tt_zero2(int n);
std::string tt_zero4(int n);

std::string tt_to_string(int n);
std::string tt_to_string(float n);
std::string tt_to_string(double n);
std::string tt_to_string(const std::type_info& info, const std::string& format = "short");
std::string tt_bytes_to_string(unsigned long long bytes);

#endif  // tt_string_h

