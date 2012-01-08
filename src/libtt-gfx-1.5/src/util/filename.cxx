#include <stdlib.h>     // getenv
#ifndef WIN32
#include <pwd.h>        // getpwnam
#include <sys/types.h>  // getpwnam
#endif // !WIN32
#include <algorithm>    // count
#include "filename.h"

std::string filename_expand_homedir(std::string path);
std::string filename_expand_env(std::string path);
std::string filename_erase_redundant_slash(std::string path);
std::string filename_expand_path(std::string path);
std::string::size_type filename_find_fname(const std::string& path);
std::string::size_type filename_find_ext(const std::string& path);

// expand '~', "~usr_name"
std::string
filename_expand_homedir(std::string path) {
	if (path[0] == '~') {
		if (path.size() == 1 || path[1] == '/') {
			// expand '~'
			path.replace(0, 1, getenv("HOME"));
		} else {
			// expand "~usr_name"
			std::string::size_type i = path.find('/');
			const std::string& usr_name = path.substr(1, i-1);
#ifndef WIN32
			struct passwd* pw = getpwnam(usr_name.c_str());
			if (pw) path.replace(0, i, pw->pw_dir);
#endif // !WIN32
		}
	}

	return path;
}

// expand ${...}
std::string
filename_expand_env(std::string path) {
	std::string::size_type i1, i2;

	while (1) {
		i1 = path.find("${"); if (i1 == std::string::npos) break;
		i2 = path.find("}");  if (i2 == std::string::npos) break;

		const std::string& s = path.substr(i1+2, i2-i1-2);
		char* env = getenv(s.c_str());
		if (env) path.replace(i1, i2-i1+1, env);
		else     path.replace(i1, i2-i1+1, "");
	}

	return path;
}

std::string
filename_erase_redundant_slash(std::string path) {
	std::string::size_type i;

	while (1) {
		i = path.find("//");
		if (i == std::string::npos) break;

		path.erase(i, 1);
	}

	return path;
}

std::string
filename_expand_path(std::string path) {
	path = filename_expand_homedir(path);
	path = filename_expand_env(path);
	path = filename_erase_redundant_slash(path);
	return path;
}

std::string::size_type
filename_find_fname(const std::string& path) {
	std::string::size_type i = path.rfind('/');
	if (i == std::string::npos)
		return 0;
	else
		return i+1;
}

std::string::size_type
filename_find_ext(const std::string& path) {
	std::string::size_type i = path.rfind('.');
	if (i != std::string::npos) {
		const std::string& s = path.substr(i);
		if (s == ".Z" || s == ".gz" || s == ".bz2") {
			std::string::size_type i2 = path.rfind('.', i-1);
			if (i2 != std::string::npos) i = i2;
		}
	}
	return i;
}

//=========================================================================
// class TtFileName
//=========================================================================

void
TtFileName::path(const std::string& path) {
	m_orig      = path;
	m_path      = filename_expand_path(path);
	m_idx_fname = filename_find_fname(m_path);
	m_idx_ext   = filename_find_ext(m_path);
}

std::string
TtFileName::dirname() const {
	return m_path.substr(0, m_idx_fname);
}

std::string
TtFileName::fname() const {
	return m_path.substr(m_idx_fname);
}

std::string
TtFileName::name() const {
	return m_path.substr(m_idx_fname, m_idx_ext - m_idx_fname);
}

std::string
TtFileName::ext() const {
	if (m_idx_ext == std::string::npos) return "";
	return m_path.substr(m_idx_ext);
}

std::string
TtFileName::prefix(int n) const {
	std::string::size_type i1 = 0;
	std::string::size_type i2 = m_idx_fname;
	while (n >= 0) {
		i1 = i2;
		i2 = m_path.find('.', i1);
		if (i2 == std::string::npos) return "";
		++i2;
		--n;
	}
	return m_path.substr(i1, i2-i1);
}

std::string
TtFileName::suffix(int n) const {
	std::string::size_type i1 = 0;
	std::string::size_type i2 = std::string::npos;
	while (n >= 0) {
		i1 = i2;
		i2 = m_path.rfind('.', i1);
		if (i2 == std::string::npos) return "";
		--i2;
		--n;
	}
	return m_path.substr(i2+1, i1-i2);
}

bool
TtFileName::isAbsPath() const {
	return (m_path[0] == '/');
}

std::istream& operator >>(std::istream& is, TtFileName& fname) {
	std::string s;
	is >> s;
	fname.path(s);
	return is;
}

std::ostream& operator <<(std::ostream& os, const TtFileName& fname) {
	os << fname.path();
	return os;
}

//=========================================================================
// class TtFrameName
//=========================================================================

TtFrameName::TtFrameName(const std::string& fname)
: TtFileName(fname) {
}

std::string
TtFrameName::path(int nr) const {
	return dirname() + name(nr) + ext();
}

std::string
TtFrameName::fname(int nr) const {
	return name(nr) + ext();
}

std::string
TtFrameName::name(int nr) const {
	char buf[256];
	std::string orgname = TtFileName::name();
	sprintf(buf, orgname.c_str(), nr);
	return std::string(buf);
}

