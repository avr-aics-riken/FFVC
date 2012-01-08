#ifndef tt_filename_h
#define tt_filename_h

#include <string>
#include <iostream>

//=========================================================================
// class TtFileName
//=========================================================================

class TtFileName
{
public:
	TtFileName() {}
	TtFileName(const char* _path) { path(_path); }
	TtFileName(const std::string& _path) { path(_path); }

	void path(const std::string& path);

	const char* c_str()           const { return m_path.c_str(); }
	std::string str()             const { return m_path; }

	// path()  == dirname() + fname()
	// fname() == name() + ext()
	std::string orig()            const { return m_orig; }
	std::string path()            const { return m_path; }
	std::string dirname()         const;
	std::string fname()           const;
	std::string name()            const;
	std::string ext()             const;
	std::string prefix(int n = 0) const;
	std::string suffix(int n = 0) const;

	bool isAbsPath() const;

private:
	std::string m_orig, m_path;
	std::string::size_type m_idx_fname, m_idx_ext;
};

std::istream& operator >>(std::istream& is, TtFileName& fname);
std::ostream& operator <<(std::ostream& os, const TtFileName& fname);

//=========================================================================
// class TtFrameName
//
// EXAMPLE:
//   TtFrameName frname("frame.%04d.ppm");
//=========================================================================

class TtFrameName : public TtFileName
{
public:
	TtFrameName(const std::string& fname);

	std::string path (int nr) const;
	std::string fname(int nr) const;
	std::string name (int nr) const;
};

#endif  // tt_filename_h

