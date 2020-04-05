//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################
//
#include "FileSystemUtil.h"

#include <vector>
#include <string>

#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h> 
#include <unistd.h>

#include <string.h>
#include "mydebug.h"

#include <iostream>
#include <fstream>
#include <errno.h>

#include <dirent.h>

#include "cdm_DFI.h"

namespace  {
	std::string parseFilename(const char* fpath, const char* splitChar)
	{
		std::string f(fpath);
		size_t p = f.rfind(splitChar);
		if (p != std::string::npos)
			f.erase(f.begin(), f.begin() + p + 1);
		return f;
	}
	/////////////////////////////////////////////////////////////////////////////

	class Dir
	{
	public:
		Dir(const char* dirname);
	
		unsigned int GetFileCount() const;
		unsigned int GetDirCount() const;
		const std::string& GetFile(unsigned int i) const;
		const std::string& GetDir(unsigned int i) const;
		std::string GetFilePath(const char* filename) const;
		bool IsOpened() const;
		const std::string& GetPath() const;
		const std::string& GetName() const;
		bool CreateDir(const char* dirname);
		
	private:
		std::string m_dirname;
		std::string m_dirpath;
		bool m_opened;
		std::vector<std::string> m_files;
		std::vector<std::string> m_dirs;
	};
	
	Dir::Dir(const char* dirname)
	{
		m_opened = false;
		m_dirpath = std::string(dirname);
		
		m_dirname = parseFilename(dirname,  "/");
		DIR* pDir = opendir(dirname);	
		if (!pDir)
			return;
		
		dirent* pEnt = readdir(pDir);
	    while (pEnt)
		{
	        if ( strcmp( pEnt->d_name, "." ) &&
	             strcmp( pEnt->d_name, ".." ) )
			{
				std::string wPathName = std::string(dirname) + std::string("/") + std::string(pEnt->d_name);
				struct stat wStat;
	            if ( stat( wPathName.c_str(), &wStat ) )
	                break;
	            
	            if ( S_ISDIR( wStat.st_mode ) )
					m_dirs.push_back(pEnt->d_name);
	            else
	                m_files.push_back(pEnt->d_name);
	        }
			pEnt = readdir(pDir);
	    }
		closedir(pDir);
		m_opened = true;
	}
		
	bool Dir::IsOpened() const
	{
		return m_opened;
	}
		
	unsigned int Dir::GetFileCount() const
	{
		return static_cast<unsigned int>(m_files.size());
	}
		
	unsigned int Dir::GetDirCount() const
	{
		return static_cast<unsigned int>(m_dirs.size());
	}
		
	const std::string& Dir::GetFile(unsigned int i) const
	{
		static std::string nullstring;
		if (i >= m_files.size())
			return nullstring;
		
		return m_files[i];
	}
		
	const std::string& Dir::GetDir(unsigned int i) const
	{
		static std::string nullstring;
		if (i >= m_dirs.size())
			return nullstring;
		
		return m_dirs[i];
	}
	
	std::string Dir::GetFilePath(const char* filename) const
	{
		return m_dirpath + "/" + filename;
	}
		
	const std::string& Dir::GetPath() const
	{
		return m_dirpath;
	}
	
	const std::string& Dir::GetName() const
	{
		return m_dirname;
	}
	
	bool Dir::CreateDir(const char* dirname)
	{
		std::string dirpath = GetFilePath(dirname);
		std::vector<std::string>::iterator it = std::find(m_files.begin(), m_files.end(), dirname);
		if (it != m_files.end())
			return false;
		it = std::find(m_dirs.begin(), m_dirs.end(), dirname);
		if (it != m_dirs.end())
			return true;
		
		const mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
		if ( mkdir(dirpath.c_str(), mode) != 0 )
		{
			usleep(10000);
			Dir check(dirpath.c_str());
			return check.IsOpened();
		}
	
		return true;
	}
	
	/////////////////////////////////////////////////////////////////////////////
	bool split(const std::string& input, const char delimiter, std::vector<std::string>& output)
	{
		output.clear();
		
		std::string istr = input;
		
		size_t p = 0;
		while( (p = istr.find(delimiter)) != std::string::npos ){
			if(p != 0){
				output.push_back( istr.substr(0, p) );
			}
			istr.erase(istr.begin(), istr.begin() + p + 1);
		}
	
		if( istr.length() != 0 ) output.push_back(istr);
	
		return true;
	}
	
} // namespace

namespace BVX_IO {
	bool CreateDirectory(const std::string& path)
  {
		std::vector<std::string> dirList;
		split(path, '/', dirList);
		
		std::string npath;
    
		if( isAbsolutePath(path) )
    {
			npath = std::string("/");
		}
    else
    {
			npath = std::string("./");
		}
	
		for(std::vector<std::string>::iterator it = dirList.begin(); it != dirList.end(); ++it){
			Dir dir(npath.c_str());
			
			if( *it == std::string(".")  ){ continue; }
			if( *it == std::string("..") ){ npath += *it + std::string("/"); continue; }
	
			if( !dir.CreateDir(it->c_str()) ){ 
				stamped_printf("can not create directory : %s\n", (npath + *it).c_str());
				return false;
			}
			
			npath += *it + std::string("/");
		}
	
		return true;
	}
  
  
  void CheckDir(std::string dirstr)
  {
      
    if( dirstr.size() == 0 ) {
      //printf("\toutput current directory\n");
      return;
    }
    
    DIR* dir;
    if( !(dir = opendir(dirstr.c_str())) ) {
      if( errno == ENOENT ) {
        mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
        //if ( !FBUtility::mkdirs(dirstr.c_str()) )
        //CDM.20131008.s
        //if ( !cdm_DFI::MakeDirectorySub(dirstr) )
        if ( cdm_DFI::MakeDirectorySub(dirstr) != 0 )
        {
          printf("\tCan't generate directory(%s).\n", dirstr.c_str());
          Exit(0);
        }
      }
      else {
        printf("Directory open error.(%s)", dirstr.c_str());
        Exit(0);
      }
    }
    else {
      if( closedir(dir) == -1 ) {
        printf("Directory close error.(%s)", dirstr.c_str());
        Exit(0);
      }
    }
    
    return;
  }
  
} // namespace BVX_IO

