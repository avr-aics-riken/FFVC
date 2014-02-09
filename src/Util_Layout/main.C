//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   main.C
 * @brief  layoutのmain関数
 * @author aics
 */

#include "layout.h"

void usage(const char *progname)
{
  std::cerr
  << "Usage: " << progname << " <option> filename <option> dirname <options> mchname <options>\n"
  << " filename: file name when -f specified\n"
  << " Options:\n"
  << "  -a              : 0 step skip option\n"
  << "  -f filename     : input file name\n"
  << "  -d dirname      : directory name\n"
  << "  -m machinename  : machine name\n"
  << "  -h              : Show usage and exit\n"
  << std::endl;
}

int main( int argc, char **argv )
{
  char *progname = argv[0];
  bool skip0 = false;
  bool out_layout = false;
  bool out_layout_dir = false;
  bool out_layout_machine = false;
  string fname;
  string dname;
  string mname;
  
  int opt;
  while ((opt = getopt(argc, argv, "af:d:m:h")) != -1) {
    switch (opt) {
      case 'a':
        skip0 = true;
        break;
      case 'f':
        fname = optarg;
        out_layout = true;
        break;
      case 'd':
        dname = optarg;
        out_layout_dir = true;
        break;
      case 'm':
        mname = optarg;
        out_layout_machine = true;
        break;
      case 'h': // Show usage and exit
        usage(progname);
        return 0;
        break;
      case ':': // not find argument
        usage(progname);
        return 0;
        break;
      case '?': // unknown option
        usage(progname);
        return 0;
        break;
    }
  }
  
  if( !(out_layout) ) {
    usage(progname);
    return 0;
  }
  
  // ##################################################################
  // LAYOUT classのインスタンス
  LAYOUT layout;
  
  // ##################################################################
  // 初期処理
  
  // 並列管理クラスのインスタンスと初期化
  // ここでMPI_Initも行う
  if ( !layout.importCPM(cpm_ParaManager::get_instance(argc, argv)) )
  {
    return CPM_ERROR_PM_INSTANCE;
  }
  
  // ##################################################################
  // 引数のセット
  layout.SetInput(skip0, fname);
  
  // ##################################################################
  // 入力ファイルの読み込み
  layout.ReadInit();
  
  // ##################################################################
  // コマンドラインからの入力を反映 ---> コマンドラインからが優先
  if( out_layout_dir )     layout.SetDirName(dname);
  if( out_layout_machine ) layout.SetMachineName(mname);
  
  // ##################################################################
  // dfiファイルの読み込み
  layout.ReadDfiFiles();
  
  // ##################################################################
  // layoutファイルの出力
  layout.OutputLayout();
  
  return 0;
}
