// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file main.C
 * @brief layoutのmain関数
 * @author kero
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
  // 初期処理

  // 並列管理クラスのインスタンスと初期化
  // ここでMPI_Initも行う
  cpm_ParaManager* m_paraMngr = cpm_ParaManager::get_instance(argc, argv);
  if ( !m_paraMngr ){
    std::cerr << "fault MPI Init" << std::endl;
    return 0;
  }

  // ##################################################################
  // LAYOUT classのインスタンス
  LAYOUT layout;
  
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
