//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   layout.C
 * @brief  LAYOUT Class
 * @author kero
 */

#include "layout.h"


// #################################################################
// コンストラクタ
LAYOUT::LAYOUT()
{
  ndfi=0;
  DI=NULL;
  
  procGrp = 0;
  myRank  = -1;
  numProc = 0;
  
  IS_DivideFunc=OFF;
}


// #################################################################
// デストラクタ
LAYOUT::~LAYOUT()
{
  if( DI ) delete [] DI;
}

// #################################################################
//
void LAYOUT::SetInput(bool m_skip0, string m_fname)
{
  skip0 = m_skip0;
  fname = m_fname;
  return;
}

// #################################################################
//
void LAYOUT::ReadInit()
{
  
  // ------------------------------------
  FILE* fp = NULL;
  
  // TPインスタンス生成
  TextParser tpCntl;

  
  //入力ファイルをセット
  int ierror = tpCntl.read(fname);
  
  //入力ファイルの読み込み--->パラメータのセット
  ReadInputFile(&tpCntl);
  
  //TextParserの破棄
  tpCntl.remove();
  
  return;
}


// #################################################################
//
void LAYOUT::ReadInputFile(TextParser* tpCntl)
{
  string str,buff;
  string label,label_base,label_leaf;
  
  // node数の取得
  int nnode=0;
  label_base = "/LayoutData";
  if ( tpCntl->chkNode(label_base) )  //nodeがあれば
  {
    nnode = tpCntl->countLabels(label_base);
  }
  
  // dfi_nameの取得
  dfi_name.clear();
  for (int i=0; i<nnode; i++) {
    
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,4).c_str(), "list") ) continue;
    label=label_base+"/"+str;
    
    if ( !(tpCntl->getInspectedValue(label, buff )) ) {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    //FList[ilist].name = str;
    dfi_name.push_back(buff.c_str());
    
  }
  
#if 0
  cout << "dfi_name.size() = " << dfi_name.size() << endl;
  vector<string>::const_iterator it;
  for (it = dfi_name.begin(); it != dfi_name.end(); it++) {
    cout << "name = " << (*it).c_str() << endl;
  }
#endif
  
  // dfi_nameの取得
  mname.clear();
  dname.clear();
  rankis.clear();
  rankie.clear();
  label_base = "/LayoutData";
  for (int i=0; i<nnode; i++) {
    
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,6).c_str(), "divide") ) continue;
    
    label=label_base+"/"+str+"/machine";
    if ( !(tpCntl->getInspectedValue(label, buff )) ) {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    mname.push_back(buff.c_str());
    
    label=label_base+"/"+str+"/rank";
    int v[2];
    for (int n=0; n<2; n++) v[n]=0;
    if ( !(tpCntl->getInspectedVector(label, v, 2)) ) {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    rankis.push_back(v[0]);
    rankie.push_back(v[1]);
    
    label=label_base+"/"+str+"/dir";
    if ( !(tpCntl->getInspectedValue(label, buff )) ) {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    dname.push_back(buff.c_str());
  }
  
  
  //出力ディレクトリの指定 ---> 実行オプションよりこちらが優先される
  label = "/LayoutData/OutputDir";
  if ( (tpCntl->getInspectedValue(label, str )) )
  {
    dirname=str;
    CheckDir(dirname);
    if( dirname.size() != 0 ) dirname=dirname+"/";
  }
  
  // DivideFunc ---> 出力を項目別にファイル分割するオプション
  label = "/LayoutData/FFVDivideFunc";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  IS_DivideFunc = ON;
    else if( !strcasecmp(str.c_str(), "off") ) IS_DivideFunc = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // FileNameGrid --- option
  label = "/LayoutData/Plot3dOptions/FileNameGrid";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    basename_g = "PLOT3DoutputGrid";
  }
  else
  {
    basename_g = str;
  }
  if ( basename_g.empty() )
  {
    basename_g = "PLOT3DoutputGrid";
  }
  
  // FileNameFunc --- option
  label = "/LayoutData/Plot3dOptions/FileNameFunc";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    basename_f = "PLOT3Doutput";
  }
  else
  {
    basename_f = str;
  }
  if ( basename_f.empty() )
  {
    basename_f = "PLOT3Doutput";
  }
  
}


// #################################################################
//
void LAYOUT::ReadDfiFiles()
{
  int ic;
  vector<string>::const_iterator it;
  
  // allocate dfi info class
  ndfi=dfi_name.size();
  DI = new DfiInfo[ndfi];
  
  // ランク情報をセット
  ic=0;
  for (it = dfi_name.begin(); it != dfi_name.end(); it++) {
    DI[ic].setRankInfo(paraMngr, procGrp);
    ic++;
  }
  
  // set dfi info class
  ic=0;
  for (it = dfi_name.begin(); it != dfi_name.end(); it++) {
    string fname=(*it).c_str();
    DI[ic].ReadDfiFile(fname);
    ic++;
  }
}


// #################################################################
//
void LAYOUT::SetDirName(string m_dname)
{
  dname.clear();
  for(int i=0; i<rankis.size(); i++ ) {
    dname.push_back(m_dname.c_str());
  }
}


// #################################################################
//
void LAYOUT::SetMachineName(string m_mname)
{
  mname.clear();
  for(int i=0; i<rankis.size(); i++ ) {
    mname.push_back(m_mname.c_str());
  }
}


// #################################################################
//
void LAYOUT::OutputLayout()
{
  string prefix;
  
  if ( IS_DivideFunc == ON )
  {
    for(int i=0;i<ndfi;i++){
      prefix=DI[i].Prefix;
      DI[i].Prefix=basename_g + "_" + prefix;
      OutLayGrid(&DI[i]);
      DI[i].Prefix=basename_f + "_" + prefix;
      OutLayFunc(&DI[i]);
    }
  }
  else
  {
    DI[0].Prefix=basename_g;
    OutLayGrid(&DI[0]);
    DI[0].Prefix=basename_f;
    OutLayFunc(&DI[0]);
  }
  
}


// #################################################################
//
void LAYOUT::OutLayGrid(DfiInfo *D)
{
  FILE* fp;
  int ifl=8;
  string layoutfile,nodefile;
  string d_name;
  string m_name;
  
  char tmp[FB_FILE_PATH_LENGTH];
  char line[FB_BUFF_LENGTH];
  
  int len ;
  int fnsize;
  int ierror;
  
  //int is=0;
  //if(skip0) is=1;
  
  // error check
  //if( D->step.size() == 0 ) {
  if( D->Sc.size() == 0 ) {
    fprintf(stderr, "error : step size == 0\n");
    return;
  }
  if( D->NodeInfoSize == 0 ) {
    fprintf(stderr, "error : node size == 0\n");
    return;
  }
  
  //xyz file
  
  //xyz.layoutファイルオープン
  //layoutfile = Generate_LayoutFileName(D->Prefix, "xyz", D->step[0]);
  layoutfile = Generate_LayoutFileName(D->Prefix, D->Sc[0]->step);
  layoutfile = dirname + layoutfile;
  if( (fp = fopen(layoutfile.c_str(), "w")) == NULL ) {
    fprintf(stderr, "Can't open file.(%s)\n", layoutfile.c_str());
    return;
  }
  
  //ファイル書き出し
  fprintf(fp,"FIELDVIEW LAYOUT 1\n");
  for(int j=0; j< D->NodeInfoSize; j++ ) {
    nodefile = Generate_FileName_Free(D->Prefix, "xyz", D->Sc[0]->step, D->Node[j].RankID, true);
    
    int found=0;
    for(int i=0; i<rankis.size(); i++ ) {
      if(rankis[i]<=D->Node[j].RankID && D->Node[j].RankID<=rankie[i]){
        m_name=mname[i];
        d_name=dname[i];
        found=1;
        break;
      }
    }
    if(!found){
      fprintf(stderr, "error : rank not found (%s)\n", nodefile.c_str());
      return;
    }
    
    fprintf(fp,"%s\n", nodefile.c_str());
    fprintf(fp,"%s\n", m_name.c_str());
    fprintf(fp,"%s\n", d_name.c_str());
  }
  
  //xyz.layoutファイルクローズ
  fclose(fp);
  
}


// #################################################################
//
void LAYOUT::OutLayFunc(DfiInfo *D)
{
  FILE* fp;
  int ifl=8;
  string layoutfile,nodefile;
  string d_name;
  string m_name;
  
  char tmp[FB_FILE_PATH_LENGTH];
  char line[FB_BUFF_LENGTH];
  
  int len ;
  int fnsize;
  int ierror;
  
  int is=0;
  if(skip0) is=1;
  
  // error check
  //if( D->step.size() == 0 ) {
  if( D->Sc.size() == 0 ) {
    fprintf(stderr, "error : step size == 0\n");
    return;
  }
  if( D->NodeInfoSize == 0 ) {
    fprintf(stderr, "error : node size == 0\n");
    return;
  }
  
  //func file
  
  //step loop
  //for(int i=is; i< D->step.size(); i++ ) {
  for(int i=is; i< D->Sc.size(); i++ ) {
    
    //func.layoutファイルオープン
    layoutfile = Generate_LayoutFileName(D->Prefix, D->Sc[i]->step);
    layoutfile = dirname + layoutfile;
    if( (fp = fopen(layoutfile.c_str(), "w")) == NULL ) {
      fprintf(stderr, "Can't open file.(%s)\n", layoutfile.c_str());
      return;
    }
    
    //ファイル書き出し
    fprintf(fp,"FIELDVIEW LAYOUT 1\n");
    for(int j=0; j< D->NodeInfoSize; j++ ) {
      nodefile = Generate_FileName_Free(D->Prefix, "func", D->Sc[i]->step, D->Node[j].RankID, true);
      
      for(int i=0; i<rankis.size(); i++ ) {
        if(rankis[i]<=D->Node[j].RankID && D->Node[j].RankID<=rankie[i]){
          m_name=mname[i];
          d_name=dname[i];
          break;
        }
      }
      fprintf(fp,"%s\n", nodefile.c_str());
      fprintf(fp,"%s\n", m_name.c_str());
      fprintf(fp,"%s\n", d_name.c_str());
    }
    
    //func.layoutファイルクローズ
    fclose(fp);
    
  }
  
}


// #################################################################
// layoutファイル名を作成する
std::string LAYOUT::Generate_LayoutFileName(const std::string prefix, const unsigned m_step)
{
  if ( prefix.empty() ) return NULL;
  
  int len = prefix.size() + 19; // step(10) + postfix(7) + 1(\0) + 1(under score)
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  sprintf(tmp, "%s_%010d.%s", prefix.c_str(), m_step, "layout");
  
  std::string filename(tmp);
  if ( tmp ) delete [] tmp;
  
  return filename;
}

// #################################################################
// ファイル名を作成する（拡張子自由）
std::string LAYOUT::Generate_FileName_Free(const std::string prefix, const std::string xxx, const unsigned m_step, const int m_id, const bool mio)
{
  if ( prefix.empty() ) return NULL;
  
  int len = prefix.size() + xxx.size() + 23; // step(10) + id(9) + 1(.拡張子) + 1(\0) + 2(under score)
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  // local出力が指定された場合、分割出力
  if ( mio )
  {
    sprintf(tmp, "%s_%06d_%010d.%s", prefix.c_str(), m_id, m_step, xxx.c_str());
    //sprintf(tmp, "%s%010d_id%06d.%s", prefix.c_str(), m_step, m_id, xxx.c_str());
  }
  else
  {
    sprintf(tmp, "%s_%010d.%s", prefix.c_str(), m_step, xxx.c_str());
  }
  
  std::string fname(tmp);
  if ( tmp ) delete [] tmp;
  
  return fname;
}


// #################################################################
//
void LAYOUT::CheckDir(string dirstr)
{
  //Hostonly_
  //{
  
#ifndef _WIN32
  
  if( dirstr.size() == 0 ) {
    printf("\toutput current directory\n");
    return;
  }
  
  DIR* dir;
  if( !(dir = opendir(dirstr.c_str())) ) {
    if( errno == ENOENT ) {
      mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
      if( mkdir(dirstr.c_str(), mode) == -1 ) {
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
  
#else // for windows
  
  if( dirstr.size() == 0 ) {
    printf("\toutput current directory\n");
    return;
  }
  
  // check to exist directory
  if (IsDirExsist(dirstr)) {
    // exist directory
    return;
  }
  
  // make directory
  if(!CreateDirectory(dirstr.c_str(), NULL)){
    printf("\tCan't generate directory(%s).\n", dirstr.c_str());
    Exit(0);
  }
  
#endif  // _WIN32
  
  //}
  
  return;
}