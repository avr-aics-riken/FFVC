/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file ParseMat.C
//@brief FlowBase ParseMat class
//@author keno, FSI Team, VCAD, RIKEN

#include "ParseMat.h"

/**
 @fn void ParseMat::chkList(FILE* mp, FILE* fp, CompoList* compo, unsigned basicEq)
 @brief 取得したCompoList[]の内容を表示する
 */
void ParseMat::chkList(FILE* mp, FILE* fp, CompoList* compo, unsigned basicEq)
{
  chkList(mp, compo, basicEq);
  chkList(fp, compo, basicEq);
}

/**
 @fn void ParseMat::chkList(FILE* fp, CompoList* compo, unsigned basicEq)
 @brief 取得したCompoList[]の内容を表示する
 */
void ParseMat::chkList(FILE* fp, CompoList* compo, unsigned basicEq)
{
  if( !fp ) Exit(0);
  
  if ( basicEq == INCMP_2PHASE ) {
    Hostonly_ fprintf(fp,"\t  No :      ID            Element      Medium   Phase                     Label : BCtype\n");
    for (unsigned i=1; i<=NoCompo; i++) {
      
      Hostonly_ fprintf(fp,"\t%4d : %7d %18d ", i, compo[i].getID(), compo[i].getElement());
      ( compo[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
      ( compo[i].getPhase() == GAS )   ? fprintf(fp, "        Gas ") : fprintf(fp, "     Liquid ") ;
      Hostonly_ fprintf(fp, " %24s : %s", (compo[i].name == NULL)?"":compo[i].name, compo[i].getBCstr().c_str() );
    }
  }
  else {
    if ( !isHeatProblem() ) {
      Hostonly_ fprintf(fp,"\t  No :      ID            Element      Medium                    Label : BCtype\n");
      for (unsigned i=1; i<=NoBC; i++) {
        Hostonly_ fprintf(fp,"\t%4d : %7d %18d ", i, compo[i].getID(), compo[i].getElement());
        ( compo[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
        Hostonly_ fprintf(fp, "%24s : %s", (compo[i].name == NULL)?"":compo[i].name, compo[i].getBCstr().c_str() );
        Hostonly_ fprintf(fp,"\n");
      }
      for (unsigned i=NoBC+1; i<=NoCompo; i++) {
        Hostonly_ fprintf(fp,"\t%4d : %7d %18d ", i, compo[i].getID(), compo[i].getElement());
        ( compo[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
        Hostonly_ fprintf(fp, "%24s : %s", (compo[i].name == NULL)?"":compo[i].name, compo[i].getBCstr().c_str() );
        Hostonly_ fprintf(fp,"\n");
      }
    }
    else {
      Hostonly_ fprintf(fp,"\t  No :      ID            Element      Medium    Init.Temp(%s)                    Label : BCtype\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C" );
      for (unsigned i=1; i<=NoBC; i++) {
        Hostonly_ fprintf(fp,"\t%4d : %7d %18d ", i, compo[i].getID(), compo[i].getElement());
        ( compo[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
        Hostonly_ fprintf(fp, "%14s %24s : %s", "-", (compo[i].name == NULL)?"":compo[i].name, compo[i].getBCstr().c_str() );
        Hostonly_ fprintf(fp,"\n");
      }
      for (unsigned i=NoBC+1; i<=NoCompo; i++) {
        Hostonly_ fprintf(fp,"\t%4d : %7d %18d ", i, compo[i].getID(), compo[i].getElement());
        ( compo[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
        Hostonly_ fprintf(fp, "%14.4e %24s : %s", FBUtility::convK2Temp(compo[i].getInitTemp(), Unit_Temp), 
                          (compo[i].name == NULL)?"":compo[i].name, compo[i].getBCstr().c_str() );
        Hostonly_ fprintf(fp,"\n");
      }
    }
  }
  Hostonly_ fprintf(fp,"\n");
}


/**
 @fn bool ParseMat::chkList4Solver(int m)
 @brief 媒質情報の内容物をチェックする
 @param m
 */
bool ParseMat::chkList4Solver(int m)
{
  unsigned c=0;
  if (  BaseMat[m].getState() == FLUID ) {
    if ( !ChkList[p_density] )                c += missingMessage(m, p_density);
    if ( !ChkList[p_kinematic_viscosity] )    c += missingMessage(m, p_kinematic_viscosity);
    if ( !ChkList[p_viscosity] )              c += missingMessage(m, p_viscosity);
    if ( !ChkList[p_thermal_conductivity] )   c += missingMessage(m, p_thermal_conductivity);
    if ( !ChkList[p_specific_heat] )          c += missingMessage(m, p_specific_heat);
    if ( !ChkList[p_sound_of_speed] )         c += missingMessage(m, p_sound_of_speed);
    if ( !ChkList[p_vol_expansion] )          c += missingMessage(m, p_vol_expansion);
  }
  else {  // solid
    if ( !ChkList[p_density] )                c += missingMessage(m, p_density);
    if ( !ChkList[p_specific_heat] )          c += missingMessage(m, p_specific_heat);
    if ( !ChkList[p_thermal_conductivity] )   c += missingMessage(m, p_thermal_conductivity);
  }
  return ( (c != 0) ? false : true ) ;
}


/**
 @fn bool ParseMat::chkStateList(CompoList* compo)
 @brief CompoList[]とMaterialList[]の対応するstateが同じであるかチェックする
 @retval エラーコード
 */
bool ParseMat::chkStateList(CompoList* compo)
{
  unsigned i, odr, id;
  
  for (i=NoBC+1; i<=NoCompo; i++) {
    odr = compo[i].getMatOdr();
    id  = compo[i].getID();
    //printf("\tcompo.odr=%d MatOdr=%d id=%d compo[%d].state=%d mat[%d].state=%d\n", i, odr, id, i, compo[i].getState(), odr, mat[odr].getState() );
    if ( compo[i].getState() != mat[odr].getState() ) {
      Hostonly_ printf("\tState between CompoList[%s] and MaterialList[%s] for cell ID[%d] is conflicting.\n",
             ( compo[i].getState() == FLUID ) ? "Fluid" : "Solid",
             ( mat[odr].getState() == FLUID ) ? "Fluid" : "Solid", id );
      Hostonly_ printf("\tCheck State between 'Model_Setting' and 'Medium_Table'.\n");
      return false;
    }
  }
  return true;
}

/**
 @fn void ParseMat::chkState_Mat_Cmp(CompoList* compo)
 @brief  cmp[]とmat[]のStateの不一致をチェックする
 */
void ParseMat::chkState_Mat_Cmp(CompoList* compo)
{
  unsigned m_id, m_state, m_odr;
  unsigned c=0;
  
  printf("\t  no :   ID  mat_odr mat_ID         mat_Name  : mat_State : specified_state\n");
  
  for (unsigned i=1; i<=NoCompo; i++) {
    m_odr = compo[i].getMatOdr();
    m_id = compo[i].getID();
    m_state = compo[i].getState();
    
    printf("\t%4d : %4d %8d %6d %16s  :     %s : %s", i, m_id, m_odr, mat[m_odr].getMatID(), mat[m_odr].getName(), 
           (mat[m_odr].getState() == FLUID ) ? "Fluid" : "Solid",
           (m_state == FLUID ) ? "Fluid" : "Solid");
    if (mat[m_odr].getState() != m_state) {
      printf(" <\n");
      c++;
    }
    else {
      printf("\n");
    }
  }
  printf("\n");
  
  if ( c>0 ) {
    printf("\tMedium ID described in 'Model_Setting is not consistent with the state of Solid/Fluid. See above lines marked '<'.'\n\n");
    Exit(0);
  }
}


/**
 @fn bool ParseMat::copyMaterials(unsigned& odr, unsigned id)
 @brief MaterialList BaseMat[]からmat[]へ属性をコピーする
 @retval エラーコード
 @param odr 登録エントリ番号
 @param id 登録対象となるmediumID
 @note
 - マテリアル数はNoMaterialで既知（ParseBC::getNoMaterial）
 - 登録対象のIDがmat[]に登録されているかを調べる
 - odrの値は現在までに登録されている個数で，更新して呼び出し元に戻る．初期値には0が入ってくる．
 - 既に登録されているかどうかの判定は，1から現在登録されている値までについて調べる．
 - もし，idが登録できなければ，関数の最後でリターン
 */
bool ParseMat::copyMaterials(unsigned& odr, unsigned id)
{
  int odr_base=0;
  
  // Base MaterialListに登録されているMatIDを調べ，登録エントリ番号を取得する
  if ( -1 == (odr_base=getOdrInMatList(id, NoBaseMat, BaseMat)) ) {
    Hostonly_ stamped_printf("\tParsing error : Medium ID [%d] is not listed in MaterialList\n", id);
    return false;
  }
  
  if ( odr == 0 ) {
    odr++;
    mat[odr].setMatID( BaseMat[odr_base].getMatID() );
    mat[odr].setState( BaseMat[odr_base].getState() );
    strcpy(mat[odr].getName(), BaseMat[odr_base].getName());
    for (int m=0; m<property_END; m++)  mat[odr].P[m] = BaseMat[odr_base].P[m];
    return true;
  }
  else {
    bool flag=true; 
    for (unsigned n=1; n<=odr; n++) { // 1から現在登録されている値までについて調べる
      if ( mat[n].getMatID() == id ) flag = false; // 既に登録されていると，falseにする
    }
    if (flag) { // 登録されていない場合，新規登録して，odrをインクリメントし，関数を抜ける
      odr++;
      mat[odr].setMatID( BaseMat[odr_base].getMatID() );
      mat[odr].setState( BaseMat[odr_base].getState() );
      strcpy(mat[odr].getName(), BaseMat[odr_base].getName());
      for (int m=0; m<property_END; m++)  mat[odr].P[m] = BaseMat[odr_base].P[m];
      return true;
    }
  }
  return true;
}

/**
 @fn void ParseMat::dbg_printRelation(FILE* mp, FILE* fp, CompoList* compo)
 @brief CompoList[]とMaterialList[]のチェック
 */
void ParseMat::dbg_printRelation(FILE* mp, FILE* fp, CompoList* compo)
{
  printRelation(mp, compo);
  printRelation(fp, compo);
}

/**
 @fn unsigned ParseMat::getMatIDinTable(unsigned id)
 @brief iTable[]にidが含まれるかどうかを調べる
 @retval 含まれれば対応するMatID，そうでなければ0を返す
 */
unsigned ParseMat::getMatIDinTable(unsigned id)
{
  for (unsigned i=1; i<=NoID; i++) {
    if ( iTable[i].getID() == id ) return iTable[i].getMatID();
  }
  return 0;
}


/**
 @fn unsigned ParseMat::getOdrInMatList(unsigned id, unsigned Msize, MaterialList* mlist)
 @brief MaterialListに含まれるidのエントリ番号を返す
 @retval 含まれれば，エントリ番号，そうでなければ0を返す
 */
unsigned ParseMat::getOdrInMatList(unsigned id, unsigned Msize, MaterialList* mlist)
{
  for (unsigned i=1; i<=Msize; i++) {
    if ( mlist[i].getMatID() == id ) return i;
  }
  return 0;
}


/**
 @fn void ParseMat::getPvalue(const CfgParam* p, REAL_TYPE &value)
 @brief paramの値を返す
 @param p param要素
 @param value 値
 */
void ParseMat::getPvalue(const CfgParam* p, REAL_TYPE &value)
{
  REAL_TYPE f;
  
  if ( !(p->GetData( &f )) ) {
    Hostonly_ stamped_printf("\tParsing error : Invalid float value for in Medium\n");
    Exit(0);
  }
  value = f;
}


/**
 @fn void ParseMat::getXMLmaterial(void)
 @brief Material情報の内容をXMLファイルをパースして，MaterialListクラスのオブジェクトBaseMatに保持する
 */
void ParseMat::getXMLmaterial(void)
{
  int id;
  const CfgElem *elmL1=NULL, *elmL2=NULL;
  const char *p=NULL, *pnt=NULL;
  
  // Check Model_Setting section
  if ( !(elmL1 = CF->GetTop(MDMTBL)) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing Medium_Table tree\n");
    Exit(0);
  }
  
  // check # of Elem
  if ( (NoBaseMat=elmL1->GetElemSize()) == 0 ) {
    Hostonly_ printf("\tNo description was found in 'Medium_Table'\n");
    Exit(0);
  }
  
  // Instance MaterialList valid only inside of this class
  BaseMat = new MaterialList[NoBaseMat+1];
  
  // load Base list
  elmL2 = elmL1->GetElemFirst();
  for (int i=1; i<=NoBaseMat; i++) {
    p = elmL2->GetName();
    if (  !strcasecmp(p, "Solid") || ( !strcasecmp(p, "Fluid") ) ) {
      
      // ID
      if ( !elmL2->isSetID() ) {
        Hostonly_ printf("\tParsing error : No ID section for Medium in 'Medium_Table'\n");
        Exit(0);
      }
      if ( -1 == (id=elmL2->GetID()) ) {
        Hostonly_ printf("\tParsing error : No valid ID for Medium in 'Medium_Table'\n");
        Exit(0);
      }
      BaseMat[i].setMatID( (unsigned)id );
      
      // state
      if      ( !strcasecmp(p, "Fluid") ) BaseMat[i].setState(FLUID);
      else if ( !strcasecmp(p, "Solid") ) BaseMat[i].setState(SOLID);
      else {
        Hostonly_ printf("\tInvalid medium keyword 'solid/fluid'\n");
        Exit(0);
      }
      
      // Material name
      pnt = NULL;
      if ( !(pnt = elmL2->GetComment()) ) {
        Hostonly_ printf("\tNo comment for Medium\n");
      }
      else strcpy(BaseMat[i].getName(), pnt);
      
      // each properties
      readMedium(elmL2, i);
    }
    else {
      Hostonly_ printf("\tInvalid keyword in Medium_Table : [%s]\n", p);
      Exit(0);
    }
    
    elmL2 = elmL1->GetElemNext(elmL2);
  }
  
  // 重複をチェック
  unsigned md;
  for (unsigned n=2; n<=NoBaseMat; n++) {
    md = BaseMat[n].getMatID();
    for (unsigned l=1; l<n; l++) {
      if ( BaseMat[l].getMatID() == md ) {
        Hostonly_ stamped_printf("\tDuplicate Medium id[%d]=%d in Base MaterialList\n", l, md);
        FILE* mp=stdout;
        
        Hostonly_ fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
        Hostonly_ fprintf(mp,"\n\t>> Base Material List\n\n");
        
        Hostonly_ printMatList(mp, NoBaseMat, BaseMat);
        Exit(0);
      }
    }
  }
}


/**
 @fn bool ParseMat::isIDinList(unsigned RefID)
 @brief BaseMatにRefIDがあるかを調べる
 @retval RefIDがあれば，trueを返す
 */
bool ParseMat::isIDinList(unsigned RefID)
{
  for (unsigned i=1; i<=NoBaseMat; i++) {
    if (BaseMat[i].getMatID() == RefID) return true;
  }
  return false;
}


/**
 @fn void ParseMat::makeLinkCmpMat(CompoList* compo)
 @brief CompoList, MaterialListの相互参照の準備を行う
 */
void ParseMat::makeLinkCmpMat(CompoList* compo)
{
  unsigned id, id_mat, odr_mat;
  
  // コンポーネントリストへ，マテリアルリストへのエントリ番号をエンコードする
  for (int i=1; i<=NoCompo; i++) {
    id = compo[i].getID();
    
    // IDに対応したMatIDをiTableから取得する
    if ( 0 == (id_mat = getMatIDinTable(id)) ) {
      Hostonly_ stamped_printf("\tSomething wrong! Medium ID = %d\n", id_mat);
      Exit(0);
    }
    
    // MatIDのMaterialListのエントリ番号
    odr_mat = getOdrInMatList(id_mat, NoMaterial, mat);

    // copy property
    compo[i].setMatOdr( odr_mat );
  }
}


/**
 @fn void ParseMat::makeMaterialList(void)
 @brief BaseListからMaterialListを作成する
 @note
 - ParseBC::getNoMaterial()で，iTableに登録されたMaterial数を取得ずみ
 */
void ParseMat::makeMaterialList(void)
{
  unsigned odr=0, matid;
  
  for (int i=1; i<=NoID; i++) {
    matid = iTable[i].getMatID();
    if ( !isIDinList(matid) ) {
      Hostonly_ stamped_printf("\tMedium ID [%d] described in 'Model_Setting' can not find in 'Medium_Table'\n", matid);
      Exit(0);
    }
    if ( !copyMaterials( odr, matid ) ) Exit(0);
    if ( odr == NoMaterial ) return;
  }
  stamped_printf("\tError\n");
  Exit(0);
}


/**
 @fn unsigned ParseMat::missingMessage(int m, unsigned key)
 @brief 警告メッセージの表示
 @param m
 @param key
 */
unsigned ParseMat::missingMessage(int m, unsigned key)
{
  Hostonly_ printf("\tMissing keyword '%s' for '%s' in %s phase\n", 
         MaterialList::getPropertyName(key).c_str(), BaseMat[m].getName(),
         ( BaseMat[m].getState() == SOLID ) ? "solid" : "fluid" );
  return 1; 
}


/**
 @fn void ParseMat::dbg_printBaseMaterialList(FILE* mp, FILE* fp)
 @brief 基本媒質情報の表示
 */
void ParseMat::dbg_printBaseMaterialList(FILE* mp, FILE* fp)
{

  Hostonly_ fprintf(mp, "DEBUG : \n");
  Hostonly_ fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
  Hostonly_ fprintf(mp,"\n\t>> Base Material List\n\n");
  Hostonly_ fprintf(fp, "DEBUG : \n");
  Hostonly_ fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  Hostonly_ fprintf(fp,"\n\t>> Base Material List\n\n");

  printMatList(mp, NoBaseMat, BaseMat);
  printMatList(fp, NoBaseMat, BaseMat);
}

/**
 @fn void ParseMat::printMaterialList(FILE* mp, FILE* fp)
 @brief MaterialListを表示する
 */
void ParseMat::printMaterialList(FILE* mp, FILE* fp)
{
  printMatList(mp, NoMaterial, mat);
  printMatList(fp, NoMaterial, mat);
}


/**
 @fn void ParseMat::printMatList(FILE* fp, unsigned Max, MaterialList* mlist)
 @brief 媒質情報の表示
 @param fp
 @param Max
 @param mlist
 */
void ParseMat::printMatList(FILE* fp, unsigned Max, MaterialList* mlist)
{  
  unsigned n, id;
  
  fprintf(fp, "\t  no :      ID           Medium\n");
  
  for(n=1; n<=Max; n++) {
    id = mlist[n].getMatID();
    fprintf(fp,"\t%4d : %7d %16s\n", n, id, mlist[n].getName());
    
    if ( mlist[n].getState() == FLUID ) {
      fprintf(fp, "\t\t\t\t\t density              %12.6e [kg/m^3]\n",   mlist[n].P[p_density]);
      fprintf(fp, "\t\t\t\t\t kinematic viscosity  %12.6e [m^2/s]\n",    mlist[n].P[p_kinematic_viscosity]);
      fprintf(fp, "\t\t\t\t\t viscosity            %12.6e [Pa s]\n",     mlist[n].P[p_viscosity]);
      fprintf(fp, "\t\t\t\t\t specific heat        %12.6e [J/(Kg K)]\n", mlist[n].P[p_specific_heat]);
      fprintf(fp, "\t\t\t\t\t thermal conductivity %12.6e [W/(m K)]\n",  mlist[n].P[p_thermal_conductivity]);
      fprintf(fp, "\t\t\t\t\t sound of speed       %12.6e [m/s]\n",      mlist[n].P[p_sound_of_speed]);
      fprintf(fp, "\t\t\t\t\t volume expansion     %12.6e [1/K]\n",      mlist[n].P[p_vol_expansion]);
    }
    else {  // solid
      fprintf(fp, "\t\t\t\t\t density              %12.6e [kg/m^3]\n",   mlist[n].P[p_density]);
      fprintf(fp, "\t\t\t\t\t specific heat        %12.6e [J/(Kg K)]\n", mlist[n].P[p_specific_heat]);
      fprintf(fp, "\t\t\t\t\t thermal conductivity %12.6e [W/(m K)]\n",  mlist[n].P[p_thermal_conductivity]);
    }
    
  }
  fprintf(fp,"\n");
  fflush(fp);
}

/**
 @fn void ParseMat::printRelation(FILE* fp, CompoList* compo)
 @brief CompoList[]とMaterialList[]のチェック
 @note debug function
 */
void ParseMat::printRelation(FILE* fp, CompoList* compo)
{  
  unsigned odr;
  
  fprintf(fp,"\n");
  fprintf(fp, "DEBUG @ ParseMat::printRelation\n\n");
  fprintf(fp, "\t  Order in CompoList :  Cell ID  Medium order     ID              Name\n");
  
  for (unsigned i=1; i<=NoCompo; i++) {
    if ( compo[i].getState() != -1 ) {  // is Medium
      odr = compo[i].getMatOdr();
      fprintf(fp,"\t\t\t%4d : %8d            %2d  %5d  %16s  : %s\n", i, compo[i].getID(), odr, mat[odr].getMatID(), mat[odr].getName(),
              (mat[odr].getState() == FLUID ) ? "Fluid" : "Solid" );
    }
    else {
      stamped_printf("\terror\n");
    }
  }
  
  fprintf(fp,"\n");
  fflush(fp);
}


/**
 @fn void ParseMat::readMedium(const CfgElem *elmL2, int m)
 @brief 媒質情報を取得する
 @param elmL2
 @param m
 */
void ParseMat::readMedium(const CfgElem *elmL2, int m)
{
  const CfgParam *param=NULL;
  const char* Pname=NULL;
  int NoParam=0;
  int key;
  
  NoParam = elmL2->GetParamSize();
  
  // clear for each medium
  for (int i=0; i<property_END; i++) ChkList[i] = false;
  
  // check for solver
  param   = elmL2->GetParamFirst();
  for (int i=0; i<NoParam; i++) {
    Pname = param->GetName();
    if ( -1 == (key = MaterialList::getKey(Pname)) ) {
      Hostonly_ printf("\tParameter error : Invalid keyword in 'Medium_Table' : [%s]\n", Pname);
      Exit(0);
    }
    ChkList[key] = true;
    
    param = elmL2->GetParamNext(param);
  }
  
  // Check
  if ( !chkList4Solver(m) ) Exit(0);
  
  // load value
  param   = elmL2->GetParamFirst();
  for (int i=0; i<NoParam; i++) {
    Pname = param->GetName();
    key = MaterialList::getKey(Pname);  // return value was already checked
    
    getPvalue(param, BaseMat[m].P[key]);
    
    param = elmL2->GetParamNext(param);
  }
}


/**
 @fn void ParseMat::setControlVars(Control* Cref, IDtable* itbl, MaterialList* m_mat, SklSolverConfig* cfg)
 @brief 媒質情報のハンドリング
 @param Cref 
 @param itbl 
 @param m_mat
 @param cfg 
 */
void ParseMat::setControlVars(Control* Cref, IDtable* itbl, MaterialList* m_mat, SklSolverConfig* cfg)
{
  if ( !itbl ) Exit(0);
  if ( !cfg )  Exit(0);
  
  NoCompo        = Cref->NoCompo;
  NoBC           = Cref->NoBC;
  NoID           = Cref->NoID;
  NoMaterial     = Cref->NoMaterial;
  guide          = Cref->guide;
  imax = size[0] = Cref->imax;
  jmax = size[1] = Cref->jmax;
  kmax = size[2] = Cref->kmax;
  Unit_Temp      = Cref->Unit.Temp;
  BaseTemp       = Cref->BaseTemp;
  DiffTemp       = Cref->DiffTemp;
  KOS            = Cref->KindOfSolver;
  iTable = itbl;
  CF     = cfg;
  
  mat = m_mat;
}
