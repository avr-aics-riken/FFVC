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
    if ( KOS != FLOW_ONLY ) {
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
 @fn bool ParseMat::chkStateList(CompoList* compo)
 @brief CompoList[]とMediumList[]の対応するstateが同じであるかチェックする
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
      Hostonly_ printf("\tState between CompoList[%s] and MediumList[%s] for cell ID[%d] is conflicting.\n",
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
void ParseMat::chkState_Mat_Cmp(CompoList* compo, FILE* fp)
{
  unsigned m_id, m_state, m_odr;
  unsigned c=0;
  
  printf("\t  no :   ID  mat_odr mat_ID         mat_Name  : mat_State : specified_state\n");
  fprintf(fp, "\t  no :   ID  mat_odr mat_ID         mat_Name  : mat_State : specified_state\n");
  
  for (unsigned i=1; i<=NoCompo; i++) {
    m_odr = compo[i].getMatOdr();
    m_id = compo[i].getID();
    m_state = compo[i].getState();
    
    printf("\t%4d : %4d %8d %6d %16s  :     %s : %s", i, m_id, m_odr, m_odr, mat[m_odr].getLabel(), 
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
 @fn void ParseMat::dbg_printRelation(FILE* mp, FILE* fp, CompoList* compo)
 @brief CompoList[]とMediumList[]のチェック
 */
void ParseMat::dbg_printRelation(FILE* mp, FILE* fp, CompoList* compo)
{
  printRelation(mp, compo);
  printRelation(fp, compo);
}



/**
 @fn void ParseMat::makeLinkCmpMat(CompoList* compo)
 @brief CompoList, MediumListの相互参照の準備を行う
 
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

    // copy property
    compo[i].setMatOdr( i );
  }
}*/



/**
 @fn void ParseMat::printMediumList(FILE* mp, FILE* fp)
 @brief MediumListを表示する
 */
void ParseMat::printMediumList(FILE* mp, FILE* fp)
{
  printMatList(mp, NoMedium, mat);
  printMatList(fp, NoMedium, mat);
}


/**
 @fn void ParseMat::printMatList(FILE* fp, int Max, MediumList* mlist)
 @brief 媒質情報の表示
 @param fp
 @param Max
 @param mlist
 */
void ParseMat::printMatList(FILE* fp, int Max, MediumList* mlist)
{
  fprintf(fp, "\t  no :           Medium  Properties\n");
  
  for (int n=1; n<=Max; n++) {
    fprintf(fp,"\t%4d : %16s\n", n, mlist[n].getLabel());
    
    if ( mlist[n].getState() == FLUID ) {
      fprintf(fp, "\t\t\t\t density              %12.6e [kg/m^3]\n",   mlist[n].P[p_density]);
      fprintf(fp, "\t\t\t\t kinematic viscosity  %12.6e [m^2/s]\n",    mlist[n].P[p_kinematic_viscosity]);
      fprintf(fp, "\t\t\t\t viscosity            %12.6e [Pa s]\n",     mlist[n].P[p_viscosity]);
      fprintf(fp, "\t\t\t\t specific heat        %12.6e [J/(Kg K)]\n", mlist[n].P[p_specific_heat]);
      fprintf(fp, "\t\t\t\t thermal conductivity %12.6e [W/(m K)]\n",  mlist[n].P[p_thermal_conductivity]);
      fprintf(fp, "\t\t\t\t sound of speed       %12.6e [m/s]\n",      mlist[n].P[p_sound_of_speed]);
      fprintf(fp, "\t\t\t\t volume expansion     %12.6e [1/K]\n",      mlist[n].P[p_vol_expansion]);
    }
    else {  // solid
      fprintf(fp, "\t\t\t\t density              %12.6e [kg/m^3]\n",   mlist[n].P[p_density]);
      fprintf(fp, "\t\t\t\t specific heat        %12.6e [J/(Kg K)]\n", mlist[n].P[p_specific_heat]);
      fprintf(fp, "\t\t\t\t thermal conductivity %12.6e [W/(m K)]\n",  mlist[n].P[p_thermal_conductivity]);
    }
    
  }
  fprintf(fp,"\n");
  fflush(fp);
}

/**
 @fn void ParseMat::printRelation(FILE* fp, CompoList* compo)
 @brief CompoList[]とMediumList[]のチェック
 @note debug function
 */
void ParseMat::printRelation(FILE* fp, CompoList* compo)
{  
  int odr;
  
  fprintf(fp,"\n");
  fprintf(fp, "DEBUG @ ParseMat::printRelation\n\n");
  fprintf(fp, "\t  Order in CompoList :  Cell ID  Medium order     ID              Name\n");
  
  for (int i=1; i<=(int)NoCompo; i++) {
    if ( compo[i].getState() != -1 ) {  // is Medium
      odr = compo[i].getMatOdr();
      fprintf(fp,"\t\t\t%4d : %8d            %2d  %5d  %16s  : %s\n", i, compo[i].getID(), odr, odr, mat[odr].getLabel(),
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
 @fn void ParseMat::setControlVars(Control* Cref, MediumList* m_mat)
 @brief 媒質情報のハンドリング
 @param Cref 
 @param m_mat 
 */
void ParseMat::setControlVars(Control* Cref)
{
  NoCompo        = Cref->NoCompo;
  NoBC           = Cref->NoBC;
  Unit_Temp      = Cref->Unit.Temp;
  KOS            = Cref->KindOfSolver;
}


/**
 @fn void ParseMat::receive_TP_Ptr(TPControl* tp)
 @brief TPのポインタを受け取る
 */
bool ParseMat::receive_TP_Ptr(TPControl* tp) 
{ 
  if ( !tp ) return false;
  tpCntl = tp;
  return true;
}



/**
 @fn int ParseMat::get_MediumTable(void)
 @brief Medium_Tableを読んでMediumTableInfoクラスオブジェクトに格納
 */
int ParseMat::get_MediumTable(void)
{
  std::string str, label;
  std::string label_base, label_m, label_leaf;
  REAL_TYPE fval;
  int n1=0, n2=0, n3=0;
  
  NoMedium = 0;
  
  label_base = "/Medium_Table";
  
  // 媒質の個数を取得
  n1 = tpCntl->countLabels(label_base);
  if ( n1 < 0) {
    stamped_printf("\tcountLabels --- %s\n", label_base.c_str());
    Exit(0);
  }
  NoMedium = n1;
  
  MTITP = new MediumTableInfo[NoMedium+1];
  
  for (int i1=1; i1<=NoMedium; i1++) {
    if ( !tpCntl->GetNodeStr(label_base, i1, &str) ) {
      stamped_printf("\tParsing error : No Leaf Node \n");
      Exit(0);
    }
    
    label_m = label_base + "/" + str;
    n2 = tpCntl->countLabels(label_m);
    
    for (int i2=1; i2<=n2; i2++) {
      if ( !tpCntl->GetNodeStr(label_m, i2, &str) ) {
        stamped_printf("\tParsing error : No Leaf Node in Medium_Table[%d]\n", i1);
        Exit(0);
      }
      
      label_leaf = label_m + "/" + str;
      
      if ( !strcasecmp(str.c_str(), "type") ) {
        if ( !(tpCntl->GetValue(label_leaf, &label)) ){
          ;
        }
        else {
          if     ( !strcasecmp(label.c_str(), "Fluid") ) MTITP[i1].type = FLUID;
          else if( !strcasecmp(label.c_str(), "Solid") ) MTITP[i1].type = SOLID;
          else {
            stamped_printf("\tParsing error : unknown type at Medium_Table \n");
            Exit(0);
          }
        }
      }
      else if( !strcasecmp(str.c_str(), "label") ){
        if ( !(tpCntl->GetValue(label_leaf, &label)) ){
          ;
        }
        else {
          MTITP[i1].label = label;
        }
      }
      else {
        if ( !(tpCntl->GetValue(label_leaf, &fval)) ){
          ;
        }
        else{
          // データをmapに追加
          MTITP[i1].m_fval.insert( std::pair<std::string, REAL_TYPE>(str, fval) );
        }
      }
    }
  }
  
  return NoMedium; 
}


/**
 @fn void ParseMat::makeMediumList(MediumList* m_mat)
 @brief MediumListを作成する
 */
void ParseMat::makeMediumList(MediumList* m_mat)
{
  
  if ( !m_mat ) return;
  
  mat = m_mat;
  
  int type;
  std::string label;
  
  for (int i=1; i<=NoMedium; i++) {
    
    type  = MTITP[i].type;
    label = MTITP[i].label;
    
    // すでに登録されているかどうか調べる
    if ( !chkDuplicateLabel(i, label) ) {
      Hostonly_ printf("\tError : Medium label '%s' is already used\n", label.c_str());
      Exit(0);
    }
    
    // state
    mat[i].setState(type);

    // Medium name
    strcpy( mat[i].getLabel(), label.c_str() );
    
    // set medium value
    copyProperty(i);
  }

}



//@brief ラベルの重複を調べる
bool ParseMat::chkDuplicateLabel(const int n, std::string m_label)
{
	for (int i=0; i<n; i++){
    if ( mat[i].getLabel() == m_label.c_str() ) return false;
	}
	return true;
}


/**
 @fn void ParseMat::copyProperty(const int n)
 @brief matの変数値を格納する
 */
void ParseMat::copyProperty(const int n)
{
	int nfval = MTITP[n].m_fval.size();
  
	if ( nfval > property_END ){
    Hostonly_ printf("\tParameter error : too big property size 'Medium_Table'\n");
    Exit(0);
	};
  
  // clear for each medium
  for (int i=0; i<property_END; i++) ChkList[i] = false;
	
	// イテレータを生成
  std::map<string, REAL_TYPE>::iterator itr;
  
	for (itr = MTITP[n].m_fval.begin(); itr != MTITP[n].m_fval.end(); itr++)
	{
    std::string a1 = itr->first;   // キー取得
		REAL_TYPE   a2 = itr->second;  // 値取得
    
    int key = MediumList::getKey(a1.c_str());
    
		mat[n].P[key]=a2;
    ChkList[key] = true;
	}
	
	// Check
	if ( !chkList4Solver(n) ){
		printf("\tParameter error : chkList4Solver \n");
		Exit(0);
	}
  
	return;
}



/**
 @fn bool ParseMat::chkList4Solver(int m)
 @brief 媒質情報の内容物をチェックする
 @param m
 */
bool ParseMat::chkList4Solver(const int m)
{
  int c=0;
  if ( mat[m].getState() == FLUID ) {
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
 @fn int ParseMat::missingMessage(const int m, const int key)
 @brief 警告メッセージの表示
 @param m
 @param key
 */
int ParseMat::missingMessage(const int m, const int key)
{
  printf("\tMissing keyword '%s' for '%s' in %s phase\n", 
         MediumList::getPropertyName(key).c_str(), mat[m].getLabel(),
         ( mat[m].getState() == SOLID ) ? "solid" : "fluid" );
  return 1; 
}
