//##################################################################################
//
// Flow Base class
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
 * @file   ParseMat.C
 * @brief  FlowBase ParseMat class
 * @author kero
 */

#include "ParseMat.h"


// #################################################################
// ラベルの重複を調べる
bool ParseMat::chkDuplicateLabel(MediumList* mat, const int n, const std::string m_label)
{
	for (int i=0; i<n; i++){
    if ( mat[i].getAlias() == m_label ) return false;
	}
	return true;
}



// #################################################################
// 取得したCompoList[]の内容を表示する
void ParseMat::chkList(FILE* fp, CompoList* compo, const int basicEq)
{
  if( !fp ) Exit(0);
  
  if ( basicEq == INCMP_2PHASE ) {
    fprintf(fp,"\t  No :     Mat            Element      Medium   Phase                     Label : BCtype\n");
    for (int i=1; i<=NoCompo; i++) {
      
      fprintf(fp,"\t%4d : %7d %18ld ", i, compo[i].getMatOdr(), compo[i].getElement());
      ( compo[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
      ( compo[i].getPhase() == GAS )   ? fprintf(fp, "        Gas ") : fprintf(fp, "     Liquid ") ;
      fprintf(fp, " %24s : %s", (compo[i].getLabel().empty()) ? "" : compo[i].getLabel().c_str(), compo[i].getBCstr().c_str() );
    }
  }
  else {
    if ( KOS == FLOW_ONLY ) {
      fprintf(fp,"\t  No :     Mat            Element      Medium                    Label : BCtype\n");

      for (int i=1; i<=NoBC; i++) {
        fprintf(fp,"\t%4d : %7d %18ld ", i, compo[i].getMatOdr(), compo[i].getElement());
        ( compo[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
        fprintf(fp, "%24s : %s", (compo[i].getLabel().empty()) ? "" : compo[i].getLabel().c_str(), compo[i].getBCstr().c_str() );
        fprintf(fp,"\n");
      }

      for (int i=NoBC+1; i<=NoCompo; i++) {
        fprintf(fp,"\t%4d : %7d %18ld ", i, compo[i].getMatOdr(), compo[i].getElement());
        ( compo[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
        fprintf(fp, "%24s : %s", (compo[i].getLabel().empty()) ? "" : compo[i].getLabel().c_str(), compo[i].getBCstr().c_str() );
        fprintf(fp,"\n");
      }
    }
    else {
      fprintf(fp,"\t  No :     Mat            Element      Medium    Init.Temp(%s)                    Label : BCtype\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C" );
      for (int i=1; i<=NoBC; i++) {
        fprintf(fp,"\t%4d : %7d %18ld ", i, compo[i].getMatOdr(), compo[i].getElement());
        ( compo[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
        fprintf(fp, "%14s %24s : %s", "-", (compo[i].getLabel().empty()) ? "" : compo[i].getLabel().c_str(), compo[i].getBCstr().c_str() );
        fprintf(fp,"\n");
      }
      for (int i=NoBC+1; i<=NoCompo; i++) {
        fprintf(fp,"\t%4d : %7d %18ld ", i, compo[i].getMatOdr(), compo[i].getElement());
        ( compo[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
        fprintf(fp, "%14.4e %24s : %s", FBUtility::convK2Temp(compo[i].getInitTemp(), Unit_Temp), 
                          (compo[i].getLabel().empty()) ? "" : compo[i].getLabel().c_str(), compo[i].getBCstr().c_str() );
        fprintf(fp,"\n");
      }
    }
  }
  fprintf(fp,"\n");
}


// #################################################################
// 媒質情報の内容物をチェックする
bool ParseMat::chkList4Solver(MediumList* mat, const int m)
{
  int c=0;
  if ( mat[m].getState() == FLUID ) {
    if ( !ChkList[p_density] )                c += missingMessage(mat, m, p_density);
    if ( !ChkList[p_kinematic_viscosity] )    c += missingMessage(mat, m, p_kinematic_viscosity);
    if ( !ChkList[p_viscosity] )              c += missingMessage(mat, m, p_viscosity);
    if ( !ChkList[p_thermal_conductivity] )   c += missingMessage(mat, m, p_thermal_conductivity);
    if ( !ChkList[p_specific_heat] )          c += missingMessage(mat, m, p_specific_heat);
    if ( !ChkList[p_speed_of_sound] )         c += missingMessage(mat, m, p_speed_of_sound);
    if ( !ChkList[p_vol_expansion] )          c += missingMessage(mat, m, p_vol_expansion);
  }
  else {  // solid
    if ( !ChkList[p_density] )                c += missingMessage(mat, m, p_density);
    if ( !ChkList[p_specific_heat] )          c += missingMessage(mat, m, p_specific_heat);
    if ( !ChkList[p_thermal_conductivity] )   c += missingMessage(mat, m, p_thermal_conductivity);
  }
  return ( (c != 0) ? false : true ) ;
}



// #################################################################
// matの変数値を格納する
void ParseMat::copyProperty(MediumList* mat, const int n)
{
	int nfval = MTITP[n].m_fval.size();

	if ( nfval > property_END ) Exit(0);
  
  // clear for each medium
  for (int i=0; i<property_END; i++) ChkList[i] = false;
	
	// イテレータを生成
  std::map<std::string, REAL_TYPE>::iterator itr;
  
	for (itr = MTITP[n].m_fval.begin(); itr != MTITP[n].m_fval.end(); itr++)
	{
    std::string a1 = itr->first;   // キー取得
		REAL_TYPE   a2 = itr->second;  // 値取得
    
    int key = MediumList::getKey(a1.c_str());
    if (key<0) {
      printf("Invalid keyword [%s]\n", a1.c_str());
    }

		mat[n].P[key]=a2;
    ChkList[key] = true;
	}
	
	// Check
	if ( !chkList4Solver(mat, n) ) Exit(0);
}



// #################################################################
// MediumTableを読んでMediumTableInfoクラスオブジェクトに格納
int ParseMat::get_MediumTable()
{
  std::string str, label;
  std::string label_base, label_m, label_leaf;
  REAL_TYPE fval;
  int n1=0, n2=0, n3=0;
  
  int NoMedium = 0;
  
  label_base = "/MediumTable";
  
  // 媒質の個数を取得
  n1 = tpCntl->countLabels(label_base);
  if ( n1 < 0) Exit(0);

  NoMedium = n1;
  
  MTITP = new MediumTableInfo[NoMedium+1];
  
  for (int i1=1; i1<=NoMedium; i1++) {
    
    if ( !tpCntl->GetNodeStr(label_base, i1, &str) )
    {
      Exit(0);
    }
    
    label_m = label_base + "/" + str;
    n2 = tpCntl->countLabels(label_m);
    
    for (int i2=1; i2<=n2; i2++) {
      
      if ( !tpCntl->GetNodeStr(label_m, i2, &str) )
      {
        Exit(0);
      }
      
      label_leaf = label_m + "/" + str;
      
      if ( !strcasecmp(str.c_str(), "Type") )
      {
        if ( !(tpCntl->GetValue(label_leaf, &label)) )
        {
          ;
        }
        else
        {
          if     ( !strcasecmp(label.c_str(), "Fluid") ) MTITP[i1].type = FLUID;
          else if( !strcasecmp(label.c_str(), "Solid") ) MTITP[i1].type = SOLID;
          else
          {
            Exit(0);
          }
        }
      }
      else if( !strcasecmp(str.c_str(), "alias") )
      {
        if ( !(tpCntl->GetValue(label_leaf, &label)) )
        {
          ;
        }
        else
        {
          MTITP[i1].label = label;
        }
      }
      else if( !strcasecmp(str.c_str(), "fxgencolor") )
      {
        ; // FXgenが吐き出す情報でffvcでは不使用
      }
      else if( !strcasecmp(str.c_str(), "fxgenid") )
      {
        ; // FXgenが吐き出す情報でffvcでは不使用
      }
      else // その他のパラメータは浮動小数点データとしてロード
      {
        if ( !(tpCntl->GetValue(label_leaf, &fval)) )
        {
          ;
        }
        else
        {
          // データをmapに追加
          MTITP[i1].m_fval.insert( std::pair<std::string, REAL_TYPE>(str, fval) );
        }
      }
    }
  }
  
  return NoMedium; 
}


// #################################################################
// MediumListを作成する
bool ParseMat::makeMediumList(MediumList* mat, const int NoMedium)
{
  int type;
  std::string label;
  
  for (int i=1; i<=NoMedium; i++) {
    
    type  = MTITP[i].type;
    label = MTITP[i].label;

    // すでに登録されているかどうか調べる
    if ( !chkDuplicateLabel(mat, i, label) ) {
      return false;
    }

    // state
    mat[i].setState(type);

    // Medium name
    mat[i].setAlias(label);

    // set medium value
    copyProperty(mat, i);

  }
  
  return true;
}



// #################################################################
// 警告メッセージの表示
int ParseMat::missingMessage(MediumList* mat, const int m, const int key)
{
  printf("\tMissing keyword '%s' for '%s' in %s phase\n", 
         MediumList::getPropertyName(key).c_str(), mat[m].getAlias().c_str(),
         ( mat[m].getState() == SOLID ) ? "solid" : "fluid" );
  return 1; 
}



// #################################################################
// 媒質情報の表示
// Hostonly
void ParseMat::printMatList(FILE* fp, MediumList* mat, const int NoMedium)
{
  fprintf(fp, "\t  no :           Medium  Properties\n");
  
  for (int n=1; n<=NoMedium; n++) {
    fprintf(fp,"\t%4d : %16s\n", n, mat[n].getAlias().c_str());
    
    if ( mat[n].getState() == FLUID ) {
      fprintf(fp, "\t\t\t\t mass density         %12.6e [kg/m^3]\n",   mat[n].P[p_density]);
      fprintf(fp, "\t\t\t\t kinematic viscosity  %12.6e [m^2/s]\n",    mat[n].P[p_kinematic_viscosity]);
      fprintf(fp, "\t\t\t\t viscosity            %12.6e [Pa s]\n",     mat[n].P[p_viscosity]);
      fprintf(fp, "\t\t\t\t specific heat        %12.6e [J/(Kg K)]\n", mat[n].P[p_specific_heat]);
      fprintf(fp, "\t\t\t\t thermal conductivity %12.6e [W/(m K)]\n",  mat[n].P[p_thermal_conductivity]);
      fprintf(fp, "\t\t\t\t speed of sound       %12.6e [m/s]\n",      mat[n].P[p_speed_of_sound]);
      fprintf(fp, "\t\t\t\t volume expansion     %12.6e [1/K]\n",      mat[n].P[p_vol_expansion]);
    }
    else {  // solid
      fprintf(fp, "\t\t\t\t mass density         %12.6e [kg/m^3]\n",   mat[n].P[p_density]);
      fprintf(fp, "\t\t\t\t specific heat        %12.6e [J/(Kg K)]\n", mat[n].P[p_specific_heat]);
      fprintf(fp, "\t\t\t\t thermal conductivity %12.6e [W/(m K)]\n",  mat[n].P[p_thermal_conductivity]);
    }
    
  }
  fprintf(fp,"\n");
  fflush(fp);
}


// #################################################################
// CompoList[]とMediumList[]のチェック
void ParseMat::printRelation(FILE* fp, CompoList* compo, MediumList* mat)
{  
  int odr;
  
  fprintf(fp,"\n");
  fprintf(fp, "DEBUG @ ParseMat::printRelation\n\n");
  fprintf(fp, "\t  Order in CompoList :  Cell ID  Medium order     ID              Name\n");
  
  for (int i=1; i<=NoCompo; i++) {
    if ( compo[i].getState() != -1 ) {  // is Medium
      odr = compo[i].getMatOdr();
      fprintf(fp,"\t\t\t%4d : %8d            %2d  %5d  %16s  : %s\n", i, compo[i].getMatOdr(), odr, odr, mat[odr].getAlias().c_str(),
              (mat[odr].getState() == FLUID ) ? "Fluid" : "Solid" );
    }
    else {
      stamped_printf("\terror\n");
    }
  }
  
  fprintf(fp,"\n");
  fflush(fp);
}


// #################################################################
//TPのポインタを受け取る
bool ParseMat::importTP(TPControl* tp) 
{ 
  if ( !tp ) return false;
  tpCntl = tp;
  
  return true;
}


// #################################################################
// 媒質情報の初期化
void ParseMat::setControlVars(const int m_NoCompo,
                              const int m_NoBC,
                              const int m_Unit_Temp,
                              const int m_KOS)
{
  NoCompo        = m_NoCompo;
  NoBC           = m_NoBC;
  Unit_Temp      = m_Unit_Temp;
  KOS            = m_KOS;
}
