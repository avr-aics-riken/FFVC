//##################################################################################
//
// Flow Base class
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
//##################################################################################

/**
 * @file   ParseMat.C
 * @brief  FlowBase ParseMat class
 * @author aics
 */

#include "ParseMat.h"


// #################################################################
/**
 * @brief mat[]に値を格納する
 * @param [in] mat     MediumList
 * @param [in] m       格納順
 * @param [in] key     Propertyのキー
 * @param [in] fval    値
 */
void ParseMat::addVaues(MediumList* mat, const int m, const property_list key, const REAL_TYPE fval)
{
	mat[m].P[key] = fval;
  ChkList[key] = true;
}


// #################################################################
/**
 * @brief ラベルの重複を調べる
 * @param [in] mat     MediumList
 * @param [in] n       現在までに登録した媒質リストの格納番号
 * @param [in] m_label 検査するラベル名
 * @retval 同じラベル名があったらfalseを返す
 */
bool ParseMat::chkDuplicateLabel(const MediumList* mat, const int n, const std::string m_label)
{
	for (int i=1; i<n; i++)
  {
    if ( FBUtility::compare(mat[i].alias, m_label) ) return false;
	}
	return true;
}


// #################################################################
// MediumListを作成する
bool ParseMat::check(const MediumList* mat)
{
  
  for (int i=1; i<=NoMedium; i++)
  {
    if ( !chkList4Solver(mat, i) ) return false;
  }
  
  return true;
}


// #################################################################
/**
 * @brief 媒質情報の内容物をチェックする
 * @param [in] mat MediumList
 * @param [in] m   媒質リストの格納番号
 */
bool ParseMat::chkList4Solver(const MediumList* mat, const int m)
{
  int c=0;
  
  if ( mat[m].getState() == FLUID )
  {
    if ( !ChkList[p_density] )                c += missingMessage(mat, m, p_density);
    if ( !ChkList[p_kinematic_viscosity] )    c += missingMessage(mat, m, p_kinematic_viscosity);
    if ( !ChkList[p_viscosity] )              c += missingMessage(mat, m, p_viscosity);
    if ( !ChkList[p_thermal_conductivity] )   c += missingMessage(mat, m, p_thermal_conductivity);
    if ( !ChkList[p_specific_heat] )          c += missingMessage(mat, m, p_specific_heat);
    if ( !ChkList[p_speed_of_sound] )         c += missingMessage(mat, m, p_speed_of_sound);
    if ( !ChkList[p_vol_expansion] )          c += missingMessage(mat, m, p_vol_expansion);
  }
  else  // solid
  {
    if ( !ChkList[p_density] )                c += missingMessage(mat, m, p_density);
    if ( !ChkList[p_specific_heat] )          c += missingMessage(mat, m, p_specific_heat);
    if ( !ChkList[p_thermal_conductivity] )   c += missingMessage(mat, m, p_thermal_conductivity);
  }
  return ( (c != 0) ? false : true ) ;
}


// #################################################################
// MediumTableを読んでMediumTableInfoクラスオブジェクトに格納
void ParseMat::getMediumTable(const int m_NoMedium, MediumList* mat)
{
  NoMedium   = m_NoMedium;
  
  std::string str, label;
  std::string label_base, label_m, label_leaf;
  REAL_TYPE fval;
  
  label_base = "/MediumTable";
  
  vector<string> nodes;
  int m=1;
  
  // /MediumTable直下のラベルを取得
  tpCntl->getLabelVector(label_base, nodes);
  
  
  // ラベルの重複チェックとセット
  for (vector<string>::iterator it = nodes.begin(); it != nodes.end(); it++)
  {
    if ( !chkDuplicateLabel(mat, m, *it) ) Exit(0);
    mat[m].alias = (*it);
    m++;
  }
  
  
  for (m=1; m<=NoMedium; m++)
  {
    if ( !(tpCntl->getNodeStr(label_base, m, str)) )
    {
      Exit(0);
    }
    
    label_m = label_base + "/" + str;
    
    int n = tpCntl->countLabels(label_m);
    
    for (int i=1; i<=n; i++)
    {
      if ( !(tpCntl->getNodeStr(label_m, i, str)) )
      {
        Exit(0);
      }
      
      
      label_leaf = label_m + "/" + str;
      
      if ( !strcasecmp(str.c_str(), "State") )
      {
        if ( !(tpCntl->getInspectedValue(label_leaf, label)) )
        {
          Exit(0);
        }
        else
        {
          if     ( !strcasecmp(label.c_str(), "Fluid") ) mat[m].setState(FLUID);
          else if( !strcasecmp(label.c_str(), "Solid") ) mat[m].setState(SOLID);
          else
          {
            Exit(0);
          }
        }
      }
      else if ( !strcasecmp(str.c_str(), "MassDensity") )
      {
        if ( !(tpCntl->getInspectedValue(label_leaf, fval )) ) Exit(0);
        addVaues(mat, m, p_density, fval);
      }
      else if ( !strcasecmp(str.c_str(), "SpecificHeat") )
      {
        if ( !(tpCntl->getInspectedValue(label_leaf, fval )) ) Exit(0);
        addVaues(mat, m, p_specific_heat, fval);
      }
      else if ( !strcasecmp(str.c_str(), "ThermalConductivity") )
      {
        if ( !(tpCntl->getInspectedValue(label_leaf, fval )) ) Exit(0);
        addVaues(mat, m, p_thermal_conductivity, fval);
      }
      else if ( !strcasecmp(str.c_str(), "KinematicViscosity") )
      {
        if ( !(tpCntl->getInspectedValue(label_leaf, fval )) ) Exit(0);
        addVaues(mat, m, p_kinematic_viscosity, fval);
      }
      else if ( !strcasecmp(str.c_str(), "Viscosity") )
      {
        if ( !(tpCntl->getInspectedValue(label_leaf, fval )) ) Exit(0);
        addVaues(mat, m, p_viscosity, fval);
      }
      else if ( !strcasecmp(str.c_str(), "SpeedOfSound") )
      {
        if ( !(tpCntl->getInspectedValue(label_leaf, fval )) ) Exit(0);
        addVaues(mat, m, p_speed_of_sound, fval);
      }
      else if ( !strcasecmp(str.c_str(), "VolumeExpansion") )
      {
        if ( !(tpCntl->getInspectedValue(label_leaf, fval )) ) Exit(0);
        addVaues(mat, m, p_vol_expansion, fval);
      }
      else if( !strcasecmp(str.c_str(), "color") )
      {
        ; // FXgenが吐き出す情報でffvcでは不使用
      }
      else
      {
        Exit(0); // 予定していないキーワードの場合エラー
      }
    }
  }
  
}


// #################################################################
//TPのポインタを受け取る
bool ParseMat::importTP(TextParser* tp)
{
  if ( !tp ) return false;
  tpCntl = tp;
  
  return true;
}



// #################################################################
/**
 * @brief 警告メッセージの表示
 * @param [in] mat MediumList
 * @param [in] m   媒質リストの格納番号
 * @param [in] key キーワードの登録番号
 */
int ParseMat::missingMessage(const MediumList* mat, const int m, const int key)
{
  printf("\tMissing keyword '%s'\t\tfor '%s'\t in %s phase\n", 
         MediumList::getPropertyName(key).c_str(),
         mat[m].alias.c_str(),
         ( mat[m].getState() == SOLID ) ? "SOLID" : "FLUID" );
  return 1; 
}



// #################################################################
// 媒質情報の表示
// Hostonly
void ParseMat::printMatList(FILE* fp, const MediumList* mat)
{
  fprintf(fp, "\t  no :            Alias  Physical constant\n");
  fprintf(fp, "\t ----------------------------------------------------------------------\n");
  
  for (int n=1; n<=NoMedium; n++)
  {
    fprintf(fp,"\t%4d : %16s\n", n, mat[n].alias.c_str() );
    
    if ( mat[n].getState() == FLUID )
    {
      fprintf(fp, "\t\t\t\t Mass density         %12.6e [kg/m^3]\n",   mat[n].P[p_density]);
      fprintf(fp, "\t\t\t\t Kinematic viscosity  %12.6e [m^2/s]\n",    mat[n].P[p_kinematic_viscosity]);
      fprintf(fp, "\t\t\t\t Viscosity            %12.6e [Pa s]\n",     mat[n].P[p_viscosity]);
      fprintf(fp, "\t\t\t\t Specific heat        %12.6e [J/(Kg K)]\n", mat[n].P[p_specific_heat]);
      fprintf(fp, "\t\t\t\t Thermal conductivity %12.6e [W/(m K)]\n",  mat[n].P[p_thermal_conductivity]);
      fprintf(fp, "\t\t\t\t Speed of sound       %12.6e [m/s]\n",      mat[n].P[p_speed_of_sound]);
      fprintf(fp, "\t\t\t\t Volume expansion     %12.6e [1/K]\n",      mat[n].P[p_vol_expansion]);
    }
    else // solid
    {
      fprintf(fp, "\t\t\t\t Mass density         %12.6e [kg/m^3]\n",   mat[n].P[p_density]);
      fprintf(fp, "\t\t\t\t Specific heat        %12.6e [J/(Kg K)]\n", mat[n].P[p_specific_heat]);
      fprintf(fp, "\t\t\t\t Thermal conductivity %12.6e [W/(m K)]\n",  mat[n].P[p_thermal_conductivity]);
    }
    
  }
  fprintf(fp,"\n");
  fflush(fp);
}
