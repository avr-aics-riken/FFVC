/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file Material.C
//@brief FlowBase MaterialList class
//@author keno, FSI Team, VCAD, RIKEN

#include "Material.h"

/**
 @fn int MaterialList::getKey(const char* p) 
 @brief プロパティ文字列に対応するキー番号を返す
 @retval キー番号
 @param p 文字列
 */
int MaterialList::getKey(const char* p) 
{
  int key=-1;
  if      ( !(strcasecmp(p, "density")) )              key = (unsigned)p_density;
  else if ( !(strcasecmp(p, "kinematic_viscosity")) )  key = (unsigned)p_kinematic_viscosity;
  else if ( !(strcasecmp(p, "viscosity")) )            key = (unsigned)p_viscosity;
  else if ( !(strcasecmp(p, "thermal_conductivity")) ) key = (unsigned)p_thermal_conductivity;
  else if ( !(strcasecmp(p, "thermal_diffusivity")) )  key = (unsigned)p_thermal_diffusivity;
  else if ( !(strcasecmp(p, "specific_heat")) )        key = (unsigned)p_specific_heat;
  else if ( !(strcasecmp(p, "sound_of_speed")) )       key = (unsigned)p_sound_of_speed;
  else if ( !(strcasecmp(p, "volume_expansion")) )     key = (unsigned)p_vol_expansion;
  return key;
}

/**
 @fn string MaterialList::getPropertyName(unsigned key)
 @brief プロパティの文字列を返す
 @retval プロパティの文字列
 @param key キー番号
 */
string MaterialList::getPropertyName(unsigned key)
{
  string name;
  
  switch (key) {
    case p_density:
      name = "density";
      break;
    case p_kinematic_viscosity:
      name = "kinematic_viscosity";
      break;
    case p_viscosity:
      name = "viscosity";
      break;
    case p_thermal_conductivity:
      name = "thermal_conductivity";
      break;
    case p_thermal_diffusivity:
      name = "thermal_diffusivity";
      break;
    case p_specific_heat:
      name = "specific_heat";
      break;
    case p_sound_of_speed:
      name = "sound_of_speed";
      break;
    case p_vol_expansion:
      name = "volume_expansion";
      break;
  }
  return name;
}
