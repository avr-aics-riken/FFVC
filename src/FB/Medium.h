#ifndef _FB_MEDIUM_H_
#define _FB_MEDIUM_H_

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
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   Medium.h
 * @brief  FlowBase Medium class Header
 * @author aics
 */

#include "cpm_Define.h"
#include <string>
#include "FB_Define.h"
#include <strings.h>


/** 属性リストのキー */
enum property_list 
{
  p_density,              // 0
  p_kinematic_viscosity,  // 1
  p_viscosity,            // 2
  p_thermal_conductivity, // 3
  p_thermal_diffusivity,  // 4
  p_specific_heat,        // 5
  p_speed_of_sound,       // 6
  p_vol_expansion,        // 7
  property_END
};

class MediumList {
private:
  int state;          ///< solid or fluid


public:
  REAL_TYPE P[property_END]; ///< プロパティリスト
  std::string alias;         ///< 媒質ラベル
  
  MediumList() {
    state = -1;
    for (int i=0; i<property_END; i++) P[i] = 0.0;
  }
  ~MediumList() {}
  
public:
  
  
  /**
   * @brief 文字列に対応するキー番号を返す
   * @param p[in] 文字列
   * @return キー
   */
  static int getKey(const char* p)
  {
    int key=-1;
    
    if      ( !(strcasecmp(p, "MassDensity")) )         key = p_density;
    else if ( !(strcasecmp(p, "KinematicViscosity")) )  key = p_kinematic_viscosity;
    else if ( !(strcasecmp(p, "Viscosity")) )           key = p_viscosity;
    else if ( !(strcasecmp(p, "ThermalConductivity")) ) key = p_thermal_conductivity;
    else if ( !(strcasecmp(p, "ThermalDiffusivity")) )  key = p_thermal_diffusivity;
    else if ( !(strcasecmp(p, "SpecificHeat")) )        key = p_specific_heat;
    else if ( !(strcasecmp(p, "SpeedOfSound")) )        key = p_speed_of_sound;
    else if ( !(strcasecmp(p, "VolumeExpansion")) )     key = p_vol_expansion;
    
    return key;
  }
  
  
  /**
   * @brief プロパティの文字列を返す
   * @param[in] key キー番号
   * @return ラベル
   */
  static std::string getPropertyName(const int key)
  {
    std::string name;
    
    switch (key) {
      case p_density:
        name = "MassDensity";
        break;
      case p_kinematic_viscosity:
        name = "KinematicViscosity";
        break;
      case p_viscosity:
        name = "Viscosity";
        break;
      case p_thermal_conductivity:
        name = "ThermalConductivity";
        break;
      case p_thermal_diffusivity:
        name = "ThermalDiffusivity";
        break;
      case p_specific_heat:
        name = "SpecificHeat";
        break;
      case p_speed_of_sound:
        name = "SpeedOfSound";
        break;
      case p_vol_expansion:
        name = "VolumeExpansion";
        break;
    }
    return name;
  }
  
  
  /**
   * @brief 媒質の属性を取得
   * @return Fluid or Solid
   */
  int getState() const
  {
    return state;
  }
  
  
  /**
   * @brief 状態をセット
   * @param [in] key fluid(1) or solid(0)
   */
  void setState(const int key)
  {
    state = key;
  }
  
};

#endif // _FB_MEDIUM_H_
