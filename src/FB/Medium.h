#ifndef _FB_MEDIUM_H_
#define _FB_MEDIUM_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file Medium.h
//@brief FlowBase Medium class Header
//@author kero

#include <string>
#include "FB_Define.h"

enum property_list {
  p_density,              // 0
  p_kinematic_viscosity,  // 1
  p_viscosity,            // 2
  p_thermal_conductivity, // 3
  p_thermal_diffusivity,  // 4
  p_specific_heat,        // 5
  p_sound_of_speed,       // 6
  p_vol_expansion,        // 7
  property_END
};

class MediumList {
private:
  int  state;         // solid or fluid
  std::string name;   // ラベル

public:
  REAL_TYPE P[property_END];
  
  MediumList() {
    state  = -1;
    for (int i=0; i<property_END; i++) P[i] = 0.0;
  }
  ~MediumList() {}
  
public:

  int getState(void) { return state; }
  
  std::string getLabel(void) { return name; }
  
  void setState(const int key) { 
    state = key; 
  }
  
  void setLabel(const std::string key) {
    name = key;
  }
  
  
  //@brief プロパティ文字列に対応するキー番号を返す
  //@param p 文字列
  static int getKey(const char* p) {
    int key=-1;
    
    if      ( !(strcasecmp(p, "density")) )              key = p_density;
    else if ( !(strcasecmp(p, "kinematic_viscosity")) )  key = p_kinematic_viscosity;
    else if ( !(strcasecmp(p, "viscosity")) )            key = p_viscosity;
    else if ( !(strcasecmp(p, "thermal_conductivity")) ) key = p_thermal_conductivity;
    else if ( !(strcasecmp(p, "thermal_diffusivity")) )  key = p_thermal_diffusivity;
    else if ( !(strcasecmp(p, "specific_heat")) )        key = p_specific_heat;
    else if ( !(strcasecmp(p, "sound_of_speed")) )       key = p_sound_of_speed;
    else if ( !(strcasecmp(p, "volume_expansion")) )     key = p_vol_expansion;
    
    return key;
  }
  

  //@brief プロパティの文字列を返す
  //@param key キー番号
  static std::string getPropertyName(const int key) {
    std::string name;
    
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
  
};

#endif // _FB_MEDIUM_H_
