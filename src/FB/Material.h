#ifndef _SKL_FB_MATERIAL_H_
#define _SKL_FB_MATERIAL_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Material.h
//@brief FlowBase Material class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <string>
#include <string.h>
#include "FBDefine.h"

using namespace std;

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

class MaterialList {
private:
  unsigned  material_id;   // color id in voxel model
  int       state;         // solid or fluid
  char      name[LABEL];

public:
  REAL_TYPE  P[property_END];
  
  MaterialList() {
    material_id = 0;
    state  = -1;
    for (int i=0; i<property_END; i++) P[i] = 0.0;
    for (int n=0; n<LABEL; n++) name[n]='\0';
  }
  ~MaterialList() {}
  
public:
  static int getKey             (const char* p);
  static string getPropertyName (unsigned key);

  /**
   @fn int getState(void) const
   @brief Solid or Fluidの状態を返す
   @retval 媒質の状態
   */
  int getState(void) const { return state; }
  
  /**
   @fn unsigned getMatID(void) const
   @brief 媒質IDを返す
   @retval 媒質ID
   */
  unsigned getMatID(void) const { return material_id; }
  
  /**
   @fn char* getName(void)
   @brief 媒質名を返す
   @retval 媒質名
   */
  char* getName(void) { return name; }
  
  /**
   @fn void setMatID(const unsigned key)
   @brief 媒質IDを設定する
   @param key 媒質ID
   */
  void setMatID(const unsigned key) { material_id = key; }
  
  /**
   @fn void setMatID(const int key)
   @brief 状態を設定する
   @param key 状態
   */
  void setState(const int key) { state = key; }
};

#endif // _SKL_FB_MATERIAL_H_
