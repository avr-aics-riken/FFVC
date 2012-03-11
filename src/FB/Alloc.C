/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Alloc.C
//@brief FlowBase Aloc class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Alloc.h"
extern SklParaComponent* ParaCmpo;

/**
 @fn bool Alloc::alloc_Real_S3D(SklSolverBase* obj, SklScalar3D<float>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
 REAL_TYPE init, unsigned long &mc, int para_key, int m_procGrp)
 @brief SklScalar3D<float>データクラスをアロケートする
 @retval エラーコード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param dc_var データクラス
 @param label データクラスに付与するキーワード
 @param sz 計算内部領域のサイズ
 @param gc ガイドセルサイズ
 @param init 初期値
 @param mc メモリ使用量のカウンタ
 @param para_key コンポーネントの識別ID
 @param m_procGrp プロセスグループ
 @note
 - SklParaComponent()はSklSolverBaseクラスのメソッドなので，SklSolverBaseから派生したクラスのオブジェクトポインタを受け取る
 */
bool Alloc::alloc_Float_S3D(SklSolverBase* obj, SklScalar3D<float>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                           float init, unsigned long &mc, int para_key, int m_procGrp)
{
  if ( !obj )   return false;
  if ( dc_var ) return false; // アロケート済みの場合
  if ( !label ) return false;
  if ( !sz )    return false;
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager(para_key);
  
  size_t dims[3];
  dims[0] = sz[0]; 
  dims[1] = sz[1]; 
  dims[2] = sz[2];
  size_t nx = 0;
  
  if ( !obj->SklIsCheckMode() ) {
    if( !(dc_var = dynamic_cast<SklScalar3D<float>*>(
                                                         obj->SklAllocateArray(
                                                                               para_mng,
                                                                               label,
                                                                               SKL_CLASS_SCALAR3D,
                                                                               SKL_ARRAY_DTYPE_FLOAT,
                                                                               dims,
                                                                               gc,
                                                                               m_procGrp))) ) {
                                                           stamped_printf("Error : Alloc::alloc_Float_S3D(%s)\n", label);
                                                           return false;
                                                         }
    nx = dc_var->GetArrayLength();
    if ( nx > LONG_MAX ) {
      stamped_printf("Error : Allocation index overflow %ld\n", nx);
      return false;
    }
    if ( nx != 0 ) {
      fb_set_real_s_(dc_var->GetData(), (int*)sz, (int*)&gc, &init);
    }
  }
  mc += (long)( nx*sizeof(float) );
  
  return true;
}

/**
 @fn bool Alloc::alloc_Real_S3D(SklSolverBase* obj, SklScalar3D<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                                REAL_TYPE init, unsigned long &mc, int para_key, int m_procGrp)
 @brief SklScalar3D<REAL_TYPE>データクラスをアロケートする
 @retval エラーコード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param dc_var データクラス
 @param label データクラスに付与するキーワード
 @param sz 計算内部領域のサイズ
 @param gc ガイドセルサイズ
 @param init 初期値
 @param mc メモリ使用量のカウンタ
 @param para_key コンポーネントの識別ID
 @param m_procGrp プロセスグループ
 @note
    - SklParaComponent()はSklSolverBaseクラスのメソッドなので，SklSolverBaseから派生したクラスのオブジェクトポインタを受け取る
 */
bool Alloc::alloc_Real_S3D(SklSolverBase* obj, SklScalar3D<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                           REAL_TYPE init, unsigned long &mc, int para_key, int m_procGrp)
{
  if ( !obj )   return false;
  if ( dc_var ) return false; // アロケート済みの場合
  if ( !label ) return false;
  if ( !sz )    return false;
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager(para_key);
  
  size_t dims[3];
  dims[0] = sz[0]; 
  dims[1] = sz[1]; 
  dims[2] = sz[2];
  size_t nx = 0;
  
  if ( !obj->SklIsCheckMode() ) {
    if( !(dc_var = dynamic_cast<SklScalar3D<REAL_TYPE>*>(
                                                        obj->SklAllocateArray(
                                                                              para_mng,
                                                                              label,
                                                                              SKL_CLASS_SCALAR3D,
                                                                              SKL_ARRAY_DTYPE_REAL,
                                                                              dims,
                                                                              gc,
                                                                              m_procGrp))) ) {
			stamped_printf("Error : Alloc::alloc_Real_S3D(%s)\n", label);
			return false;
		}
    nx = dc_var->GetArrayLength();
    if ( nx > LONG_MAX ) {
      stamped_printf("Error : Allocation index overflow %ld\n", nx);
      return false;
    }
    if ( nx != 0 ) {
      fb_set_real_s_(dc_var->GetData(), (int*)sz, (int*)&gc, &init);
    }
  }
  mc += (long)( nx*sizeof(REAL_TYPE) );
  
  return true;
}

/**
 @fn bool Alloc::alloc_Int_S3D(SklSolverBase* obj, SklScalar3D<int>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                               int init, unsigned long &mc, int para_key, int m_procGrp)
 @brief SklScalar3D<int>データクラスをアロケートする
 @retval エラーコード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param dc_var データクラス
 @param label データクラスに付与するキーワード
 @param sz 計算内部領域のサイズ
 @param gc ガイドセルサイズ
 @param init 初期値
 @param mc メモリ使用量のカウンタ
 @param para_key コンポーネントの識別ID
 @param m_procGrp プロセスグループ
 @note
    - SklParaComponent()はSklSolverBaseクラスのメソッドなので，SklSolverBaseから派生したクラスのオブジェクトポインタを受け取る
 */
bool Alloc::alloc_Int_S3D(SklSolverBase* obj, SklScalar3D<int>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                          int init, unsigned long &mc, int para_key, int m_procGrp)
{
  if ( !obj )    return false;
  if ( dc_var ) return false; // アロケート済みの場合
  if ( !label )  return false;
  if ( !sz )     return false;
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager(para_key);
  
  size_t dims[3];
  dims[0] = sz[0]; 
  dims[1] = sz[1]; 
  dims[2] = sz[2];
  size_t nx = 0;
  
  if ( !obj->SklIsCheckMode() ) {
    if( !(dc_var = dynamic_cast<SklScalar3D<int>*>(
                                                   obj->SklAllocateArray(
                                                                         para_mng,
                                                                         label,
                                                                         SKL_CLASS_SCALAR3D,
                                                                         SKL_ARRAY_DTYPE_INT,
                                                                         dims,
                                                                         gc,
                                                                         m_procGrp))) ) {
			stamped_printf("Error : Alloc::alloc_Int_S3D(%s)\n", label);
			return false;
		}
    nx = dc_var->GetArrayLength();
    if ( nx > LONG_MAX ) {
      stamped_printf("Error : Allocation index overflow %ld\n", nx);
      return false;
    }
    fb_set_int_s_(dc_var->GetData(), (int*)sz, (int*)&gc, &init);
  }
  mc = (long)( nx*sizeof(int) );
  
  return true;
}

/**
 @fn bool Alloc::alloc_Uint_S3D(SklSolverBase* obj, SklScalar3D<unsigned>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                                unsigned init, unsigned long &mc, int para_key, int m_procGrp)
 @brief SklScalar3D<unsigned>データクラスをアロケートする
 @retval エラーコード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param dc_var データクラス
 @param label データクラスに付与するキーワード
 @param sz 計算内部領域のサイズ
 @param gc ガイドセルサイズ
 @param init 初期値
 @param mc メモリ使用量のカウンタ
 @param para_key コンポーネントの識別ID
 @param m_procGrp プロセスグループ
 @note
    - SklParaComponent()はSklSolverBaseクラスのメソッドなので，SklSolverBaseから派生したクラスのオブジェクトポインタを受け取る
 */
bool Alloc::alloc_Uint_S3D(SklSolverBase* obj, SklScalar3D<unsigned>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                           unsigned init, unsigned long &mc, int para_key, int m_procGrp)
{
  if ( !obj )    return false;
  if ( dc_var ) return false; // アロケート済みの場合
  if ( !label )  return false;
  if ( !sz )     return false;
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager(para_key);
  
  size_t dims[3];
  dims[0] = sz[0]; 
  dims[1] = sz[1]; 
  dims[2] = sz[2];
  size_t nx = 0;
  
  if ( !obj->SklIsCheckMode() ) {
    if( !(dc_var = dynamic_cast<SklScalar3D<unsigned>*>(
                                                        obj->SklAllocateArray(
                                                                              para_mng,
                                                                              label,
                                                                              SKL_CLASS_SCALAR3D,
                                                                              SKL_ARRAY_DTYPE_UINT,
                                                                              dims,
                                                                              gc,
                                                                              m_procGrp))) ) {
			stamped_printf("Error : Alloc::alloc_Uint_S3D(%s)\n", label);
			return false;
		}
    nx = dc_var->GetArrayLength();
    if ( nx > LONG_MAX ) {
      stamped_printf("Error : Allocation index overflow %ld\n", nx);
      return false;
    }
    fb_set_int_s_((int*)dc_var->GetData(), (int*)sz, (int*)&gc, (int*)&init);
  }
  mc = (long)( nx*sizeof(unsigned) );
  
  return true;
}

/**
 @fn bool Alloc::alloc_Real_V3DEx(SklSolverBase* obj, SklVector3DEx<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                                  REAL_TYPE init, unsigned long &mc, int para_key, int m_procGrp)
 @brief SklVector3DEx<REAL_TYPE>データクラスをアロケートする
 @retval エラーコード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param dc_var データクラス
 @param label データクラスに付与するキーワード
 @param sz 計算内部領域のサイズ
 @param gc ガイドセルサイズ
 @param init 初期値
 @param mc メモリ使用量のカウンタ
 @param para_key コンポーネントの識別ID
 @param m_procGrp プロセスグループ
 @note
    - SklParaComponent()はSklSolverBaseクラスのメソッドなので，SklSolverBaseから派生したクラスのオブジェクトポインタを受け取る
 */
bool Alloc::alloc_Real_V3DEx(SklSolverBase* obj, SklVector3DEx<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                             REAL_TYPE init, unsigned long &mc, int para_key, int m_procGrp)
{
  if ( !obj )    return false;
  if ( dc_var ) return false; // アロケート済みの場合
  if ( !label )  return false;
  if ( !sz )     return false;
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager(para_key);
  
  size_t dims[3];
  dims[0] = sz[0]; 
  dims[1] = sz[1]; 
  dims[2] = sz[2];
  size_t nx = 0;
  
  if ( !obj->SklIsCheckMode() ) {
    if( !(dc_var = dynamic_cast<SklVector3DEx<REAL_TYPE>*>(
                                                          obj->SklAllocateArray(
                                                                                para_mng,
                                                                                label,
                                                                                SKL_CLASS_VECTOR3DEX,
                                                                                SKL_ARRAY_DTYPE_REAL,
                                                                                dims,
                                                                                gc,
                                                                                m_procGrp))) ) {
			stamped_printf("Error : Alloc::alloc_Real_V3DEx(%s)\n", label);
			return false;
		}
    nx = dc_var->GetArrayLength();
    if ( nx > LONG_MAX ) {
      stamped_printf("Error : Allocation index overflow %ld\n", nx);
      return false;
    }
    fb_set_real_v_(dc_var->GetData(), (int*)sz, (int*)&gc, &init);
  }
  mc = (long)( nx*sizeof(REAL_TYPE) );
  
  return true;
}

/**
 @fn bool Alloc::alloc_Int_S1D(SklSolverBase* obj, SklScalar<int>* &dc_var, const char* label, unsigned sz, unsigned gc, 
                               int init, unsigned long &mc, int para_key, int m_procGrp)
 @brief SklScalar3D<int>データクラスをアロケートする
 @retval エラーコード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param dc_var データクラス
 @param label データクラスに付与するキーワード
 @param sz 計算内部領域のサイズ
 @param gc ガイドセルサイズ
 @param init 初期値
 @param mc メモリ使用量のカウンタ
 @param para_key コンポーネントの識別ID
 @param m_procGrp プロセスグループ
 @note
 - SklParaComponent()はSklSolverBaseクラスのメソッドなので，SklSolverBaseから派生したクラスのオブジェクトポインタを受け取る
 */
bool Alloc::alloc_Int_S1D(SklSolverBase* obj, SklScalar<int>* &dc_var, const char* label, unsigned sz, unsigned gc, 
                          int init, unsigned long &mc, int para_key, int m_procGrp)
{
  if ( !obj )    return false;
  if ( dc_var ) return false; // アロケート済みの場合
  if ( !label )  return false;
  if ( !sz )     return false;
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager(para_key);
  
  size_t dims = sz; 
  size_t nx = 0;
  
  if ( !obj->SklIsCheckMode() ) {
    if( !(dc_var = dynamic_cast<SklScalar<int>*>(
                                                   obj->SklAllocateArray(
                                                                         para_mng,
                                                                         label,
                                                                         SKL_CLASS_SCALAR,
                                                                         SKL_ARRAY_DTYPE_INT,
                                                                         &dims,
                                                                         gc,
                                                                         m_procGrp))) ) {
                                                     stamped_printf("Error : Alloc::alloc_Int_S3D(%s)\n", label);
                                                     return false;
                                                   }
    nx = dc_var->GetArrayLength();
    if ( nx > LONG_MAX ) {
      stamped_printf("Error : Allocation index overflow %ld\n", nx);
      return false;
    }
    fb_set_int_s_(dc_var->GetData(), (int*)sz, (int*)&gc, &init);
  }
  mc = (long)( nx*sizeof(int) );
  
  return true;
}

/**
 @fn bool Alloc::alloc_Uint_S1D(SklSolverBase* obj, SklScalar<unsigned>* &dc_var, const char* label, unsigned sz, unsigned gc, 
                                int init, unsigned long &mc, int para_key, int m_procGrp)
 @brief SklScalar3D<int>データクラスをアロケートする
 @retval エラーコード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param dc_var データクラス
 @param label データクラスに付与するキーワード
 @param sz 計算内部領域のサイズ
 @param gc ガイドセルサイズ
 @param init 初期値
 @param mc メモリ使用量のカウンタ
 @param para_key コンポーネントの識別ID
 @param m_procGrp プロセスグループ
 @note
 - SklParaComponent()はSklSolverBaseクラスのメソッドなので，SklSolverBaseから派生したクラスのオブジェクトポインタを受け取る
 */
bool Alloc::alloc_Uint_S1D(SklSolverBase* obj, SklScalar<unsigned>* &dc_var, const char* label, unsigned sz, unsigned gc, 
                           unsigned init, unsigned long &mc, int para_key, int m_procGrp)
{
  if ( !obj )    return false;
  if ( dc_var ) return false; // アロケート済みの場合
  if ( !label )  return false;
  if ( !sz )     return false;
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager(para_key);
  
  size_t dims = sz; 
  size_t nx = 0;
  
  if ( !obj->SklIsCheckMode() ) {
    if( !(dc_var = dynamic_cast<SklScalar<unsigned>*>(
                                                 obj->SklAllocateArray(
                                                                       para_mng,
                                                                       label,
                                                                       SKL_CLASS_SCALAR,
                                                                       SKL_ARRAY_DTYPE_UINT,
                                                                       &dims,
                                                                       gc,
                                                                       m_procGrp))) ) {
                                                   stamped_printf("Error : Alloc::alloc_Uint_S3D(%s)\n", label);
                                                   return false;
                                                 }
    nx = dc_var->GetArrayLength();
    if ( nx > LONG_MAX ) {
      stamped_printf("Error : Allocation index overflow %ld\n", nx);
      return false;
    }
    fb_set_int_s_((int*)dc_var->GetData(), (int*)sz, (int*)&gc, (int*)&init);
  }
  mc = (long)( nx*sizeof(unsigned) );
  
  return true;
}
