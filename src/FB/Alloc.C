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
 @fn bool Alloc::alloc_Real_S4DEx(SklSolverBase* obj, SklScalar4D<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                                  REAL_TYPE init, unsigned long &mc, unsigned dlen, int para_key, int m_procGrp)
 @brief SklScalar4DEx<REAL_TYPE>データクラスをアロケートする
 @retval エラーコード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param dc_var データクラス
 @param label データクラスに付与するキーワード
 @param sz 計算内部領域のサイズ
 @param gc ガイドセルサイズ
 @param init 初期値
 @param mc メモリ使用量のカウンタ
 @param dlen 内部に保持するデータ長
 @param para_key コンポーネントの識別ID
 @param m_procGrp プロセスグループ
 @note
    - SklParaComponent()はSklSolverBaseクラスのメソッドなので，SklSolverBaseから派生したクラスのオブジェクトポインタを受け取る
 */
bool Alloc::alloc_Real_S4DEx(SklSolverBase* obj, SklScalar4DEx<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                             REAL_TYPE init, unsigned long &mc, unsigned dlen, int para_key, int m_procGrp)
{
  if ( !obj )   return false;
  if ( dc_var ) return false; // アロケート済みの場合
  if ( !label ) return false;
  if ( !sz )    return false;
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager(para_key);
  
  size_t dims[4];
  dims[0] = dlen;
  dims[1] = sz[0]; 
  dims[2] = sz[1]; 
  dims[3] = sz[2];
  size_t nx = 0;
  
  if ( !obj->SklIsCheckMode() ) {
    if( !(dc_var = dynamic_cast<SklScalar4DEx<REAL_TYPE>*>(
                                                          obj->SklAllocateArray(
                                                                                para_mng,
                                                                                label,
                                                                                SKL_CLASS_SCALAR4DEX,
                                                                                SKL_ARRAY_DTYPE_REAL,
                                                                                dims,
                                                                                gc,
                                                                                m_procGrp))) ) {
			stamped_printf("Error : Alloc::alloc_Real_S4D(%s)\n", label);
			return false;
		}
    nx = dc_var->GetArrayLength();
    if ( nx > LONG_MAX ) {
      stamped_printf("Error : Allocation index overflow %ld\n", nx);
      return false;
    }
    fb_set_value_real_(dc_var->GetData(), (int*)&nx, &init);
  }
  mc += (long)( nx*sizeof(REAL_TYPE) );
  
  return true;
}

/**
 @fn bool Alloc::alloc_Real_S4D(SklSolverBase* obj, SklScalar4D<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                                REAL_TYPE init, unsigned long &mc, unsigned dlen, int para_key, int m_procGrp)
 @brief SklScalar4D<REAL_TYPE>データクラスをアロケートする
 @retval エラーコード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param dc_var データクラス
 @param label データクラスに付与するキーワード
 @param sz 計算内部領域のサイズ
 @param gc ガイドセルサイズ
 @param init 初期値
 @param mc メモリ使用量のカウンタ
 @param dlen 内部に保持するデータ長
 @param para_key コンポーネントの識別ID
 @param m_procGrp プロセスグループ
 @note
    - SklParaComponent()はSklSolverBaseクラスのメソッドなので，SklSolverBaseから派生したクラスのオブジェクトポインタを受け取る
 */
bool Alloc::alloc_Real_S4D(SklSolverBase* obj, SklScalar4D<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                           REAL_TYPE init, unsigned long &mc, unsigned dlen, int para_key, int m_procGrp)
{
  if ( !obj )   return false;
  if ( dc_var ) return false; // アロケート済みの場合
  if ( !label ) return false;
  if ( !sz )    return false;
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager(para_key);
  
  size_t dims[4];
  dims[0] = sz[0]; 
  dims[1] = sz[1]; 
  dims[2] = sz[2];
  dims[3] = dlen;
  size_t nx = 0;
  
  if ( !obj->SklIsCheckMode() ) {
    if( !(dc_var = dynamic_cast<SklScalar4D<REAL_TYPE>*>(
                                                        obj->SklAllocateArray(
                                                                              para_mng,
                                                                              label,
                                                                              SKL_CLASS_SCALAR4D,
                                                                              SKL_ARRAY_DTYPE_REAL,
                                                                              dims,
                                                                              gc,
                                                                              m_procGrp))) ) {
			stamped_printf("Error : Alloc::alloc_Real_S4D(%s)\n", label);
			return false;
		}
    nx = dc_var->GetArrayLength();
    if ( nx > LONG_MAX ) {
      stamped_printf("Error : Allocation index overflow %ld\n", nx);
      return false;
    }
    fb_set_value_real_(dc_var->GetData(), (int*)&nx, &init);
  }
  mc += (long)( nx*sizeof(REAL_TYPE) );
  
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
      fb_set_value_real_(dc_var->GetData(), (int*)&nx, &init);
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
    fb_set_value_int_(dc_var->GetData(), (int*)&nx, &init);
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
    fb_set_value_int_((int*)dc_var->GetData(), (int*)&nx, (int*)&init);
  }
  mc = (long)( nx*sizeof(unsigned) );
  
  return true;
}

/**
 @fn bool Alloc::alloc_Real_V3D(SklSolverBase* obj, SklVector3D<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                                REAL_TYPE init, unsigned long &mc, int para_key, int m_procGrp)
 @brief SklVector3D<REAL_TYPE>データクラスをアロケートする
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
bool Alloc::alloc_Real_V3D(SklSolverBase* obj, SklVector3D<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
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
    if( !(dc_var = dynamic_cast<SklVector3D<REAL_TYPE>*>(
                                                        obj->SklAllocateArray(
                                                                              para_mng,
                                                                              label,
                                                                              SKL_CLASS_VECTOR3D,
                                                                              SKL_ARRAY_DTYPE_REAL,
                                                                              dims,
                                                                              gc,
                                                                              m_procGrp))) ) {
			stamped_printf("Error : Alloc::alloc_Real_V3D(%s)\n", label);
			return false;
		}
    nx = dc_var->GetArrayLength();
    if ( nx > LONG_MAX ) {
      stamped_printf("Error : Allocation index overflow %ld\n", nx);
      return false;
    }
    fb_set_value_real_(dc_var->GetData(), (int*)&nx, &init);
  }
  mc = (long)( nx*sizeof(REAL_TYPE) );
  
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
    fb_set_value_real_(dc_var->GetData(), (int*)&nx, &init);
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
    fb_set_value_int_(dc_var->GetData(), (int*)&nx, &init);
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
    fb_set_value_int_((int*)dc_var->GetData(), (int*)&nx, (int*)&init);
  }
  mc = (long)( nx*sizeof(unsigned) );
  
  return true;
}
