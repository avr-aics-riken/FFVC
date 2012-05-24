/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file FileIO.C
//@brief FlowBase FileIO class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "FileIO.h"
extern SklParaComponent* ParaCmpo;

/**
 @fn void FileIO::cnv_Div(SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE coef, REAL_TYPE& flop)
 @brief ファイル出力時，発散値を計算する
 @param dst 単位変換のデータクラス
 @param src 単位変換前のデータクラス
 @param coef 係数
 @param flop 浮動小数演算数
 @see SklUtil::cpyS3D()
 */
void FileIO::cnv_Div(SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE coef, REAL_TYPE& flop)
{
  if( !dst || !src ) Exit(0);
  
  const unsigned* dst_sz = dst->GetSize();
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) Exit(0);
  
  REAL_TYPE* dst_data = dst->GetData();
  const REAL_TYPE* src_data = src->GetData();
  if( !dst_data || !src_data ) Exit(0);
  
  dst_sz = dst->_GetSize();
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  register unsigned i, j, k;
  
  unsigned long dst_lsz[3];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  // compare whole size
  if(   (src_sz[0] == dst_sz[0])
     && (src_sz[1] == dst_sz[1])
     && (src_sz[2] == dst_sz[2]) ) {
    unsigned long idx;
    
    for(k=0; k<dst_sz[2]; k++){
      for(j=0; j<dst_sz[1]; j++){
        for(i=0; i<dst_sz[0]; i++){
          idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          dst_data[idx] = src_data[idx]*coef;
        }
      }
    }
    flop += (REAL_TYPE)(dst_sz[0]*dst_sz[1]*dst_sz[2]);
  }
  else {
    Exit(0);
  }
}

/**
 @fn void FileIO::cnv_TP_ND2D(SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, REAL_TYPE& flop)
 @brief 全圧データについて，無次元から有次元単位に変換する
 @param dst 単位変換のデータクラス
 @param src 単位変換前のデータクラス
 @param Ref_rho 代表密度(kg/m^3)
 @param Ref_v 代表速度(m/s)
 @param flop 浮動小数演算数
 @see SklUtil::cpyS3D()
 */
void FileIO::cnv_TP_ND2D(SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, REAL_TYPE& flop)
{
  if( !dst || !src ) Exit(0);
  
  const unsigned* dst_sz = dst->GetSize();
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) Exit(0);
  
  REAL_TYPE* dst_data = dst->GetData();
  const REAL_TYPE* src_data = src->GetData();
  if( !dst_data || !src_data ) Exit(0);
  
  dst_sz = dst->_GetSize();
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  register unsigned i, j, k;
  
  REAL_TYPE c = Ref_rho*Ref_v*Ref_v;
  flop += (REAL_TYPE)(dst_sz[0]*dst_sz[1]*dst_sz[2]*1);
  
  unsigned long src_idx, dst_idx, idx, src_lsz[3], dst_lsz[3];
  src_lsz[0] = src_sz[0];
  src_lsz[1] = src_sz[1];
  src_lsz[2] = src_sz[2];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  // compare whole size
  if(   (src_sz[0] == dst_sz[0])
     && (src_sz[1] == dst_sz[1])
     && (src_sz[2] == dst_sz[2]) ) {
    unsigned idx;
    
    for(k=0; k<dst_sz[2]; k++){
      for(j=0; j<dst_sz[1]; j++){
        for(i=0; i<dst_sz[0]; i++){
          idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          dst_data[idx] = src_data[idx]*c;
        }
      }
    }
  }
  else {
    CalcIndex(dst_sz[0], dst_sz[1], dst_sz[2], dst_gc, src_gc, ix, jx, kx, diff, sta);
    
    for(k=sta; k<kx; k++){
      unsigned kk = k+diff;
      for(j=sta; j<jx; j++){
        unsigned jj = j+diff;
        for(i=sta; i<ix; i++){
          src_idx = src_lsz[0]*src_lsz[1]*kk + src_lsz[0]*jj + i+diff;
          dst_idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          dst_data[dst_idx] = src_data[src_idx]*c;
        }
      }
    }
  }
}


/**
 @fn void FileIO::readSVX(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide,
                          SklScalar3D<int>* dc_mid, const bool vf_mode, SklScalar3D<REAL_TYPE>* dc_ws)
 @brief svxファイルをロードする
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param fname InFileのattrラベル名
 @param size グローバルな分割数
 @param guide ガイドセルサイズ
 @param dc_mid IDを保持するデータクラス
 @param vf_mode Volume Fractionを読み込むモード(trueのとき，デフォルト false）
 @param dc_ws  体積率を保持するデータクラス
 @todo 
    - cpyS3DCenter(dc_ws, dc_tmp_svx) が倍精度のときの挙動
 */
void FileIO::readSVX(SklSolverBase* obj, FILE* fp, const char* fname, unsigned* size, unsigned guide, 
                     SklScalar3D<int>* dc_mid, bool vf_mode, SklScalar3D<REAL_TYPE>* dc_ws)
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  if ( !obj )   Exit(0);
  if (!dc_mid)  Exit(0);
  // if ( !fp )    Exit(0); 並列時にはfpをオープンしない．ホストのみ
  if ( !fname ) Exit(0);
  
  SklVoxDataSet* svxSet = NULL;

  if ( !(svxSet=obj->LoadFile(para_mng, fname, guide, pn.procGrp, size)) ) { // < Global size
    Hostonly_ {
      printf     ("\tError: during SVX file reading\n");
      fprintf(fp, "\tError: during SVX file reading\n");
    }
    Exit(0);
  }

  Hostonly_ {
    printf     ("\n\treading %s\n", svxSet->GetFileName());
    printf     ("\tThis file includes : ");
    fprintf(fp, "\n\treading %s\n", svxSet->GetFileName());
    fprintf(fp, "\tThis file includes : ");
  }
  
  // Medium ID
  if ( !(svxSet->GetData(4)) ) {
    Hostonly_ {
      printf     ("\tThis SVX file does not include Medium ID\n");
      fprintf(fp, "\tThis SVX file does not include Medium ID\n");
    }
    Exit(0);
  }
  else {
    Hostonly_ {
      printf     ("\tMedium ID");
      fprintf(fp, "\tMedium ID");
    }
  }
  SklScalar3D<int>* dc_tmp_mid;
  if ( !(dc_tmp_mid = dynamic_cast<SklScalar3D<int>*>(svxSet->GetData(SklVoxDataSet::SVX_MID))) ) Exit(0);
  if ( !SklUtil::cpyS3DCenter(dc_mid, dc_tmp_mid) ) Exit(0);
  if (dc_tmp_mid) { 
    delete dc_tmp_mid; dc_tmp_mid=NULL; }
  
  // Volume Fraction　オプション
  if (vf_mode == true) {
    if (!dc_ws)   Exit(0);
    if ( svxSet->GetData(0) ) { // 含まれている場合
      if ( sizeof(REAL_TYPE) == sizeof(double) ) {
        Hostonly_ printf("double is not supported at this moment\n");
        if (svxSet) { delete svxSet; svxSet=NULL; }
        return;
      }
      Hostonly_ {
        printf     ("\tVolume Fraction");
        fprintf(fp, "\tVolume Fraction");
      }
      SklScalar3D<REAL_TYPE>* dc_tmp_svx;
      if ( !(dc_tmp_svx = dynamic_cast<SklScalar3D<REAL_TYPE>*>(svxSet->GetData(SklVoxDataSet::SVX_VOLUME))) ) Exit(0);
      if ( !SklUtil::cpyS3DCenter(dc_ws, dc_tmp_svx) ) Exit(0);
      if (dc_tmp_svx) { delete dc_tmp_svx; dc_tmp_svx=NULL; }
    }
    else {
      Hostonly_ {
        printf     ("\tThis SVX file does not include Volume Fraction\n");
        fprintf(fp, "\tThis SVX file does not include Volume Fraction\n");
      }
    }
  }
  
  if (svxSet) { delete svxSet; svxSet=NULL; }
  
  Hostonly_ {
    printf     ("\n");
    fprintf(fp, "\n");
  }
}

/**
 @fn void FileIO::readSBX(SklSolverBase* obj, FILE* fp, const char* mid_str, unsigned* size, unsigned guide, 
                          SklScalar3D<int>* dc_mid, SklScalar3D<REAL_TYPE>* dc_vol)
 @brief SBXファイルの読み込み
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param mid_str コンフィギュレーションファイルにおけるIDの入力情報ラベル
 @param size グローバルな分割数
 @param guide ガイドセルサイズ
 @param dc_mid IDを保持するデータクラス
 @param dc_vol volume rateを保持するデータクラス
 */
void FileIO::readSBX(SklSolverBase* obj, FILE* fp, const char* mid_str, unsigned* size, unsigned guide, 
                     SklScalar3D<int>* dc_mid, SklScalar3D<REAL_TYPE>* dc_vol)
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  if ( !obj )   Exit(0);
  if (!dc_mid)  Exit(0);
  if (!mid_str) Exit(0);
  // if ( !fp )    Exit(0); 並列時にはfpをオープンしない．ホストのみ
  
  const unsigned int* mid_sz = dc_mid->_GetSize();
  
  unsigned long mid_lsz[3];
  mid_lsz[0] = mid_sz[0];
  mid_lsz[1] = mid_sz[1];
  mid_lsz[2] = mid_sz[2];
  
  memset(dc_mid->GetData(), 0x00, sizeof(unsigned int)*mid_lsz[0]*mid_lsz[1]*mid_lsz[2]);
  if (dc_vol != NULL) {
    const unsigned int* vol_sz = dc_vol->_GetSize();
    unsigned long vol_lsz[3];
    vol_lsz[0] = vol_sz[0];
    vol_lsz[1] = vol_sz[1];
    vol_lsz[2] = vol_sz[2];
    memset(dc_vol->GetData(), 0x00, sizeof(unsigned int)*vol_lsz[0]*vol_lsz[1]*vol_lsz[2]);
  }

  // load sbx file for MediumID
  loadSBXfile(obj, fp, mid_str, size, guide, dc_mid, dc_vol);
}

/**
 @fn void FileIO::readSBX(SklSolverBase* obj, FILE* fp, const char* file_attr, unsigned* size, unsigned guide, 
                  SklScalar3D<unsigned char>* dc_mid)
 @brief SBXファイルの読み込み
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param file_str コンフィギュレーションファイルにおける体積率の入力情報ラベル
 @param size グローバルな分割数
 @param guide ガイドセルサイズ
 @param dc_mid IDを保持するデータクラス
 */
void FileIO::readSBX(SklSolverBase* obj, FILE* fp, const char* file_attr, unsigned* size, unsigned guide, 
                     SklScalar3D<unsigned char>* dc_mid)
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  if ( !obj )   Exit(0);
  if (!dc_mid)  Exit(0);
  if (!file_attr) Exit(0);
  // if ( !fp )    Exit(0); 並列時にはfpをオープンしない．ホストのみ
  
  unsigned char* dataP = dc_mid->GetData();
  const unsigned int* mid_sz = dc_mid->_GetSize();
  unsigned long mid_lsz[3];
  mid_lsz[0] = mid_sz[0];
  mid_lsz[1] = mid_sz[1];
  mid_lsz[2] = mid_sz[2];
  memset(dataP, NULL, sizeof(unsigned int)*mid_lsz[0]*mid_lsz[1]*mid_lsz[2]);
  
  // load sbx file for MediumID
  loadSBXfile(obj, fp, file_attr, size, guide, dc_mid);
}

/**
 @fn void FileIO::loadSBXfile(SklSolverBase* obj, FILE* fp, const char* name, unsigned* size, unsigned guide, 
                              SklScalar3D<int>* mid_data, SklScalar3D<REAL_TYPE>* vol_data)
 @brief SBXデータのロード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param name コンフィギュレーションファイルにおける入力情報ラベル
 @param size グローバルなサイズ
 @param guide ガイドセルサイズ
 @param mid_data IDを保持するデータクラス
 @param vol_data volume rateを保持するデータクラス
 */
void FileIO::loadSBXfile(SklSolverBase* obj, FILE* fp, const char* name, unsigned* size, unsigned guide, 
                         SklScalar3D<int>* mid_data, SklScalar3D<REAL_TYPE>* vol_data)
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  if ( !obj )   Exit(0);
  if (!name)   Exit(0);
  if (!mid_data)   Exit(0);
  // if ( !fp )   Exit(0); 並列時にはfpをオープンしない．ホストのみ
  
  SklSbxDataSet* sbxSet = NULL;
  
  // read SBX file
  if ( !(sbxSet=dynamic_cast<SklSbxDataSet*>(obj->LoadFile(para_mng, name, guide, pn.procGrp))) ) {
    Hostonly_ printf     ("\tError: during SBX reading : %s\n", name);
    Hostonly_ fprintf(fp, "\tError: during SBX reading : %s\n", name);
    Exit(0);
  }
  const char* sbx_fname = sbxSet->GetFileName();
  if (sbx_fname) {
    Hostonly_ printf     ("\n\treading %s\n", sbxSet->GetFileName());
    Hostonly_ fprintf(fp, "\n\treading %s\n", sbxSet->GetFileName());
  }
  Hostonly_ printf     ("\n\tThis file includes : ");
  Hostonly_ fprintf(fp, "\n\tThis file includes : ");
  
  // check size/org/pitch between in XML and sbx
  unsigned ix=0, jx=0, kx=0;
  if ( !(sbxSet->GetSize(ix, jx, kx)) ) {
    Hostonly_ printf     ("\tSize in SBX file is out of order : %s\n", name);
    Hostonly_ fprintf(fp, "\tSize in SBX file is out of order : %s\n", name);
    Exit(0);
  }
  
  // Get SBX Data
  SklScalar3D<int>* dc_tmp_sbx = NULL;
  SklVoxDataSet::VszType vsz = sbxSet->GetVszType();
  const unsigned int sz[3] = {ix, jx, kx};
  if( !(dc_tmp_sbx = dynamic_cast<SklScalar3D<int>*>(
                                  obj->SklAllocateArray(
                                    para_mng,
                                    "sbx",
                                    SKL_CLASS_SCALAR3D,
                                    SKL_ARRAY_DTYPE_INT,
                                    sz,
                                    guide,
                                    pn.procGrp))) ) {
                                      Hostonly_ printf     ("Error : SBX AllocateArray(sbx) fault : %s\n", name);
                                      Hostonly_ fprintf(fp, "Error : SBX AllocateArray(sbx) fault : %s\n", name);
                                      Exit(0);
                                    }
  
  if (vsz == SklVoxDataSet::TYPE_CHAR) {
    SklScalar3D<unsigned char>* dc_tmp_char;
    if ( !(dc_tmp_char = dynamic_cast<SklScalar3D<unsigned char>*>(sbxSet->GetData(SklVoxDataSet::SBX_DATA))) ) {
      Hostonly_ printf     ("Error : SBX Data GetData fault : %s\n", name);
      Hostonly_ fprintf(fp, "Error : SBX Data GetData fault : %s\n", name);
      Exit(0);
    }
    const unsigned char* src_ds = dc_tmp_char->GetData();
    int* dest_ds = dc_tmp_sbx->GetData();
    
    // AUX data
    int lowerBits = sbxSet->GetLowerBits();
    int upperBits = sbxSet->GetUpperBits();
    int bitsFlag = sbxSet->GetValidBitsFlag();
    if (lowerBits < 0 || lowerBits > 8) {
      Hostonly_ printf     ("Error : SBX LowerBits[=%d] fault : %s\n", lowerBits, name);
      Hostonly_ fprintf(fp, "Error : SBX LowerBits[=%d] fault : %s\n", lowerBits, name);
      Exit(0);
    }
    if (upperBits < 0 || upperBits > 8) {
      Hostonly_ printf     ("Error : SBX UpperBits[=%d] fault : %s\n", upperBits, name);
      Hostonly_ fprintf(fp, "Error : SBX UpperBits[=%d] fault : %s\n", upperBits, name);
      Exit(0);
    }
    
    Hostonly_ printf     ("\tSBX Data type : 1 byte ");
    Hostonly_ fprintf(fp, "\tSBX Data type : 1 byte ");

    switch (bitsFlag) {
      case 2:
        Hostonly_ printf     ("-> ID only\n");
        Hostonly_ fprintf(fp, "->ID only\n");
        break;
        
      case 0:
        Hostonly_ printf     ("-> 6bits(ID) + 2bits(VF)\n");
        Hostonly_ fprintf(fp, "-> 6bits(ID) + 2bits(VF)\n");
        break;
        
      case 1:
        Hostonly_ printf     ("-> VF only\n");
        Hostonly_ fprintf(fp, "-> VF only\n");
        break;
        
      default:
        Hostonly_ printf     ("-> bit encoding error\n");
        Hostonly_ fprintf(fp, "-> bit encoding error\n");
        Exit(0);
    }
    
    for (register unsigned i=0; i<dc_tmp_char->GetArrayLength(); i++) {
      if (bitsFlag == 0 || bitsFlag == 2) {
        dest_ds[i] = (int)SklUtil::shiftLowerBits(src_ds[i], lowerBits);
      }
    }
    
    if (vol_data != NULL) {
      SklScalar3D<REAL_TYPE>* sbx_vol =
      dynamic_cast<SklScalar3D<REAL_TYPE>*>(
                                           obj->SklAllocateArray(
                                             para_mng,
                                             "sbx_vol",
                                             SKL_CLASS_SCALAR3D,
                                             SKL_ARRAY_DTYPE_REAL,
                                             sz,
                                             guide,
                                             pn.procGrp) );
      if (!sbx_vol) {
        Hostonly_ printf     ("Error : SBX AllocateArray(sbx_vol) fault : %s\n", name);
        Hostonly_ fprintf(fp, "Error : SBX AllocateArray(sbx_vol) fault : %s\n", name);
        Exit(0);
      }
      REAL_TYPE* dest_vol = sbx_vol->GetData();
      for (register unsigned i=0; i<dc_tmp_char->GetArrayLength(); i++) {
        if (bitsFlag == 0 || bitsFlag == 1) {
          int bits = SklUtil::maskUpperBits(src_ds[i], upperBits);
          if (bits == 0x00) dest_vol[i] = 0.0;
          else if (bits == 0x01) dest_vol[i] = 0.5;
          else if (bits == 0x02) dest_vol[i] = 1.0;
        }
      }
      if (!SklUtil::cpyS3DOverlap(vol_data, sbx_vol)) {
        Hostonly_ printf     ("Error : SBX Data cpyS3DOverlap[REAL_TYPE] fault : %s\n", name);
        Hostonly_ fprintf(fp, "Error : SBX Data cpyS3DOverlap[REAL_TYPE] fault : %s\n", name);
        Exit(0);
      }
    }
  }
  else if (vsz == SklVoxDataSet::TYPE_SHORT) {
    Hostonly_ stamped_printf     ("\tSBX Data type : 2 bytes\n");
    Hostonly_ stamped_fprintf(fp, "\tSBX Data type : 2 bytes\n");
    SklScalar3D<unsigned short>* dc_tmp_short;
    if ( !(dc_tmp_short = dynamic_cast<SklScalar3D<unsigned short>*>(sbxSet->GetData(SklVoxDataSet::SBX_DATA))) ) {
      Hostonly_ printf     ("Error : SBX Data GetData fault : %s\n", name);
      Hostonly_ fprintf(fp, "Error : SBX Data GetData fault : %s\n", name);
      Exit(0);
    }
    const unsigned short* src_ds = dc_tmp_short->GetData();
    int* dest_ds = dc_tmp_sbx->GetData();
    for (register unsigned i=0; i<dc_tmp_short->GetArrayLength(); i++) {
      dest_ds[i] = (int)src_ds[i];
    }
  }
  else if (vsz == SklVoxDataSet::TYPE_INT) {
    Hostonly_ stamped_printf     ("\tSBX Data type : 4 bytes\n");
    Hostonly_ stamped_fprintf(fp, "\tSBX Data type : 4 bytes\n");
    SklScalar3D<unsigned int>* dc_tmp_int;
    if ( !(dc_tmp_int = dynamic_cast<SklScalar3D<unsigned int>*>(sbxSet->GetData(SklVoxDataSet::SBX_DATA))) ) {
      Hostonly_ printf     ("Error : SBX Data GetData fault : %s\n", name);
      Hostonly_ fprintf(fp, "Error : SBX Data GetData fault : %s\n", name);
      Exit(0);
    }
    const unsigned int* src_ds = dc_tmp_int->GetData();
    int* dest_ds = dc_tmp_sbx->GetData();
    memcpy(dest_ds, src_ds, dc_tmp_int->GetArrayLength()*sizeof(int));
  }
  
  if (!dc_tmp_sbx) {
    Hostonly_ printf     ("Error : can't set dc_tmp_sbx : %s\n", name);
    Hostonly_ fprintf(fp, "Error : can't set dc_tmp_sbx : %s\n", name);
    Exit(0);
  }
  
  // copy data
  if (!SklUtil::cpyS3DOverlap(mid_data, dc_tmp_sbx)) {
    Hostonly_ printf     ("Error : SBX Data cpyS3DOverlap[int] fault : %s\n", name);
    Hostonly_ fprintf(fp, "Error : SBX Data cpyS3DOverlap[int] fault : %s\n", name);
    Exit(0);
  }
  
  Hostonly_ printf     ("\tget sbx data : %s\n", name);
  Hostonly_ fprintf(fp, "\tget sbx data : %s\n", name);

  if (dc_tmp_sbx) {
    obj->SklGetDataManager()->DeleteDataObj("sbx");
    dc_tmp_sbx=NULL;
  }
  if (sbxSet)     { delete sbxSet; sbxSet=NULL; }
  
  Hostonly_ printf     ("\t\n");
  Hostonly_ fprintf(fp, "\t\n");
}

/**
 @fn void FileIO::loadSBXfile(SklSolverBase* obj, FILE* fp, const char* name, unsigned* size, unsigned guide, SklScalar3D<unsigned char>* data)
 @brief SBXデータのロード
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param name コンフィギュレーションファイルにおける入力情報ラベル
 @param size グローバルなサイズ
 @param guide ガイドセルサイズ
 @param data データクラス
 */
void FileIO::loadSBXfile( SklSolverBase* obj,
                         FILE* fp,
                         const char* name,
                         unsigned* size,
                         unsigned guide,
                         SklScalar3D<unsigned char>* data)
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  if ( !obj )  Exit(0);
  if (!name)   Exit(0);
  if (!data)   Exit(0);
  // if ( !fp )   Exit(0); 並列時にはfpをオープンしない．ホストのみ
  
  SklSbxDataSet* sbxSet = NULL;
  
  // read SBX file
  if ( !(sbxSet=dynamic_cast<SklSbxDataSet*>(obj->LoadFile(para_mng, name, guide, pn.procGrp))) ) {
    Hostonly_ printf     ("\tError: during SBX reading : %s\n", name);
    Hostonly_ fprintf(fp, "\tError: during SBX reading : %s\n", name);
    Exit(0);
  }
  
  const char* sbx_fname = sbxSet->GetFileName();
  if (sbx_fname) {
    Hostonly_ printf     ("\n\treading %s\n", sbxSet->GetFileName());
    Hostonly_ fprintf(fp, "\n\treading %s\n", sbxSet->GetFileName());
  }
  
  Hostonly_ printf     ("\n\tThis file includes : ");
  Hostonly_ fprintf(fp, "\n\tThis file includes : ");
  
  // Get SBX Data
  SklVoxDataSet::VszType vsz = sbxSet->GetVszType();
  
  register int i, j, k;
  if (vsz == SklVoxDataSet::TYPE_CHAR) {
    SklScalar3D<unsigned char>* dc_src_char;
    // get data
    if ( !(dc_src_char = dynamic_cast<SklScalar3D<unsigned char>*>(sbxSet->GetData(SklVoxDataSet::SBX_DATA))) ) {
      Hostonly_ printf     ("Error : SBX Data GetData fault : %s\n", name);
      Hostonly_ fprintf(fp, "Error : SBX Data GetData fault : %s\n", name);
      Exit(0);
    }
    
    // copy data
    if (!SklUtil::cpyS3DOverlap(data, dc_src_char)) {
      Hostonly_ printf     ("Error : SBX Data cpyS3DOverlap fault : %s\n", name);
      Hostonly_ fprintf(fp, "Error : SBX Data cpyS3DOverlap fault : %s\n", name);
      Exit(0);
    }
  }
  else if (vsz == SklVoxDataSet::TYPE_SHORT) {
    Hostonly_ printf     ("Error : SBX Data Type[short] fault : %s\n", name);
    Hostonly_ fprintf(fp, "Error : SBX Data Type[short] fault : %s\n", name);
    Exit(0);
  }
  else if (vsz == SklVoxDataSet::TYPE_INT) {
    Hostonly_ printf     ("Error : SBX Data Type[int] fault : %s\n", name);
    Hostonly_ fprintf(fp, "Error : SBX Data Type[int] fault : %s\n", name);
    Exit(0);
  }
  
  Hostonly_ printf     ("\tget sbx data : %s\n", name);
  Hostonly_ fprintf(fp, "\tget sbx data : %s\n", name);
  
  if (sbxSet)     { delete sbxSet; sbxSet=NULL; }
  
  Hostonly_ printf     ("\t\n");
  Hostonly_ fprintf(fp, "\t\n");
}

/**
 @fn void FileIO::writeRawSPH(const REAL_TYPE *vf, const unsigned* size, const unsigned gc, const REAL_TYPE* org, const REAL_TYPE* ddx, const unsigned m_ModePrecision)
 @brief sphファイルの書き出し（内部領域のみ）
 @param vf スカラデータ
 @param size 配列サイズ
 @param gc ガイドセル
 @param org 基点
 @param ddx ピッチ
 @param m_ModePrecision 浮動小数点の精度
 @note 標記上，long 対応になっているが，ファイルフォーマットとの対応を確認のこと
 */
void FileIO::writeRawSPH(const REAL_TYPE *vf, const unsigned* size, const unsigned gc, const REAL_TYPE* org, const REAL_TYPE* ddx, const unsigned m_ModePrecision)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  int sz, dType, stp, svType;
  int ix, jx, kx, i, j, k;
  unsigned long l, nx;
  REAL_TYPE ox, oy, oz, dx, dy, dz, tm;
  long long szl[3], stpl;
  unsigned m;
  
  char sph_fname[512];
  
  if( para_mng->IsParallel() ){
    sprintf( sph_fname, "field%010d.sph", pn.ID );
  } else {
    sprintf( sph_fname, "field.sph" );
  }
  
  ofstream ofs(sph_fname, ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << sph_fname << " file" << endl;
    Exit(0);
  }
  
  ix = size[0]; //+2*gc;
  jx = size[1]; //+2*gc;
  kx = size[2]; //+2*gc;
  nx = ix * jx * kx;
  ox = org[0]; //-ddx[0]*(REAL_TYPE)gc;
  oy = org[1]; //-ddx[1]*(REAL_TYPE)gc;
  oz = org[2]; //-ddx[2]*(REAL_TYPE)gc;
  dx = ddx[0];
  dy = ddx[1];
  dz = ddx[2];
  //printf("org: %f %f %f\n", ox, oy, oz);
  //printf("dx : %f %f %f\n", dx, dy, dz);
  
  svType = kind_scalar;
  if ( sizeof(REAL_TYPE) == sizeof(double) ) {
    for (i=0; i<3; i++)   szl[i] = (long long)size[i];
  }
  
  REAL_TYPE *f = new REAL_TYPE[nx];
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        l = ix*jx*(k-1) + ix*(j-1) + i-1;
        m = FBUtility::getFindexS3D(size, gc, i, j, k);
        f[l] = (REAL_TYPE)vf[m];
      }
    }
  }
  
  // data property
  ( m_ModePrecision == SPH_SINGLE ) ? dType=1 : dType=2;
  sz = sizeof(unsigned)*2;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&svType, sizeof(int) );
  ofs.write( (char*)&dType, sizeof(int) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // voxel size
  if (dType == 1) {
    sz = sizeof(unsigned)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&ix, sizeof(unsigned) );
    ofs.write( (char*)&jx, sizeof(unsigned) );
    ofs.write( (char*)&kx, sizeof(unsigned) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    sz = sizeof(long long)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&szl[0], sizeof(long long) );
    ofs.write( (char*)&szl[1], sizeof(long long) );
    ofs.write( (char*)&szl[2], sizeof(long long) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  // original point of domain
  if (dType == 1) {
    sz = sizeof(float)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(float) );
    ofs.write( (char*)&oy, sizeof(float) );
    ofs.write( (char*)&oz, sizeof(float) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    sz = sizeof(double)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(double) );
    ofs.write( (char*)&oy, sizeof(double) );
    ofs.write( (char*)&oz, sizeof(double) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  // pitch of voxel
  if (dType == 1) {
    sz = sizeof(float)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(float) );
    ofs.write( (char*)&dy, sizeof(float) );
    ofs.write( (char*)&dz, sizeof(float) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    sz = sizeof(double)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(double) );
    ofs.write( (char*)&dy, sizeof(double) );
    ofs.write( (char*)&dz, sizeof(double) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  // time stamp
  if (dType == 1) {
    stp = 0;
    tm = 0.0;
    sz = sizeof(int)+sizeof(float);
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&stp, sizeof(int) );
    ofs.write( (char*)&tm, sizeof(float) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    stpl =0;
    tm = 0.0;
    sz = sizeof(long long)+sizeof(double);
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&stpl, sizeof(long long) );
    ofs.write( (char*)&tm, sizeof(double) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  if (svType == kind_scalar) {
    sz = (m_ModePrecision == SPH_SINGLE) ? nx * sizeof(float) : nx * sizeof(double);
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)f,   sz );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  ofs.close();
  
  if (f) { delete [] f; f=NULL; }
}



//@fn void FileIO::readPressure()
//@brief 圧力ファイルをロードする
void FileIO::readPressure(FILE* fp,                    /// @param fp ファイルポインタ（ファイル出力）
                          const std::string fname,     /// @param fname ファイル名
                          const unsigned* size,        /// @param size サイズ
                          const unsigned gc,           /// @param guide ガイドセルサイズ
                          REAL_TYPE* p,                /// @param p 圧力データ
                          int& step,                   /// @param step[out] ステップ
                          REAL_TYPE& time,             /// @param time[out] 時刻
                          const unsigned Dmode,        /// @param Dmode 次元（無次元-0 / 有次元-1）
                          const REAL_TYPE BasePrs,     /// @param BasePrs 基準圧力
                          const REAL_TYPE RefDensity,  /// @param RefDensity　代表密度
                          const REAL_TYPE RefVelocity, /// @param RefVelocity 代表速度
                          REAL_TYPE& flop,             /// @param flop
                          const int guide_out,         /// @param guide_out 出力ガイドセル数
                          const bool mode,             /// @param mode 平均値出力指示（瞬時値のときtrue，平均値のときfalse）
                          int& step_avr,               /// @param step_avr 平均操作したステップ数
                          REAL_TYPE& time_avr          /// @param time_avr 平均操作した時間
                          )
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  

  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int g = guide_out;
  int avs = (mode == true) ? 1 : 0;
  
  fb_read_sph_s_ (p, (int*)size, (int*)&gc, tmp, &step, &time, &g, &avs, &step_avr, &time_avr);
  
  if ( !mode ) {
    if ( (step_avr == 0) || (time_avr <= 0.0) ) {
      Hostonly_ printf ("Error : restarted step[%d] or time[%e] is invalid\n", step_avr, time_avr);
      Exit(0);
    }
  }

  
  // 有次元ファイルの場合，無次元に変換する
  if ( Dmode == DIMENSIONAL ) {
    int d_length = (size[0]+2*gc) * (size[1]+2*gc) * (size[2]+2*gc);
    REAL_TYPE scale = (mode == true) ? 1.0 : (REAL_TYPE)step_avr; // 瞬時値の時スケールは1.0、平均値の時は平均数
    REAL_TYPE basep = BasePrs;
    REAL_TYPE ref_d = RefDensity;
    REAL_TYPE ref_v = RefVelocity;
  
    fb_prs_d2nd_(p, &d_length, &basep, &ref_d, &ref_v, &scale, &flop);
  }
  
  if ( mode ) {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }
  else {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n",
                          tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n", 
                          tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
  }

}

/**
 @fn void FileIO::readVelocity(FILE* fp, 
 const std::string fname,
 const unsigned* size, 
 const unsigned gc, 
 REAL_TYPE* v, 
 int& step, 
 REAL_TYPE& time, 
 const REAL_TYPE *v00, 
 const unsigned Dmode, 
 const REAL_TYPE RefVelocity, 
 REAL_TYPE& flop, 
 const bool mode)
 @brief 速度をロードする
 @param fp ファイルポインタ（ファイル出力）
 @param fname ファイル名
 @param size サイズ
 @param guide ガイドセルサイズ
 @param block ブロック数（粗い格子のロードのときに指定、通常は1）
 @param v  結果を保持するデータ
 @param step ステップ
 @param time 時刻
 @param v00[4]
 @param Dmode 次元（無次元-0 / 有次元-1）
 @param RefVelocity 代表速度
 @param flop
 @param guide_out 出力ガイドセル数
 @param mode 平均値出力指示（瞬時値のときtrue，平均値のときfalse）
 @param step_avr 平均操作したステップ数
 @param time_avr 平均操作した時間
 */
void FileIO::readVelocity(FILE* fp, 
                          const std::string fname,
                          const unsigned* size, 
                          const unsigned gc, 
                          REAL_TYPE* v, 
                          int& step, 
                          REAL_TYPE& time, 
                          const REAL_TYPE *v00, 
                          const unsigned Dmode, 
                          const REAL_TYPE RefVelocity, 
                          REAL_TYPE& flop, 
                          const int guide_out,
                          const bool mode,
                          int& step_avr,
                          REAL_TYPE& time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int g = guide_out;
  int avs = (mode == true) ? 1 : 0;
  
  fb_read_sph_v_ (v, (int*)size, (int*)&gc, tmp, &step, &time, &g, &avs, &step_avr, &time_avr);
  
  if ( !mode ) {
    if ( (step_avr == 0) || (time_avr <= 0.0) ) {
      Hostonly_ printf ("Error : restarted step[%d] or time[%e] is invalid\n", step_avr, time_avr);
      Exit(0);
    }
  }

  REAL_TYPE refv = (Dmode == DIMENSIONAL) ? RefVelocity : 1.0;
  REAL_TYPE scale = (mode == true) ? 1.0 : (REAL_TYPE)step_avr; // 瞬時値の時スケールは1.0、平均値の時は平均数
  REAL_TYPE u0[4];
  u0[0] = v00[0];
  u0[1] = v00[1];
  u0[2] = v00[2];
  u0[3] = v00[3];

  fb_shift_refv_in_(v, (int*)size, (int*)&gc, u0, &scale, &refv, &flop);

  if ( mode ) {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }
  else {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n",
                          tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n", 
                      tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
  }

}

/**
 @fn void FileIO::readTemperature(FILE* fp, 
 const std::string fname,
 const unsigned* size, 
 const unsigned gc, 
 REAL_TYPE* t, 
 int& step, 
 REAL_TYPE& time, 
 const unsigned Dmode, 
 const REAL_TYPE Base_tmp, 
 const REAL_TYPE Diff_tmp, 
 const REAL_TYPE Kelvin, 
 REAL_TYPE& flop, 
 const bool mode)
 @brief 温度をロードする
 @param fp ファイルポインタ（ファイル出力）
 @param fname ファイル名
 @param size グローバルなサイズ
 @param gc ガイドセルサイズ
 @param block ブロック数（粗い格子のロードのときに指定、通常は1）
 @param t  結果を保持するデータクラス
 @param step[out] ステップ
 @param time[out] 時刻
 @param Dmode 次元（無次元-0 / 有次元-1）
 @param Base_tmp 基準温度
 @param Diff_tmp　代表温度差
 @param Kelvin 定数
 @param flop
 @param guide_out 出力ガイドセル数
 @param mode 平均値出力指示（瞬時値のときtrue，平均値のときfalse）
 @param step_avr 平均操作したステップ数
 @param time_avr 平均操作した時間
 */
void FileIO::readTemperature(FILE* fp, 
                             const std::string fname,
                             const unsigned* size, 
                             const unsigned gc, 
                             REAL_TYPE* t, 
                             int& step, 
                             REAL_TYPE& time, 
                             const unsigned Dmode, 
                             const REAL_TYPE Base_tmp, 
                             const REAL_TYPE Diff_tmp, 
                             const REAL_TYPE Kelvin, 
                             REAL_TYPE& flop, 
                             const int guide_out,
                             const bool mode,
                             int& step_avr,
                             REAL_TYPE& time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());
  
  int g = guide_out;
  int avs = (mode == true) ? 1 : 0;
  
  fb_read_sph_s_ (t, (int*)size, (int*)&gc, tmp, &step, &time, &g, &avs, &step_avr, &time_avr);
  if ( !mode ) {
    if ( (step_avr == 0) || (time_avr <= 0.0) ) {
      Hostonly_ printf ("Error : restarted step[%d] or time[%e] is invalid\n", step_avr, time_avr);
      Exit(0);
    }
  }
  
  // 有次元ファイルの場合，無次元に変換する
  int d_length = (size[0]+2*gc) * (size[1]+2*gc) * (size[2]+2*gc);
  REAL_TYPE scale = (mode == true) ? 1.0 : (REAL_TYPE)step_avr; // 瞬時値の時スケールは1.0、平均値の時は平均数
  REAL_TYPE base_t = Base_tmp;
  REAL_TYPE diff_t = Diff_tmp;
  REAL_TYPE klv    = Kelvin;
  
  if ( Dmode == DIMENSIONAL ) {
    fb_tmp_d2nd_(t, &d_length, &base_t, &diff_t, &klv, &scale, &flop);
  }
  
  if ( mode ) {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }
  else {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n",
                          tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n", 
                      tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
  }

}

/**
 @fn void FileIO::writeScalar(const std::string fname, 
 const unsigned* size, 
 const unsigned gc,
 REAL_TYPE* s, 
 const int step, 
 const REAL_TYPE time, 
 const REAL_TYPE* org, 
 const REAL_TYPE* pit, 
 const int guide_out)
 @brief スカラー場を出力する
 @param fname ファイル名
 @param size
 @param gc
 @param s スカラー場
 @param step ステップ
 @param time 時刻
 @param org
 @param pit
 @param guide_out ガイドセル数
 @param mode 平均値出力指示（瞬時値のときtrue，平均値のときfalse）
 @param step_avr 平均操作したステップ数
 @param time_avr 平均操作した時間
 */
void FileIO::writeScalar(const std::string fname, 
                         const unsigned* size, 
                         const unsigned gc,
                         REAL_TYPE* s, 
                         const int step, 
                         const REAL_TYPE time, 
                         const REAL_TYPE* org, 
                         const REAL_TYPE* pit, 
                         const int guide_out,
                         const bool mode,
                         const int step_avr,
                         const REAL_TYPE time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int stp = step;
  REAL_TYPE tm = time;
  int g = guide_out;
  REAL_TYPE o[3], p[3];
  o[0] = org[0];
  o[1] = org[1];
  o[2] = org[2];
  p[0] = pit[0];
  p[1] = pit[1];
  p[2] = pit[2];
  
  int avs = (mode == true) ? 1 : 0;
  int stp_a = step_avr;
  REAL_TYPE tm_a = time_avr;
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  
  fb_write_sph_s_ (s, (int*)size, (int*)&gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);
  
}

/**
 @fn void FileIO::writeVector(const std::string fname, 
 const unsigned* size, 
 const unsigned gc, 
 REAL_TYPE* v, 
 const int step, 
 const REAL_TYPE time, 
 const REAL_TYPE* org, 
 const REAL_TYPE* pit, 
 const int guide_out)
 @brief ベクトル場を出力する
 @param fname ファイル名
 @param size
 @param gc
 @param v ベクトル場
 @param step ステップ
 @param time 時刻
 @param org
 @param pit
 @param guide_out ガイドセル数
 @param mode 平均値出力指示（瞬時値のときtrue，平均値のときfalse）
 @param step_avr 平均操作したステップ数
 @param time_avr 平均操作した時間
 */
void FileIO::writeVector(const std::string fname, 
                         const unsigned* size, 
                         const unsigned gc, 
                         REAL_TYPE* v, 
                         const int step, 
                         const REAL_TYPE time, 
                         const REAL_TYPE* org, 
                         const REAL_TYPE* pit, 
                         const int guide_out,
                         const bool mode,
                         const int step_avr,
                         const REAL_TYPE time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());
  
  int stp = step;
  REAL_TYPE tm = time;
  int g = guide_out;
  REAL_TYPE o[3], p[3];
  o[0] = org[0];
  o[1] = org[1];
  o[2] = org[2];
  p[0] = pit[0];
  p[1] = pit[1];
  p[2] = pit[2];
  
  int avs = (mode == true) ? 1 : 0;
  int stp_a = step_avr;
  REAL_TYPE tm_a = time_avr;
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  
  fb_write_sph_v_ (v, (int*)size, (int*)&gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);

}

