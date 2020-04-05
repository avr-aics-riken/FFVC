#ifndef _FFV_DEFINE_H_
#define _FFV_DEFINE_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/*
 * @file   ffv_Define.h
 * @brief  FFV Class Definition Header
 * @author aics
 */

#define TM_LABEL_MAX 64


// PLOT3DのときのIOバッファのブロックサイズ
// FALLOC::allocArray_Main()
#define IO_BLOCK_SIZE_FLOW 9
#define IO_BLOCK_SIZE_HEAT 16

// PMlibの登録ラベル個数
#define PM_NUM_MAX 200


#endif // _FFV_DEFINE_H_
