#ifndef _CIO_DFI_F_FUNC_H_
#define _CIO_DFI_F_FUNC_H_

/* #################################################################
 *
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All right reserved.
 *
 * #################################################################
 */

/**
 * @file   cio_DFI_Ffunc.h
 * @brief  Fortran function Header
 * @author kero
 */

extern "C" {

void cio_dfi_coarse_ijkn_(REAL_TYPE* dst, int* dhead, int* dtail, int* dg,
                     REAL_TYPE* src, int* shead, int* stail, int* sg,
                     int* ncomp, int* ista, int* iend, int *irank);

void cio_dfi_coarse_nijk_(REAL_TYPE* dst, int* dhead, int* dtail, int* dg,
                     REAL_TYPE* src, int* shead, int* stail, int* sg,
                     int* ncomp, int* ista, int* iend, int *irank);

void cio_dfi_copy_ijkn_(REAL_TYPE* dst, int* sz, int* dg, REAL_TYPE* src, int *sg, int* ncomp);

void cio_dfi_copy_nijk_(REAL_TYPE* dst, int* sz, int* dg, REAL_TYPE* src, int *sg, int* ncomp);

}

#endif //_CIO_DFI_F_FUNC_H_
