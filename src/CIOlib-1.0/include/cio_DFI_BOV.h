#ifndef _CIO_DFI_BOV_H_
#define _CIO_DFI_BOV_H_

/* #################################################################
 *
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) AICS, RIKEN. All right reserved. 2013
 *
 * #################################################################
 */

/** 
 * @file   cio_DFI_BOV.h
 * @brief  cio_DFI_BOV Class Header
 * @author kero    
 */

#include "cio_DFI.h"

using namespace std;


class cio_DFI_BOV : public cio_DFI {

public:
  /** コンストラクタ */
  cio_DFI_BOV();
  
  /**　デストラクタ */
  ~cio_DFI_BOV();
};

#endif // _cio_DFI_BOV_H_
