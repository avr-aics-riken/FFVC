#ifndef _CPM_TESTPROG_DEBUG_H_
#define _CPM_TESTPROG_DEBUG_H_

#include "cpm_ParaManager.h"
#include "cpm_TextParserDomain.h"

//読み込んだ領域情報のデバッグライト
void printDomainInfo( cpm_GlobalDomainInfo* dInfo )
{
  const REAL_TYPE *org = dInfo->GetOrigin();
  const REAL_TYPE *pch = dInfo->GetPitch();
  const REAL_TYPE *rgn = dInfo->GetRegion();
  const int       *vox = dInfo->GetVoxNum();
  const int       *div = dInfo->GetDivNum();

  std::cout << "####### read parameters ########" << std::endl;
  std::cout << " G_org      = " << org[0] << "," << org[1] << "," << org[2] << std::endl;
  std::cout << " G_voxel    = " << vox[0] << "," << vox[1] << "," << vox[2] << std::endl;
  std::cout << " G_pitch    = " << pch[0] << "," << pch[1] << "," << pch[2] << std::endl;
  std::cout << " G_region   = " << rgn[0] << "," << rgn[1] << "," << rgn[2] << std::endl;
  std::cout << " G_div      = " << div[0] << "," << div[1] << "," << div[2] << std::endl;
  std::cout << " #subdomain = " << dInfo->GetSubDomainNum() << std::endl;
  for( size_t i=0;i<dInfo->GetSubDomainNum();i++ )
  {
    const cpm_ActiveSubDomainInfo *dom = dInfo->GetSubDomainInfo(i);
    const int *pos  = dom->GetPos();
    const int *bcid = dom->GetBCID();
    std::cout << "  domain" << i
              << "  pos="   << pos[0]  << "," << pos[1]  << "," << pos[2]
              << "  bcid="  << bcid[0] << "," << bcid[1] << "," << bcid[2]
              <<       ","  << bcid[3] << "," << bcid[4] << "," << bcid[5] << std::endl;
  }
}

#endif /* _CPM_TESTPROG_DEBUG_H_ */
