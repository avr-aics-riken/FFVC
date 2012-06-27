  ///////////// MPI test /////////////
  {
    int rank = paraMngr->GetMyRankID();
    char fname[512];
    sprintf ( fname, "mpiTest_%04d.log", rank );
    std::ofstream ofs( fname );
  
    // Broardcast
#if 1
    {
      // process group 0
      int buf = 0;
      int rank = paraMngr->GetMyRankID();
      if( rank==0 ) buf = 123;
      if( rank==0 ) cout << endl << "***** Bcast test pg[" << 0 << "]" << endl;
      ofs << endl << "***** Bcast test pg[" << 0 << "]" << endl;
      ofs << " [" << rank << "] before Bcast = " << buf << endl;
      paraMngr->Bcast( &buf, 1, 0 );
      ofs << " [" << rank << "] after  Bcast = " << buf << endl;
    }
    {
      // process group 1
      int pg = 1;
      REAL_TYPE buf = 0.0;
      int rank0 = paraMngr->GetMyRankID();
      int rank1 = paraMngr->GetMyRankID(pg);
      if( rank1==0 ) buf = 456.789;
      if( rank1==0 ) cout << endl << "***** Bcast test pg[" << pg << "]" << endl;
      ofs << endl << "***** Bcast test pg[" << pg << "]" << endl;
      ofs << " [" << rank0 << "/" << rank1 << "] before Bcast = " << buf << endl;
      paraMngr->Bcast( &buf, 1, 0, pg );
      ofs << " [" << rank0 << "/" << rank1 << "] after  Bcast = " << buf << endl;
    }
#endif

    // Send/Recv
#if 1
    {
      // process group 0
      int buf = 0;
      int rank = paraMngr->GetMyRankID();
      if( rank==0 ) buf = 123;
      if( rank==0 ) cout << endl << "***** Send/Recv test pg[" << 0 << "]" << endl;
      ofs << endl << "***** Send/Recv test pg[" << 0 << "]" << endl;
      ofs << " [" << rank << "] before Send/Recv = " << buf << endl;
      if( rank==0 )//send
        paraMngr->Send( &buf, 1, 1 );
      if( rank==1 )//recv
        paraMngr->Recv( &buf, 1, 0 );
      ofs << " [" << rank << "] after  Send/Recv = " << buf << endl;
    }
    {
      // process group 1
      int pg = 1;
      double buf = 0;
      int rank0 = paraMngr->GetMyRankID();
      int rank1 = paraMngr->GetMyRankID(pg);
      if( rank1==0 ) buf = 789.123;
      if( rank1==0 ) cout << endl << "***** Send/Recv test pg[" << 0 << "]" << endl;
      ofs << endl << "***** Send/Recv test pg[" << 0 << "]" << endl;
      ofs << " [" << rank0 << "/" << rank1 << "] before Send/Recv = " << buf << endl;
      if( rank1==0 )//send
        paraMngr->Send( &buf, 1, 1, pg );
      if( rank1==1 )//recv
        paraMngr->Recv( &buf, 1, 0, pg );
      ofs << " [" << rank0 << "/" << rank1 << "] after  Send/Recv = " << buf << endl;
    }
#endif

    // Isend/Irecv
#if 1
    MPI_Request req;
    {
      // process group 0
      int buf = 0;
      int rank = paraMngr->GetMyRankID();
      if( rank==0 ) buf = 123;
      if( rank==0 ) cout << endl << "***** Isend/Irecv test pg[" << 0 << "]" << endl;
      ofs << endl << "***** Isend/Irecv test pg[" << 0 << "]" << endl;
      ofs << " [" << rank << "] before Isend/Irecv = " << buf << endl;
      if( rank==0 )//Isend
        paraMngr->Isend( &buf, 1, 1, &req );
      if( rank==1 )//Irecv
        paraMngr->Irecv( &buf, 1, 0, &req );

      // process group 0 wait
      if( rank==0 || rank==1 )
        paraMngr->Wait( &req );
      ofs << " [" << rank << "] after  Isend/Irecv = " << buf << endl;
    }
    {
      // process group 1
      int pg = 1;
      double buf = 0;
      int rank0 = paraMngr->GetMyRankID();
      int rank1 = paraMngr->GetMyRankID(pg);
      if( rank1==0 ) buf = 789.123;
      if( rank1==0 ) cout << endl << "***** Isend/Irecv test pg[" << 0 << "]" << endl;
      ofs << endl << "***** Isend/Irecv test pg[" << 0 << "]" << endl;
      ofs << " [" << rank0 << "/" << rank1 << "] before Isend/Irecv = " << buf << endl;
      if( rank1==0 )//Isend
        paraMngr->Isend( &buf, 1, 1, &req, pg );
      if( rank1==1 )//Irecv
        paraMngr->Irecv( &buf, 1, 0, &req, pg );
      // process group 1 wait
      if( rank1==0 || rank1==1 )
        paraMngr->Wait( &req );
      ofs << " [" << rank0 << "/" << rank1 << "] after  Isend/Irecv = " << buf << endl;
    }
#endif

    // Allreduce
#if 1
    {
      // process group 0
      int rank = paraMngr->GetMyRankID();
      int sendbuf = rank+1;
      int recvbuf = rank+1;
      if( rank==0 ) cout << endl << "***** Allreduce test pg[" << 0 << "]" << endl;
      ofs << endl << "***** Allreduce test pg[" << 0 << "]" << endl;
      ofs << " [" << rank << "] before Allreduce(SUM) = " << sendbuf << endl;
      paraMngr->Allreduce( &sendbuf, &recvbuf, 1, MPI_SUM );
      ofs << " [" << rank << "] after  Allreduce(SUM) = " << recvbuf << endl;
    }
    {
      // process group 1
      int pg = 1;
      int rank0 = paraMngr->GetMyRankID();
      int rank1 = paraMngr->GetMyRankID(pg);
      REAL_TYPE sendbuf = (rank1+1) * 1e-2;
      REAL_TYPE recvbuf = (rank1+1) * 1e-2;
      if( rank1==0 ) cout << endl << "***** Allreduce test pg[" << pg << "]" << endl;
      ofs << endl << "***** Allreduce test pg[" << pg << "]" << endl;
      ofs << " [" << rank0 << "/" << rank1 << "] before Allreduce(MAX) = " << sendbuf << endl;
      paraMngr->Allreduce( &sendbuf, &recvbuf, 1, MPI_MAX, pg );
      ofs << " [" << rank0 << "/" << rank1 << "] after  Allreduce(MAX) = " << recvbuf << endl;
    }
#endif

    // Gather / Allgather
#if 1
    {
      // process group 0
      int rank = paraMngr->GetMyRankID();
      int nrank = paraMngr->GetNumRank();
      int sendbuf = rank+1;
      int *recvbuf = new int[nrank];
      for( int i=0;i<nrank;i++ ) recvbuf[i] = 0;
      if( rank==0 ) cout << endl << "***** Gather test pg[" << 0 << "]" << endl;
      ofs << endl << "***** Gather test pg[" << 0 << "]" << endl;
      ofs << " [" << rank << "] before Gather = " << sendbuf << endl;
      paraMngr->Gather( &sendbuf, 1, recvbuf, 1, 1 );
      string str ="";
      for( int i=0;i<nrank;i++ )
      {
        char buf[100];
        sprintf( buf, "%d", recvbuf[i] );
        if( i!=nrank-1 ) strcat( buf, "," );
        str += buf;
      }
      ofs << " [" << rank << "] after  Gather = " << str << endl;
      delete [] recvbuf;
    }
    {
      // process group 1
      int pg = 1;
      int rank0 = paraMngr->GetMyRankID();
      int nrank0 = paraMngr->GetNumRank();
      int rank1 = paraMngr->GetMyRankID(pg);
      REAL_TYPE sendbuf = (rank1+1) * 1e-2;
      REAL_TYPE *recvbuf = new REAL_TYPE[nrank0];
      for( int i=0;i<nrank0;i++ ) recvbuf[i] = 0.0;
      if( rank1==0 ) cout << endl << "***** Allgather test pg[" << pg << "]" << endl;
      ofs << endl << "***** Allgather test pg[" << pg << "]" << endl;
      ofs << " [" << rank0 << "/" << rank1 << "] before Allgather = " << sendbuf << endl;
      paraMngr->Allgather( &sendbuf, 1, recvbuf, 1, pg );
      string str ="";
      for( int i=0;i<nrank0;i++ )
      {
        char buf[100];
        sprintf( buf, "%g", recvbuf[i] );
        if( i!=nrank0-1 ) strcat( buf, "," );
        str += buf;
      }
      ofs << " [" << rank0 << "/" << rank1 << "] after  Allgather = " << str << endl;
      delete [] recvbuf;
    }
#endif

    // Gatherv / Allgatherv
#if 1
    {
      // process group 0
      int rank = paraMngr->GetMyRankID();
      int np   = paraMngr->GetNumRank();
      int sendcnt = rank+1;
      int *sendbuf = new int[sendcnt];
      for( int i=0;i<sendcnt;i++ ) sendbuf[i] = rank+1;
      int *displs = new int[np];
      int *recvcnts = new int[np];
      int recvcnt = 0;
      for( int i=0;i<np;i++ )
      {
        recvcnts[i] = i+1;
        displs[i] = recvcnt;
        recvcnt += recvcnts[i];
      }
      int *recvbuf = new int[recvcnt+1];
      for( int i=0;i<recvcnt+1;i++ ) recvbuf[i] = 0;
      if( rank==0 ) cout << endl << "***** Gatherv test pg[" << 0 << "]" << endl;
      ofs << endl << "***** Gatherv test pg[" << 0 << "]" << endl;
      string strsend ="";
      for( int i=0;i<sendcnt;i++ )
      {
        char buf[100];
        sprintf( buf, "%d", sendbuf[i] );
        if( i!=sendcnt-1 ) strcat( buf, "," );
        strsend += buf;
      }
      ofs << " [" << rank << "] before Gatherv = " << strsend << endl;
      paraMngr->Gatherv( sendbuf, sendcnt, recvbuf, recvcnts, displs, 1 );
      string str ="";
      for( int i=0;i<recvcnt+1;i++ )
      {
        char buf[10];
        sprintf( buf, "%d", recvbuf[i] );
        if( i!=recvcnt ) strcat( buf, "," );
        str += buf;
      }
      ofs << " [" << rank << "] after  Gatherv = " << str << endl;
      delete [] sendbuf;
      delete [] displs;
      delete [] recvcnts;
      delete [] recvbuf;
    }
    {
      // process group 1
      int pg = 1;
      int rank0 = paraMngr->GetMyRankID();
      int rank  = paraMngr->GetMyRankID(pg);
      int np0   = paraMngr->GetNumRank();
      int np    = paraMngr->GetNumRank(pg);
      if( rank >= 0 )
      {
        int sendcnt = rank+1;
        int *sendbuf = new int[sendcnt];
        for( int i=0;i<sendcnt;i++ ) sendbuf[i] = rank+1;
        int *displs = new int[np];
        int *recvcnts = new int[np];
        int recvcnt = 0;
        for( int i=0;i<np;i++ )
        {
          recvcnts[i] = i+1;
          displs[i] = recvcnt;
          recvcnt += recvcnts[i];
        }
        int *recvbuf = new int[recvcnt+1];
        for( int i=0;i<recvcnt+1;i++ ) recvbuf[i] = 0;
        if( rank==0 ) cout << endl << "***** Allgatherv test pg[" << 0 << "]" << endl;
        ofs << endl << "***** Allgatherv test pg[" << 0 << "]" << endl;
        string strsend ="";
        for( int i=0;i<sendcnt;i++ )
        {
          char buf[100];
          sprintf( buf, "%d", sendbuf[i] );
          if( i!=sendcnt-1 ) strcat( buf, "," );
          strsend += buf;
        }
        ofs << " [" << rank0 << "/" << rank << "] before Allgatherv = " << strsend << endl;
        paraMngr->Allgatherv( sendbuf, sendcnt, recvbuf, recvcnts, displs, pg );
        string str ="";
        for( int i=0;i<recvcnt+1;i++ )
        {
          char buf[10];
          sprintf( buf, "%d", recvbuf[i] );
          if( i!=recvcnt ) strcat( buf, "," );
          str += buf;
        }
        ofs << " [" << rank0 << "/" << rank << "] after  Allgatherv = " << str << endl;
        delete [] sendbuf;
        delete [] displs;
        delete [] recvcnts;
        delete [] recvbuf;

        if( rank==0 ) cout << endl;
      }
    }
#endif
  }

//#define _NOWAIT_TEST_

#define _PRINT_S4D(_OFS,_TITLE,_IMAX,_JMAX,_KMAX,_NMAX,_VC,_N) \
{ \
    _OFS << endl; \
    _OFS << "#### " << _TITLE << " I-J ####" << endl; \
    { \
      _OFS << "J |" << endl; \
      int k=_KMAX/2; \
      for( int j=_JMAX+_VC-1;j>=0-_VC;j-- ){ \
      stringstream ss; \
      for( int i=0-_VC;i<_IMAX+_VC;i++ ){ \
        size_t idx = _IDX_S4D(i,j,k,_N,_IMAX,_JMAX,_KMAX,_VC); \
        ss.width(www); \
        ss << pp[idx]; \
        if( i<_IMAX+_VC-1 ) ss << ", "; \
      } \
      _OFS << ss.str(); \
      if( j==0-_VC ) _OFS << "  -> I"; \
      _OFS << endl; \
      } \
    } \
\
    _OFS << endl; \
    _OFS << "#### " << _TITLE << " J-K ####" << endl; \
    { \
      _OFS << "K |" << endl; \
      int i=_IMAX/2; \
      for( int k=_KMAX+_VC-1;k>=0-_VC;k-- ){ \
      stringstream ss; \
      for( int j=0-_VC;j<_JMAX+_VC;j++ ){ \
        size_t idx = _IDX_S4D(i,j,k,_N,_IMAX,_JMAX,_KMAX,_VC); \
        ss.width(www); \
        ss << pp[idx]; \
        if( j<_JMAX+_VC-1 ) ss << ", "; \
      } \
      _OFS << ss.str(); \
      if( k==0-_VC ) _OFS << "  -> J"; \
      _OFS << endl; \
      } \
    } \
}

#define _PRINT_S4DEX(_OFS,_TITLE,_NMAX,_IMAX,_JMAX,_KMAX,_VC,_N) \
{ \
    _OFS << endl; \
    _OFS << "#### " << _TITLE << " I-J ####" << endl; \
    { \
      _OFS << "J |" << endl; \
      int k=_KMAX/2; \
      for( int j=_JMAX+_VC-1;j>=0-_VC;j-- ){ \
      stringstream ss; \
      for( int i=0-_VC;i<_IMAX+_VC;i++ ){ \
        size_t idx = _IDX_S4DEX(_N,i,j,k,_NMAX,_IMAX,_JMAX,_KMAX,_VC); \
        ss.width(www); \
        ss << pp[idx]; \
        if( i<_IMAX+_VC-1 ) ss << ", "; \
      } \
      _OFS << ss.str(); \
      if( j==0-_VC ) _OFS << "  -> I"; \
      _OFS << endl; \
      } \
    } \
\
    _OFS << endl; \
    _OFS << "#### " << _TITLE << " J-K ####" << endl; \
    { \
      _OFS << "K |" << endl; \
      int i=_IMAX/2; \
      for( int k=_KMAX+_VC-1;k>=0-_VC;k-- ){ \
      stringstream ss; \
      for( int j=0-_VC;j<_JMAX+_VC;j++ ){ \
        size_t idx = _IDX_S4DEX(_N,i,j,k,_NMAX,_IMAX,_JMAX,_KMAX,_VC); \
        ss.width(www); \
        ss << pp[idx]; \
        if( j<_JMAX+_VC-1 ) ss << ", "; \
      } \
      _OFS << ss.str(); \
      if( k==0-_VC ) _OFS << "  -> J"; \
      _OFS << endl; \
      } \
    } \
}

  // BndComm, PeriodicComm(S4D)
#if 1
  {
    int rank = paraMngr->GetMyRankID();
    const int *sz = paraMngr->GetLocalVoxelSize();
    int imax = sz[0];
    int jmax = sz[1];
    int kmax = sz[2];
    int vc  = 3;
    size_t nw = size_t(imax+2*vc) * size_t(jmax+2*vc) * size_t(kmax+2*vc);

    REAL_TYPE *pp = paraMngr->AllocRealS3D(vc);
    for( int k=0-vc;k<kmax+vc;k++ ){
    for( int j=0-vc;j<jmax+vc;j++ ){
    for( int i=0-vc;i<imax+vc;i++ ){
      pp[_IDX_S3D(i,j,k,imax,jmax,kmax,vc)] = -(rank+1);
    }}}
    int cnt = 0;
    for( int k=0;k<kmax;k++ ){
    for( int j=0;j<jmax;j++ ){
    for( int i=0;i<imax;i++ ){
      cnt++;
      pp[_IDX_S3D(i,j,k,imax,jmax,kmax,vc)] = cnt + (rank+1)*1000;
    }}}

    char fname[512];
    sprintf ( fname, "commS3D_%04d.log", rank );
    std::ofstream ofs( fname );
#ifdef _NOWAIT_TEST_
    ofs << "#### commS3D nowait test ####" << endl << endl;
    if( rank==0 ) cout << "#### commS3D nowait test ####" << endl << endl;
#else
    ofs << "#### commS3D test ####" << endl << endl;
    if( rank==0 ) cout << "#### commS3D test ####" << endl << endl;
#endif

    size_t mem = paraMngr->GetBndCommBufferSize(-1);
    ofs << "buff mem = " << cpm_Base::GetMemString(mem) << endl;


    int n=0;
    int www = 6;
    ofs << "n=" << n << endl;

    _PRINT_S4D(ofs,"before comm",imax,jmax,kmax,1,vc,n);

#ifndef _NOWAIT_TEST_
    paraMngr->BndCommS3D( pp, imax, jmax, kmax, vc, 2 );
#else
    MPI_Request req[12];
    paraMngr->BndCommS3D_nowait( pp, imax, jmax, kmax, vc, 2, req );
    paraMngr->wait_BndCommS3D( pp, imax, jmax, kmax, vc, 2, req );
#endif

    _PRINT_S4D(ofs,"after BndComm",imax,jmax,kmax,1,vc,n);

    // PeriodicCommX
    paraMngr->PeriodicCommS3D( pp, imax, jmax, kmax, vc, 2, X_DIR, BOTH );

    _PRINT_S4D(ofs,"after PeriodicCommX",imax,jmax,kmax,1,vc,n);

    // PeriodicCommY
    paraMngr->PeriodicCommS3D( pp, imax, jmax, kmax, vc, 2, Y_DIR, BOTH );

    _PRINT_S4D(ofs,"after PeriodicCommY",imax,jmax,kmax,1,vc,n);

    // PeriodicCommZ
    paraMngr->PeriodicCommS3D( pp, imax, jmax, kmax, vc, 2, Z_DIR, BOTH );

    _PRINT_S4D(ofs,"after PeriodicCommZ",imax,jmax,kmax,1,vc,n);
  }
#endif

  // BndComm, PeriodicComm(V3D)
#if 1
  {
    int rank = paraMngr->GetMyRankID();
    const int *sz = paraMngr->GetLocalVoxelSize();
    int imax = sz[0];
    int jmax = sz[1];
    int kmax = sz[2];
    int nmax = 3;
    int vc  = 3;
    size_t nw = size_t(imax+2*vc) * size_t(jmax+2*vc) * size_t(kmax+2*vc) * size_t(nmax);

    REAL_TYPE *pp = paraMngr->AllocRealV3D(vc);
    for( int n=0;n<nmax;n++ ){
    for( int k=0-vc;k<kmax+vc;k++ ){
    for( int j=0-vc;j<jmax+vc;j++ ){
    for( int i=0-vc;i<imax+vc;i++ ){
      pp[_IDX_V3D(i,j,k,n,imax,jmax,kmax,vc)] = -(rank+1);
    }}}}
    int cnt = 0;
    for( int n=0;n<nmax;n++ ){
    for( int k=0;k<kmax;k++ ){
    for( int j=0;j<jmax;j++ ){
    for( int i=0;i<imax;i++ ){
      cnt++;
      pp[_IDX_V3D(i,j,k,n,imax,jmax,kmax,vc)] = n + cnt * 10 + (rank+1)*10000;
    }}}}

    char fname[512];
    sprintf ( fname, "commV3D_%04d.log", rank );
    std::ofstream ofs( fname );
#ifdef _NOWAIT_TEST_
    ofs << "#### commV3D nowait test ####" << endl << endl;
    if( rank==0 ) cout << "#### commV3D nowait test ####" << endl << endl;
#else
    ofs << "#### commV3D test ####" << endl << endl;
    if( rank==0 ) cout << "#### commV3D test ####" << endl << endl;
#endif

    size_t mem = paraMngr->GetBndCommBufferSize(-1);
    ofs << "buff mem = " << cpm_Base::GetMemString(mem) << endl;

    int n=1;
    int www = 6;
    ofs << "n=" << n << endl;

    _PRINT_S4D(ofs,"before comm",imax,jmax,kmax,nmax,vc,n);

#ifndef _NOWAIT_TEST_
    paraMngr->BndCommV3D( pp, imax, jmax, kmax, vc, 2 );
#else
    MPI_Request req[12];
    paraMngr->BndCommV3D_nowait( pp, imax, jmax, kmax, vc, 2, req );
    paraMngr->wait_BndCommV3D( pp, imax, jmax, kmax, vc, 2, req );
#endif

    _PRINT_S4D(ofs,"after BndComm",nmax,imax,jmax,kmax,vc,n);

    // PeriodicCommX
    paraMngr->PeriodicCommV3D( pp, imax, jmax, kmax, vc, 2, X_DIR, BOTH );

    _PRINT_S4D(ofs,"after PeriodicCommX",imax,jmax,kmax,nmax,vc,n);

    // PeriodicCommY
    paraMngr->PeriodicCommV3D( pp, imax, jmax, kmax, vc, 2, Y_DIR, BOTH );

    _PRINT_S4D(ofs,"after PeriodicCommY",imax,jmax,kmax,nmax,vc,n);

    // PeriodicCommZ
    paraMngr->PeriodicCommV3D( pp, imax, jmax, kmax, vc, 2, Z_DIR, BOTH );

    _PRINT_S4D(ofs,"after PeriodicCommZ",imax,jmax,kmax,nmax,vc,n);
  }
#endif

  // BndComm, PeriodicComm(S4D)
#if 1
  {
    int rank = paraMngr->GetMyRankID();
    const int *sz = paraMngr->GetLocalVoxelSize();
    int imax = sz[0];
    int jmax = sz[1];
    int kmax = sz[2];
    int nmax = 5;
    int vc  = 3;
    size_t nw = size_t(imax+2*vc) * size_t(jmax+2*vc) * size_t(kmax+2*vc) * size_t(nmax);

    REAL_TYPE *pp = paraMngr->AllocRealS4D(nmax,vc);
    for( int n=0;n<nmax;n++ ){
    for( int k=0-vc;k<kmax+vc;k++ ){
    for( int j=0-vc;j<jmax+vc;j++ ){
    for( int i=0-vc;i<imax+vc;i++ ){
      pp[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = -(rank+1);
    }}}}
    int cnt = 0;
    for( int n=0;n<nmax;n++ ){
    for( int k=0;k<kmax;k++ ){
    for( int j=0;j<jmax;j++ ){
    for( int i=0;i<imax;i++ ){
      cnt++;
      pp[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = n + cnt * 10 + (rank+1)*10000;
    }}}}

    char fname[512];
    sprintf ( fname, "commS4D_%04d.log", rank );
    std::ofstream ofs( fname );
#ifdef _NOWAIT_TEST_
    ofs << "#### commS4D nowait test ####" << endl << endl;
    if( rank==0 ) cout << "#### commS4D nowait test ####" << endl << endl;
#else
    ofs << "#### commS4D test ####" << endl << endl;
    if( rank==0 ) cout << "#### commS4D test ####" << endl << endl;
#endif

    size_t mem = paraMngr->GetBndCommBufferSize(-1);
    ofs << "buff mem = " << cpm_Base::GetMemString(mem) << endl;

    int n=3;
    int www = 6;
    ofs << "n=" << n << endl;

    _PRINT_S4D(ofs,"before comm",imax,jmax,kmax,nmax,vc,n);

#ifndef _NOWAIT_TEST_
    paraMngr->BndCommS4D( pp, imax, jmax, kmax, nmax, vc, 2 );
#else
    MPI_Request req[12];
    paraMngr->BndCommS4D_nowait( pp, imax, jmax, kmax, nmax, vc, 2, req );
    paraMngr->wait_BndCommS4D( pp, imax, jmax, kmax, nmax, vc, 2, req );
#endif

    _PRINT_S4D(ofs,"after BndComm",nmax,imax,jmax,kmax,vc,n);

    // PeriodicCommX
    paraMngr->PeriodicCommS4D( pp, imax, jmax, kmax, nmax, vc, 2, X_DIR, BOTH );

    _PRINT_S4D(ofs,"after PeriodicCommX",imax,jmax,kmax,nmax,vc,n);

    // PeriodicCommY
    paraMngr->PeriodicCommS4D( pp, imax, jmax, kmax, nmax, vc, 2, Y_DIR, BOTH );

    _PRINT_S4D(ofs,"after PeriodicCommY",imax,jmax,kmax,nmax,vc,n);

    // PeriodicCommZ
    paraMngr->PeriodicCommS4D( pp, imax, jmax, kmax, nmax, vc, 2, Z_DIR, BOTH );

    _PRINT_S4D(ofs,"after PeriodicCommZ",imax,jmax,kmax,nmax,vc,n);
  }
#endif

  // BndComm, PeriodicComm(V3DEx)
#if 1
  {
    int rank = paraMngr->GetMyRankID();
    const int *sz = paraMngr->GetLocalVoxelSize();
    int imax = sz[0];
    int jmax = sz[1];
    int kmax = sz[2];
    int nmax = 3;
    int vc  = 3;
    size_t nw = size_t(imax+2*vc) * size_t(jmax+2*vc) * size_t(kmax+2*vc) * size_t(nmax);

    REAL_TYPE *pp = paraMngr->AllocRealV3DEx(vc);
    for( int k=0-vc;k<kmax+vc;k++ ){
    for( int j=0-vc;j<jmax+vc;j++ ){
    for( int i=0-vc;i<imax+vc;i++ ){
    for( int n=0;n<nmax;n++ ){
      pp[_IDX_V3DEX(n,i,j,k,imax,jmax,kmax,vc)] = -(rank+1);
    }}}}
    int cnt = 0;
    for( int k=0;k<kmax;k++ ){
    for( int j=0;j<jmax;j++ ){
    for( int i=0;i<imax;i++ ){
    for( int n=0;n<nmax;n++ ){
      cnt++;
      pp[_IDX_V3DEX(n,i,j,k,imax,jmax,kmax,vc)] = n + cnt * 10 + (rank+1)*10000;
    }}}}

    char fname[512];
    sprintf ( fname, "commV3DEx_%04d.log", rank );
    std::ofstream ofs( fname );
#ifdef _NOWAIT_TEST_
    ofs << "#### commV3DEx nowait test ####" << endl << endl;
    if( rank==0 ) cout << "#### commV3DEx nowait test ####" << endl << endl;
#else
    ofs << "#### commV3DEx test ####" << endl << endl;
    if( rank==0 ) cout << "#### commV3DEx test ####" << endl << endl;
#endif

    size_t mem = paraMngr->GetBndCommBufferSize(-1);
    ofs << "buff mem = " << cpm_Base::GetMemString(mem) << endl;

    int n=2;
    int www = 6;
    ofs << "n=" << n << endl;

    _PRINT_S4DEX(ofs,"before comm",nmax,imax,jmax,kmax,vc,n);

#ifndef _NOWAIT_TEST_
    paraMngr->BndCommV3DEx( pp, imax, jmax, kmax, vc, 2 );
#else
    MPI_Request req[12];
    paraMngr->BndCommV3DEx_nowait( pp, imax, jmax, kmax, vc, 2, req );
    paraMngr->wait_BndCommV3DEx( pp, imax, jmax, kmax, vc, 2, req );
#endif

    _PRINT_S4DEX(ofs,"after BndComm",nmax,imax,jmax,kmax,vc,n);

    // PeriodicCommX
    paraMngr->PeriodicCommV3DEx( pp, imax, jmax, kmax, vc, 2, X_DIR, BOTH );

    _PRINT_S4DEX(ofs,"after PeriodicCommX",nmax,imax,jmax,kmax,vc,n);

    // PeriodicCommY
    paraMngr->PeriodicCommV3DEx( pp, imax, jmax, kmax, vc, 2, Y_DIR, BOTH );

    _PRINT_S4DEX(ofs,"after PeriodicCommY",nmax,imax,jmax,kmax,vc,n);

    // PeriodicCommZ
    paraMngr->PeriodicCommV3DEx( pp, imax, jmax, kmax, vc, 2, Z_DIR, BOTH );

    _PRINT_S4DEX(ofs,"after PeriodicCommZ",nmax,imax,jmax,kmax,vc,n);
  }
#endif

  // BndComm, PeriodicComm(S4DEx)
#if 1
  {
    int rank = paraMngr->GetMyRankID();
    const int *sz = paraMngr->GetLocalVoxelSize();
    int imax = sz[0];
    int jmax = sz[1];
    int kmax = sz[2];
    int nmax = 6;
    int vc  = 3;
    size_t nw = size_t(imax+2*vc) * size_t(jmax+2*vc) * size_t(kmax+2*vc) * size_t(nmax);

    REAL_TYPE *pp = paraMngr->AllocRealS4DEx(nmax,vc);
    for( int k=0-vc;k<kmax+vc;k++ ){
    for( int j=0-vc;j<jmax+vc;j++ ){
    for( int i=0-vc;i<imax+vc;i++ ){
    for( int n=0;n<nmax;n++ ){
      pp[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)] = -(rank+1);
    }}}}
    int cnt = 0;
    for( int k=0;k<kmax;k++ ){
    for( int j=0;j<jmax;j++ ){
    for( int i=0;i<imax;i++ ){
    for( int n=0;n<nmax;n++ ){
      cnt++;
      pp[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)] = n + cnt * 10 + (rank+1)*10000;
    }}}}

    char fname[512];
    sprintf ( fname, "commS4DEx_%04d.log", rank );
    std::ofstream ofs( fname );
#ifdef _NOWAIT_TEST_
    ofs << "#### commS4DEx nowait test ####" << endl << endl;
    if( rank==0 ) cout << "#### commS4DEx nowait test ####" << endl << endl;
#else
    ofs << "#### commS4DEx test ####" << endl << endl;
    if( rank==0 ) cout << "#### commS4DEx test ####" << endl << endl;
#endif

    size_t mem = paraMngr->GetBndCommBufferSize(-1);
    ofs << "buff mem = " << cpm_Base::GetMemString(mem) << endl;

    int n=4;
    int www = 6;
    ofs << "n=" << n << endl;

    _PRINT_S4DEX(ofs,"before comm",nmax,imax,jmax,kmax,vc,n);

#ifndef _NOWAIT_TEST_
    paraMngr->BndCommS4DEx( pp, nmax, imax, jmax, kmax, vc, 2 );
#else
    MPI_Request req[12];
    paraMngr->BndCommS4DEx_nowait( pp, nmax, imax, jmax, kmax, vc, 2, req );
    paraMngr->wait_BndCommS4DEx( pp, nmax, imax, jmax, kmax, vc, 2, req );
#endif

    _PRINT_S4DEX(ofs,"after BndComm",nmax,imax,jmax,kmax,vc,n);

    // PeriodicCommX
    paraMngr->PeriodicCommS4DEx( pp, nmax, imax, jmax, kmax, vc, 2, X_DIR, BOTH );

    _PRINT_S4DEX(ofs,"after PeriodicCommX",nmax,imax,jmax,kmax,vc,n);

    // PeriodicCommY
    paraMngr->PeriodicCommS4DEx( pp, nmax, imax, jmax, kmax, vc, 2, Y_DIR, BOTH );

    _PRINT_S4DEX(ofs,"after PeriodicCommY",nmax,imax,jmax,kmax,vc,n);

    // PeriodicCommZ
    paraMngr->PeriodicCommS4DEx( pp, nmax, imax, jmax, kmax, vc, 2, Z_DIR, BOTH );

    _PRINT_S4DEX(ofs,"after PeriodicCommZ",nmax,imax,jmax,kmax,vc,n);
  }
#endif

