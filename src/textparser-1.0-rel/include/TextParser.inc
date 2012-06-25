! -*- mode: Fortran -*-
!****************************************************************************
!**
!** Copyright (C) 2012 Tokyo University.
!**
!****************************************************************************
!** @file TextParser.inc
! * ここには TextParser ライブラリのFortran用インターフェースが定義されています。
! *
! *

integer TP_READ
integer TP_WRITE
integer TP_REMOVE
integer TP_GET_NUMBER_OF_LEAVES
integer TP_GET_LABEL
integer TP_GET_VALUE
integer TP_GET_TYPE
integer*1 TP_CONVERT_CHAR
integer*2 TP_CONVERT_SHORT
integer*4 TP_CONVERT_INT
integer*4 TP_CONVERT_LONG
integer*1 TP_CONVERT_INT1
integer*2 TP_CONVERT_INT2
integer*4 TP_CONVERT_INT4
integer*8 TP_CONVERT_INT8
integer*8 TP_CONVERT_LONG_LONG
real TP_CONVERT_FLOAT
real*8 TP_CONVERT_DOUBLE
logical TP_CONVERT_LOGICAL
integer TP_GET_NUMBER_OF_ELEMENTS
integer TP_GET_ITH_ELEMENT
integer TP_CURRENT_NODE
integer TP_CHANGE_NODE
integer TP_GET_NUMBER_OF_CNODES
integer TP_GET_NUMBER_OF_CLEAVES
integer TP_GET_ITH_NODE
integer TP_GET_ITH_LEAF
integer TP_GET_ITH_NODE_ORDER
integer TP_GET_ITH_LEAF_ORDER

external TP_READ
external TP_write
external TP_REMOVE
external TP_GET_NUMBER_OF_LEAVES
external TP_GET_LABEL
external TP_GET_VALUE
external TP_GET_TYPE
external TP_CONVERT_CHAR
external TP_CONVERT_SHORT
external TP_CONVERT_INT
external TP_CONVERT_LONG
external TP_CONVERT_LONG_LONG
external TP_CONVERT_FLOAT
external TP_CONVERT_DOUBLE
external TP_CONVERT_LOGICAL
external TP_GET_NUMBER_OF_ELEMENTS
external TP_GET_ITH_ELEMENT
external TP_CURRENT_NODE
external TP_CHANGE_NODE
external TP_GET_NUMBER_OF_CNODES
external TP_GET_NUMBER_OF_CLEAVES
external TP_GET_ITH_NODE
external TP_GET_ITH_LEAF
external TP_GET_ITH_NODE_ORDER
external TP_GET_ITH_LEAF_ORDER

external TP_CONVERT_INT1
external TP_CONVERT_INT2
external TP_CONVERT_INT4
external TP_CONVERT_INT8
