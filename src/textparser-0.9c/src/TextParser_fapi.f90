! -*- mode: Fortran -*-
!***************************************************************************!**
!** Copyright (C) 2012 Tokyo University.
!**
!***************************************************************************
!** @file TextParser.f90
!* ここには TextParser ライブラリのFortran用インターフェースが実装されています。


      integer function TP_READ(file)
        integer TP_READ_FORT
        external TP_READ_FORT
        character(LEN=*),INTENT(IN) :: file
        CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
        integer length
!        length=trim_len(file)
        TP_READ = TP_READ_FORT(trim(file)//NULL)
      
       return

      end function 

      integer function TP_WRITE(file)
        integer TP_WRITE_FORT
        external TP_WRITE_FORT
        character(LEN=*),INTENT(IN) :: file
        CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
        integer length
!        length=trim_len(file)
        TP_WRITE = TP_WRITE_FORT(trim(file)//NULL)
      
       return
       end function 


       integer function TP_REMOVE()
       
        integer TP_REMOVE_FORT
        external TP_REMOVE_FORT

        TP_REMOVE = TP_REMOVE_FORT()
       
        return
       
       end function 

      integer function TP_GET_NUMBER_OF_LEAVES(nleaves)
        integer,intent(out) :: nleaves

        integer TP_GET_NUMBER_OF_LEAVES_FORT
        external TP_GET_NUMBER_OF_LEAVES_FORT
!        write(6,*) 'TP_GET_NUMBER_OF_LEAVES',nleaves
        TP_GET_NUMBER_OF_LEAVES=TP_GET_NUMBER_OF_LEAVES_FORT(nleaves)
!        write(6,*) 'TP_GET_NUMBER_OF_LEAVES',nleaves
        return

      end function TP_GET_NUMBER_OF_LEAVES

      integer function TP_GET_LABEL(ileaves,label)
        integer,intent(in) :: ileaves
        character(len=*),intent(out) ::label
        integer TP_GET_LABEL_FORT
        external TP_GET_LABEL_FORT
        integer ileaves_shift
       
        ileaves_shift=ileaves-1
        TP_GET_LABEL=TP_GET_LABEL_FORT(ileaves_shift,label)
        
        return

      end function TP_GET_LABEL

      integer function TP_GET_VALUE(label,value)
        character(len=*),intent(in) ::label
        character(len=*),intent(out) ::value
        integer TP_GET_VALUE_FORT
        external TP_GET_VALUE_FORT
        CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
!        write(6,*) 'TP_GET_VALUE label|', trim(label),'|'
!        write(6,*) 'TP_GET_VALUE value|', trim(value),'|'
        TP_GET_VALUE=TP_GET_VALUE_FORT(trim(label)//NULL,value)
!        write(6,*) 'tp_get_value',TP_GET_VALUE
        
        return

      end function TP_GET_VALUE

      integer function TP_GET_TYPE(label,value)
        character(len=*),intent(in) ::label
        character(len=*),intent(out) ::value
        integer TP_GET_TYPE_FORT
        external TP_GET_TYPE_FORT
        CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)       

        TP_GET_TYPE=TP_GET_TYPE_FORT(trim(label)//NULL,value)
        
        return

      end function TP_GET_TYPE

      integer*1 function TP_CONVERT_CHAR(value,ierror)
       character(len=*),intent(in) :: value
       integer ierror
       integer*1 TP_CONVERT_CHAR_FORT
       external  TP_CONVERT_CHAR_FORT
        CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)

       TP_CONVERT_CHAR=TP_CONVERT_CHAR_FORT(trim(value)//NULL,ierror)
       return
      end function TP_CONVERT_CHAR

      integer*2 function TP_CONVERT_SHORT(value,ierror)
       character(len=*),intent(in) :: value
        CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
       integer ierror
       integer*2 TP_CONVERT_SHORT_FORT
       external  TP_CONVERT_SHORT_FORT

       TP_CONVERT_SHORT=TP_CONVERT_SHORT_FORT(trim(value)//NULL,ierror)
       return
      end function TP_CONVERT_SHORT

      integer function TP_CONVERT_INT(value,ierror)
       character(len=*),intent(in) :: value
        CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
       integer ierror
       integer TP_CONVERT_INT_FORT
       external  TP_CONVERT_INT_FORT

       TP_CONVERT_INT=TP_CONVERT_INT_FORT(trim(value)//NULL,ierror)
       return
      end function TP_CONVERT_INT

      integer*8 function TP_CONVERT_LONG(value,ierror)
       character(len=*),intent(in) :: value
        CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
       integer*8 ierror
       integer TP_CONVERT_LONG_FORT
       external  TP_CONVERT_LONG_FORT

       TP_CONVERT_LONG=TP_CONVERT_LONG_FORT(trim(value)//NULL,ierror)
       return
      end function TP_CONVERT_LONG

      integer*8 function TP_CONVERT_LONG_LONG(value,ierror)
       character(len=*),intent(in) :: value
        CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
       integer ierror
       integer*8 TP_CONVERT_LONG_LONG_FORT
       external  TP_CONVERT_LONG_LONG_FORT

       TP_CONVERT_LONG_LONG=TP_CONVERT_LONG_LONG_FORT(trim(value)//NULL,ierror)
       return
      end function TP_CONVERT_LONG_LONG

      real function TP_CONVERT_FLOAT(value,ierror)
       character(len=*),intent(in) :: value
       CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
       integer ierror
       real TP_CONVERT_FLOAT_FORT
       external  TP_CONVERT_FLOAT_FORT

       TP_CONVERT_INT=TP_CONVERT_FLOAT_FORT(trim(value)//NULL,ierror)
       return
      end function TP_CONVERT_FLOAT

      real*8 function TP_CONVERT_DOUBLE(value,ierror)
       character(len=*),intent(in) :: value
       CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
       integer ierror
       real*8 TP_CONVERT_DOUBLE_FORT
       external  TP_CONVERT_DOUBLE_FORT

       TP_CONVERT_INT=TP_CONVERT_DOUBLE_FORT(trim(value)//NULL,ierror)
       return
      end function TP_CONVERT_DOUBLE
      

      logical function TP_CONVERT_LOGICAL(value,ierror)
       character(len=*),intent(in) :: value
       CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
       integer ierror
       integer TP_CONVERT_LOGICAL_FORT
       external  TP_CONVERT_LOGICAL_FORT
       integer tmp
       tmp=TP_CONVERT_FLOAT_FORT(trim(value)//NULL,ierror)
       TP_CONVERT_LOGICAL=.false.
       if(tmp.eq.1) TP_CONVERT_LOGICAL=.true.
       return
      end function TP_CONVERT_LOGICAL

      integer function TP_GET_NUMBER_OF_ELEMENTS(vector_value,nvec)
      CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
      character(len=*),intent(in)::vector_value

      integer TP_GET_NUMBER_OF_ELEMENTS_FORT
      external TP_GET_NUMBER_OF_ELEMENTS_FORT

      TP_GET_NUMBER_OF_ELEMENTS=TP_GET_NUMBER_OF_ELEMENTS_FORT(trim(vector_value)//NULL,nvec)

      return

      end function TP_GET_NUMBER_OF_ELEMENTS

      integer function TP_GET_ITH_ELEMENT(vector_value,ivec,velem)

      CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
      character(len=*),intent(in)::vector_value
      character(len=*),intent(out)::velem
      integer TP_GET_ITH_ELEMENT_FORT
      external TP_GET_ITH_ELEMENT_FORT
      integer ivec_shift
      ivec_shift=ivec-1
      TP_GET_ITH_ELEMENT=TP_GET_ITH_ELEMENT_FORT(trim(vector_value)//NULL,ivec_shift,velem)
      return

      end function TP_GET_ITH_ELEMENT

      integer function TP_CURRENT_NODE(label)
      character(len=*),intent(in)::label
      CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
      integer TP_CURRENT_NODE_FORT
      external TP_CURRENT_NODE_FORT

      TP_CURRENT_NODE=TP_CURRENT_NODE_FORT(trim(label)//NULL)

      return 
      end function TP_CURRENT_NODE

      integer function TP_CHANGE_NODE(label)
      character(len=*),intent(in)::label
      CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
      integer TP_CHANGE_NODE_FORT
      external TP_CHANGE_NODE_FORT
      TP_CHANGE_NODE=TP_CHANGE_NODE_FORT(trim(label)//NULL)
      return 
      end function TP_CHANGE_NODE

      integer function TP_GET_NUMBER_OF_CNODES(inode)
      integer,intent(in) :: inode
      integer TP_GET_NUMBER_OF_CNODES_FORT
      external TP_GET_NUMBER_OF_CNODES_FORT

      TP_GET_NUMBER_OF_CNODES=TP_GET_NUMBER_OF_CNODES_FORT(inode)
      return

      end function TP_GET_NUMBER_OF_CNODES

      integer function TP_GET_NUMBER_OF_CLEAVES(ileaf)
      integer,intent(OUT) :: ileaf
      integer TP_GET_NUMBER_OF_CLEAVES_FORT
      external TP_GET_NUMBER_OF_CLEAVES_FORT

      TP_GET_NUMBER_OF_CLEAVES=TP_GET_NUMBER_OF_CLEAVES_FORT(ileaf)
      return

      end function TP_GET_NUMBER_OF_CLEAVES
      

      integer function TP_GET_ITH_NODE(inode,node)
      integer,intent(in) :: inode
      character(len=*),intent(OUT) :: node
      integer TP_GET_ITH_NODE_FORT
      external TP_GET_ITH_NODE_FORT
      integer iinode
      iinode=inode-1
!      write(6,*) 'TP_GET_ITH_NODE',iinode,inode,trim(node)
      TP_GET_ITH_NODE=TP_GET_ITH_NODE_FORT(iinode,node)
!      write(6,*) 'TP_GET_ITH_NODE',iinode,inode,trim(node)
      return
      end function TP_GET_ITH_NODE


      integer function TP_GET_ITH_NODE_ORDER(inode,node,order)
      integer,intent(in) :: inode,order
      character(len=*),intent(OUT) :: node
      integer TP_GET_ITH_NODE_FORT
      external TP_GET_ITH_NODE_FORT
      integer iinode
      iinode=inode-1
!      write(6,*) 'TP_GET_ITH_NODE',iinode,inode,trim(node)
      TP_GET_ITH_NODE=TP_GET_ITH_NODE_ORDER_FORT(iinode,node,order)
!      write(6,*) 'TP_GET_ITH_NODE',iinode,inode,trim(node)
      return
      end function TP_GET_ITH_NODE_ORDER


      integer function TP_GET_ITH_LEAF(ileaf,leaf)
      integer,intent(in) :: ileaf
      character(len=*),intent(OUT) :: leaf
      integer TP_GET_ITH_LEAF_FORT
      external TP_GET_ITH_LEAF_FORT
      integer iileaf
      iileaf=ileaf-1
      TP_GET_ITH_LEAF=TP_GET_ITH_LEAF_FORT(iileaf,leaf)
      return
      end function TP_GET_ITH_LEAF

      integer function TP_GET_ITH_LEAF_ORDER(ileaf,leaf,order)
      integer,intent(in) :: ileaf,order
      character(len=*),intent(OUT) :: leaf
      integer TP_GET_ITH_LEAF_FORT
      external TP_GET_ITH_LEAF_FORT
      integer iileaf
      iileaf=ileaf-1
      TP_GET_ITH_LEAF=TP_GET_ITH_LEAF_ORDER_FORT(iileaf,leaf,order)
      return
      end function TP_GET_ITH_LEAF_ORDER

      integer function TP_CONVERT_INT1(value,ierror)
      integer*1 TP_CONVERT_CHAR
      character(len=*),intent(in)::value
      integer ierror
      TP_CONVERT_INT1=TP_CONVERT_CHAR(value,ierror)
      return 
      end function TP_CONVERT_INT1

      integer function TP_CONVERT_INT2(value,ierror)
      integer*2 TP_CONVERT_SHORT
      character(len=*),intent(in)::value
      integer ierror
      TP_CONVERT_INT2=TP_CONVERT_SHORT(value,ierror)
      return 
      end function TP_CONVERT_INT2

      integer function TP_CONVERT_INT4(value,ierror)
      integer TP_CONVERT_INT
      character(len=*),intent(in)::value
      integer ierror
      TP_CONVERT_INT4=TP_CONVERT_INT(value,ierror)
      return 
      end function TP_CONVERT_INT4

      integer function TP_CONVERT_INT8(value,ierror)
      integer*8 TP_CONVERT_LONG_LONG
      character(len=*),intent(in)::value
      integer ierror
      TP_CONVERT_INT8=TP_CONVERT_LONG_LONG(value,ierror)
      return 
      end function TP_CONVERT_INT8
