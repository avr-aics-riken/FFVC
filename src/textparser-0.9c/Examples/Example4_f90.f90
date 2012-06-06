      program Example4

      implicit none 



      integer move_and_get_parameters
      integer status
  
      status = move_and_get_parameters('Input0-1.txt')
      status = move_and_get_parameters('Input4-1.txt')
  
      end program Example4


!! move_and_get_parameters
!      integer function move_and_get_parameters(filename)
      recursive function move_and_get_parameters(filename) RESULT(status)
      character(*) filename
      character*5 path
      integer status

      include 'TextParser.inc'

      write(6,*) 'input filename: ', trim(filename)
      status = TP_READ(filename)
	  if (status /= 0) then
          write(6,*) 'TP_READ   status: ', status
	  end if

      path = '/'
      call get_directory_parameters(filename, len_trim(path), trim(path))

      status = TP_REMOVE()  
      if (status /= 0) then
         write(6,*) 'TP_REMOVE status: ', status
      end if

      end function


!! get_directory_parameters
      recursive subroutine get_directory_parameters(filename, l_label, label)
      character(*) filename, label
	  integer l_label
	  integer dir_number
	  integer parm_number
	  integer l_dir
	  character*512 dir_label
	  integer l_parm
	  character*512 parm_label
      character(80) value
	  integer l_value
	  integer value_type
      integer status

      include 'TextParser.inc'

      write(6,*) filename,l_label,label


      status = TP_CHANGE_NODE(label)
	  if (status /= 0) then
          write(6,*) 'TP_CHANGE_NODE  status: ', status
	  end if

      status = TP_CURRENT_NODE(label)
	  if (status /= 0) then
          write(6,*) 'TP_CURRENT_NODE   status: ', status
	  end if
      write(6,*) 'Current directory: ', trim(label)

      status = TP_GET_NUMBER_OF_CNODES(dir_number)
	  if (status /= 0) then
          write(6,*) 'TP_GET_NUMBER_OF_CNODES   status: ', status
	  end if
          write(6,*) 'Number of (daughter) Nodes in Current node.',dir_number
          if(dir_number.gt.0) then
             do i = 1, dir_number
                dir_label=''
                write(6,*) 'initialized ',trim(dir_label)
                status = TP_GET_ITH_NODE(i, dir_label)
                write(6,*) 'initialized 2',trim(dir_label)
                if (status /= 0) then
                   write(6,*) 'TP_GET_ITH_NODE   status: ', status
                end if
                call get_directory_parameters(filename, len_trim(dir_label),trim(dir_label))
             end do
	  endif

             status = TP_GET_NUMBER_OF_CLEAVES(parm_number)
	  if (status /= 0) then
          write(6,*) 'TP_GET_NUMBER_OF_CLEAVES   status: ', status
	  end if
          write(6,*) 'Number of Leaves in Current node.',parm_number
          if(parm_number.gt.0) then
	   do i = 1, parm_number 
           parm_label=''
           value=''
           status = TP_GET_ITH_LEAF(i, parm_label)
           if (status /= 0) then
              write(6,*) 'TP_GET_ITH_LEAF   status: ', status
	       end if
               write(6,*) i, ' ', trim(parm_label)
	 	  
          status = TP_GET_VALUE(parm_label, value)
		  if (status /= 0) then
              write(6,*) 'TP_GET_VALUE   status: ', status
		  end if
          write(6,*) i ,' ', trim(value)

          status = TP_GET_TYPE(parm_label, value_type)
		  if (status /= 0) then
              write(6,*) 'TP_GET_TYPE   status: ', status
		  end if
           write(6,*) i ,' type: ', value_type
	   end do
          endif
          status = TP_CURRENT_NODE(label)
	  if (status /= 0) then
          write(6,*) 'TP_CURRENT_NODE   status: ', status
	  end if
	  if (trim(label) /= '/') then
          label = '..'
          status = TP_CHANGE_NODE(label)
	      if (status /= 0) then
              write(6,*) 'TP_CHANGE_NODE   status: ', status
	      end if

          status = TP_Current_NODE(label)
	      if (status /= 0) then
              write(6,*) 'TP_CURRENT_NODE   status: ', status
	      end if
          write(6,*) 'Current directory: ', trim(label)
	  end if
		  
      end subroutine get_directory_parameters
