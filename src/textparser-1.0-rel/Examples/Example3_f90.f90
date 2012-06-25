      program Example3

      implicit none 


      integer scan_all_parameters
      integer status
  
      status = scan_all_parameters('Input0-1.txt')
      status = scan_all_parameters('Input1-1.txt')
      status = scan_all_parameters('Input4-1.txt')
      status = scan_all_parameters('Input4-2.txt')
  
      end program Example3

      integer function scan_all_parameters(filename)
      character(*) filename
      integer status
	  integer number
      character*512 label
	  integer l_label
      character*80 value
	  integer l_value
      character*80 val
	  integer l_val
	  integer value_type
	  integer*1 cval
	  integer*2 sval
	  integer ival
	  integer*8 llval
	  real fval
	  real*8 dval
	  integer vector_number
   include 'TextParser.inc'

      write(6,*) 'input filename: ', trim(filename)
      status = TP_READ(filename)
	  if (status /= 0) then
          write(6,*) 'TP_READ   status: ', status
	  end if

      status = TP_GET_NUMBER_OF_LEAVES(number)
	  if (status /= 0) then
          write(6,*) 'TP_GET_NUMBER_OF_LEAVES   status: ', status
	  end if
      write(6,*) 'number: ', number
	  
	  do i = 1, number 
           label=''
           value=''
!           write(6,*)'initilaized label|',trim(label),'|'
!           write(6,*)'initilaized value|',trim(value),'|'

          status = TP_GET_LABEL(i, label)
		  if (status /= 0) then
              write(6,*) 'TP_GET_LABEL   status: ', status
		  end if
          write(6,*) i ,' ', trim(label)
!          write(6,*) i ,' ', label(1:len_trim(label)-1)

          status = TP_GET_VALUE(label, value)
          write(6,*) 'TP_GET_VALUE  status: ', status
		  if (status /= 0) then
              write(6,*) 'TP_GET_VALUE  status: ', status
		  end if
          write(6,*) i ,' ', trim(value)

          status = TP_GET_TYPE(label, value_type)
		  if (status /= 0) then
              write(6,*) 'TP_GET_TYPE   status: ', status
		  end if
          write(6,*) i ,' value type: ', value_type

          if (value_type == 1) then
              ival = TP_CONVERT_INT(value,status);
		      if (status /= 0) then
                  write(6,*) 'TP_CONVERT_INT   status: ', status
		      end if
              write(6,*) i ,' convert to int: ', ival
              fval = TP_CONVERT_FLOAT(value,status);
		      if (status /= 0) then
           write(6,*) 'TP_CONVERT_FLOAT   status: ', status
		      end if
              write(6,*) i ,' convert to float: ', fval
              dval = TP_CONVERT_DOUBLE(value,status);
		      if (status /= 0) then
                  write(6,*) 'TP_CONVERT_DOUBLE   status: ', status
		      end if
              write(6,*) i ,' convert to double: ', dval
          else if (value_type == 5) then
              status  = TP_GET_NUMBER_OF_ELEMENTS(value, vector_number);
		      if (status /= 0) then
                  write(6,*) 'TP_GET_NUMBER_OF_ELEMENTS   status: ', status
		      end if
        do j = 1, vector_number 
                  val=''
                  status = TP_GET_ITH_ELEMENT(value, j, val)
		          if (status /= 0) then
                      write(6,*) 'TP_GET_ITH_ELEMENTS   status: ', status
		          end if
                    ival = TP_CONVERT_INT(val,status);
		          if (status /= 0) then
                      write(6,*) 'TO_CONVERT_INT   status: ', status
		          end if
                    fval = TP_CONVERT_FLOAT(val,status);
		          if (status /= 0) then
                      write(6,*) 'TP_CONVERT_FLOAT   status: ', status
		          end if
                    dval = TP_CONVERT_DOUBLE(val,status);
		          if (status /= 0) then
                      write(6,*) 'TP_CONVERT_DOUBLE  status: ', status
		          end if
                  write(6,*) j, ' ', trim(val), ' int ', ival, ' float ', fval, ' double ', dval
	          end do
          else if (value_type == 6) then
              status  = TP_GET_NUMBER_OF_ELEMENTS(value, vector_number)
		      if (status /= 0) then
                  write(6,*) 'TP_GET_NUMBER_OF_ELEMENTS   status: ', status
		      end if
              write(6,*) 'vector_number: ', vector_number
	          do j = 1, vector_number 
                  val=''              
                  status = TP_GET_ITH_ELEMENT(value, j, val);
		          if (status /= 0) then
                      write(6,*) 'TP_GET_ITH_ELEMENT   status: ', status
		          end if
                  write(6,*) j, ' ', trim(val)
	          end do
		  end if
	  end do

          status = TP_REMOVE()
	  if (status /= 0) then
          write(6,*) 'Remove status: ', status
	  end if

      end function
