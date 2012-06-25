      program Example2

      implicit none 


      integer read_parameters
      integer status



    status= read_parameters("incorrect_basic_1.txt","incorrect_basic_1_out.txt")
    status= read_parameters("incorrect_basic_2.txt","incorrect_basic_2_out.txt")
    status= read_parameters("incorrect_basic_3.txt","incorrect_basic_3_out.txt")
    status= read_parameters("incorrect_basic_4.txt","incorrect_basic_4_out.txt")
    status= read_parameters("incorrect_basic_5.txt","incorrect_basic_5_out.txt")
    status= read_parameters("incorrect_basic_6.txt","incorrect_basic_6_out.txt")
    status= read_parameters("incorrect_basic_7.txt","incorrect_basic_7_out.txt")
    status= read_parameters("incorrect_basic_8.txt","incorrect_basic_8_out.txt")
    status= read_parameters("incorrect_basic_9.txt","incorrect_basic_9_out.txt")
    status= read_parameters("incorrect_basic_10.txt","incorrect_basic_10_out.txt")
    status= read_parameters("incorrect_basic_11.txt","incorrect_basic_11_out.txt")
    status= read_parameters("incorrect_basic_12.txt","incorrect_basic_12_out.txt")
    status= read_parameters("incorrect_basic_13.txt","incorrect_basic_13_out.txt")
    status= read_parameters("incorrect_basic_14.txt","incorrect_basic_14_out.txt")
    status= read_parameters("incorrect_basic_15.txt","incorrect_basic_15_out.txt")
    status= read_parameters("incorrect_basic_16.txt","incorrect_basic_16_out.txt")
    status= read_parameters("incorrect_basic_17.txt","incorrect_basic_17_out.txt")
    status= read_parameters("incorrect_basic_18.txt","incorrect_basic_18_out.txt")
    status= read_parameters("incorrect_basic_19.txt","incorrect_basic_19_out.txt")
    status= read_parameters("incorrect_cond_1.txt","incorrect_cond_1_out.txt")
    status= read_parameters("incorrect_cond_2.txt","incorrect_cond_2_out.txt")
    status= read_parameters("incorrect_cond_3.txt","incorrect_cond_3_out.txt")
    status= read_parameters("incorrect_cond_4.txt","incorrect_cond_4_out.txt")
    status= read_parameters("incorrect_cond_5.txt","incorrect_cond_5_out.txt")
    status= read_parameters("incorrect_label_1.txt","incorrect_label_1_out.txt")
    status= read_parameters("incorrect_label_2.txt","incorrect_label_2_out.txt")
!    status= read_parameters("incorrect_label_3.txt","incorrect_label_3_out.txt")
!    status= read_parameters("incorrect_label_4.txt","incorrect_label_4_out.txt")
!    status= read_parameters("incorrect_label_5.txt","incorrect_label_5_out.txt")
    status= read_parameters("incorrect_label_6.txt","incorrect_label_6_out.txt")
    status= read_parameters("incorrect_label_7.txt","incorrect_label_7_out.txt")
    status= read_parameters("incorrect_label_8.txt","incorrect_label_8_out.txt")
    status= read_parameters("incorrect_label_9.txt","incorrect_label_9_out.txt")
    status= read_parameters("incorrect_label_10.txt","incorrect_label_10_out.txt")
    status= read_parameters("incorrect_label_11.txt","incorrect_label_11_out.txt")
    status= read_parameters("incorrect_label_12.txt","incorrect_label_12_out.txt")
    status= read_parameters("incorrect_label_13.txt","incorrect_label_13_out.txt")
    status= read_parameters("incorrect_label_14.txt","incorrect_label_14_out.txt")
    status= read_parameters("incorrect_label_15.txt","incorrect_label_15_out.txt")
    status= read_parameters("incorrect_label_16.txt","incorrect_label_16_out.txt")
!    status= read_parameters("incorrect_label_17.txt","incorrect_label_17_out.txt")
    status= read_parameters("incorrect_labelarray_1.txt","incorrect_labelarray_1_out.txt")
    status= read_parameters("incorrect_labelarray_2.txt","incorrect_labelarray_2_out.txt")
    status= read_parameters("incorrect_labelarray_3.txt","incorrect_labelarray_3_out.txt")

  
      end program Example2

      integer function read_parameters(ifilename, ofilename)
      character(*) ifilename, ofilename
      integer status
      include 'TextParser.inc'

      write(6,*) 'input filename: ', trim(ifilename)
      status = TP_READ(ifilename)
      write(6,*) 'Read   status: ', status

      status = TP_REMOVE()  
      write(6,*) 'Remove status: ', status

      end function
