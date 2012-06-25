      program Example1

      implicit none 


      integer read_write_parameters
      integer status
  
    status= read_write_parameters("correct_basic_1.txt","correct_basic_1_out.txt")
    status= read_write_parameters("correct_basic_2.txt","correct_basic_2_out.txt")
    status= read_write_parameters("correct_basic_3.txt","correct_basic_3_out.txt")
    status= read_write_parameters("correct_basic_4.txt","correct_basic_4_out.txt")
    status= read_write_parameters("correct_basic_5.txt","correct_basic_5_out.txt")
    status= read_write_parameters("correct_basic_6.txt","correct_basic_6_out.txt")
    status= read_write_parameters("correct_basic_7.txt","correct_basic_7_out.txt")
    status= read_write_parameters("correct_basic_8.txt","correct_basic_8_out.txt")
    status= read_write_parameters("correct_basic_9.txt","correct_basic_9_out.txt")
    status= read_write_parameters("correct_basic_10.txt","correct_basic_10_out.txt")
    status= read_write_parameters("correct_basic_11.txt","correct_basic_11_out.txt")
    status= read_write_parameters("correct_basic_12.txt","correct_basic_12_out.txt")
    status= read_write_parameters("correct_basic_13.txt","correct_basic_13_out.txt")

    status= read_write_parameters("correct_string_1.txt","correct_string_1_out.txt")
    status= read_write_parameters("correct_string_2.txt","correct_string_2_out.txt")
    status= read_write_parameters("correct_string_3.txt","correct_string_3_out.txt")
    status= read_write_parameters("correct_string_4.txt","correct_string_4_out.txt")

    status= read_write_parameters("correct_label_1.txt","correct_label_1_out.txt")
    status= read_write_parameters("correct_label_2.txt","correct_label_2_out.txt")
    status= read_write_parameters("correct_label_3.txt","correct_label_3_out.txt")
    status= read_write_parameters("correct_label_4.txt","correct_label_4_out.txt")
    status= read_write_parameters("correct_label_5.txt","correct_label_5_out.txt")
    status= read_write_parameters("correct_label_6.txt","correct_label_6_out.txt")
    status= read_write_parameters("correct_label_7.txt","correct_label_7_out.txt")
    status= read_write_parameters("correct_label_8.txt","correct_label_8_out.txt")

    status= read_write_parameters("correct_cond_1.txt","correct_cond_1_out.txt")
    status= read_write_parameters("correct_cond_2.txt","correct_cond_2_out.txt")
    status= read_write_parameters("correct_cond_3.txt","correct_cond_3_out.txt")
    status= read_write_parameters("correct_cond_4.txt","correct_cond_4_out.txt")
    status= read_write_parameters("correct_cond_5.txt","correct_cond_5_out.txt")
    status= read_write_parameters("correct_cond_6.txt","correct_cond_6_out.txt")
    status= read_write_parameters("correct_cond_7.txt","correct_cond_7_out.txt")
    status= read_write_parameters("correct_cond_8.txt","correct_cond_8_out.txt")
    status= read_write_parameters("correct_cond_9.txt","correct_cond_9_out.txt")
    status= read_write_parameters("correct_cond_10.txt","correct_cond_10_out.txt")

    status= read_write_parameters("correct_cond_11.txt","correct_cond_11_out.txt")
    status= read_write_parameters("correct_cond_12.txt","correct_cond_12_out.txt")
    status= read_write_parameters("correct_cond_13.txt","correct_cond_13_out.txt")
    status= read_write_parameters("correct_cond_14.txt","correct_cond_14_out.txt")
    status= read_write_parameters("correct_cond_15.txt","correct_cond_15_out.txt")
    status= read_write_parameters("correct_cond_16.txt","correct_cond_16_out.txt")
    status= read_write_parameters("correct_cond_17.txt","correct_cond_17_out.txt")
    status= read_write_parameters("correct_cond_18.txt","correct_cond_18_out.txt")
    status= read_write_parameters("correct_cond_19.txt","correct_cond_19_out.txt")
    status= read_write_parameters("correct_cond_20.txt","correct_cond_20_out.txt")

    status= read_write_parameters("correct_labelarray_1.txt","correct_labelarray_1_out.txt")
    status= read_write_parameters("correct_labelarray_2.txt","correct_labelarray_2_out.txt")
    status= read_write_parameters("correct_labelarray_3.txt","correct_labelarray_3_out.txt")
    status= read_write_parameters("correct_labelarray_4.txt","correct_labelarray_4_out.txt")
  
      end program Example1

      integer function read_write_parameters(ifilename, ofilename)
      character(*) ifilename, ofilename
      integer status
      include 'TextParser.inc'

      write(6,*) 'input filename: ', trim(ifilename)
      status = TP_READ(ifilename)
      write(6,*) 'Read   status: ', status

      status = TP_WRITE(ofilename)  
      write(6,*) 'Write  status: ', status

      status = TP_REMOVE()  
      write(6,*) 'Remove status: ', status

      end function
