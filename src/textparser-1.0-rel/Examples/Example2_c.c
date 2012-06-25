/** @file Example2_c.c
 * サンプルプログラム
 */

#include <stdio.h>
#include <stdlib.h>
#include "TextParser.h"

void read_parameters(char *ifilename, char *ofilename)
{
    int stat = 0;

    // ファイルの読み込み
    printf(">>> ファイルの読み込み(%s)\n", ifilename);
    stat=tp_read(ifilename);
    if (stat == 0) {
        printf("# ファイルの読み込みが正常に終了しました : %s\n", ifilename);
    } else {
        printf("# ファイルの読み込み中にエラーが発生しました : C%d\n", stat);
    }


    // 
    stat=tp_remove();
    if (stat == 0) {
        printf("# データベースの削除が正常に終了しました\n");
    } else {
        printf("# データベースの削除中にエラーが発生しました : C%d\n", stat);
    }
    printf("\n");
}

int main(int argc, char* argv[])
{

    read_parameters("incorrect_basic_1.txt","incorrect_basic_1_out.txt");
    read_parameters("incorrect_basic_2.txt","incorrect_basic_2_out.txt");
    read_parameters("incorrect_basic_3.txt","incorrect_basic_3_out.txt");
    read_parameters("incorrect_basic_4.txt","incorrect_basic_4_out.txt");
    read_parameters("incorrect_basic_5.txt","incorrect_basic_5_out.txt");
    read_parameters("incorrect_basic_6.txt","incorrect_basic_6_out.txt");
    read_parameters("incorrect_basic_7.txt","incorrect_basic_7_out.txt");
    read_parameters("incorrect_basic_8.txt","incorrect_basic_8_out.txt");
    read_parameters("incorrect_basic_9.txt","incorrect_basic_9_out.txt");
    read_parameters("incorrect_basic_10.txt","incorrect_basic_10_out.txt");
    read_parameters("incorrect_basic_11.txt","incorrect_basic_11_out.txt");
    read_parameters("incorrect_basic_12.txt","incorrect_basic_12_out.txt");
    read_parameters("incorrect_basic_13.txt","incorrect_basic_13_out.txt");
    read_parameters("incorrect_basic_14.txt","incorrect_basic_14_out.txt");
    read_parameters("incorrect_basic_15.txt","incorrect_basic_15_out.txt");
    read_parameters("incorrect_basic_16.txt","incorrect_basic_16_out.txt");
    read_parameters("incorrect_basic_17.txt","incorrect_basic_17_out.txt");
    read_parameters("incorrect_basic_18.txt","incorrect_basic_18_out.txt");
    read_parameters("incorrect_basic_19.txt","incorrect_basic_19_out.txt");
    read_parameters("incorrect_cond_1.txt","incorrect_cond_1_out.txt");
    read_parameters("incorrect_cond_2.txt","incorrect_cond_2_out.txt");
    read_parameters("incorrect_cond_3.txt","incorrect_cond_3_out.txt");
    read_parameters("incorrect_cond_4.txt","incorrect_cond_4_out.txt");
    read_parameters("incorrect_cond_5.txt","incorrect_cond_5_out.txt");
    read_parameters("incorrect_label_1.txt","incorrect_label_1_out.txt");
    read_parameters("incorrect_label_2.txt","incorrect_label_2_out.txt");
    //    read_parameters("incorrect_label_3.txt","incorrect_label_3_out.txt");
    //read_parameters("incorrect_label_4.txt","incorrect_label_4_out.txt");
    //    read_parameters("incorrect_label_5.txt","incorrect_label_5_out.txt");
    read_parameters("incorrect_label_6.txt","incorrect_label_6_out.txt");
    read_parameters("incorrect_label_7.txt","incorrect_label_7_out.txt");
    read_parameters("incorrect_label_8.txt","incorrect_label_8_out.txt");
    read_parameters("incorrect_label_9.txt","incorrect_label_9_out.txt");
    read_parameters("incorrect_label_10.txt","incorrect_label_10_out.txt");
    read_parameters("incorrect_label_11.txt","incorrect_label_11_out.txt");
    read_parameters("incorrect_label_12.txt","incorrect_label_12_out.txt");
    read_parameters("incorrect_label_13.txt","incorrect_label_13_out.txt");
    read_parameters("incorrect_label_14.txt","incorrect_label_14_out.txt");
    read_parameters("incorrect_label_15.txt","incorrect_label_15_out.txt");
    read_parameters("incorrect_label_16.txt","incorrect_label_16_out.txt");
    //    read_parameters("incorrect_label_17.txt","incorrect_label_17_out.txt");
    read_parameters("incorrect_labelarray_1.txt","incorrect_labelarray_1_out.txt");
    read_parameters("incorrect_labelarray_2.txt","incorrect_labelarray_2_out.txt");
    read_parameters("incorrect_labelarray_3.txt","incorrect_labelarray_3_out.txt");


    return 0;
}

