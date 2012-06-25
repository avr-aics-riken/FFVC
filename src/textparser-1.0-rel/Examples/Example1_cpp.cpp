/** @file Example1_cpp.cpp
 * サンプルプログラム
 */

#include <string>
#include <vector>
#include <iostream>
#include "TextParser.h"

void read_write_parameters(std::string ifilename, std::string ofilename)
{
    int stat = 0;

    TextParser* tp=TextParser::get_instance();

    // ファイルの読み込み
    std::cout << ">>> ファイルの読み込み(" << ifilename << ") <<<" << std::endl;
    stat=tp->read(ifilename);
    if (stat == 0) {
        std::cout << "# ファイルの読み込みが正常に終了しました : " << ifilename << std::endl;
        // パラメータ書き出し
	//stat=tp->write(ofilename,1);
	stat=tp->write(ofilename);
        if (stat == 0) {
            std::cout << "# ファイルの書き出しが正常に終了しました : " << ofilename << std::endl;
        } else {
            std::cout << "# ファイルの書き出し中にエラーが発生しました : C" << stat << std::endl;
        }
    } else {
        std::cout << "# ファイルの読み込み中にエラーが発生しました : C" << stat << std::endl;
    }

    // 

    stat=tp->remove();
    if (stat == 0) {
        std::cout << "# データベースの削除が正常に終了しました" << std::endl;
    } else {
        std::cout << "# データベースの削除中にエラーが発生しました : C" << stat << std::endl;
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[])
{
    //MGPPDebugWrite(true);
    // 正常


  read_write_parameters("correct_basic_1.txt","correct_basic_1_out.txt");  
  read_write_parameters("correct_basic_2.txt","correct_basic_2_out.txt");
  read_write_parameters("correct_basic_3.txt","correct_basic_3_out.txt");
  read_write_parameters("correct_basic_4.txt","correct_basic_4_out.txt");
  read_write_parameters("correct_basic_5.txt","correct_basic_5_out.txt");
  read_write_parameters("correct_basic_6.txt","correct_basic_6_out.txt");
  read_write_parameters("correct_basic_7.txt","correct_basic_7_out.txt");
  read_write_parameters("correct_basic_8.txt","correct_basic_8_out.txt");
  read_write_parameters("correct_basic_9.txt","correct_basic_9_out.txt");
  read_write_parameters("correct_basic_10.txt","correct_basic_10_out.txt");
  read_write_parameters("correct_basic_11.txt","correct_basic_11_out.txt");
  read_write_parameters("correct_basic_12.txt","correct_basic_12_out.txt");
  read_write_parameters("correct_basic_13.txt","correct_basic_13_out.txt");

  read_write_parameters("correct_string_1.txt","correct_string_1_out.txt");
  read_write_parameters("correct_string_2.txt","correct_string_2_out.txt");
  read_write_parameters("correct_string_3.txt","correct_string_3_out.txt");
  read_write_parameters("correct_string_4.txt","correct_string_4_out.txt");

  read_write_parameters("correct_label_1.txt","correct_label_1_out.txt");
  read_write_parameters("correct_label_2.txt","correct_label_2_out.txt");
  read_write_parameters("correct_label_3.txt","correct_label_3_out.txt");
  read_write_parameters("correct_label_4.txt","correct_label_4_out.txt");
  read_write_parameters("correct_label_5.txt","correct_label_5_out.txt");
  read_write_parameters("correct_label_6.txt","correct_label_6_out.txt");
  read_write_parameters("correct_label_7.txt","correct_label_7_out.txt");
  read_write_parameters("correct_label_8.txt","correct_label_8_out.txt");

  read_write_parameters("correct_cond_1.txt","correct_cond_1_out.txt");
  read_write_parameters("correct_cond_2.txt","correct_cond_2_out.txt");
  read_write_parameters("correct_cond_3.txt","correct_cond_3_out.txt");
  read_write_parameters("correct_cond_4.txt","correct_cond_4_out.txt");
  read_write_parameters("correct_cond_5.txt","correct_cond_5_out.txt");
  read_write_parameters("correct_cond_6.txt","correct_cond_6_out.txt");
  read_write_parameters("correct_cond_7.txt","correct_cond_7_out.txt");
  read_write_parameters("correct_cond_8.txt","correct_cond_8_out.txt");
  read_write_parameters("correct_cond_9.txt","correct_cond_9_out.txt");
  read_write_parameters("correct_cond_10.txt","correct_cond_10_out.txt");

  read_write_parameters("correct_cond_11.txt","correct_cond_11_out.txt");
  read_write_parameters("correct_cond_12.txt","correct_cond_12_out.txt");
  read_write_parameters("correct_cond_13.txt","correct_cond_13_out.txt");
  read_write_parameters("correct_cond_14.txt","correct_cond_14_out.txt");
  read_write_parameters("correct_cond_15.txt","correct_cond_15_out.txt");
  read_write_parameters("correct_cond_16.txt","correct_cond_16_out.txt");
  read_write_parameters("correct_cond_17.txt","correct_cond_17_out.txt");
  read_write_parameters("correct_cond_18.txt","correct_cond_18_out.txt");
  read_write_parameters("correct_cond_19.txt","correct_cond_19_out.txt");
  read_write_parameters("correct_cond_20.txt","correct_cond_20_out.txt");

  read_write_parameters("correct_labelarray_1.txt",
			"correct_labelarray_1_out.txt");
  read_write_parameters("correct_labelarray_2.txt",
			"correct_labelarray_2_out.txt");
  read_write_parameters("correct_labelarray_3.txt",
			"correct_labelarray_3_out.txt");
  read_write_parameters("correct_labelarray_4.txt",
			"correct_labelarray_4_out.txt");
  return 0;
}

