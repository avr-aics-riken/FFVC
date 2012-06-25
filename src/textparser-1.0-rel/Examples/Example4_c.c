/** @file Example4_c.c
 * サンプルプログラム
 */

#include <stdio.h>
#include <stdlib.h>
#include "TextParser.h"
#define TP_STRING_MAX_SIZE 512


int get_directory_parameters(char *filename, char *label)
{
    int ierror;
    int i;
    char dir_label[TP_STRING_MAX_SIZE];
    char parm_label[TP_STRING_MAX_SIZE];
    int dir_number;
    int parm_number;
    char value[TP_STRING_MAX_SIZE];
    TextParserValueType type;

    //    MGPPChangeDirectory(label, &ierror);
    ierror=tp_changeNode(label);
    if (ierror != 0) {
        printf("ERROR in tp_changeNode file: %s ERROR CODE %d\n", filename, ierror);
        return ierror;
    }

    //label = MGPPCurrentDirectory(&ierror);
    ierror=tp_currentNode(label);
    if (ierror != 0) {
        printf("ERROR in tp_currentNode file: %s ERROR CODE %d\n", filename, ierror);
        return ierror;
    }
    printf("checking Current node: %s\n", label);
    //  dir_number = tp_getDirectoryNumber(&ierror);
    ierror = tp_getNumberOfCNodes(&dir_number);
    if (ierror != 0) {
        printf("ERROR in tp_getNumberOfNodes file: %s ERROR CODE %d\n", filename, ierror);
        return ierror;
    }

    for (i = 0; i < dir_number; i++) {
        ierror = tp_getIthNode(i,dir_label);
        if (ierror != 0) {
            printf("ERROR in tp_getIthNode file: %s ERROR CODE %d\n", filename, ierror);
            return ierror;
        }
        ierror = get_directory_parameters(filename, dir_label);
    }

    //    parm_number = MGPPGetLabelNumber(&ierror);
    ierror = tp_getNumberOfCLeaves(&parm_number);
    if (ierror != 0) {
        printf("ERROR in tp_getNumberOfCLeaves file: %s ERROR CODE %d\n", filename, ierror);
        return ierror;
    }

    for (i = 0; i < parm_number; i++) {
        ierror = tp_getIthLeaf(i,parm_label);
        if (ierror != 0) {
            printf("ERROR in getIthLeaf file: %s ERROR CODE %d\n", filename, ierror);
            return ierror;
        }
        printf("%d %s\n", i, parm_label);

	//        value = MGPPGetValue(parm_label, &ierror);
	ierror = tp_getValue(parm_label,value);
        if (ierror != 0) {
            printf("ERROR in tp_getValue file: %s ERROR CODE %d\n", filename, ierror);
            return ierror;
        }
        printf("%d %s\n", i, value);

	//  type = MGPPGetType(parm_label, &ierror);
	ierror=tp_getType(parm_label, (int*)&type);
        if (ierror != 0) {
            printf("ERROR in tp_getType file: %s ERROR CODE %d\n", filename, ierror);
            return ierror;
        }
        printf("%d value type: %d\n", i, type);
    }

    //label = MGPPCurrentDirectory(&ierror);
    ierror=tp_currentNode(label);
    if (ierror != 0) {
        printf("ERROR in tp_currentNode file: %s ERROR CODE %d\n", filename, ierror);
        return ierror;
    }
    if (strcmp(label, "/")) {
      ierror=tp_changeNode("..");
        if (ierror != 0) {
            printf("ERROR in tp_changeNode file: %s ERROR CODE %d\n", filename, ierror);
            return ierror;
        }

        //label = MGPPCurrentDirectory(&ierror);
	//        MGPPCurrentDirectory(&label, &ierror);
	ierror = tp_currentNode(label);
        if (ierror != 0) {
            printf("ERROR in tp_currentNode file: %s ERROR CODE %d\n", filename, ierror);
            return ierror;
        }
		printf("Current directory: %s\n", label);
    }
}

int move_and_get_parameters(char *filename)
{
    int ierror;
    char root[TP_STRING_MAX_SIZE];
    // ファイルの読み込み
    printf("filename: %s\n", filename);
    ierror=tp_read(filename);
    if (ierror != 0) {
        printf("ERROR in MgppReadParameters file: %s ERROR CODE %d\n", filename, ierror);
        return ierror;
    }
    strcpy(root,"/");
    get_directory_parameters(filename, root);

    // パラメータの削除
    //tp_remove(&ierror);
    tp_remove();
    if (ierror != 0) {
        printf("ERROR in MGPPRemove file: %s ERROR CODE %d\n", filename, ierror);
        return ierror;
    }
    printf("\n");

    return 0;
}

int main(int argc, char* argv[])
{
    move_and_get_parameters ("Input0-1.txt");
    move_and_get_parameters ("Input4-1.txt");

    return 0;
}

