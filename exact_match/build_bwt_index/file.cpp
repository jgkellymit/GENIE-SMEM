#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "file.h"

/*long long int read_multi_fasta_file_batch(FILE *input, int64_t first, int64_t last, int64_t cur, char* str_db)
{
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    if(input==NULL){
        return -1;
    }
    int64_t i = cur;
    char name[1000]="", content[1000]="";

    name[0]='\0',content[0]='\0';

    while( (read=getline(&line,&len,input))!=-1 ){
        if( strlen(line)==0 || line[0] == '>' ){ // Identifier marker
            if( strlen(name)!=0 ){

                if(i >= first)
                {
                    memcpy(str_db,content,(strlen(content)));
                    str_db+=(strlen(content)-1);
                }
                content[0]='\0';
                i++;
                if(i == last)
                    break;
                name[0]='\0';
            }
            if(( strlen(line)!=0) ){
                memcpy(name,&line[1],strlen(line));
            }
            content[0]='\0';
        } else if( strlen(name)!=0 ){
            if( strstr(line," ") != NULL ){ // Invalid sequence--no spaces allowed
                name[0]='\0';
                content[0]='\0';
            } else {
                strcat(content,line);
            }
        }
    }
    if(strlen(content) > 0)
    {
        memcpy(str_db,content,strlen(content));
        i++;
    }
    return i - first;
}
*/
long long int read_multi_fasta_file(char *file_path, char* str_db)
{
    FILE *input;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    input=fopen(file_path,"r");

    if(input==NULL){
        return -1;
    }
    int i=0;
    char name[1000]="", content[1000]="";

    name[0]='\0',content[0]='\0';

    while( (read=getline(&line,&len,input))!=-1 ){
        if( strlen(line)==0 || line[0] == '>' ){ /* Identifier marker*/
            if( strlen(name)!=0 ){

                memcpy(str_db,content,(strlen(content)));
                str_db+=(strlen(content)-1);
                i++;
                name[0]='\0';
            }
            if(( strlen(line)!=0) ){
                memcpy(name,&line[1],strlen(line));
            }
            content[0]='\0';
        } else if( strlen(name)!=0 ){
            if( strstr(line," ") != NULL ){ /* Invalid sequence--no spaces allowed*/
                name[0]='\0';
                content[0]='\0';
            } else {
                strcat(content,line);
            }
        }
    }
    memcpy(str_db,content,strlen(content));
    i++;
    return i;
}

long long int count_fasta_file(char *file_path) {
    /*assume file has just one sequence*/
    FILE *input;
    char *line = NULL;
    size_t len = 0;

    input = fopen(file_path, "r");

    if (input == NULL) {
        return -1;
    }
    long long int i = 0;
    char name[1000], content[1000];
    name[0] = '\0', content[0] = '\0';

    while (getline(&line, &len, input) != -1) {
        if (strlen(line) == 0 || line[0] == '>') { /* Identifier marker*/
            if (strlen(name) != 0) {

                i += strlen(content) - 1;
                name[0] = '\0';
            }
            if ((strlen(line) != 0)) {
                memcpy(name, &line[1], strlen(line));

            }
            content[0] = '\0';
        } else if (strlen(name) != 0) {
            if (strstr(line, " ") != NULL) { /* Invalid sequence--no spaces allowed*/
                name[0] = '\0';
                content[0] = '\0';
            } else {
                i += strlen(line) - 1;
            }
        }
    }

    //i++;
    return i;
}
