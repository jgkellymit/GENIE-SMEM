#include <stdlib.h>
#include <stdio.h>

#ifndef _FILE_H_INCLUDED
#define _FILE_H_INCLUDED

//long long int read_multi_fasta_file_batch(FILE *input, int64_t first, int64_t last, int64_t cur, char* str_db);

long long int read_multi_fasta_file(char *file_path, char* str_db);

long long int count_fasta_file(char *file_path);

#endif