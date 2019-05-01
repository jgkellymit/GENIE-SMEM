#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <cstring>
#include "../build_bwt_index/file.h"
#define QUERY_DB_SIZE 1250000000
#define LUT_KEYSIZE 8

int read_fasta_file(char *file_path, std::string &genome_seq)
{
    /*assume file has just one sequence*/
    genome_seq="";
    std::ifstream input(file_path);
    if(!input.good()){
        return -1;
    }

    std::string line, name, content;

    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ /* Identifier marker*/
            if( !name.empty() ){ /* Print out what we read from the last entry
                                    std::cout << name << " : " << content << std::endl;
                                    ref_name=name;*/
                genome_seq=content;
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ /* Invalid sequence--no spaces allowed*/
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){/*  Print out what we read from the last entry*/
        //std::cout << name << " : " << content << std::endl;
    }

    /*name is name of sequence and content is content of query sequence*/


    /*query_name=name;*/
    genome_seq=content;
    return 0;
}


int main(int argc, char** argv) {
    if (argc != 7) {
        printf("need 6 argument: lutfile indexfile reffile query_set size outputfile");
        exit(1);
    }

    FILE* lutFile = fopen(argv[1], "rb");
    std::fstream outputFile (argv[6] , std::ios::out | std::ios::binary);

    char* const ref_file_path = argv[3];
    // read in reference sequence length
    long reference_seq_len = count_fasta_file(ref_file_path);
    char *ref_seq = (char*)malloc(reference_seq_len * sizeof(char));
    std::string ref_seq_str(ref_seq);
    read_fasta_file(ref_file_path, ref_seq_str);
    if (reference_seq_len == -1) {
        printf("Error opening database '%s'. Bailing out.", ref_file_path);
        return -1;
    }

    // read in suffix array
    char* const index_file_path = argv[2];
    // open the index file
    FILE *stream = fopen(index_file_path, "rb");
    unsigned int *sa_bwt = (unsigned int *) malloc(reference_seq_len * sizeof(unsigned int));
    fread(sa_bwt, sizeof(unsigned int), reference_seq_len, stream);
    fclose(stream);

    uint32_t *lutshit = (uint32_t*)malloc((1U << (2*LUT_KEYSIZE)) * 2 * sizeof(uint32_t));
    std::fread(lutshit, sizeof(uint32_t), (1U << (2*LUT_KEYSIZE)) * 2, lutFile);

    uint32_t querySize = atoi(argv[5]);

    char* query_seq = (char*)malloc(QUERY_DB_SIZE*sizeof(char));
    long numQueries = read_multi_fasta_file(argv[4], query_seq);

    long *enc_qdb = (long*)malloc(numQueries * sizeof(long));
    int keySize = querySize < LUT_KEYSIZE ? querySize : LUT_KEYSIZE;
    for (long st = 0; st < numQueries; st++) {
        long temp = 0;
        long cid = st * querySize;
        for (int r = 0; r < keySize; r++) {
            int index = cid + r + querySize - keySize;
            switch (query_seq[index]) {
                case 'A': temp = (temp << 2) | 0; break;
                case 'C': temp = (temp << 2) | 1; break;
                case 'G': temp = (temp << 2) | 2; break;
                case 'T': temp = (temp << 2) | 3; break;
                default:
                    printf("Illegal character in query: %c at index: %d\n", query_seq[index], index);
                    return -1;
            }
        }
        if (keySize < LUT_KEYSIZE) {
            temp = temp << (2 * (LUT_KEYSIZE - keySize));
        }
        enc_qdb[st] = temp;
    }

    uint32_t* lutput = (uint32_t*)malloc(2 * numQueries * sizeof(uint32_t));

    for (int i = 0; i < numQueries; i++) {
        uint32_t sp = lutshit[enc_qdb[i] * 2];
        long key = sp;
        for (key = sp; key <= sp + LUT_KEYSIZE; key++) {
            bool isShit = false;
            for (int d = 0; d < keySize; d++) {
                if (sa_bwt[key] + d >= reference_seq_len) {
                  printf("fucking hell\n");
                  isShit = true;
                  break;
                }
                if (ref_seq_str[sa_bwt[key] + d] != query_seq[i * querySize + querySize - keySize + d]) isShit = true;
            }
            if (!isShit) break;
        }
        sp = key;
        key = lutshit[enc_qdb[i] * 2 + 1];
        for (key = lutshit[enc_qdb[i] * 2 + 1]; key < reference_seq_len; key++) {
            bool isFucked = false;
            for (int d = 0; d < querySize; d++) {
                if (ref_seq_str[sa_bwt[key] + d] != query_seq[i * querySize + d]) {
                    isFucked = true;
                }
            }
            if (isFucked) break;
        }
        uint32_t ep = key;
        lutput[2 * i] = sp;
        lutput[2 * i + 1] = ep;
    }

    outputFile.write((char*)lutput, 2 * numQueries * sizeof(uint32_t));

    outputFile.close();
    fclose(lutFile);
}
