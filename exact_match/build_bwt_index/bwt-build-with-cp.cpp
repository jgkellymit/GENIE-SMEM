/*============================================================================
 *  Name        : bwt-build-with-cp.cpp
 *  Author      : Sanchit Misra
 *  Version     : 1
 *  Description : Print COUNT, CP_OCC, SA.
 *============================================================================
 */
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include<cstring>
#include<vector>
#include<set>
#include<seqan/index.h>
#include <ctime>
#include "file.h"

#define DUMMY_CHAR 6

using namespace seqan;
using namespace std;

#if defined (__2BIT_LEVEL__)

#if (CP_BLOCK_SIZE == 32) || (CP_BLOCK_SIZE == 64)

#if (CP_BLOCK_SIZE == 64)

#define CP_MASK 63
#define CP_SHIFT 6
#define BIT_DATA_TYPE uint64_t
#define PADDING 24

#else
#define CP_MASK 31
#define CP_SHIFT 5
#define BIT_DATA_TYPE uint32_t
#define PADDING 4

#endif


typedef struct checkpoint_occ
{
    BIT_DATA_TYPE bwt_str_bit0;
    BIT_DATA_TYPE bwt_str_bit1;
    BIT_DATA_TYPE dollar_mask;
    uint32_t cp_count[4];
    uint8_t  pad[PADDING];
}CP_OCC;

#elif (CP_BLOCK_SIZE == 128)

#define CP_MASK 127
#define CP_SHIFT 7
#define BIT_DATA_TYPE uint64_t
#define PADDING 0
#define ITER 2
#define _MM_COUNTBITS _mm_countbits_64



typedef struct checkpoint_occ
{
    uint64_t bwt_str_bit0[ITER];
    uint64_t bwt_str_bit1[ITER];
    uint64_t dollar_mask[ITER];
    uint32_t cp_count[4];
}CP_OCC;

#elif (CP_BLOCK_SIZE == 256)

#define CP_MASK 255
#define CP_SHIFT 8
#define BIT_DATA_TYPE uint64_t
#define PADDING 16
#define ITER 4
#define _MM_COUNTBITS _mm_countbits_64

typedef struct checkpoint_occ
{
    uint64_t bwt_str_bit0[ITER];
    uint64_t bwt_str_bit1[ITER];
    uint64_t dollar_mask[ITER];
    uint32_t cp_count[4];
    uint8_t  pad[PADDING];
}CP_OCC;

#endif
#else

#define CP_BLOCK_SIZE 32
#define CP_MASK 31
#define CP_SHIFT 5
#define PADDING 16

typedef struct checkpoint_occ
{
    uint8_t  bwt_str[CP_BLOCK_SIZE];
    uint32_t cp_count[4];
    uint8_t  pad[PADDING];
}CP_OCC;

#endif


int main(int argc, char **argv) {
    if(argc != 3)
    {
        printf("Need three arguments : bwt-file ref_seq.fa\n");
        return 1;
    }
    char *reference_seq = NULL;
    int64_t ref_seq_len;
    ref_seq_len = count_fasta_file(argv[2]);

    if(ref_seq_len == -1){
        printf("Error opening database '%s'. Bailing out.", argv[3]);
        return -1;
    }
    ref_seq_len++;

    FILE *instream = fopen(argv[1], "rb");
    char outname[200];
#ifdef __2BIT_LEVEL__
    sprintf(outname, "%s.2bit.%d", argv[1], CP_BLOCK_SIZE);
#else
    sprintf(outname, "%s.8bit.%d", argv[1], CP_BLOCK_SIZE);
#endif
    std::fstream outstream (outname, ios::out | ios::binary);
    outstream.seekg(0);	
    uint32_t count[5];

    uint32_t *sa_bwt;
    uint32_t *occ;
    uint8_t *bwt;
    int64_t size_occ = (ref_seq_len + 1);

    sa_bwt = (uint32_t *)_mm_malloc(ref_seq_len * sizeof(uint32_t), 64);
    uint32_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE - 1) / CP_BLOCK_SIZE) * CP_BLOCK_SIZE;
    bwt = (uint8_t *)_mm_malloc(ref_seq_len_aligned * sizeof(uint8_t), 64);

    fread(sa_bwt, sizeof(uint32_t), ref_seq_len, instream);
    fread(&count[0], sizeof(uint32_t), 5, instream);
    outstream.write((char*)count, 5 * sizeof(uint32_t));
    printf("count = %u, %u, %u, %u, %u\n", count[0], count[1], count[2], count[3], count[4]);
    occ = (uint32_t *)_mm_malloc(size_occ * 4 * sizeof(uint32_t), 64);
    fread(occ, sizeof(uint32_t), 4 * (ref_seq_len + 1), instream);
    fread(bwt, sizeof(uint8_t), ref_seq_len, instream);

    int64_t i;
    for(i = 0; i < ref_seq_len; i++)
    {
        switch(bwt[i])
        {
            case 'A': bwt[i] = 0;
                      break;
            case 'C': bwt[i] = 1;
                      break;
            case 'G': bwt[i] = 2;
                      break;
            case 'T': bwt[i] = 3;
                      break;
            default:
                      break;
        }
    }
    for(i = ref_seq_len; i < ref_seq_len_aligned; i++)
        bwt[i] = DUMMY_CHAR;
    fclose(instream);

    printf("CP_SHIFT = %d, CP_MASK = %d\n", CP_SHIFT, CP_MASK);
    printf("sizeof CP_OCC = %ld\n", sizeof(CP_OCC));
    // create checkpointed occ
    uint32_t cp_occ_size = (ref_seq_len >> CP_SHIFT) + 1;
    CP_OCC *cp_occ = NULL;

    cp_occ = (CP_OCC *)_mm_malloc(cp_occ_size * sizeof(CP_OCC), 64);
    for(i = 0; i < ref_seq_len; i += CP_BLOCK_SIZE)
    {
        CP_OCC cpo;
        cpo.cp_count[0] = occ[i * 4 + 0];
        cpo.cp_count[1] = occ[i * 4 + 1];
        cpo.cp_count[2] = occ[i * 4 + 2];
        cpo.cp_count[3] = occ[i * 4 + 3];
#if defined (__2BIT_LEVEL__)
#if (CP_BLOCK_SIZE >= 128)
        int64_t x;
        for(x = 0; x < ITER; x++)
        {
            BIT_DATA_TYPE bwt_str_bit0 = 0;
            BIT_DATA_TYPE bwt_str_bit1 = 0;
            BIT_DATA_TYPE dollar_mask = 0;
            int64_t j;
            for(j = 0; j < 64; j++)
            {
                uint8_t c = bwt[i + x * 64 + j];
                if((c == '$') || (c == DUMMY_CHAR))
                {
                    dollar_mask <<= 1;
                    dollar_mask += 1;
                    c = 0;
                }
                else if(c > 3)
                {
                    printf("ERROR! [%d, %d] c = %u\n", i, j, c);
                    exit(0);
                }
                else
                {
                    dollar_mask <<= 1;
                    dollar_mask += 0;
                }
                bwt_str_bit0 = bwt_str_bit0 << 1;
                bwt_str_bit0 += (c & 1);
                bwt_str_bit1 = bwt_str_bit1 << 1;
                bwt_str_bit1 += ((c >> 1) & 1);
            }
            cpo.bwt_str_bit0[x] = bwt_str_bit0;
            cpo.bwt_str_bit1[x] = bwt_str_bit1;
            cpo.dollar_mask[x]  = dollar_mask;
        }
#else
        BIT_DATA_TYPE bwt_str_bit0 = 0;
        BIT_DATA_TYPE bwt_str_bit1 = 0;
        BIT_DATA_TYPE dollar_mask = 0;
        int32_t j;
        for(j = 0; j < CP_BLOCK_SIZE; j++)
        {
            uint8_t c = bwt[i + j];
            if((c == '$') || (c == DUMMY_CHAR))
            {
                dollar_mask <<= 1;
                dollar_mask += 1;
                c = 0;
            }
            else if(c > 3)
            {
                printf("ERROR! [%d, %d] c = %u\n", i, j, c);
                exit(0);
            }
            else
            {
                dollar_mask <<= 1;
                dollar_mask += 0;
            }
            bwt_str_bit0 = bwt_str_bit0 << 1;
            bwt_str_bit0 += (c & 1);
            bwt_str_bit1 = bwt_str_bit1 << 1;
            bwt_str_bit1 += ((c >> 1) & 1);
        }
        cpo.bwt_str_bit0 = bwt_str_bit0;
        cpo.bwt_str_bit1 = bwt_str_bit1;
        cpo.dollar_mask  = dollar_mask;
#endif
#else
        memcpy(cpo.bwt_str, bwt + i, CP_BLOCK_SIZE * sizeof(uint8_t));
#endif
        cp_occ[i >> CP_SHIFT] = cpo;
    }
    _mm_free(occ);
	outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC));
	outstream.write((char*)sa_bwt, ref_seq_len * sizeof(uint32_t));
    //outstream.write((char*)bwt, ref_seq_len * sizeof(uint8_t));
    outstream.close();
    printf("max_occ_ind = %lld\n", i >> CP_SHIFT);    

    _mm_free(cp_occ);
    _mm_free(sa_bwt);
    _mm_free(bwt);
    return 0;
}

