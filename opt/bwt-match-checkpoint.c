//#include "sampling.h"
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include "../build_bwt_index/file.h"
#include <omp.h>
#include <string.h>
#include <immintrin.h>

/*static __inline__ unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}*/

#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif

#ifdef IACA_ANALYSIS
#include "iacaMarks.h"
#endif

#define QUERY_DB_SIZE 12500000000L

#define DUMMY_CHAR 6

#ifdef __AVX512PF__
#include <hbwmalloc.h>

#define TID_SECOND_NUMA 999999

#else

#ifdef __AVX512BW__
#define TID_SECOND_NUMA 56
#else
#define TID_SECOND_NUMA 18
#endif

#endif

#define MAX_FIFO_SIZE 1048576

#if defined (__2BIT_LEVEL__)

#if (CP_BLOCK_SIZE == 32) || (CP_BLOCK_SIZE == 64)

#if (CP_BLOCK_SIZE == 64)

#define CP_MASK 63
#define CP_SHIFT 6
#define BIT_DATA_TYPE uint64_t
#define PADDING 24
#define _MM_COUNTBITS __builtin_popcountl

#else
#define CP_MASK 31
#define CP_SHIFT 5
#define BIT_DATA_TYPE uint32_t
#define PADDING 4
#define _MM_COUNTBITS __builtin_popcount

#endif

#define \
GET_OCC(pp, c, occ_id_pp, y_pp, occ_pp, bwt_str_bit0_pp, bwt_str_bit1_pp, bit0_cmp_pp, bit1_cmp_pp, mismatch_mask_pp, num_match_pp) \
                uint32_t occ_id_pp = pp >> CP_SHIFT; \
                uint32_t y_pp = pp & CP_MASK; \
                uint32_t occ_pp = my_cp_occ[occ_id_pp].cp_count[c]; \
                if(y_pp > 0) \
                { \
                BIT_DATA_TYPE bwt_str_bit0_pp = my_cp_occ[occ_id_pp].bwt_str_bit0; \
                BIT_DATA_TYPE bwt_str_bit1_pp = my_cp_occ[occ_id_pp].bwt_str_bit1; \
                BIT_DATA_TYPE bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                BIT_DATA_TYPE bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                uint64_t mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | my_cp_occ[occ_id_pp].dollar_mask; \
                mismatch_mask_pp = mismatch_mask_pp >> (CP_BLOCK_SIZE - y_pp); \
                uint32_t num_match_pp = 0; \
                occ_pp += y_pp - _MM_COUNTBITS(mismatch_mask_pp); \
                }

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
#define _MM_COUNTBITS __builtin_popcountl


#define \
GET_OCC(pp, c, occ_id_pp, y_pp, y0_pp, y1_pp, occ_pp, bwt_str_bit0_pp, bwt_str_bit1_pp, bit0_cmp_pp, bit1_cmp_pp, mismatch_mask_pp, num_match_pp) \
                uint32_t occ_id_pp = pp >> CP_SHIFT; \
                uint32_t y_pp = pp & CP_MASK; \
                uint32_t occ_pp = my_cp_occ[occ_id_pp].cp_count[c]; \
                \
                uint64_t bwt_str_bit0_pp, bwt_str_bit1_pp; \
                uint64_t bit0_cmp_pp, bit1_cmp_pp; \
                uint64_t mismatch_mask_pp; \
                uint32_t y0_pp = y_pp; \
                if(y0_pp > 64) y0_pp = 64;\
                if(y0_pp > 0) \
                { \
                bwt_str_bit0_pp = my_cp_occ[occ_id_pp].bwt_str_bit0[0]; \
                bwt_str_bit1_pp = my_cp_occ[occ_id_pp].bwt_str_bit1[0]; \
                bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | my_cp_occ[occ_id_pp].dollar_mask[0]; \
                mismatch_mask_pp = mismatch_mask_pp >> (64 - y0_pp); \
                occ_pp += y0_pp - _MM_COUNTBITS(mismatch_mask_pp); \
                \
                uint32_t y1_pp = y_pp - y0_pp; \
                if(y1_pp > 0) \
                {\
                bwt_str_bit0_pp = my_cp_occ[occ_id_pp].bwt_str_bit0[1]; \
                bwt_str_bit1_pp = my_cp_occ[occ_id_pp].bwt_str_bit1[1]; \
                bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | my_cp_occ[occ_id_pp].dollar_mask[1]; \
                mismatch_mask_pp = mismatch_mask_pp >> (64 - y1_pp); \
                occ_pp += (y1_pp - _MM_COUNTBITS(mismatch_mask_pp)); \
                }\
                }\

typedef struct checkpoint_occ
{
    uint64_t bwt_str_bit0[ITER];
    uint64_t bwt_str_bit1[ITER];
    uint64_t dollar_mask[ITER];
    uint32_t cp_count[4];
    //uint8_t  pad[0];
}CP_OCC;

#elif (CP_BLOCK_SIZE == 256)

#define CP_MASK 255
#define CP_SHIFT 8
#define BIT_DATA_TYPE uint64_t
#define PADDING 16
#define ITER 4
#define _MM_COUNTBITS __builtin_popcountl


#define \
GET_OCC(pp, c, occ_id_pp, y_pp, y0_pp, y1_pp, occ_pp, bwt_str_bit0_pp, bwt_str_bit1_pp, bit0_cmp_pp, bit1_cmp_pp, mismatch_mask_pp, num_match_pp) \
                uint32_t occ_id_pp = pp >> CP_SHIFT; \
                uint32_t y_pp = pp & CP_MASK; \
                uint32_t occ_pp = my_cp_occ[occ_id_pp].cp_count[c]; \
                \
                uint64_t bwt_str_bit0_pp, bwt_str_bit1_pp; \
                uint64_t bit0_cmp_pp, bit1_cmp_pp; \
                uint64_t mismatch_mask_pp; \
                x = 0; \
                uint32_t y0_pp = y_pp; \
                while(y0_pp > 64) \
                {\
                bwt_str_bit0_pp = my_cp_occ[occ_id_pp].bwt_str_bit0[x]; \
                bwt_str_bit1_pp = my_cp_occ[occ_id_pp].bwt_str_bit1[x]; \
                bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | my_cp_occ[occ_id_pp].dollar_mask[x]; \
                occ_pp += 64 - _MM_COUNTBITS(mismatch_mask_pp); \
                x++;\
                y0_pp -= 64;\
                }\
                \
                if(y0_pp > 0) \
                {\
                bwt_str_bit0_pp = my_cp_occ[occ_id_pp].bwt_str_bit0[x]; \
                bwt_str_bit1_pp = my_cp_occ[occ_id_pp].bwt_str_bit1[x]; \
                bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | my_cp_occ[occ_id_pp].dollar_mask[x]; \
                mismatch_mask_pp = mismatch_mask_pp >> (64 - y0_pp); \
                occ_pp += (y0_pp - _MM_COUNTBITS(mismatch_mask_pp)); \
                }\

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

#define \
GET_OCC(pp, c, c256, occ_id_pp, y_pp, occ_pp, bwt_str_pp, bwt_pp_vec, mask_pp_vec, mask_pp) \
                uint32_t occ_id_pp = pp >> CP_SHIFT; \
                uint32_t y_pp = pp & CP_MASK; \
                uint32_t occ_pp = my_cp_occ[occ_id_pp].cp_count[c]; \
                uint8_t *bwt_str_pp = my_cp_occ[occ_id_pp].bwt_str; \
                __m256i bwt_pp_vec = _mm256_load_si256((const __m256i *)(bwt_str_pp)); \
                __m256i mask_pp_vec = _mm256_cmpeq_epi8(bwt_pp_vec, c256); \
                uint64_t mask_pp = _mm256_movemask_epi8(mask_pp_vec); \
                mask_pp = mask_pp << (32 - y_pp); \
                occ_pp += __builtin_popcount(mask_pp);

typedef struct checkpoint_occ
{
    uint8_t  bwt_str[CP_BLOCK_SIZE];
    uint32_t cp_count[4];
    uint8_t  pad[16];
}CP_OCC;

#endif

typedef struct batch_meta
{
    uint32_t k, l;
    int32_t rid;
    int32_t cur_ind;
}BatchMetadata;

typedef struct inexact_batch_meta
{
    int32_t rid;
    int32_t cur_fifo_id;
    int32_t fifo_size;
    int32_t num_matches;
}InexactBatchMetadata;

typedef struct align_state
{
    uint32_t k, l;
    int32_t i;
    int32_t z;
}AlignState;



int compare_uint(const void *a, const void *b)
{
    uint32_t *pa = (uint32_t *)a;
    uint32_t *pb = (uint32_t *)b;

    uint32_t va = *pa;
    uint32_t vb = *pb;

    if(va < vb) return -1;
    if(va > vb) return 1;
    return 0;
}


int64_t loopTicks;
int64_t loopCount;

void exact_match(CP_OCC *cp_occ, CP_OCC *cp_occ2,
                 uint32_t *count,
                 uint8_t *enc_qdb0, uint8_t *enc_qdb1,
                 int64_t query_db_size, uint32_t readlength,
                 uint32_t *k_l_range0,
                 int32_t batch_size, int32_t numthreads, uint32_t *lutput, int lutkeysize)
{
    uint8_t c_bcast_array[256] __attribute__((aligned(64)));
    int64_t i;
    for(i = 0; i < 4; i++)
    {
        int32_t j;
        for(j = 0; j < 64; j++)
        {
            c_bcast_array[i * 64 + j] = i;
        }
    }

#if defined (__2BIT_LEVEL__)
#if (CP_BLOCK_SIZE >= 64)
    BIT_DATA_TYPE base_mask[4][2];
    base_mask[0][0] = 0;
    base_mask[0][1] = 0;
    base_mask[1][0] = 0xffffffffffffffffL;
    base_mask[1][1] = 0;
    base_mask[2][0] = 0;
    base_mask[2][1] = 0xffffffffffffffffL;
    base_mask[3][0] = 0xffffffffffffffffL;
    base_mask[3][1] = 0xffffffffffffffffL;
#else
    BIT_DATA_TYPE base_mask[4][2];
    base_mask[0][0] = 0;
    base_mask[0][1] = 0;
    base_mask[1][0] = 0xffffffff;
    base_mask[1][1] = 0;
    base_mask[2][0] = 0;
    base_mask[2][1] = 0xffffffff;
    base_mask[3][0] = 0xffffffff;
    base_mask[3][1] = 0xffffffff;
#endif
#endif

#pragma omp parallel num_threads(numthreads)
    {
        int32_t tid = omp_get_thread_num();
        CP_OCC *my_cp_occ;
        uint8_t *my_enc_qdb;
        uint32_t *my_k_l_range;
        if(tid < TID_SECOND_NUMA)
        {
            my_cp_occ = cp_occ;
            my_enc_qdb = enc_qdb0;
            my_k_l_range = k_l_range0;
        }
        else
        {
            my_cp_occ = cp_occ2;
            my_enc_qdb = enc_qdb1;
            my_k_l_range = k_l_range0;
        }
        uint32_t perThreadQuota = (query_db_size + numthreads - 1) / numthreads;
        uint32_t first = tid * perThreadQuota;
        uint32_t last  = (tid + 1) * perThreadQuota;
        if(last > query_db_size) last = query_db_size;

        BatchMetadata metadata[batch_size];

        int64_t ii;

        for(ii = 0; ii < batch_size; ii++)
        {
            BatchMetadata m;
            int64_t rid = ii + first;
            m.k = lutput[2 * rid];
            m.l = lutput[2 * rid + 1];
            m.rid = rid;
            m.cur_ind = readlength - lutkeysize - 1;
            metadata[ii] = m;
        }
        int64_t next_read = first + batch_size;
        int32_t num_active = batch_size;

        while(num_active > 0)
        {
            int32_t new_num_active = num_active;
            for(ii = 0; ii < num_active; ii++)
            {
                int64_t sp, ep;
                int64_t rid;
                uint32_t cur_ind;
                sp = metadata[ii].k;
                ep = metadata[ii].l;
                rid = metadata[ii].rid;
                int j = metadata[ii].cur_ind;
                if (j < 0) {
                    my_k_l_range[2 * rid] = metadata[ii].k;
                    my_k_l_range[2 * rid + 1] = metadata[ii].l;
                    if(next_read < last)
                    {
                        BatchMetadata m;
                        m.k = lutput[2 * next_read];
                        m.l = lutput[2 * next_read + 1];
                        m.rid = next_read;
                        m.cur_ind = readlength - lutkeysize - 1;
                        metadata[ii] = m;
                        next_read++;
                    }
                    else
                    {
                        metadata[ii].rid = -1;
                        new_num_active--;
                    }
                    continue;
                }
                int64_t ind = my_enc_qdb[rid * readlength + j];
                uint8_t *bwt_str;
                uint32_t y;
                int32_t x;
                uint64_t mask;
#if defined (__2BIT_LEVEL__)
#if (CP_BLOCK_SIZE <= 64)
                GET_OCC(sp, ind, occ_id_sp, y_sp, occ_sp, bwt_str_bit0_sp, bwt_str_bit1_sp, bit0_cmp_sp, bit1_cmp_sp, mismatch_mask_sp, num_match_sp);
#else
                GET_OCC(sp, ind, occ_id_sp, y_sp, y0_sp, y1_sp, occ_sp, bwt_str_bit0_sp, bwt_str_bit1_sp, bit0_cmp_sp, bit1_cmp_sp, mismatch_mask_sp, num_match_sp);
#endif
#else
                __m256i c256;
                c256 = _mm256_load_si256((const __m256i *)(c_bcast_array + ind * 64));
                GET_OCC(sp, ind, c256, occ_id_sp, y_sp, occ_sp, bwt_str_sp, bwt_sp_vec, mask_sp_vec, mask_sp);
#endif
                sp  = count[ind] + occ_sp;

                metadata[ii].k = sp;

#if defined (__2BIT_LEVEL__)
#if (CP_BLOCK_SIZE <= 64)
                GET_OCC(ep, ind, occ_id_ep, y_ep, occ_ep, bwt_str_bit0_ep, bwt_str_bit1_ep, bit0_cmp_ep, bit1_cmp_ep, mismatch_mask_ep, num_match_ep);
#else
                GET_OCC(ep, ind, occ_id_ep, y_ep, y0_ep, y1_ep, occ_ep, bwt_str_bit0_ep, bwt_str_bit1_ep, bit0_cmp_ep, bit1_cmp_ep, mismatch_mask_ep, num_match_ep);
#endif
#else
                GET_OCC(ep, ind, c256, occ_id_ep, y_ep, occ_ep, bwt_str_ep, bwt_ep_vec, mask_ep_vec, mask_ep);
#endif
                ep  = count[ind] + occ_ep;

                metadata[ii].l = ep;
                metadata[ii].cur_ind--;
                if((ep <= sp) || (metadata[ii].cur_ind < 0))
                {
                    my_k_l_range[2 * rid] = metadata[ii].k;
                    my_k_l_range[2 * rid + 1] = metadata[ii].l;
                    if(next_read < last)
                    {
                        BatchMetadata m;
                        m.k = lutput[2 * next_read];
                        m.l = lutput[2 * next_read + 1];
                        m.rid = next_read;
                        m.cur_ind = readlength - lutkeysize - 1;
                        metadata[ii] = m;
                        next_read++;
                    }
                    else
                    {
                        metadata[ii].rid = -1;
                        new_num_active--;
                    }
                }
                _mm_prefetch((const char *)(&my_cp_occ[(metadata[ii].k) >> CP_SHIFT]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&my_cp_occ[(metadata[ii].l) >> CP_SHIFT]),_MM_HINT_T0);
#if defined (__2BIT_LEVEL__) && (CP_BLOCK_SIZE == 256)
                _mm_prefetch((const char *)(&my_cp_occ[(metadata[ii].k) >> CP_SHIFT]) + 64, _MM_HINT_T0);
                _mm_prefetch((const char *)(&my_cp_occ[(metadata[ii].l) >> CP_SHIFT]) + 64,_MM_HINT_T0);
#endif
            }
            if(num_active > new_num_active)
            {
                int p1 = 0;
                for(ii = 0; ii < num_active; ii++)
                {
                    if(metadata[ii].rid >= 0)
                    {
                        metadata[p1] = metadata[ii];
                        p1++;
                    }
                }
            }
            num_active = new_num_active;
        }
    }

}

int main(int argc, char **argv) {
    loopTicks = 0;
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    {
        printf("Running:\n");
        int i;
        for(i = 0; i < argc; i++)
        {
            printf("%s ", argv[i]);
        }
        printf("\n");
    }
    printf("TID_SECOND_NUMA = %d\n", TID_SECOND_NUMA);
    if(argc != 10)
    {
        printf("Need 9 arguments : indexfile query_set referencesequence.fa batch_size readlength z n_threads lutputfile lutkeysize\n");
        return 1;
    }

    int lutkeysize = atoi(argv[9]);

    int64_t ref_seq_len;
    ref_seq_len = count_fasta_file(argv[3]);
    if(ref_seq_len == -1){
        printf("Error opening database '%s'. Bailing out.", argv[3]);
        return -1;
    }
    //ref_seq_len--;

    char cp_file_name[1000];
#ifdef __2BIT_LEVEL__
    sprintf(cp_file_name, "%s.2bit.%d", argv[1], CP_BLOCK_SIZE);
#else
    sprintf(cp_file_name, "%s.8bit.%d", argv[1], CP_BLOCK_SIZE);
#endif
    printf("Index file: %s\n", cp_file_name);
    FILE *cpstream = fopen(cp_file_name, "rb");
    FILE *lutputfile = fopen(argv[8], "rb");

    uint32_t count[5];

    uint32_t *sa_bwt = NULL;
    // create checkpointed occ
    uint32_t cp_occ_size = (ref_seq_len >> CP_SHIFT) + 1;
    CP_OCC *cp_occ = NULL;
    CP_OCC *cp_occ2 = NULL;
    uint32_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE - 1) / CP_BLOCK_SIZE) * CP_BLOCK_SIZE;
    int64_t size_occ = (ref_seq_len + 1);
    int32_t readlength = atoi(argv[5]);
    int32_t z = atoi(argv[6]);
    int numthreads = atoi(argv[7]);

    if((z != 0) && (z != 1)) {
        printf("ERROR! z can only be zero or one.\n");
        exit(-1);
    }

    int32_t max_threads = omp_get_max_threads();
#pragma omp parallel num_threads(max_threads)
    {
        int tid = omp_get_thread_num();
#ifdef __AVX512PF__
        if(tid == 0)
        {
            fread(&count[0], sizeof(uint32_t), 5, cpstream);

            cp_occ = (CP_OCC *)hbw_malloc(cp_occ_size * sizeof(CP_OCC));
            fread(cp_occ, sizeof(CP_OCC), cp_occ_size, cpstream);
        }
#else
        if(tid == 0)
        {
            fread(&count[0], sizeof(uint32_t), 5, cpstream);
            cp_occ = (CP_OCC *)_mm_malloc(cp_occ_size * sizeof(CP_OCC), 64);

            fread(cp_occ, sizeof(CP_OCC), cp_occ_size, cpstream);
        }
#endif
    }
    int64_t i = 0;
    for(i = 0; i < 5; i++)// update read count structure
    {
        count[i] = count[i] + 1;
    }


    printf("sizeof CP_OCC = %ld\n", sizeof(CP_OCC));
    uint8_t *enc_qdb0 = NULL;
    uint8_t *enc_qdb1 = NULL;
    uint8_t *query_seq = NULL;
    int64_t query_db_size;
    int64_t orig_query_count;
    int64_t numReads = 0;
    int64_t numSeq = 0;
    uint32_t *k_l_range0 = NULL;

#pragma omp parallel num_threads(numthreads)
    {
        int tid = omp_get_thread_num();
        if(tid == 0)
        {
            printf("Running %d threads\n", omp_get_num_threads());
            query_seq = (uint8_t *)_mm_malloc(QUERY_DB_SIZE * sizeof(uint8_t), 64);
            orig_query_count = query_db_size = read_multi_fasta_file(argv[2], (char*)query_seq);
            enc_qdb0 = query_seq;
            int64_t st;
            for (st=0; st < query_db_size; st++) {
                int64_t cindQ = st * readlength;
                int64_t cind = numReads * readlength;

                int32_t nflag = 0;
                int32_t r;
                for(r = 0; r < readlength; ++r) {
                    switch(query_seq[r + cindQ])

                    {
                        case 'A': enc_qdb0[r + cind] = 0;
                                  break;
                        case 'C': enc_qdb0[r + cind] = 1;
                                  break;
                        case 'G': enc_qdb0[r + cind] = 2;
                                  break;
                        case 'T': enc_qdb0[r + cind] = 3;
                                  break;
                        default: nflag = 1;
                    }
                }
                if(nflag == 0)
                    numReads++;
            }

            query_db_size = numReads;
            printf("query_db_size = %lld\n", query_db_size);
            int64_t numSeq = numReads * (readlength * 3 * z + 1);
            printf("k_l_range size = %llu\n", numSeq * 2 * sizeof(uint32_t));
            k_l_range0 = (uint32_t *)_mm_malloc(numSeq * 2 * sizeof(uint32_t), 64);
            memset(k_l_range0, 0, numSeq * 2 * sizeof(uint32_t));
        }
#pragma omp barrier
        if(tid == TID_SECOND_NUMA)
        {
            cp_occ2 = (CP_OCC *)_mm_malloc(cp_occ_size * sizeof(CP_OCC), 64);
            memcpy(cp_occ2, cp_occ, cp_occ_size * sizeof(CP_OCC));
            enc_qdb1 = (uint8_t *)_mm_malloc(query_db_size * readlength * sizeof(uint8_t), 64);
            memcpy(enc_qdb1, enc_qdb0, query_db_size * readlength * sizeof(uint8_t));
        }
    }

    if(query_db_size == -1)
    {
        printf("Error opening query'%s'. Bailing out.",argv[2]);
        return -1;
    }
    printf("max_occ_ind = %lld\n", i >> CP_SHIFT);

    int32_t batch_size = 0;
    batch_size = atoi(argv[4]);


    int qind;

    if(query_db_size < (batch_size * numthreads))
    {
        printf("very few reads for %d threads\n", numthreads);
        printf("query_db_size = %lld, batch_size = %d, numthreads = %d\n", query_db_size, batch_size, numthreads);
        exit(0);
    }

    uint32_t* lutput = (uint32_t*)malloc(2 * query_db_size * sizeof(uint32_t));
    fread(lutput, sizeof(uint32_t), 2 * query_db_size, lutputfile);

#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
    int64_t startTick = __rdtsc();
    if(z > 0)
    {
        printf("inexact match unsupported\n");
    }
    else
    {
        exact_match(cp_occ, cp_occ2,
                 count,
                 enc_qdb0, enc_qdb1,
                 query_db_size, readlength,
                 k_l_range0,
                 batch_size, numthreads, lutput, lutkeysize);
    }
    int64_t endTick = __rdtsc();
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    int64_t totalTicks = endTick - startTick;

    int64_t numMatches = 0;
    int64_t numMapped  = 0;
    uint32_t buf[1048576];

#ifdef __AVX512PF__
    hbw_free(cp_occ);
#else
    _mm_free(cp_occ);
#endif
    if(cp_occ2 != NULL)
        _mm_free(cp_occ2);
    _mm_free(enc_qdb0);
    if(enc_qdb1 != NULL)
        _mm_free(enc_qdb1);

    sa_bwt = (uint32_t *)_mm_malloc(ref_seq_len * sizeof(uint32_t), 64);
    fread(sa_bwt, sizeof(uint32_t), ref_seq_len, cpstream);

    for(i = 0; i < query_db_size; i++)
    {
        int32_t j;
        int32_t max_map_per_read = z * readlength * 3 + 1;
        int32_t read_num_matches = 0;
        for(j = 0; j < max_map_per_read; j++)
        {
            long u1, u2;
            u1 = k_l_range0[(i * max_map_per_read + j) * 2];
            u2 = k_l_range0[(i * max_map_per_read + j) * 2 + 1];
#ifdef PRINT_OUTPUT
            //printf("%ld] %ld, %ld\n", i, u1, u2);
            memcpy(buf + read_num_matches, sa_bwt + u1, (u2 - u1) * sizeof(uint32_t));
#endif
            read_num_matches += (u2 - u1);
        }
#ifdef PRINT_OUTPUT
            qsort(buf, read_num_matches, sizeof(uint32_t), compare_uint);
            for(j = 0; j < read_num_matches; j++)
            {
                printf("%u ", buf[j]);
            }
            printf("\n");
#endif
        if(read_num_matches > 0) numMapped++;
        numMatches += read_num_matches;
    }
    printf("numReads = %lld\n",   numReads);
    printf("numMatches = %lld\n", numMatches);
    printf("numMapped = %lld\n",  numMapped);
    printf("Consumed %ld cycles\n", totalTicks);

    fclose(cpstream);

    _mm_free(sa_bwt);
    _mm_free(k_l_range0);
    return 0;
}
