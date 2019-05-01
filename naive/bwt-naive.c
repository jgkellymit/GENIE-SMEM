#include <stdio.h>
#include <stdlib.h>
#include "../build_bwt_index/file.h"
#include <omp.h>
#include <string.h>
#include <xmmintrin.h>
#define QUERY_DB_SIZE 12500000000L

int main(int argc, char **argv) {

	if(argc!=5)
	{
		printf("Need four arguments : indexfile query_set referencesequence.fa readlength\n");
		return 1;
	}
	char *query_seq=(char *)malloc(QUERY_DB_SIZE*sizeof(char));
	char *reference_seq=NULL;
	long reference_seq_len;
	long query_db_size;
	query_db_size=read_multi_fasta_file(argv[2],query_seq);

	if(query_db_size==-1)
	{
		printf("Error opening query'%s'. Bailing out.",argv[2]);
		return -1;
	}
	reference_seq_len=count_fasta_file(argv[3]);

	if(reference_seq_len==-1){
		printf("Error opening database '%s'. Bailing out.",argv[3]);
		return -1;
	}
	reference_seq_len++;

	FILE *stream=fopen(argv[1],"rb");

	int count[5];
	unsigned int* sa_bwt=(unsigned int *)malloc(reference_seq_len*sizeof(unsigned int));
	unsigned int *occ;
	int k=0,l=0,r;
	long i=0;
	long size_occ=(reference_seq_len+1);
	occ = (unsigned int *)_mm_malloc(size_occ*4*sizeof(unsigned int),64);

	double elapsed_secs;
	int readlength=atoi(argv[4]);

	int *enc_qdb=(int *)malloc(query_db_size*readlength*sizeof(int));

	fread(sa_bwt,sizeof(unsigned int),reference_seq_len,stream);
	fread(&count[0],sizeof(unsigned int),5,stream);

	for(i=0;i<5;i++)// update read count structure
	{
		count[i]=count[i]+1;
	}


	fread(occ,sizeof(unsigned int),4*(reference_seq_len+1),stream);


	fclose(stream);
	long cind,st;
	for (st=0; st < query_db_size; st++) {
		cind=st*readlength;
		for(r = 0; r < readlength; ++r) {
			switch(query_seq[r+cind])

			{
				case 'A': enc_qdb[r+cind]=0;
					  break;
				case 'C': enc_qdb[r+cind]=1;
					  break;
				case 'G': enc_qdb[r+cind]=2;
					  break;
				case 'T': enc_qdb[r+cind]=3;
					  break;

			}
		}
	}


	long int st1,en1;

	unsigned int *q_result=(unsigned int *)malloc(query_db_size*sizeof(unsigned int));
	unsigned int *k_l_range = (unsigned int *)_mm_malloc(query_db_size*2*sizeof(unsigned int),64);

	int qind;


	for(i=0;i<query_db_size;i++)
	{
		int j=readlength-1;
		long long int sp,ep;
		int ind;
		cind=i*readlength;
		ind=enc_qdb[cind+j];
		sp=count[ind];
		ep=count[ind+1];
		j--;
		while(sp<ep&&j>=0)
		{
			ind=enc_qdb[cind+j];
			sp=count[ind]+occ[sp*4+ind];
			ep=count[ind]+occ[ep*4+ind];
			j--;
		}
		k_l_range[i*2]=sp;
		k_l_range[i*2+1]=ep;
	}

    unsigned int numMatches = 0;
	long u1,u2,u3;
	for(i=0;i<query_db_size;i++)
	{
		u1=k_l_range[i*2];
		u2=k_l_range[i*2+1];
        numMatches += (u2 - u1);
#ifdef PRINT_OUTPUT
		for(u3=u1;u3<u2;u3++)
		{
			printf("%u ",sa_bwt[u3]);
		}
		printf("\n");
#endif
	}
    printf("numMatches = %u\n", numMatches);
	return 0;
}
