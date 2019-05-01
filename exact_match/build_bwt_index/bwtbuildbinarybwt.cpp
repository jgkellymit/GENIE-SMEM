/*============================================================================
 *  Name        : bwtbuildbinarybwt.cpp
 *  Author      : Kanak M
 *  Version     : 1
 *  Description : Print  SA, COUNT, OCC.
 *                Written during internship at Intel Labs.
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
#include<fstream>
#define LUT_KEYSIZE 8
using namespace seqan;
using namespace std;

int read_fasta_file(char *file_path,std::string &a);
int build_bowtie_index(std::string ref_db,String<unsigned> sa,String<char> bwt,unsigned int *count,unsigned int **occ,int,char* index, char* lut_name);
int read_multi_fasta_file(char *file_path, std::vector<std::string> &str_db);



int main(int argc, char **argv) {

    if(argc!=5)
    {
        std::cout<<"Need four arguments: reference_sequence checkpoint_frequency index_name lut_file";
        return 1;
    }
    std::vector<std::string> ref_db;
    int freq=atoi(argv[2]);
    std::string reference_seq;
    int status;
    status=read_fasta_file(argv[1],reference_seq);
    if(status==-1){
        std::cerr << "Error opening database '"<<argv[1]<<"'. Bailing out." << std::endl;
        return -1;
    }
    reference_seq+='$';
    unsigned int count[16];
    String<unsigned> sa;
    String<char> bwt;
    unsigned int **occ;
    occ =(unsigned int **)malloc((reference_seq.length()+1)*sizeof(unsigned int *));
    occ[0]=(unsigned int *)malloc(4*(reference_seq.length()+1)*sizeof(unsigned int));
    for(long i=1;i<reference_seq.length()+1;i++)
        occ[i]=occ[i-1]+4;

    //occ = new unsigned int *[reference_seq.length()+1];
    //for(long i = 0; i <reference_seq.length()+1; i++)
    //	occ[i] = new unsigned int[4];	

    status=build_bowtie_index(reference_seq,sa,bwt,count,occ,freq,argv[3],argv[4]);
    return 0;
}
int read_multi_fasta_file(char *file_path, std::vector<std::string> &str_db)
{
    std::ifstream input(file_path);
    if(!input.good()){
        return -1;
    }

    std::string line, name, content;
    std::string gen_seq;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ /* Identifier marker*/
            if( !name.empty() ){ /* Print out what we read from the last entry
                                    std::cout << name << " : " << content << std::endl;
                                    ref_name=name;*/
                gen_seq=content;
                str_db.push_back(gen_seq);
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
        //        }
        //
}
gen_seq=content;
str_db.push_back(gen_seq);
return 0;
}
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
int build_bowtie_index(std::string ref_db,String<unsigned> sa,String<char> bwt,unsigned int *count,unsigned int **occ,int freq,char *index_name, char *lut_name)
{

    String<char> text;
    text=ref_db;
    unsigned int j;
    std::fstream myFile (index_name, ios::out | ios::binary);
    myFile.seekg(0);	
    std::fstream lutFile (lut_name, ios::out | ios::binary);
    lutFile.seekg(0);
    unsigned int *suffix_array=(unsigned int *)malloc(length(text)*sizeof(unsigned int));
    resize(bwt,length(text));
    resize(sa, length(text));
    createSuffixArray(sa, text, Skew7());
    for(unsigned i=0;i<length(text);i++)
    {
        if(sa[i]==0)
            bwt[i]='$';
        else
            bwt[i]=text[sa[i]-1];

    }

    for(unsigned i=0;i<length(text);i++)
    {
        suffix_array[i]=sa[i];
    }
    myFile.write((char*)suffix_array,length(text)*sizeof(unsigned int));

    std::fill_n(count, 16, 0);
    for (unsigned i=0; i<length(text); i++)  
    {
        switch(text[i])
        {
            case 'A': ++count[0];
                      break;
            case 'C': ++count[1];
                      break;
            case 'G': ++count[2];
                      break;
            case 'T': ++count[3];
                      break;
            default: break;
        }
    }	
    count[4]=count[0]+count[1]+count[2]+count[3];
    count[3]=count[0]+count[1]+count[2];
    count[2]=count[0]+count[1];
    count[1]=count[0];
    count[0]=0;	
    myFile.write((char*)count,5*sizeof(unsigned int));

    /*std::cout<<"C strucuture is:\n"
      for(unsigned i=0;i<5;i++)
      {
      std::cout<<i<<",\t"<<count[i]<<"\n";
      }*/

    for(unsigned i=0;i<length(text)+1;i++)
    {
        for(j=0;j<4;j++)
        {
            occ[i][j]=0;
        }
    }
    for(unsigned i=0;i<length(text);i++)
    {

        j=i+1;
        switch(bwt[i])
        {
            case 'A':occ[j][0]=occ[j-1][0]+1;
                     occ[j][1]=occ[j-1][1];
                     occ[j][2]=occ[j-1][2];
                     occ[j][3]=occ[j-1][3];
                     break;
            case 'C':
                     occ[j][0]=occ[j-1][0];
                     occ[j][1]=occ[j-1][1]+1;
                     occ[j][2]=occ[j-1][2];
                     occ[j][3]=occ[j-1][3];
                     break;
            case 'G':
                     occ[j][0]=occ[j-1][0];
                     occ[j][1]=occ[j-1][1];
                     occ[j][2]=occ[j-1][2]+1;
                     occ[j][3]=occ[j-1][3];
                     break;
            case 'T':
                     occ[j][0]=occ[j-1][0];
                     occ[j][1]=occ[j-1][1];
                     occ[j][2]=occ[j-1][2];
                     occ[j][3]=occ[j-1][3]+1;
                     break;
            default:
                     occ[j][0]=occ[j-1][0];
                     occ[j][1]=occ[j-1][1];
                     occ[j][2]=occ[j-1][2];
                     occ[j][3]=occ[j-1][3];
                     break;
        }
    }

    myFile.write((char*)&occ[0][0],((length(text)+1)*4*sizeof(unsigned int)));

    unsigned char *bwt_str = (unsigned char *)malloc(length(text) * sizeof(unsigned char));
    for(uint32_t i = 0; i < length(text); i++)
    {
        bwt_str[i] = bwt[i];
    }
    myFile.write((char*)bwt_str, length(text) * sizeof(unsigned char));

    uint64_t lastKey = 0;
    uint32_t oldInd = 1;
    uint32_t *lut = (uint32_t*)malloc((1 << (2*LUT_KEYSIZE)) * 2 * sizeof(uint32_t));
    for (uint32_t i = 0; i < length(text); i++) {
        uint64_t key = 0;
        bool resetIndex = false;
        if (sa[i] + LUT_KEYSIZE >= length(text)) continue;
        for (int j = sa[i]; j < sa[i] + LUT_KEYSIZE; j++) {
            switch (text[j]) {
                case 'A':
                    key = (key << 2) | 0;
                    break;
                case 'C':
                    key = (key << 2) | 1;
                    break;
                case 'G':
                    key = (key << 2) | 2;
                    break;
                case 'T':
                    key = (key << 2) | 3;
                    break;
            }
        }
        if (key != lastKey) resetIndex = true;

        if (resetIndex) {
            for (uint64_t j = lastKey; j < key; j++) {
                lut[2*j] = oldInd;
                lut[2*j + 1] = i;
            }
            oldInd = i;
            lastKey = key;
        }
    }

    for (uint64_t j = lastKey; j < (1 << 2*LUT_KEYSIZE); j++) {
        lut[2*j] = oldInd;
        lut[2*j + 1] = length(text);
    }

    lutFile.write((char*)lut, (1 << (2*LUT_KEYSIZE)) * 2 * sizeof(uint32_t));
    myFile.close();
    lutFile.close();
    return 0;
}
