//**********************************************************************************
//* Title: SEED: Efficient Clustering Software of Next Generation Sequences
//* Platform: 32-Bit/64-Bit Windows/Linux/Mac
//* Author: Ergude Bao
//* Affliation: Department of Computer Science & Engineering
//* University of California, Riverside
//* Date: 08/08/2011
//* Version: 1.4.1
//* Copy Right: For Purpose of Study Only
//**********************************************************************************
//In this updated version with pre-sorting, N base detection in function Hash::build() and read sequence ID recording 
//in fucntion Hash::seqInsert should have been rolled back. However, for possible usage in future versions, they stay
//the same.

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <cstring>
using namespace std;

int QV = 0;
int reversed = 0;
int seedsCount = 10;
int seedsWeight = 16 * 1024;

#define OFFSET 33
#define RANGE 94 

static int seeds[10][30] = 
{
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
	1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
	0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

static int fastSeeds[4][52] =
{
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

static int shortSeeds[10][15] =
{
	1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
	1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
	1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
	0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
	0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1,
	0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1	
};

class Hash
{
	ifstream in;
	int lowerSizeInBit, lowerSizeInChar, upperSizeInBit, upperSizeInChar, num, mismatchAllowed;
	char * seqHead;
	unsigned int ** indexHead;
//	unsigned int * offsetCount;
public:
	unsigned int * offsetCount;
	Hash(char [], int, int, int);
	void build();
	void seqInsert(char [], int, char *, unsigned int);
	void indexInsert(char [], unsigned int);
	char change(char);
	char changeBack(char);
	unsigned int calOffset(int, char []);
	int searchByIndex(int, unsigned int, int, char []);
	int searchBySeq(unsigned int, char []);
	void deleteByIndex(int, unsigned int, int);
	void deleteBySeq(unsigned int);
	int calSeqID(int, unsigned int, int);
	void QVInsert(char [], int, char *, unsigned int);
	int searchByIndex(int, unsigned int, int, char [], char [], unsigned int &);
	int searchBySeq(unsigned int, char [], char [], unsigned int &);
	void seqInsert(char [], int, char *, unsigned int, unsigned int);
	int searchByIndex(int, unsigned int, int, char [], unsigned int &);
	int searchBySeq(unsigned int, char [], unsigned int &);
	void adjust();
	void tmpDeleteByIndex(int, unsigned int, int);
	void recoverByIndex(int, unsigned int, int);
	unsigned int getOffsetCount(unsigned int, int);
};

class FastqGenerator
{
	ifstream in;
	char addiInput[100];
	ifstream addiIn;
	char outputq[100];
	ofstream out;
	int num;
	char * seq;
public:
	FastqGenerator(char [], char [], int);
	void record();
	void generateFastq();
};

class Cluster
{
	ofstream out;
//	ofstream dis;
	int lowerSizeInChar;
	int upperSizeInChar;
	int lowerSizeInBit;
	int upperSizeInBit;
	int num;
	int mismatchAllowed;
	int shiftAllowed;
	int CLID;
//	int numInCL;
//	Hash * h;
	int lowerQV;
	int upperQV;
	ofstream addiOut;
	char addiOutput[100];
	unsigned long seqNum;
	unsigned long adjustNum;
	unsigned int ** mappingTable;
	unsigned int * mappingNum;
	char midInput[100];
public:
	Hash * h;
	Cluster(char [], char [], int, int, int, int, int, unsigned int **, unsigned int *);
	void cluster();
	int compare(char [], char [], int, int, int &);
	void clusterWithMismatches(char []);
	void clusterWithShifts(char []);
	char max(int, int, int, int);
	void clusterByConsensus();
	void calConsensus(char [], unsigned int, int &);
	void preprocess();
	Cluster(char [], char [], int, int, int, int, int, int, int, unsigned int **, unsigned int *);
	int compare(char [], char [], char [], char [], int, int, int &);
	void clusterWithMismatches(char [], char []);
	void clusterWithShifts(char [], char []);
	char reverseChange(char);
};

class FileAnalyzer
{
public:
	void inputAnalyze(char [], int &, int &, int &, int &);
	void outputAnalyze(int, int);
};

class Sorter
{
	ifstream in;
	ofstream midOut;
	char midOutput[100];
	int num;
	unsigned int ** mappingTable;
	unsigned int * mappingNum;
	int lowerSizeInChar;
	int realNum;
	typedef struct 
	{
		int ID;
		int realID;
	} Order;
public:
	Sorter(char [], int, int);
	void sort();
	void suffixSort(int, int, int, char [], Order []);
	int getRealNum();
	unsigned int ** getMappingTable();
	unsigned int * getMappingNum();
};

Hash::Hash(char input[], int num, int lowerSizeInChar, int upperSizeInChar)
{
	long int i;

	if(upperSizeInChar % 4)
		this->upperSizeInBit = upperSizeInChar / 4 + 1;
	else
		this->upperSizeInBit = upperSizeInChar / 4;
	if(lowerSizeInChar % 4)
		this->lowerSizeInBit = lowerSizeInChar / 4 + 1;
	else
		this->lowerSizeInBit = lowerSizeInChar / 4;
	this->lowerSizeInChar = lowerSizeInChar;
	this->upperSizeInChar = upperSizeInChar;
	this->num = num;
	if(QV)
	{
		seqHead = new char[(long int)num * (upperSizeInBit + 5 + upperSizeInChar)];
		for(i = 0; i < (long int)num * (upperSizeInBit + 5 + upperSizeInChar); i ++)
			seqHead[i] = 0;
	}
	else
	{
		seqHead = new char[(long int)num * (upperSizeInBit + 5)];
		for(i = 0; i < (long int)num * (upperSizeInBit + 5); i ++)
			seqHead[i] = 0;
	}
//	indexHead = new unsigned int * [1024 * 1024 * 16 * 10];
//	indexHead = new unsigned int * [1024 * 1024 * 64 * 4];
	indexHead = new unsigned int * [1024 * seedsWeight * seedsCount];
	offsetCount = new unsigned int [1024 * seedsWeight * seedsCount];
	for(i = 0; i < (long int)1024 * seedsWeight * seedsCount; i ++)
		offsetCount[i] = 0;
	in.open(input);
}

void Hash::build()
{
	char buf[201];
	int i, count = 0, tag = 1;
	unsigned int seqOffset = 0, seqID;

//	int seqNext;
//	char base[4];
//	ofstream out;

	int totalLength, j;

	if(in.is_open())
	{
		seqID = 0;
cont:
		while(seqID < num * 4)
		{
			in.getline(buf, 201);
			if(seqID % 4 == 1)
			{
				for(i = 0; i < in.gcount() - 1; i ++)
				{
					if(buf[i] == 'N')
					{
						tag = 0;
						seqID ++;
						goto cont;
					}//to deal with N
					else
					{
						tag = 1;
						buf[i] = change(buf[i]);
					}
				}
				seqInsert(buf, in.gcount() - 1, seqHead, seqOffset, seqID/4);
				indexInsert(buf, seqOffset);
				seqOffset = seqOffset + (upperSizeInBit + 5);
			}
			if(QV)
			{
				if(seqID % 4 == 3 && tag == 1)
				{
					if(buf[in.gcount()] == '\n')
						QVInsert(buf, in.gcount() - 1, seqHead, seqOffset);
					else
						QVInsert(buf, in.gcount(), seqHead, seqOffset);//for the last line in Windows OS
					seqOffset = seqOffset + upperSizeInChar;
				}
				for(i = 0; i < 201; i ++)
					buf[i] = 0;
			}
			seqID ++;
		}
	}
	else
	{
		cout << "CANNOT OPEN INPUT FILE!" << endl;
		exit(-1);
	}

//verification
//	cout << endl;
//#ifdef QV
//#ifdef REALID
//	for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 5 + upperSizeInChar); seqOffset = seqOffset + upperSizeInBit + 5 + upperSizeInChar)
//#else
//	for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 1 + upperSizeInChar); seqOffset = seqOffset + upperSizeInBit + 1 + upperSizeInChar)
//#endif
//#else
//#ifdef REALID
//	for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 5); seqOffset = seqOffset + upperSizeInBit + 5)
//#else
//	for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 1); seqOffset = seqOffset + upperSizeInBit + 1)
//#endif
//#endif
//	{
//		cout << (int)*(seqHead + seqOffset) << ": " << endl;
//		for(seqNext = 1; seqNext < 5; seqNext ++)
//		{
//			cout << (int)*(seqHead + seqOffset + seqNext) << " ";
//		}
//		cout << endl;
//#ifdef REALID
//		for(seqNext = 5; seqNext < upperSizeInBit + 5; seqNext ++)
//#else
//		for(seqNext = 1; seqNext < upperSizeInBit + 1; seqNext ++)
//#endif
//		{
//			base[0] = (*(seqHead + seqOffset + seqNext) >> 6) & 0x03;
//			base[1] = (*(seqHead + seqOffset + seqNext) >> 4) & 0x03;
//			base[2] = (*(seqHead + seqOffset + seqNext) >> 2) & 0x03;
//			base[3] = (*(seqHead + seqOffset + seqNext)) & 0x03;
//			cout << (unsigned int)base[0] << " " << (unsigned int)base[1] << " " << (unsigned int)base[2] << " " << (unsigned int)base[3] << " ";
//		}
//#ifdef QV
//		cout << endl;
//#ifdef REALID
//		for(; seqNext < upperSizeInBit + 1 + upperSizeInChar; seqNext ++)
//#else
//		for(; seqNext < upperSizeInBit + 5 + upperSizeInChar; seqNext ++)
//#endif
//			cout << *(seqHead + seqOffset + seqNext) << " ";
//#endif
//		cout << endl;
//	}
//verification
/*
        totalLength = 0;
        count = 0;
        for(i = 0; i < seedsCount; i ++)
        {
                for(j = 0; j < 1024 * seedsWeight; j ++)
                        if(c.h->offsetCount[j * seedsCount + i] > 1000)
                        {
                                totalLength = totalLength + c.h->offsetCount[j * seedsCount + i];
                                count ++;
                        }
        }
        cout << "#buckets longer than 1000: " << count << endl;
	if(count != 0)
	        cout << "average length of the buckets: " << totalLength / count << endl;
*/
}

void Hash::adjust()
{
	int i, j, k ,p;

	for(i = 0; i < seedsCount; i ++)
		for(j = 0; j < 1024 * seedsWeight; j ++)
			if(offsetCount[j * seedsCount + i] != 0)
			{
				p = offsetCount[j * seedsCount + i] - 1;
				for(k = 0; k < offsetCount[j * seedsCount + i] && p > k; k ++)
				{
					if(seqHead[indexHead[j * seedsCount + i][k]] == 0)
					{
						while(seqHead[indexHead[j * seedsCount + i][p]] == 0 && p > k)
						{
							p --;
							offsetCount[j * seedsCount + i] --;
						}
						if(p > k)
						{
							indexHead[j * seedsCount + i][k] = indexHead[j * seedsCount + i][p];
							p --;
							offsetCount[j * seedsCount + i] --;
						}
					}
				}
				if(seqHead[indexHead[j * seedsCount + i][k]] == 0)
					offsetCount[j * seedsCount + i] --;
			}
}

void Hash::seqInsert(char buf[], int realSize, char * seqHead, unsigned int seqOffset, unsigned int seqID)
{
	int i, seqNext;
	char bitBuf = 0x00;

	*(seqHead + seqOffset) = (char)realSize;
	*(seqHead + seqOffset + 1) = (char) ((seqID & 0xff000000) >> 24);
	*(seqHead + seqOffset + 2) = (char) ((seqID & 0x00ff0000) >> 16);
	*(seqHead + seqOffset + 3) = (char) ((seqID & 0x0000ff00) >> 8);
	*(seqHead + seqOffset + 4) = (char) (seqID & 0x000000ff);
	for(i = 0, seqNext = 5; i < realSize; i ++)
	{
		if((i + 1) % 4 == 0)
		{
			bitBuf = (bitBuf | buf[i]);
			*(seqHead + seqOffset + seqNext) = bitBuf;
			seqNext ++;
			bitBuf = 0x00;
		}
		else if(i == realSize - 1)
		{
			bitBuf = (bitBuf | buf[i]) << (4 - realSize % 4) * 2;
			*(seqHead + seqOffset + seqNext) = bitBuf;
			seqNext ++;
			bitBuf = 0x00;
		}
		else
			bitBuf = (bitBuf | buf[i]) << 2;
	}
}

void Hash::QVInsert(char buf[], int realSize, char * seqHead, unsigned int seqOffset)
{
	int i, seqNext = 0;

	for(i = 0; i < realSize; i ++)
	{
		*(seqHead + seqOffset + seqNext) = buf[i];
		seqNext ++;
	}
}

void Hash::indexInsert(char buf[], unsigned int seqOffset)
{
	int indexNext;
	unsigned int indexOffset = 0;

	for(indexNext = 0; indexNext < seedsCount; indexNext ++)
	{
		indexOffset = calOffset(indexNext, buf);
		if(offsetCount[indexOffset + indexNext] == 0)
		{
			indexHead[indexOffset + indexNext] = (unsigned int *) malloc((++ offsetCount[indexOffset + indexNext]) * sizeof(unsigned int));
			if(indexHead[indexOffset + indexNext] == NULL)
			{
				cout << "CANNOT ALLOCATE MEMORY!" << endl;
				exit(-1);
			}
		}
		else
		{
			indexHead[indexOffset + indexNext] = (unsigned int *) realloc(indexHead[indexOffset + indexNext], (++ offsetCount[indexOffset + indexNext]) * sizeof(unsigned int));
			if(indexHead[indexOffset + indexNext] == NULL)
			{
				cout << "CANNOT ALLOCATE MEMORY!" << endl;
				exit(-1);
			}
		}
		indexHead[indexOffset + indexNext][offsetCount[indexOffset + indexNext] - 1] = seqOffset;
	}
}

char Hash::change(char base)
{
	switch(base)
	{
		case 'A': return 0x00;
		case 'C': return 0x01;
		case 'G': return 0x02;
		case 'T': return 0x03;
		default: cout << "INPUT ERROR!" << endl; exit(-1);
	}
}

char Hash::changeBack(char base)
{
	switch(base)
	{
		case 0x00: return 'A';
		case 0x01: return 'C';
		case 0x02: return 'G';
		case 0x03: return 'T';
		default: cout << "MEMORY ERROR!" << endl; exit(-1);
	}
}

unsigned int Hash::calOffset(int indexNext, char buf[])
{
	unsigned int indexOffset = 0;
	int i, j = 0;

	if(seedsWeight == 1024 * 16)
	{
		for(i = 0; i < 30; i ++)
			if(seeds[indexNext][i] == 1)
				indexOffset = indexOffset + buf[i + lowerSizeInChar - 33] * (unsigned int)pow(4, j ++);
	}
	else if(seedsWeight == 1024 * 64)
	{
		for(i = 0; i < 52; i ++)
			if(fastSeeds[indexNext][i] == 1)
				indexOffset = indexOffset + buf[i + lowerSizeInChar - 55] * (unsigned int)pow(4, j ++);
	}
	else
	{
		for(i = 0; i < 15; i ++)
                        if(shortSeeds[indexNext][i] == 1)
                                indexOffset = indexOffset + buf[i + lowerSizeInChar - 18] * (unsigned int)pow(4, j ++);
	}
	return indexOffset * seedsCount;
}

int Hash::searchByIndex(int indexNext, unsigned int indexOffset, int no, char buf[], unsigned int & seqID)
{
	int seqNext, p = 0;
	int i;

	seqID = 0;
	for(i = 1; i < 5; i ++)
	{
		seqID = seqID | (((unsigned int) seqHead[indexHead[indexOffset + indexNext][no] + i]) & 0x000000ff);
		if(i < 4)
			seqID = seqID << 8;
	}
	for(seqNext = 5; seqNext < upperSizeInBit + 5; seqNext ++)
	{
		buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 6) & 0x03;
		buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 4) & 0x03;
		buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 2) & 0x03;
		buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext]) & 0x03;
	}
	return (int)seqHead[indexHead[indexOffset + indexNext][no]];
}

int Hash::searchBySeq(unsigned int seqOffset, char buf[], unsigned int & seqID)
{
	int seqNext, p = 0;
	int i;

	seqID = 0;
	for(i = 1; i < 5; i ++)
	{
		seqID = seqID | (((unsigned int) seqHead[seqOffset + i]) & 0x000000ff);
		if(i < 4)
			seqID = seqID << 8;
	}
	for(seqNext = 5; seqNext < upperSizeInBit + 5; seqNext ++)
	{
		buf[p ++] = (seqHead[seqOffset + seqNext] >> 6) & 0x03;
		buf[p ++] = (seqHead[seqOffset + seqNext] >> 4) & 0x03;
		buf[p ++] = (seqHead[seqOffset + seqNext] >> 2) & 0x03;
		buf[p ++] = (seqHead[seqOffset + seqNext]) & 0x03;
	}
	return (int)seqHead[seqOffset];
}

int Hash::searchByIndex(int indexNext, unsigned int indexOffset, int no, char buf[], char QVBuf[], unsigned int & seqID)
{
	int seqNext, p = 0;
	int i;

	seqID = 0;
	for(i = 1; i < 5; i ++)
	{
		seqID = seqID | (((unsigned int) seqHead[indexHead[indexOffset + indexNext][no] + i]) & 0x000000ff);
		if(i < 4)
			seqID = seqID << 8;
	}
	for(seqNext = 5; seqNext < upperSizeInBit + 5; seqNext ++)
	{
		buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 6) & 0x03;
		buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 4) & 0x03;
		buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 2) & 0x03;
		buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext]) & 0x03;
	}
	p = 0;
	for(; seqNext < upperSizeInBit + 5 + upperSizeInChar; seqNext ++)
		QVBuf[p ++] = seqHead[indexHead[indexOffset + indexNext][no] + seqNext];
	return (int)seqHead[indexHead[indexOffset + indexNext][no]];
}

int Hash::searchBySeq(unsigned int seqOffset, char buf[], char QVBuf[], unsigned int & seqID)
{
	int seqNext, p = 0;
	int i;

	seqID = 0;
	for(i = 1; i < 5; i ++)
	{
		seqID = seqID | (((unsigned int) *(seqHead + seqOffset + i)) & 0x000000ff);
		if(i < 4)
			seqID = seqID << 8;
	}
	for(seqNext = 5; seqNext < upperSizeInBit + 5; seqNext ++)
	{
		buf[p ++] = (seqHead[seqOffset + seqNext] >> 6) & 0x03;
		buf[p ++] = (seqHead[seqOffset + seqNext] >> 4) & 0x03;
		buf[p ++] = (seqHead[seqOffset + seqNext] >> 2) & 0x03;
		buf[p ++] = (seqHead[seqOffset + seqNext]) & 0x03;
	}
	p = 0;
	for(; seqNext < upperSizeInBit + 5 + upperSizeInChar; seqNext ++)
		QVBuf[p ++] = seqHead[seqOffset + seqNext];
	return (int)seqHead[seqOffset];
}

void Hash::deleteByIndex(int indexNext, unsigned int indexOffset, int no)
{
	seqHead[indexHead[indexOffset + indexNext][no]] = 0;
}

void Hash::tmpDeleteByIndex(int indexNext, unsigned int indexOffset, int no)
{
	seqHead[indexHead[indexOffset + indexNext][no]] = seqHead[indexHead[indexOffset + indexNext][no]] | 0x80;
}

void Hash::recoverByIndex(int indexNext, unsigned int indexOffset, int no)
{
	seqHead[indexHead[indexOffset + indexNext][no]] = seqHead[indexHead[indexOffset + indexNext][no]] & 0x7f;
}

void Hash::deleteBySeq(unsigned int seqOffset)
{
	seqHead[seqOffset] = 0;
}

int Hash::calSeqID(int indexNext, unsigned int indexOffset, int no)
{
	if(QV)
		return indexHead[indexOffset + indexNext][no]/(upperSizeInBit + 1 + upperSizeInChar);
	else
		return indexHead[indexOffset + indexNext][no]/(upperSizeInBit + 1);
}

unsigned int Hash::getOffsetCount(unsigned int indexOffset, int indexNext)
{
	return offsetCount[indexOffset + indexNext];
}

Cluster::Cluster(char input[], char output[], int num, int lowerSizeInChar, int upperSizeInChar, int mismatchAllowed, int shiftAllowed, int lowerQV, int upperQV, unsigned int ** mappingTable, unsigned int * mappingNum)
{
	if(lowerSizeInChar % 4)
		this->lowerSizeInBit = lowerSizeInChar / 4 + 1;
	else
		this->lowerSizeInBit = lowerSizeInChar / 4;
	if(upperSizeInChar % 4)
		this->upperSizeInBit = upperSizeInChar / 4 + 1;
	else
		this->upperSizeInBit = upperSizeInChar / 4;
	this->lowerSizeInChar = lowerSizeInChar;
	this->upperSizeInChar = upperSizeInChar;
	this->num = num;
	this->mismatchAllowed = mismatchAllowed;
	this->shiftAllowed = shiftAllowed;
	this->CLID = 0;
	this->lowerQV = lowerQV;
	this->upperQV = upperQV;
	this->seqNum = 0;
	this->adjustNum = num / 20;
	this->mappingTable = mappingTable;
	this->mappingNum = mappingNum;
	strcpy(midInput, input);
	strcat(midInput, ".mid.fastq");
	out.open(output);
	strcpy(addiOutput, output);
	strcat(addiOutput, ".fasta");
	addiOut.open(addiOutput);
//	dis.open("distribution.txt");
	out << "CLID	SeqID" << endl;
//	dis << "CLID	No" << endl;
	h = new Hash(midInput, num, lowerSizeInChar, upperSizeInChar);
	h->build();
}

Cluster::Cluster(char input[], char output[], int num, int lowerSizeInChar, int upperSizeInChar, int mismatchAllowed, int shiftAllowed, unsigned int ** mappingTable, unsigned int * mappingNum)
{
	if(lowerSizeInChar % 4)
		this->lowerSizeInBit = lowerSizeInChar / 4 + 1;
	else
		this->lowerSizeInBit = lowerSizeInChar / 4;
	if(upperSizeInChar % 4)
		this->upperSizeInBit = upperSizeInChar / 4 + 1;
	else
		this->upperSizeInBit = upperSizeInChar / 4;
	this->lowerSizeInChar = lowerSizeInChar;
	this->upperSizeInChar = upperSizeInChar;
	this->num = num;
	this->mismatchAllowed = mismatchAllowed;
	this->shiftAllowed = shiftAllowed;
	this->CLID = 0;
	this->seqNum = 0;
	this->adjustNum = num / 20;
	this->mappingTable = mappingTable;
	this->mappingNum = mappingNum;
	strcpy(midInput, input);
	strcat(midInput, ".mid.fastq");
	out.open(output);
	strcpy(addiOutput, output);
	strcat(addiOutput, ".fasta");
	addiOut.open(addiOutput);
//	dis.open("distribution.txt");
	out << "CLID	SeqID" << endl;
//	dis << "CLID	No" << endl;
	h = new Hash(midInput, num, lowerSizeInChar, upperSizeInChar);
	h->build();
}

char Cluster::max(int a, int b, int c, int d)
{
	if(b >= a && b >= c && b >= d) return 0x01;
	if(c >= a && c >= b && c >= d) return 0x02;
	if(d >= a && d >= b && d >= c) return 0x03;
	if(a >= b && a >= c && a >= d) return 0x00;
	cout << "UNKNOWN ERROR!" << endl; exit(-1);
}

void Cluster::calConsensus(char sBuf[], unsigned int sSeqID, int & tagReverse)
{
	long int A[200] = {0}, C[200] = {0}, G[200] = {0}, T[200] = {0};
	int indexNext, no, i, realSize, similarity = 200, diff, j;
	char sBufBak[200], tBuf[200], buf[200];
	unsigned int indexOffset;
	unsigned int seqID, centerSeqID;

	for(indexNext = 0; indexNext < seedsCount; indexNext ++)
	{
		indexOffset = h->calOffset(indexNext, sBuf);
		for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
		{
			realSize = h->searchByIndex(indexNext, indexOffset, no, tBuf, seqID);
			if(realSize > 0 && compare(sBuf, tBuf, 0, lowerSizeInChar, tagReverse) <= mismatchAllowed * 2)
			{
				if(reversed && tagReverse)
					for(i = 0, j = lowerSizeInChar - 1; i < lowerSizeInChar; i ++, j --)
						switch(tBuf[j])
						{
							case 0x00: T[i] = T[i] + mappingNum[seqID]; break;
							case 0x01: G[i] = G[i] + mappingNum[seqID]; break;
							case 0x02: C[i] = C[i] + mappingNum[seqID]; break;
							case 0x03: A[i] = A[i] + mappingNum[seqID]; break;
							default: cout << "MEMORY ERROR!" << endl; exit(-1);
						}	
				else
					for(i = 0; i < lowerSizeInChar; i ++)
						switch(tBuf[i])
						{
							case 0x00: A[i] = A[i] + mappingNum[seqID]; break;
							case 0x01: C[i] = C[i] + mappingNum[seqID]; break;
							case 0x02: G[i] = G[i] + mappingNum[seqID]; break;
							case 0x03: T[i] = T[i] + mappingNum[seqID]; break;
							default: cout << "MEMORY ERROR!" << endl; exit(-1);
						}
				h->tmpDeleteByIndex(indexNext, indexOffset, no);
			}
		}
	}
	for(i = 0; i < lowerSizeInChar; i ++)
	{
		sBufBak[i] = sBuf[i];
		sBuf[i] = max(A[i], C[i], G[i], T[i]);
	}
//find the most similar seq to the consensus and write to the output
	for(indexNext = 0; indexNext < seedsCount; indexNext ++)
	{
		indexOffset = h->calOffset(indexNext, sBufBak);
		for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
		{
			realSize = h->searchByIndex(indexNext, indexOffset, no, tBuf, seqID);
			if(realSize < 0)
			{
				h->recoverByIndex(indexNext, indexOffset, no);
				diff = compare(sBuf, tBuf, 0, lowerSizeInChar, tagReverse);
				if(diff < similarity)
				{
					similarity = diff;
					for(i = 0; i < lowerSizeInChar; i ++)
						buf[i] = tBuf[i];
					centerSeqID = seqID;
				}
			}
		}
	}
//if the virtual center and the center are too far
	if(compare(sBuf, buf, 0, lowerSizeInChar, tagReverse) <= mismatchAllowed)
	{
		addiOut << ">" << mappingTable[centerSeqID][0] << endl;
		for(i = 0; i < lowerSizeInChar; i ++)
		{
			out << h->changeBack(buf[i]);
			addiOut << h->changeBack(buf[i]);
		}
		out << endl;
		addiOut << endl;
	}
	else
	{
		addiOut << ">" << mappingTable[sSeqID][0] << endl;
		for(i = 0; i < lowerSizeInChar; i ++)
		{
			out << h->changeBack(sBufBak[i]); //out << h->changeBack(sBuf[i]);
			addiOut << h->changeBack(sBufBak[i]); //addiOut << h->changeBack(sBuf[i]);
		}
		out << endl;
		addiOut << endl;
//still need to decide if the source sequence is reverse complementary to the virtual center
		compare(sBuf, sBufBak, 0, lowerSizeInChar, tagReverse);
	}
}

void Cluster::cluster()
{
	char sBuf[200];
	int realSize, tagReverse, i;
	unsigned int seqOffset;
	char sQVBuf[200];
	unsigned int seqID;
//	long int big = 0, small = 0;

	if(QV)
	{
		for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 5 + upperSizeInChar); seqOffset = seqOffset + (upperSizeInBit + 5 + upperSizeInChar))
		{
			realSize = h->searchBySeq(seqOffset, sBuf, sQVBuf, seqID);
			if(realSize == 0) continue;
			else realSize = 0;
			calConsensus(sBuf, seqID, tagReverse);
			if(out.is_open())
				if(reversed)
					for(i = 0; i < mappingNum[seqID]; i ++)
						out << CLID << "	" << mappingTable[seqID][i] << "	" << tagReverse << endl;
				else
					for(i = 0; i < mappingNum[seqID]; i ++)
						out << CLID << "	" << mappingTable[seqID][i] << endl;
			else
			{
				cout << "CANNOT OPEN OUTPUT FILE!" << endl;
				exit(-1);
			}
			h->deleteBySeq(seqOffset);
			seqNum ++;
			if(seqNum > adjustNum)
			{
				h->adjust();
				adjustNum = adjustNum + num / 20;
			}
//forcefully write this seq to avoid lost of it
//			numInCL = 1;
			clusterWithMismatches(sBuf, sQVBuf);
			clusterWithShifts(sBuf, sQVBuf);
//			if(numInCL > 100)
//			{
//				dis << CLID << " " << numInCL << endl;
//				big ++;
//			}
//			else if(numInCL == 1)
//			{
//				small ++;
//			}
			CLID ++;
		}
//		cout << "#clusters of more than 100 seqs is " << big <<	endl;
//		cout << "#singleton clusters is " << small << endl;
	}
	else
	{
		for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 5); seqOffset = seqOffset + (upperSizeInBit + 5))
		{
			realSize = h->searchBySeq(seqOffset, sBuf, seqID);
			if(realSize == 0) continue;
			else realSize = 0;
			calConsensus(sBuf, seqID, tagReverse);
			if(out.is_open())
				if(reversed)
					for(i = 0; i < mappingNum[seqID]; i ++)
						out << CLID << "	" << mappingTable[seqID][i] << "	" << tagReverse << endl;
				else
					for(i = 0; i < mappingNum[seqID]; i ++)
						out << CLID << "	" << mappingTable[seqID][i] << endl;
			else
			{
				cout << "CANNOT OPEN OUTPUT FILE!" << endl;
				exit(-1);
			}
			h->deleteBySeq(seqOffset);
			seqNum ++;
			if(seqNum > adjustNum)
			{
				h->adjust();
				adjustNum = adjustNum + num / 20;
			}
			clusterWithMismatches(sBuf);
			clusterWithShifts(sBuf);
			CLID ++;
		}
	}
	delete h;
	addiOut.close();
}

void Cluster::clusterWithMismatches(char sBuf[])
{
	char tBuf[200];
	int indexNext, no, realSize, tagReverse, i;
	unsigned int indexOffset;
	unsigned int seqID;

	for(indexNext = 0; indexNext < seedsCount; indexNext ++)
	{
		indexOffset = h->calOffset(indexNext, sBuf);
		for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
		{
			realSize = h->searchByIndex(indexNext, indexOffset, no, tBuf, seqID);
			if(realSize && compare(sBuf, tBuf, 0, lowerSizeInChar, tagReverse) <= mismatchAllowed)
			{
				h->deleteByIndex(indexNext, indexOffset, no);
				if(out.is_open())
					if(reversed)
						for(i = 0; i < mappingNum[seqID]; i ++)
							out << CLID << "	" << mappingTable[seqID][i] << "	" << tagReverse << endl;
					else
						for(i = 0; i < mappingNum[seqID]; i ++)
							out << CLID << "	" << mappingTable[seqID][i] << endl;
				else
				{
					cout << "CANNOT OPEN OUTPUT FILE!" << endl;
					exit(-1);
				}
				seqNum ++;
//				numInCL ++;
			}
		}
	}
}

void Cluster::clusterWithShifts(char sBuf[])
{
	char tBuf[200], buf[200];
	int no, i, lShift, rShift, realSize, indexNext, tagReverse;
	unsigned int indexOffset;
	unsigned int seqID;

	for(lShift = 1; lShift <= shiftAllowed; lShift ++)
	{
		for(i = 0; i < lowerSizeInChar - lShift; i ++)
			buf[i] = sBuf[i + lShift];
		for(indexNext = 0; indexNext < seedsCount; indexNext ++)
		{
			indexOffset = h->calOffset(indexNext, buf);
			for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
			{
				realSize = (int)h->searchByIndex(indexNext, indexOffset, no, tBuf, seqID);
				if(realSize && compare(buf, tBuf, 0, lowerSizeInChar - lShift, tagReverse) <= mismatchAllowed)//from 0 to 3 - lShift
				{
					h->deleteByIndex(indexNext, indexOffset, no);
					if(out.is_open())
						if(reversed)
							for(i = 0; i < mappingNum[seqID]; i ++)
								out << CLID << "	" << mappingTable[seqID][i] << "	" << tagReverse << endl;
						else
							for(i = 0; i < mappingNum[seqID]; i ++)
								out << CLID << "	" << mappingTable[seqID][i] << endl;
					else
					{
						cout << "CANNOT OPEN OUTPUT FILE!" << endl;
						exit(-1);
					}
//					numInCL ++;
					seqNum ++;
				}
			}
		}
	}

	for(rShift = 1; rShift <= shiftAllowed; rShift ++)
	{
		for(i = rShift; i < lowerSizeInChar; i ++)
			buf[i] = sBuf[i - rShift];
		for(indexNext = 0; indexNext < seedsCount; indexNext ++)
		{
			indexOffset = h->calOffset(indexNext, buf);
			for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
			{
				realSize = (int)h->searchByIndex(indexNext, indexOffset, no, tBuf, seqID);
				if(realSize && compare(buf, tBuf, rShift, lowerSizeInChar, tagReverse) <= mismatchAllowed)//from rShift to 3
				{
					h->deleteByIndex(indexNext, indexOffset, no);
					if(out.is_open())
						if(reversed)
							for(i = 0; i < mappingNum[seqID]; i ++)
								out << CLID << "	" << mappingTable[seqID][i] << "	" << tagReverse << endl;
						else
							for(i = 0; i < mappingNum[seqID]; i ++)
								out << CLID << "	" << mappingTable[seqID][i] << endl;
					else
					{
						cout << "CANNOT OPEN OUTPUT FILE!" << endl;
						exit(-1);
					}
//					numInCL ++;
					seqNum ++;
				}
			}
		}
	}
}

int Cluster::compare(char sBuf[], char tBuf[], int start, int end, int & tagReverse)
{
	int i, j, count = 0, reverseCount = 0;

	for(i = start; i < end; i ++)
		if(sBuf[i] != tBuf[i])
			count ++;

	if(reversed)
	{
		for(i = start, j = end - 1; i < end; i ++, j --)
			if(sBuf[i] != reverseChange(tBuf[j]))
				reverseCount ++;
		if(count <= reverseCount) 
		{
			tagReverse = 0;
			return count;
		}
		else
		{
			tagReverse = 1;
			return reverseCount;
		}
	}
	else
		return count;
}

char Cluster::reverseChange(char base)
{
	switch(base)
	{
		case 0x00: return 0x03;
		case 0x01: return 0x02;
		case 0x02: return 0x01;
		case 0x03: return 0x00;
		default: cout << "MEMORY ERROR!!" << endl; exit(-1); 
	}
}

void Cluster::clusterWithMismatches(char sBuf[], char sQVBuf[])
{
	char tBuf[200], tQVBuf[200];
	int indexNext, no, realSize, tagReverse, i;
	unsigned int indexOffset;
	unsigned int seqID;

	for(indexNext = 0; indexNext < seedsCount; indexNext ++)
	{
		indexOffset = h->calOffset(indexNext, sBuf);
		for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
		{
			realSize = h->searchByIndex(indexNext, indexOffset, no, tBuf, tQVBuf, seqID);
			if(realSize && compare(sBuf, sQVBuf, tBuf, tQVBuf, 0, lowerSizeInChar, tagReverse) <= mismatchAllowed)
			{
				h->deleteByIndex(indexNext, indexOffset, no);
				if(out.is_open())
					if(reversed)
						for(i = 0; i < mappingNum[seqID]; i ++)
							out << CLID << "	" << mappingTable[seqID][i] << "	" << tagReverse << endl;
					else
						for(i = 0; i < mappingNum[seqID]; i ++)
							out << CLID << "	" << mappingTable[seqID][i] << endl;
				else
				{
					cout << "CANNOT OPEN OUTPUT FILE!" << endl;
					exit(-1);
				}
//				numInCL ++;
				seqNum ++;
			}
		}
	}
}

void Cluster::clusterWithShifts(char sBuf[], char sQVBuf[])
{
	char tBuf[200], buf[200], tQVBuf[200], QVBuf[200];
	int no, i, lShift, rShift, realSize, indexNext, tagReverse;
	unsigned int indexOffset;
	unsigned int seqID;

	for(lShift = 1; lShift <= shiftAllowed; lShift ++)
	{
		for(i = 0; i < lowerSizeInChar - lShift; i ++)
		{
			buf[i] = sBuf[i + lShift];
			QVBuf[i] = sQVBuf[i + lShift];
		}
		for(indexNext = 0; indexNext < seedsCount; indexNext ++)
		{
			indexOffset = h->calOffset(indexNext, buf);
			for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
			{
				realSize = (int)h->searchByIndex(indexNext, indexOffset, no, tBuf, tQVBuf, seqID);
				if(realSize && compare(buf, QVBuf, tBuf, tQVBuf, 0, lowerSizeInChar - lShift, tagReverse) <= mismatchAllowed)//from 0 to 3 - lShift
				{
					h->deleteByIndex(indexNext, indexOffset, no);
					if(out.is_open())
						if(reversed)
							for(i = 0; i < mappingNum[seqID]; i ++)
								out << CLID << "	" << mappingTable[seqID][i] << "	" << tagReverse << endl;
						else
							for(i = 0; i < mappingNum[seqID]; i ++)
								out << CLID << "	" << mappingTable[seqID][i] << endl;
					else
					{
						cout << "CANNOT OPEN OUTPUT FILE!" << endl;
						exit(-1);
					}
//					numInCL ++;
					seqNum ++;
				}
			}
		}
	}

	for(rShift = 1; rShift <= shiftAllowed; rShift ++)
	{
		for(i = rShift; i < lowerSizeInChar; i ++)
		{
			buf[i] = sBuf[i - rShift];
			QVBuf[i] = sQVBuf[i - rShift];
		}
		for(indexNext = 0; indexNext < seedsCount; indexNext ++)
		{
			indexOffset = h->calOffset(indexNext, buf);
			for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
			{
				realSize = (int)h->searchByIndex(indexNext, indexOffset, no, tBuf, tQVBuf, seqID);
				if(realSize && compare(buf, QVBuf, tBuf, tQVBuf, rShift, lowerSizeInChar, tagReverse) <= mismatchAllowed)//from rShift to 3
				{
					h->deleteByIndex(indexNext, indexOffset, no);
					if(out.is_open())
						if(reversed)
							for(i = 0; i < mappingNum[seqID]; i ++)
								out << CLID << "	" << mappingTable[seqID][i] << "	" << tagReverse << endl;
						else
							for(i = 0; i < mappingNum[seqID]; i ++)
								out << CLID << "	" << mappingTable[seqID][i] << endl;
					else
					{
						cout << "CANNOT OPEN OUTPUT FILE!" << endl;
						exit(-1);
					}
//					numInCL ++;
					seqNum ++;
				}
			}
		}
	}
}

int Cluster::compare(char sBuf[], char sQVBuf[], char tBuf[], char tQVBuf[], int start, int end, int & tagReverse)
{
	int i, j, count = 0, qv = 0, reverseCount = 0, reverseQv = 0;

	for(i = start; i < end; i ++)
		if(sBuf[i] != tBuf[i] && !((sQVBuf[i] - OFFSET) + (tQVBuf[i] - OFFSET) < lowerQV))
		{
			count ++;
			qv = qv + (sQVBuf[i] - OFFSET) + (tQVBuf[i] - OFFSET);
		}

	if(reversed)
	{
		for(i = start, j = end - 1; i < end; i ++, j --)
			if(sBuf[i] != reverseChange(tBuf[j]) && !((sQVBuf[i] - OFFSET) + (tQVBuf[j] - OFFSET) < lowerQV))
			{
				reverseCount ++;
				reverseQv = reverseQv + (sQVBuf[i] - OFFSET) + (tQVBuf[j] - OFFSET);
			}
		if(count <= reverseCount)
		{
			tagReverse = 0;
			return qv > upperQV ? 10 : count;
		}
		else
		{
			tagReverse = 1;
			return reverseQv > upperQV ? 10 : reverseCount;
		}
	}
	else
	{
		return qv > upperQV ? 10 : count;
	}
}

void FileAnalyzer::outputAnalyze(int num, int correctCluster)
{
	ifstream in;
	char buf[20] = {0}, CLID[10] = {0}, seqID[10] = {0};
	int i = 0, j = 0, clusterID = -1, correctLower, correctUpper, subCluster = 0, multiCluster = 0;

	in.open("output.txt");
	if(in.is_open())
	{
		if(in.good())
		{
			in.getline(buf, 20);
			for(i = 0; i < 20; i ++)
				buf[i] = 0;
			i = 0;
		}
		while(in.good())
		{
			in.getline(buf, 20);
			if(buf[0] == 0)
				break;
			while(buf[i] != ' ')
				CLID[i ++] = buf[i];
			i ++;
			while(buf[i] != 0)
				seqID[j ++] = buf[i ++];

			if(atoi(CLID) > clusterID)
			{
				clusterID = atoi(CLID);
				correctLower = atoi(seqID) / (num / correctCluster) * (num / correctCluster);
				correctUpper = correctLower + num / correctCluster - 1;
				if(multiCluster == 0)
					subCluster ++;
				multiCluster = 0;
			}
			else
			{
				if(atoi(seqID) < correctLower || atoi(seqID) > correctUpper)
					multiCluster = 1;
			}
			for(i = 0; i < 10; i ++)
				CLID[i] = seqID[i] = 0;
			for(i = 0; i < 20; i ++)
				buf[i] = 0;
			i = j = 0;
		}
	}
	else
	{
		cout << "CANNOT OPEN INPUT FILE!" << endl;
		exit(-1);
	}
	cout << "(4) output analysis finished" << endl;
	cout << " - " << clusterID + 1 << " clusters in total" << endl;
	cout << " - " << subCluster << " sub clusters of degree " << (double)subCluster / (clusterID + 1) << endl;
}

void FileAnalyzer::inputAnalyze(char input[], int & num, int & tNum, int & lower, int & upper)
{
	ifstream in;
	int seqID = 0, i;
	char buf[201];

	num = upper = 0;
	lower = 200;
	in.open(input);

	if(in.is_open())
	{
//cont:
		while(in.good())
		{
			in.getline(buf, 201);
			if(seqID % 4 == 1)
			{
//				for(i = 0; i < in.gcount() - 1; i ++)
//					if(buf[i] == 'N') 
//					{
//						seqID ++;
//						goto cont;
//					}
				if(in.gcount() - 1 < lower) lower = in.gcount() - 1;
				if(in.gcount() - 1 > upper) upper = in.gcount() - 1;
//				num ++;
			}
			seqID ++;
		}
	}
	else
	{
		cout << "CANNOT OPEN INPUT FILE!" << endl;
		exit(-1);
	}
//	tNum = seqID;

//Filter reads based on first "lower" bases rather than all bases to avoid underestimate the number of valid reads
	in.clear(); in.seekg(0); seqID = 0; num = 0;
	if(in.is_open())
	{
conti:
		while(in.good())
		{
			in.getline(buf, 201);
			if(seqID % 4 == 1)
			{
				for(i = 0; i < lower; i ++)
					if(buf[i] == 'N')
					{
						seqID ++;
						goto conti;
					}
				num ++;
			}
			seqID ++;
		}
	}
	else
	{
		cout << "CANNOT OPEN INPUT FILE!" << endl;
		exit(-1);
	}
	tNum = seqID;
}

FastqGenerator::FastqGenerator(char input[], char output[], int num)
{
	long i;

	in.open(input);
	strcpy(addiInput, output);
	strcat(addiInput, ".fasta");
	addiIn.open(addiInput);
	strcpy(outputq, output);
	strcat(outputq, ".fastq");
	out.open(outputq);
	this->num = num;
	seq = new char [num];
	for(i = 0; i < num; i ++)
		seq[i] = 0;
}

void FastqGenerator::record()
{
	long i;
	char buf[201];

	if(addiIn.is_open())
	{
		while(addiIn.good())
		{
			addiIn.getline(buf, 201);
			if(buf[0] == 0)
				break;
			for(i = 1; i < addiIn.gcount() - 1; i ++)
				buf[i - 1] = buf[i];
			buf[i - 1] = '\0';
			seq[atoi(buf)] = 1;
			addiIn.getline(buf, 201);
		}
	}
	else
	{
		cout << "CANNOT OPEN INPUT FILE!" << endl;
		exit(-1);
	}

//	for(i = 0; i < num; i ++)
//		if(seq[i])
//			cout << i << endl;
}

void FastqGenerator::generateFastq()
{
	unsigned long seqID, i;
	char buf[201];

	record();
	if(in.is_open())
	{
		for(seqID = 0; seqID < num * 4; seqID ++)
		{
			in.getline(buf, 201);
			if(seq[seqID / 4])
			{
				for(i = 0; i < in.gcount() - 1; i ++)
					out << buf[i];
				out << endl;
			}
		}
	}
	else
	{
		cout << "CANNOT OPEN INPUT FILE!" << endl;
		exit(-1);
	}
}

Sorter::Sorter(char input[], int num, int lower)
{
	in.open(input);
	strcpy(midOutput, input);
	strcat(midOutput, ".mid.fastq");
	midOut.open(midOutput);
	this->num = num;
	mappingTable = new unsigned int * [num];
	mappingNum = new unsigned int [num];
	this->lowerSizeInChar = lower;
	this->realNum = 0;
}

int Sorter::getRealNum()
{
	return realNum;
}

unsigned int ** Sorter::getMappingTable()
{
	return mappingTable;
}

unsigned int * Sorter::getMappingNum()
{
	return mappingNum;
}

void Sorter::sort()
{
	char * seqs;
	Order * order;
	int i, j, seqID = 0, tag = 1, localNum = 0;
	char buf[201];

	seqs = new char [(long int)num * lowerSizeInChar * 2];//2 means both bases and QVs
	order = new Order [num];

	if(in.is_open())
	{
cont:
		while(in.good())
		{
			in.getline(buf, 201);
			if(seqID % 4 == 1)
				for(i = 0; i < lowerSizeInChar; i ++)
				{
					if(buf[i] == 'N')
					{
						seqID ++;
						tag = 0;
						goto cont;
					}
					else
						tag = 1;
					seqs[(long int)localNum * lowerSizeInChar * 2 + i] = buf[i];
				}
			if(seqID % 4 == 3 && tag == 1)
			{
				for(i = lowerSizeInChar; i < lowerSizeInChar * 2; i ++)
					seqs[(long int)localNum * lowerSizeInChar * 2 + i] = buf[i - lowerSizeInChar];
				order[localNum ++].realID = seqID / 4;
			}
			seqID ++;
			if(localNum == num) break;
//Must finish to avoid crash. Otherwise:
//if there are reads with N in the end, they are not counted to initialize the seqs array but are put in seqs array
		}
	}
	else
	{
		cout << "CANNOT OPEN INPUT FILE!" << endl;
		exit(-1);
	}
	for(i = 0; i < localNum; i ++)
		order[i].ID = i;
//verification
//	cout << "-------------------------------------" << endl;
//	for(i = 0; i < localNum; i ++)
//		cout << order[i].ID << "|" << order[i].realID << " ";
//	cout << endl;
//	cout << "-------------------------------------" << endl;
//	for(i = 0; i < localNum; i ++)
//	{
//		for(j = 0; j < lowerSizeInChar * 2; j ++)
//			cout << seqs[i * lowerSizeInChar * 2 + j];
//		cout << endl;
//	}
//	cout << "-------------------------------------" << endl;
//verification

	suffixSort(0, localNum - 1, 0, seqs, order);
//verification
//	cout << "-------------------------------------" << endl;
//	for(i = 0; i < localNum; i ++)
//		cout << order[i].ID << "|" << order[i].realID << " ";
//	cout << endl;
//	cout << "-------------------------------------" << endl;
//	for(i = 0; i < realNum; i ++)
//	{
//		cout << mappingNum[i] << ": ";
//		for(j = 0; j < mappingNum[i]; j ++)
//			cout << mappingTable[i][j] << " ";
//		cout << endl;
//	}
//	cout << "-------------------------------------" << endl;
//verification
	delete seqs;
	delete order;
}

void Sorter::suffixSort(int start, int end, int depth, char seqs[], Order order[])
{
	int i, seqID = 0, tag;
	Order * buf;
	int s[4], e[4];
	int startBuf;
	int j;

	if(start == -1 && end == -1)
		return;

	if(start == end || depth == lowerSizeInChar - 1)
	{
		mappingTable[realNum] = new unsigned int [end - start + 1];
		for(i = start; i <= end; i ++)
			mappingTable[realNum][i - start] = order[i].realID;
		mappingNum[realNum] = end - start + 1;

		if(midOut.is_open())
		{
			midOut << "@" << realNum << endl;
			for(i = 0; i < lowerSizeInChar; i ++)
				midOut << seqs[(long int)order[start].ID * lowerSizeInChar * 2 + i];
			midOut << endl;
			midOut << "+" << realNum << endl;
			for(i = lowerSizeInChar; i < lowerSizeInChar * 2; i ++)
				midOut << seqs[(long int)order[start].ID * lowerSizeInChar * 2 + i];
			midOut << endl;
		}
		else
		{
			cout << "CANNOT OPEN OUTPUT FILE!" << endl;
			exit(-1);
		}

		realNum ++;
		return;
	}

	buf = new Order [end - start + 1];

	startBuf = start;
	tag = 0;
	for(i = start; i <= end; i ++)
		if(seqs[(long int)order[i].ID * lowerSizeInChar * 2 + depth] == 'A')
		{
			tag = 1;
			buf[seqID].ID = order[i].ID;
			buf[seqID ++].realID = order[i].realID;
		}
	if(tag == 1)
	{
		s[0] = startBuf;
		e[0] = start + seqID - 1;
		startBuf = start + seqID;
	}
	else
		s[0] = e[0] = -1;

	tag = 0;
	for(i = start; i <= end; i ++)
		if(seqs[(long int)order[i].ID * lowerSizeInChar * 2 + depth] == 'C')
		{
			tag = 1;
			buf[seqID].ID = order[i].ID;
			buf[seqID ++].realID = order[i].realID;
		}
	if(tag == 1)
	{
		s[1] = startBuf;
		e[1] = start + seqID - 1;
		startBuf = start + seqID;
	}
	else
		s[1] = e[1] = -1;

	tag = 0;
	for(i = start; i <= end; i ++)
		if(seqs[(long int)order[i].ID * lowerSizeInChar * 2 + depth] == 'G')
		{
			tag = 1;
			buf[seqID].ID = order[i].ID;
			buf[seqID ++].realID = order[i].realID;
		}
	if(tag == 1)
	{
		s[2] = startBuf;
		e[2] = start + seqID - 1;
		startBuf = start + seqID;
	}
	else
		s[2] = e[2] = -1;

	tag = 0;
	for(i = start; i <= end; i ++)
		if(seqs[(long int)order[i].ID * lowerSizeInChar * 2 + depth] == 'T')
		{
			tag = 1;
			buf[seqID].ID = order[i].ID;
			buf[seqID ++].realID = order[i].realID;
		}
	if(tag == 1)
	{
		s[3] = startBuf;
		e[3] = start + seqID - 1;
		startBuf = start + seqID;
	}
	else
		s[3] = e[3] = -1;

	for(i = start, seqID = 0; i <= end; i ++, seqID ++)
	{
		order[i].ID = buf[seqID].ID;
		order[i].realID = buf[seqID].realID;
	}

//	cout << s[0] << ", " << s[1] << ", " << s[2] << ", " << s[3] << endl;
//	cout << e[0] << ", " << e[1] << ", " << e[2] << ", " << e[3] << endl;

	suffixSort(s[0], e[0], depth + 1, seqs, order);
	suffixSort(s[1], e[1], depth + 1, seqs, order);
	suffixSort(s[2], e[2], depth + 1, seqs, order);
	suffixSort(s[3], e[3], depth + 1, seqs, order);
}

char change(char base)
{
	switch(base)
	{
		case 0x00: return 'A';
		case 0x01: return 'C';
		case 0x02: return 'G';
		case 0x03: return 'T';
		default: cout << "UNKNOWN ERROR!" << endl; exit(-1);
	}
}

bool within(int p, int pos[], int size)
{
	int i;
	for(i = 0; i < size; i ++)
		if(p == pos[i]) return true;
	return false;
}

void introduceMismatches(char sBuf[], char buf[], int mismatch, int size)
{
	int i, j;
	int pos[200], mBuf[200] = {0};

	for(i = 0; i < mismatch; i ++)
	{
		do
			pos[i] = rand() % size;
		while(mBuf[pos[i]] == 1);
		mBuf[pos[i]] = 1;
	}
	for(i = 0, j = 0; i < size; i ++)
	{
		if(within(i, pos, mismatch) && j < mismatch)
		{
			do
				buf[i] = change((char)(rand() % 4));
			while(sBuf[i] == buf[i]);
			j ++;
		}
		else
			buf[i] = sBuf[i];
	}
}

void introduceShifts(char sBuf[], char buf[], int shift, int size)
{
	int i;

	if(shift < 0)
	{
		for(i = 0; i < size - abs(shift); i ++)
			buf[i] = sBuf[i + abs(shift)];
		for(i = size - abs(shift); i < size; i ++)
			buf[i] = change((char)(rand() % 4));
	}
	else
	{
		for(i = 0; i < shift; i ++)
			buf[i] = change((char)(rand() % 4));
		for(i = shift; i < size; i ++)
			buf[i] = sBuf[i - shift];
	}
}

#ifdef WITHSIMILARITY
void generateClusteredSeq(int num, int lower, int upper, int mismatchAllowed, int shiftAllowed, int correctCluster, int distance)
#else
void generateClusteredSeq(int num, int lower, int upper, int mismatchAllowed, int shiftAllowed, int correctCluster)
#endif
{
	ofstream out;
	int size, mismatch, shift, i, j, k;
	char s[200], sBuf[200], mismatchBuf[200], shiftBuf[200], buf[200];
	int mBuf[200] = {0};

#ifdef WITHSIMILARITY
	if(distance < 2) 
	{
		cout << "incorrect distance" << endl;
		return;
	}
	for(i = 0; i < 200; i ++)
		s[i] = change((char)(rand() % 4));
#endif
	out.open("input.txt");

	for(k = 0; k < correctCluster; k ++)
	{
#ifdef WITHSIMILARITY
		generateCenter(s, sBuf, distance, mBuf);
#else
		for(i = 0; i < 200; i ++)
			sBuf[i] = change((char)(rand() % 4));
#endif
		for(i = 0; i < num / correctCluster; i ++)
		{
			size = lower + rand() % (upper - lower + 1);
			mismatch = rand() % (mismatchAllowed + 1);
			introduceMismatches(sBuf, mismatchBuf, mismatch, size);
			shift = rand() % (shiftAllowed * 2 + 1) - shiftAllowed;
			introduceShifts(mismatchBuf, shiftBuf, shift, size);
//			if(rand() % 2 == 1)
//				introduceMismatches(sBuf, buf, mismatch, size);
//			else
//				introduceShifts(sBuf, buf, shift, size);
			out << "@title " << k * num / correctCluster + i << " size = " << size << " mismatches = " << mismatch << " shifts = " << shift << endl;
			out.write(shiftBuf, size);
//			out.write(buf, size);
			out << endl;
			out << "+title " << k * num / correctCluster + i << " size = " << size << " mismatches = " << mismatch << " shifts = " << shift << endl;
			for(j = 0; j < size; j ++)
				shiftBuf[j] = (char)(rand() % RANGE) + OFFSET;
			out.write(shiftBuf, size);
			out << endl;
		}
	}
}

void itoa(char buf[], unsigned int v)
{
	if(v / 1000)
	{
		buf[0] = v / 1000 + 48;
		buf[1] = (v % 1000) / 100 + 48;
		buf[2] = (v % 100) / 10 + 48;
		buf[3] = v % 10 + 48;
		buf[4] = '\0';
	}
	else if(v / 100)
	{
		buf[0] = v / 100 + 48;
		buf[1] = (v % 100) / 10 + 48;
		buf[2] = v % 10 + 48;
		buf[3] = '\0';
	}
	else if(v / 10)
	{
		buf[0] = v / 10 + 48;
		buf[1] = v % 10 + 48;
		buf[2] = '\0';
	}
	else
	{
		buf[0] = v + 48;
		buf[1] = '\0';
	}
}

void print()
{
	cout << "SEED --input input.fastq --output output.txt [--mismatch M] [--shift S] [--QV1 L] [--QV2 U] [--fast/short] [--reverse]" << endl;
	cout << "--mismatch is the maximum number of mismatches allowed from the center sequence in each cluster (0 - 3, default 3)" << endl;
	cout << "--shift is the maximum number of shifts allowed from the center sequence in each cluster (0 - 6, default 3)" << endl;
	cout << "--QV1 is the threshold for the base call quality values (QV) that are provided in the FASTQ files as Phred scores. SEED ignores those mismatches where the sum of the Phred scores of the mismatching bases is lower than the specified QV1 threshold value (0 - 2 * 93). The default value for QV1 is 0" << endl;
	cout << "--QV2 is another QV threshold. It prevents co-clustering of sequences where the sum of all mismatched positions is higher than the threshold value (0 - 6 * 93). The default value for QV2 is 6 * 93" << endl;
	cout << "--fast uses a bigger spaced seed weight to save running time. It is only applicable for sequences longer than 58 bp and may need more memory" << endl;
	cout << "--short is to use a smaller spaced seeds weight for sequences as short as 21 bp. This setting often results in longer compute times" << endl;
	cout << "--reverse is to co-cluster sequences in sense and anti-sense orientation (reverse and complement)" << endl;
}

int main(int argc, char * argv[])
{
	time_t start, end;
	int num, tNum, lower, upper, mismatch = 3, shift = 3, lowerQV = 0, upperQV = 6 * 93, i, tagMismatch = 0, tagShift = 0, tagInput = 0, tagOutput = 0, tagFast = 0, tagShort = 0, tagReverse = 0;
	int tagQV1 = 0, tagQV2 = 0;
	char buf[5], input[100], output[100], midOutput[100];
	ifstream in;
	int io = 0;
	int totalLength, count, j;

	for(i = 1; i < argc; i ++)
		if(strcmp(argv[i], "--input") == 0)
		{
			if(tagInput == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			in.open(argv[++ i]);
			if(!in.is_open())
			{
				cout << "CANNOT OPEN INPUT FILE!" << endl;
				print();
				return 0;
			}
			in.close();
			strcpy(input, argv[i]);
			tagInput = 1;
		}
		else if(strcmp(argv[i], "--output") == 0)
		{
			if(tagOutput == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			strcpy(output, argv[++ i]);
			tagOutput = 1;
		}
		else if(strcmp(argv[i], "--mismatch") == 0)
		{
			if(tagMismatch == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			mismatch = atoi(argv[++ i]);
			itoa(buf, mismatch);
			if(strcmp(argv[i], buf) != 0)
			{
				print();
				return 0;
			}
			tagMismatch = 1;
		}
		else if(strcmp(argv[i], "--shift") == 0)
		{
			if(tagShift == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			shift = atoi(argv[++ i]);
			itoa(buf, shift);
			if(strcmp(argv[i], buf) != 0)
			{
				print();
				return 0;
			}
			tagShift = 1;
		}
		else if(strcmp(argv[i], "--QV1") == 0)
		{
			if(tagQV1 == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			lowerQV = atoi(argv[++ i]);
			itoa(buf, lowerQV);
			if(strcmp(argv[i], buf) != 0)
			{
				print();
				return 0;
			}
			tagQV1 = 1;
			QV = 1;
		}
		else if(strcmp(argv[i], "--QV2") == 0)
		{
			if(tagQV2 == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			upperQV = atoi(argv[++ i]);
			itoa(buf, upperQV);
			if(strcmp(argv[i], buf) != 0)
			{
				print();
				return 0;
			}
			tagQV2 = 1;
			QV = 1;
		}
		else if(strcmp(argv[i], "--fast") == 0)
		{
			if(tagFast == 1 || tagShort == 1)
			{
				print();
				return 0;
			}
			tagFast = 1;
			seedsCount = 4;
			seedsWeight = 64 * 1024;
		}
		else if(strcmp(argv[i], "--short") == 0)
                {
                        if(tagFast == 1 || tagShort == 1)
                        {
                                print();
                                return 0;
                        }
                        tagShort = 1;
                        seedsWeight = 4;
                }
		else if(strcmp(argv[i], "--reverse") == 0)
		{
			if(tagReverse == 1)
			{
				print();
				return 0;
			}
			tagReverse = 1;
			reversed = 1;
		}
		else
		{
			print();
			return 0;
		}

	if(tagInput == 0 || tagOutput == 0 || mismatch < 0 || mismatch > 3 || shift < 0 || shift > 6 || lowerQV < 0 || lowerQV > 2 * 93 || upperQV < 0 || upperQV > 6 * 93)
	{
		print();
		return 0;
	}

	if(QV)
		cout << "#mismatch = " << mismatch << "; #shift = " << shift << "; QV1 = " << lowerQV << "; QV2 = " << upperQV << endl;
	else
		cout << "#mismatch = " << mismatch << "; #shift = " << shift << endl;

//	generateClusteredSeq(10000, 195, 200, 3, 3, 100);

	FileAnalyzer fa;
	start = time(NULL);
	fa.inputAnalyze(input, num, tNum, lower, upper);
	cout << "(1) input analysis finished" << endl;

	cout << " - " << num << " seqs with lengths between " << lower << " and " << upper << endl;
	if(num == 0 || lower < 36 && seedsWeight == 1024 * 16 || lower < 58 && seedsWeight == 1024 * 64 || lower < 21 && seedsWeight == 4  || upper > 200 || upper - lower > 5)
	{
		cout << "INVALID INPUT FILE. VALID FILE SHOULD HAVE READS OF LENGTH 21-200 WITH MAX VARIATION 5!" << endl;
		return 0;
	}

//	produce realNum, mappingTable and mappingNum here, and the intermediate file is produced/opened by protocol
	Sorter s(input, num, lower);
	s.sort();
	cout << "(2) sorting finished" << endl;

	Cluster c(input, output, s.getRealNum(), lower, upper, mismatch, shift, lowerQV, upperQV, s.getMappingTable(), s.getMappingNum());
	cout << "(3) init finished" << endl;

	c.cluster();
	end = time(NULL);
	cout << "(4) clustering finished" << endl;

	FastqGenerator f(input, output, tNum);
	f.generateFastq();
	cout << "(5) fastq file generated" << endl;
	
	cout << " - " << end - start << " seconds" << endl;

	return 1;
}
