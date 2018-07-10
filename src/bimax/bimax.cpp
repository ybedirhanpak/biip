/*-------------------------------------------------------------------------------------
 * bimax.c: C implementation for the Bimax algorithm as described in the paper
 *          "A Systematic Comparison and Evaluation of Biclustering Methods for Gene
 *          Expression Data" by A. Prelic et al., Bioinformatics, 2006.
 *
 * USAGE:   Compile the program using any C compiler, e.g., gcc -O4 -o bimax bimax.c
 *          under Linux. The program assumes the data matrix to be contained in the
 *          file 'matrix.txt' in the same directory; this file is structured as follows:
 *
 *          - the first number specifies the number n of rows
 *          - the second number number specifies the number m of columns
 *          - the third number defines the minimal number of rows a bicluster must
 *            contain
 *          - the fourth number defines the minimal number of columns a bicluster must
 *            contain
 *          - the succeeding numbers, which are either 0s or 1s, represent the contents
 *            of the data matrix where the order is first columnwise and then rowwise,
 *            i.e., <row_1_column_1> <row_1_column_2> ... <row_1_column_m>
 *            <row_2_column_1> ... <row_2_column_m> ... <row_n_column_m>
 *
 *          The biclusters found are outputted on the terminal, where for each
 *          bicluster three lines are printed:
 *
 *          - the first line contains the number of the bicluster
 *          - the second line contains the selected rows (tab delimited)
 *          - the third line contains the selected columns (tab delimited)
 *
 *          Subsequent biclusters are separated by blank lines.
 *
 * (c)2005/6 Eckart Zitzler, ETH Zurich, Switzerland 
 *-------------------------------------------------------------------------------------*/

#include<limits.h>
#include<stdio.h>
#include <iostream>
#include<stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <vector>

using namespace std;

#define LINUX

typedef unsigned long int  bitvector_t;
typedef bitvector_t        *cs_t;

typedef enum {identity, complement}  cmode_t;

typedef struct row_s {
  long  originalRowNumber;
  cs_t  columnSet;
} row_t;

int          bitsPerBV;
int          noBVs;
bitvector_t  bitMaskLastBV;
int          biclusterCount;
long   noRows;
long   noColumns;
long   minNoRows;
long   minNoColumns;
long   maxLevels;
row_t  *rows;
cs_t   *consideredColumns;
cs_t   *mandatoryColumns;
cs_t   columnIntersection;
std::vector<std::vector<int> > store;

int  isSet(cs_t  columnSet, long  column)
{
  bitvector_t  bv;

  if (column >= 0L && column < noColumns) {
    bv = 1U << (column % bitsPerBV);
    return ((columnSet[column / bitsPerBV] & bv) != 0);
  }
  return 0;
} /* setColumn */

void  setColumn(cs_t  columnSet, long  column)
{
  bitvector_t  bv;

  if (column >= 0L && column < noColumns) {
    bv = 1U << (column % bitsPerBV);
    columnSet[column / bitsPerBV] |= bv;
  }
} /* setColumn */

void  unsetColumn(cs_t  columnSet, long  column)
{
  bitvector_t  bv;

  if (column >= 0L && column < noColumns) {
    bv = ~(~(columnSet[column / bitsPerBV]) | (1U << (column % bitsPerBV)));
    columnSet[column / bitsPerBV] &= bv;
  }
} /* unsetColumn */

long  columnCount(cs_t  columnSet)
{
  long         i, j;
  long         counter;
  bitvector_t  bv;
  
  columnSet[noBVs - 1] &= bitMaskLastBV;
  counter = 0L;  
  for (i = noBVs - 1; i >= 0; i--) {
    bv = columnSet[i];
    if (bv != 0U) {
      for (j = 0L; j < bitsPerBV; j++) {
	if (bv & 1U)  counter++;
	bv >>= 1;
      }
    }
  }
  return counter;
} /* columnCount */

int  compareColumns(cs_t  columnSetA, cs_t  columnSetB, cs_t  mask)
{
  int          i;
  int          contained, disjoint;
  bitvector_t  bitMask, sharedColumns;

  contained = disjoint = 1;
  bitMask = bitMaskLastBV;
  for (i = noBVs - 1; i >= 0; i--) {
    sharedColumns = ((columnSetA[i] & columnSetB[i]) & mask[i]) & bitMask;
    if ((sharedColumns | columnSetB[i]) != sharedColumns)
      contained = 0;
    if (sharedColumns != 0)
      disjoint = 0;
    bitMask = ~0U;
  }
  if (contained && disjoint)
    return -2; /* either set is empty */
  if (contained)
    return -1; /* set A contains set B */
  if (disjoint)
    return 1; /* set A and set B are disjoint */
  return 0; /* set B is larger than set A and the intersection is not empty */
} /* compareColumns */

void  copyColumnSet(cs_t  columnSet, cs_t  columnSetDest, cmode_t  copyMode)
{
  int  i;

  for (i = noBVs - 1; i >= 0; i--)
    if (copyMode == complement)
      columnSetDest[i] = ~columnSet[i];
    else
      columnSetDest[i] = columnSet[i];
} /* copyColumnSet */

void  intersectColumnSets(cs_t  columnSetA, cs_t  columnSetB, cs_t  columnSetDest)
{
  int  i;

  for (i = noBVs - 1; i >= 0; i--)
    columnSetDest[i] = (columnSetA[i] & columnSetB[i]);
} /* intersectColumnSets */

void  determineColumnsInCommon(long  firstRow, long  lastRow, cs_t  sharedColumnSet)
{
  int   i;
  long  j;

  if (firstRow >= 0L && lastRow >= firstRow && lastRow < noRows) {
    for (i = noBVs - 1; i >= 0; i--) {
      sharedColumnSet[i] = ~0U;
      for (j = firstRow; j <= lastRow; j++)
	sharedColumnSet[i] &= rows[j].columnSet[i];
    }
  }
} /* determineColumnsInCommon */

int  containsMandatoryColumns(cs_t  columnSet, int  noSets)
{
  int   contains, j;
  long  i;

  contains = 1;
  for (i = 0; i < noSets; i++) {
    if ((mandatoryColumns[i][noBVs - 1] & columnSet[noBVs - 1] & bitMaskLastBV) == 0U) {
      j = noBVs - 2;
      while (j >= 0 && (mandatoryColumns[i][j] & columnSet[j]) == 0U)
	j--;
      if (j < 0) {
	contains = 0;
	i = noSets;
      }
    }
  }
  return contains;
} /* containsMandatoryColumns */

void  swapRows(long  a, long  b)
{
  long   tempOriginalRowNumber;
  cs_t  tempColumnSet;

  if (a != b && a >= 0L && a < noRows && b >= 0L && b < noRows) {
    tempOriginalRowNumber = rows[a].originalRowNumber;
    tempColumnSet = rows[a].columnSet;
    rows[a].originalRowNumber = rows[b].originalRowNumber;
    rows[a].columnSet = rows[b].columnSet;
    rows[b].originalRowNumber = tempOriginalRowNumber;
    rows[b].columnSet = tempColumnSet;
  }
} /* swapRows */

long  chooseSplitRow(long  firstRow, long  lastRow, int  level)
{
  long  i;

  for (i = firstRow; i <= lastRow &&
	 compareColumns(rows[i].columnSet, consideredColumns[level],
			consideredColumns[0]) < 0; i++); 
  if (i <= lastRow)
    return i;
  return firstRow;
} /* chooseSplitRow */

long  selectRows(long  firstRow, long  lastRow, long  level, int  *overlapping)
{
  long  selected;

  selected = 0L;
  *overlapping = 0;
  while (firstRow <= lastRow) {
    switch (compareColumns(consideredColumns[level], rows[firstRow].columnSet,
			   consideredColumns[level - 1L])) {
    case -2:
    case 1:
      swapRows(lastRow, firstRow);      
      lastRow--;
      break;
    case 0:
      *overlapping = 1;
    default:
      selected++;
      firstRow++;
      break;
    }
  }
  return selected;
} /* selectRows */

void  writeBicluster(long  firstRow, long  lastRow, cs_t  columnSet)
{
  static long  biclusterCounter = 0;
  long  i;
  FILE *fptr;
  bool flag = false;
  long  j;
  /* to check bicluster content */
  for (i = firstRow; i <= lastRow; i++)
  {
    for (j = 0L; j < noColumns; j++) {
      if (isSet(columnSet, j)){
        if (store[rows[i].originalRowNumber][j] ){
		flag = true;
		break;
	}
      }
    }
    if ( flag == true )
	break;
  }
  flag = false;

  if ( flag == false ){
  fptr = fopen( "BIMAXResult.txt", "a" );

  biclusterCounter++;
  //cout << "Bicluster " << biclusterCount << " is calculated for BIMAX\n";
//   printf("\n%ld\n", biclusterCounter);
  int count = 0;
  for (i = 0; i < noColumns; i++)
    if (isSet(columnSet, i))
      count++;

  //cout << "K:" << lastRow - firstRow + 1 << " - " << "L:" << count << endl;
  fprintf( fptr, "%d\t%d\n", lastRow - firstRow + 1, count );
  for (i = firstRow; i <= lastRow; i++){
    fprintf( fptr, "%d\t", rows[i].originalRowNumber + 1L ); 
//     printf("%ld\t", rows[i].originalRowNumber + 1L);
  }
//   printf("\n");
  fprintf( fptr, "\n" );
  for (i = 0; i < noColumns; i++)
    if (isSet(columnSet, i)){
//       printf("%ld\t", i + 1L);
      fprintf( fptr, "%d\t", i );
    }
  fprintf( fptr, "\n" );
  fclose( fptr );
//   printf("\n");
  
// rows[i].columnSet

// Remove following block not to show bicluster

  long  j;
  for (i = firstRow; i <= lastRow; i++)
  {
    cout << i << ":";
    for (j = 0L; j < noColumns; j++) {
      if (isSet(columnSet, j))
	cout << store[rows[i].originalRowNumber][j];
    }
    cout << endl;
  }


  	biclusterCount--;
  }
} /* writeBicluster */

void  conquer(long  firstRow, long  lastRow, long  level, long noMandatorySets)
{
  int   overlapping = 0;
  long  splitRow, noSelectedRows;
  if( biclusterCount != 0 ){
	determineColumnsInCommon(firstRow, lastRow, columnIntersection);
	if (compareColumns(columnIntersection, consideredColumns[level],
			  consideredColumns[level]) == -1){
	  writeBicluster(firstRow, lastRow, columnIntersection);	}
	else {
	  splitRow = chooseSplitRow(firstRow, lastRow, level);
	  intersectColumnSets(consideredColumns[level], rows[splitRow].columnSet,
			      consideredColumns[level + 1L]);
	  if (columnCount(consideredColumns[level + 1L]) >= minNoColumns &&
	      containsMandatoryColumns(consideredColumns[level + 1L], noMandatorySets)) {
	    noSelectedRows = selectRows(firstRow, lastRow, level + 1L, &overlapping);
	    if (noSelectedRows >= minNoRows)
	      conquer(firstRow, firstRow + noSelectedRows - 1L, level + 1L, noMandatorySets);
	  }
	  copyColumnSet(consideredColumns[level + 1L], consideredColumns[level + 1L],
			complement);
	  intersectColumnSets(consideredColumns[level], consideredColumns[level + 1L],
			      consideredColumns[level + 1L]);
	  if (overlapping) {
	    copyColumnSet(consideredColumns[level + 1L], mandatoryColumns[noMandatorySets],
			  identity);
	    noMandatorySets++;
	  }
	  noSelectedRows = selectRows(firstRow, lastRow, level + 1L, &overlapping);
	  copyColumnSet(consideredColumns[level], consideredColumns[level + 1L], identity);
	  if (noSelectedRows >= minNoRows)
	    conquer(firstRow, firstRow + noSelectedRows - 1L, level + 1L, noMandatorySets);
	}
  }
} /* conquer */

int  initialize()
{
  bitvector_t  dummy;
  int          failed;
  long         i;

  /* initilization for handling bit vectors */
  dummy = 1;
  bitsPerBV = 0;
  while (dummy != 0) {
    dummy <<= 1;
    bitsPerBV++;
  }
  bitMaskLastBV = (~0U >> (bitsPerBV - (noColumns % bitsPerBV)));
  noBVs = (noColumns / bitsPerBV) + ((noColumns % bitsPerBV) == 0 ? 0 : 1);

  /* memory allocation */
  failed = 0;
  rows = (row_t*)malloc(sizeof(row_t) * noRows);
  if (rows == NULL)  failed = 1;
  for (i = 0L; i < noRows; i++) {
    rows[i].originalRowNumber = i;
    rows[i].columnSet = (bitvector_t*)calloc(sizeof(bitvector_t), noBVs);
    if (rows[i].columnSet == NULL)
      failed = 1;
  }
  maxLevels = (noRows + 2L);
  consideredColumns = (cs_t*)calloc(sizeof(cs_t), maxLevels);
  if (consideredColumns == NULL)  failed = 1;
  else {
    for (i = 0L; i < maxLevels; i++) {
      consideredColumns[i] = (bitvector_t*)calloc(sizeof(bitvector_t), noBVs);
      if (consideredColumns[i] == NULL)  failed = 1;
    }
    if (!failed)
      for (i = 0L; i < noColumns; i++)
	setColumn(consideredColumns[0], i);
  }
  mandatoryColumns = (cs_t*)calloc(sizeof(cs_t), maxLevels);
  if (mandatoryColumns == NULL)  failed = 1;
  else {
    for (i = 0L; i < maxLevels; i++) {
      mandatoryColumns[i] = (bitvector_t*)calloc(sizeof(bitvector_t), noBVs);
      if (mandatoryColumns[i] == NULL)  failed = 1;
    }
  }
  columnIntersection = (bitvector_t*)calloc(sizeof(bitvector_t), noBVs);
  if (columnIntersection == NULL)  failed = 1;

  return !failed;
} /* initializeMemory */

void  readInDataMatrix(FILE  *fp)
{
  long  i, j, cell;
// 	cout << endl;
  for (i = 0L; i < noRows; i++) {
    std::vector<int> tmp_vStoreRow;
    for (j = 0L; j < noColumns; j++) {
      fscanf(fp, "%d", &cell);
      tmp_vStoreRow.push_back(cell);
      if ( cell != 0 )
 	cout << i << " - " << j << endl;	
// 	cout << cell << "\t";
      if (cell == 0)
	unsetColumn(rows[i].columnSet, j);
      else
	setColumn(rows[i].columnSet, j);
    }
    store.push_back(tmp_vStoreRow);
// 	cout << endl;
  }
} /* readInDataMatrix */


double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*10)/CLOCKS_PER_SEC;
	return diffms;
}

long int returnMiliSeconds()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
}

int main(int argc, char *argv[])
{
  double solnTime = 0;
  int numberOfBiclusters = 0;
  int maxDim1 = 0, maxDim2 = 0;
  long int timeBefore = 0;
  long int timeAfter = 0;

  biclusterCount = 10000000;
  FILE  *fp;
  FILE *fptr, *fptr2;
  fptr = fopen( "BIMAXResult.txt", "w" );

  cout << "/**************************************************/" << endl;
  cout << "	       BIMAX'S RUN " 			 << endl;
  cout << "/**************************************************/" << endl;
  fclose( fptr );
	/*
  if( fp = fopen( "ex6", "r") )
  {
    int bipLayer1Max = 0;
	int bipLayer2Max = 0;
	int row, column;
	while( !feof(fp) )
	{
		fscanf(fp, "%d", &row);
		fscanf(fp, "%d", &column);
		if ( row > bipLayer1Max )

	}
	fclose( fp );
  }
  return 1;*/
  fp = fopen( argv[1], "r");

  if (fp == NULL)  exit(0);
  fscanf(fp, "%ld", &noRows);
  fscanf(fp, "%ld", &noColumns);
  fscanf(fp, "%ld", &minNoRows);
  fscanf(fp, "%ld", &minNoColumns);
//   cout << noRows << " - " << noColumns << " - " << minNoRows << " - " << minNoColumns << endl;
  if (minNoRows < 1L)
    minNoRows = 1L;
  if (minNoColumns < 1L)
    minNoColumns = 1L;
  if (noColumns > 0L && noRows > 0L && initialize()) {
    
	readInDataMatrix(fp);
	clock_t begin=clock();
	timeBefore = returnMiliSeconds();
    	conquer(0L, noRows - 1L, 0L, 0L);
	//cout << "Bimax finished" << endl;
        timeAfter = returnMiliSeconds();
	clock_t end=clock();
	solnTime = double(diffclock(end,begin));

	if ( fptr2 = fopen( "BIMAXResult.txt", "r" ) )
	{
		char strToRead[16] = "";
		long maxCase = 0;
		int dim1 = 0, dim2 = 0;
		cout << "Check File until the end" << endl;

		while( !feof( fptr2 ) )
		{
			fscanf(fptr2,"%d%d", &dim1, &dim2 );
			//cout << dim1 << " - " << dim2 << endl;
			for ( int count = 0; count < dim1; count++ )
			{
				fscanf(fptr2,"%s",strToRead);
			}

			for ( int count = 0; count < dim2; count++ )
			{
				fscanf(fptr2,"%s",strToRead);
			}

			if ( dim1 * dim2 > maxCase )
			{
				maxDim1 = dim1;
				maxDim2 = dim2;
				maxCase = dim1 * dim2;
			}
			numberOfBiclusters++;

		}
		numberOfBiclusters--;
		fclose(fptr2);	
	}
	else
	{
		cout << "File Could not Be opened" << endl;
	}


	fptr2 = fopen( "statusBimax.txt", "w" );

	fprintf( fptr2, "%d\t%d\t%d\t%d\t%lf\t%lf", numberOfBiclusters, maxDim1, maxDim2, 1, solnTime/*double(diffclock(end,begin))*/, (timeAfter - timeBefore) / 1000.0  );
	fclose( fptr2 );

  }

  return 0;
} /* main */
