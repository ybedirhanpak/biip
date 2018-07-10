#include "bigraph2.h"
#include "simple_timer.h"
#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <sys/time.h>

using namespace std;

typedef int LT;		// type of left column. change the type here according to the data
typedef int RT;		// type of right column.

long int returnMiliSeconds()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
}

double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*10)/CLOCKS_PER_SEC;
	return diffms;
}

int main(int argc, char ** argv)
{
	FILE *fptr2;
	long int timeBefore = 0;
	long int timeAfter = 0;
	double solnTime = 0;

  if (argc != 4) {
    cout << "Usage: " << argv[0] 
	 << " <bigraph file> <biclique file> <size dist file>" << endl;
    return 0;
  }

  // change the type here according to the input file.
  SimpleBigraph<LT, RT> SBG;

  ofstream bicliq_fs(argv[2]), size_fs(argv[3]);

  if (!bicliq_fs || !size_fs) {
    cout << "Error output files!" << endl;
    return 0;
  }

  Timer T;

  cout << "Read edges: " << SBG.read(argv[1]) << endl;
  cout << "Left: " << SBG.l_size() << endl;
  cout << "Right: " << SBG.r_size() << endl;
	
	clock_t begin=clock();
	timeBefore = returnMiliSeconds();

  cout << T.start() << endl;
  //SBG.mica(argv[2], cout);
  SBG.mica(bicliq_fs, size_fs, cout);
  cout << endl;
  cout << T.stop() << endl;
  cout << T.report() << endl;

	timeAfter = returnMiliSeconds();
	clock_t end= clock();
	solnTime = double(diffclock(end,begin));

  bicliq_fs.close();
  size_fs.close();

  	fptr2 = fopen( "statusMica.txt", "w" );

	fprintf( fptr2, "%d\t%d\t%d\t%d\t%lf\t%lf", SBG.numberOfBiclusters, SBG.maxL, SBG.maxK, 1, solnTime/*double(diffclock(end,begin))*/, (timeAfter - timeBefore) / 1000.0  );
	fclose( fptr2 );

  return 0;
}
