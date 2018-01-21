/*

  Function returning the column-wise minimum and minimum index. For easier integration with R, the index
will also be stored as a double. NA's are ignored.

*/


#include "minWhichMin.h"

void minWhichMin(double * matrix, int * nRows, int * nColumns, double * min, double * whichMin)
{
  int nrows = *nRows, ncols = *nColumns;

  for (int i=0; i<ncols; i++)
  {
    double * col = matrix + i*nrows;
    double curmin = *col;
    double curwhich = 0;
    for (int j=1; j<nrows; j++)
    {
      col++;
      if ( ISNAN(curmin) || (!ISNAN(*col) && (*col < curmin))) { curmin = *col; curwhich = (double) j; }
    }
    min[i] = curmin;
    whichMin[i] = curwhich;
  }
}

