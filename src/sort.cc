// sort routines

#include "sort.h"

// $Id: sort.cc,v 1.2 2001/11/29 00:04:43 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/sort.cc,v 1.2 2001/11/29 00:04:43 bruns Exp $
// $Log: sort.cc,v $
// Revision 1.2  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

extern "C"
{
#include <stdlib.h> // random() on linux
#include <math.h> // random() on sgi
#include <sys/time.h> // gettimeofday()
#include <unistd.h> // gettimeofday()
}
// Sun5 does not have random()
#ifdef NO_RANDOM
#define RANDOM rand
#define SRANDOM srand
#else
#define RANDOM random
#define SRANDOM srandom
#endif

// This quicksort will return a sorted array, and a similarly sorted
// array of integers, to preserve the re-ordering information.
// Both arrays must be provided by the user

// "Type" must have "=", ">", and "==" operators defined

template <class Type>
void quicksort(Type * inarray, int size, int * order = NULL)
{
  if (size <= 1) return;

  // if necessary, randomize random sequence
  static bool have_seed = false;
  if (!have_seed)
    {
      // struct timeval tv;
      // struct timezone tz;
      // gettimeofday(&tv, &tz);
      // SRANDOM(tv.tv_sec);

      time_t rawtime = time(NULL);
      unsigned int seed = (unsigned int) rawtime;
      SRANDOM(seed);

      have_seed = true;
    }
  int top = size - 1; // highest unsorted element of array
  int bottom = 0; // lowest unsorted element of array
  int equals = 0; // number of elements equal to key

  Type key = inarray[ RANDOM() % size ];
  Type temp;
  int tempint;

  // sort into three groups
  while (bottom <= top)
    {
      if (inarray[bottom] > key) // swap with top
	{
	  temp = inarray[top];
	  inarray[top] = inarray[bottom];
	  inarray[bottom] = temp;

	  if (order != NULL)
	    {
	      tempint = order[top];
	      order[top] = order[bottom];
	      order[bottom] = tempint;
	    }

	  --top;
	}
      else
	{
	  if (inarray[bottom] == key) ++equals;
	  // overwrite any elements equal to key
	  else 
	    {
	      temp = inarray[bottom - equals];
	      inarray[bottom - equals] = inarray[bottom];
	      inarray[bottom] = temp;

	      if (order != NULL)
		{
		  tempint = order[bottom - equals];
		  order[bottom - equals] = order[bottom];
		  order[bottom] = tempint;
		}
	    }
	  ++ bottom;
	}
    }

  // move bottom and top (now bottom = top + 1) to actual boundary
  --bottom; ++top;

  // call recursive quicksort on rest
  if (order != NULL)
    {
      quicksort(inarray, bottom - equals + 1, order);
      quicksort(inarray + bottom + 1, size - bottom - 1, order +
		bottom + 1);
    }
  else
    {
      quicksort(inarray, bottom - equals + 1, NULL);
      quicksort(inarray + bottom + 1, size - bottom - 1, NULL);
    }
}

// instantiations
template void quicksort(float *, int, int *);
template void quicksort(double *, int, int *);
