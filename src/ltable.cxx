/* lookup tables for slow mathematical functions */

// $Id: ltable.cxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/ltable.cxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: ltable.cxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

#include "ltable.hxx"
#include "rmat.hxx"
#include "MSVC_erfc.h"

/* probability of a z-score greater than x */
double zprob(double x)
{
  if (x < 0) x = -x;
  return erfc(x/SQR2);
}

#define ZPROB_BINS 100
/* probability of a z-score greater than x */
/* this is a bit faster than zprob on sun5 */
double table_zprob(double x)
{
  /* Use two tables, one from zero to ten, the other from ten to 110 */
  static int have_lookup_table = 0;
  double minval = 0;
  double maxval = 10.0;
  double minval2 = 10.0;
  double maxval2 = 110.0;
  static double table[ZPROB_BINS];
  static double table2[ZPROB_BINS];
  double deltaf, index;
  int index2;

  if (x < 0) x = -x;
  if (!have_lookup_table)
    {
      int i;
      for (i=0; i < ZPROB_BINS; ++i)
	table[i] = zprob( minval + (i * (maxval - minval))/(ZPROB_BINS - 1));
      for (i=0; i < ZPROB_BINS; ++i)
	table2[i] = zprob( minval2 + (i * (maxval2 - minval2))/(ZPROB_BINS - 1));
      have_lookup_table = 1;
    }
  if (x == 0) return 1.0;
  else if ((x < minval) || (x > maxval2)) return zprob(x);
  else if (x == maxval) return table[ZPROB_BINS - 1];
  else if (x == maxval2) return table2[ZPROB_BINS - 1];
  /* WARNING : very low probabilities will not be accurate! */
  /* but at least they are monotonic */
  else if (x > maxval2) return table2[ZPROB_BINS - 1] * maxval2 / x;
  /* else if (x > maxval2) return zprob(x); */
  else
    {
      if (x < minval2)
	{
	  index = (ZPROB_BINS - 1) * (x - minval) / (maxval - minval);
	  index2  = (int) index;
	  deltaf = table[index2 + 1] - table[index2];
	  return table[index2] + (index - index2) * deltaf;
	}
      else
	{
	  index = (ZPROB_BINS - 1) * (x - minval2) / (maxval2 - minval2);
	  index2  = (int) index;      
	  deltaf = table2[index2 + 1] - table2[index2];
	  return table2[index2] + (index - index2) * deltaf;
	}
    }
}

/* probability that two randomly chosen orientations are
   within kappa angle phi of one another */
double angprob(double x)
{
  if (x < 0) x = -x;
  while (x > PI) x -= PI;
  if (x == 0) return 0;
  else
    {
      return (x - sin(x)) / PI;
    }
}

#define ANGPROB_BINS 100
/* slightly faster on sun5 */
double table_angprob(double x)
{
  static int have_lookup_table = 0;
  double minval = 0.0;
  double maxval = PI;
  static double table[ANGPROB_BINS];
  double deltaf, index;
  int index2;

  if (x < 0) x = -x;
  while (x > PI) x -= PI;
  if (x == 0) return 0;

  if (!have_lookup_table)
    {
      int i;
      for (i=0; i < ANGPROB_BINS; ++i)
	table[i] = angprob( minval + (i * (maxval - minval))/(ANGPROB_BINS - 1));
      have_lookup_table = 1;
    }
  if ((x < minval) || (x > maxval)) return angprob(x);
  else if (x == maxval) return table[ANGPROB_BINS - 1];
  else
    {
      index = (ANGPROB_BINS - 1) * (x - minval) / (maxval - minval);
      index2 = (int) index; 
      deltaf = table[index2 + 1] - table[index2];
      return table[index2] + (index - index2) * deltaf;
    }
}

// This one uses the COSINE of the angle as the argument 
/* probability that two randomly chosen orientations are
   within kappa angle phi of one another */
double angprobC(double x)
{
  x = acos(x);
  if (x < 0) x = -x;
  while (x > PI) x -= PI;
  if (x == 0) return 0;
  else
    {
      return (x - sin(x)) / PI;
    }
}

#define ANGPROBC_BINS 100
double table_angprobC(double x)
{
  static int have_lookup_table = 0;
  double minval = -1.0;
  double maxval = 1.0;
  static double table[ANGPROBC_BINS];
  double deltaf, index;
  int index2;

  if (!have_lookup_table)
    {
      int i;
      for (i=0; i < ANGPROBC_BINS; ++i)
	table[i] = angprobC( minval + (i * (maxval - minval))/(ANGPROBC_BINS - 1));
      have_lookup_table = 1;
    }
  if ((x < minval) || (x > maxval)) return angprobC(x);
  else if (x == maxval) return table[ANGPROBC_BINS - 1];
  else
    {
      index = (ANGPROBC_BINS - 1) * (x - minval) / (maxval - minval);
      index2 = (int) index; 
      deltaf = table[index2 + 1] - table[index2];
      return table[index2] + (index - index2) * deltaf;
    }
}

#define LOG_BINS 100
/* this is about the same as log as sun5 */
double table_log(double x)
{
  static int have_lookup_table = 0;

  static double e16 = 8886110.521;
  static double e8 = 2980.957987;
  static double e4 = 54.59815003;
  static double e2 = 7.389056099;
  static double e1 = 2.718281828459045;
  double minval = 1.0;
  double maxval = 2.718281828459045;

  static double table[LOG_BINS];
  double deltaf, index;
  int index2;
  double answer_mult = 1.0;
  double answer_add = 0.0;
  double answer;

  if (x < 0) x = -x;
  if (x == 0) return log(x);
  if (x > 1000000000) return log(x);
  if (x < 1) {x = 1/x; answer_mult = -1.0;}
  while (x > e16) {x /= e16; answer_add += 16;}
  while (x > e8) {x /= e8; answer_add += 8;}
  while (x > e4) {x /= e4; answer_add += 4;}
  while (x > e2) {x /= e2; answer_add += 2;}
  while (x > e1) {x /= e1; answer_add += 1;}

  if (!have_lookup_table)
    {
      int i;
      for (i=0; i < LOG_BINS; ++i)
	table[i] = log( minval + (i * (maxval - minval))/(LOG_BINS - 1));
      have_lookup_table = 1;
    }
  if ((x < minval) || (x > maxval)) answer = log(x);
  else if (x == maxval) answer = table[LOG_BINS - 1];
  else
    {
      index = (LOG_BINS - 1) * (x - minval) / (maxval - minval);
      index2 = (int) index; 
      deltaf = table[index2 + 1] - table[index2];
      answer = table[index2] + (index - index2) * deltaf;
    }
  return (answer + answer_add) * answer_mult;
}

/* probability that two randomly chosen points in a unit-radius sphere are
   within distance x of one another */
double radprob(double x)
{
  if (x < 0) x = -x;
  if (x == 0) return 0;
  else if (x >= 2.0) return 1.0;
  else
    {
      double x3 = x * x * x;
      return x3 - 9.0/16.0 * x3 * x + 1.0/32.0 * x3 * x3;
    }
}

#define RADPROB_BINS 100
/* this is actually SLOWER than radprob on sun5 */
double table_radprob(double x)
{
  static int have_lookup_table = 0;
  double minval = 0.0;
  double maxval = 2.0;
  static double table[RADPROB_BINS];
  double deltaf, index;
  int index2;

  if (x < 0) x = -x;
  if (x == 0) return 0;
  else if (x >= 2.0) return 1.0;

  if (!have_lookup_table)
    {
      int i;
      for (i=0; i < RADPROB_BINS; ++i)
	table[i] = radprob( minval + (i * (maxval - minval))/(RADPROB_BINS - 1));
      have_lookup_table = 1;
    }
  if ((x < minval) || (x > maxval)) return radprob(x);
  else if (x == maxval) return table[RADPROB_BINS - 1];
  else
    {
      index = (RADPROB_BINS - 1) * (x - minval) / (maxval - minval);
      index2 = (int) index; 
      deltaf = table[index2 + 1] - table[index2];
      return table[index2] + (index - index2) * deltaf;
    }
}

#define SQRT_BINS 100
double table_sqrt(double x)
{
  static int have_lookup_table = 0;
  double minval = 1.0/SQRT_BINS;
  double maxval = 1.0;
  static double table[SQRT_BINS];
  double deltaf, index;
  int index2;

  if (!have_lookup_table)
    {
      int i;
      for (i=0; i < SQRT_BINS; ++i)
	table[i] = sqrt( minval + (i * (maxval - minval))/(SQRT_BINS - 1));
      have_lookup_table = 1;
    }
  if (x == 0) return 0;
  if ((x < minval) || (x > maxval)) return sqrt(x);
  else if (x == maxval) return table[SQRT_BINS - 1];
  else
    {
      index = (SQRT_BINS - 1) * (x - minval) / (maxval - minval);
      index2 = (int) index; 
      deltaf = table[index2 + 1] - table[index2];
      return table[index2] + (index - index2) * deltaf;
    }
}

#define ASIN_BINS 100
/* this is faster than asin on sun5 */
double table_asin(double x)
{
  static int have_lookup_table = 0;
  double minval = -1.0;
  double maxval = 1.0;
  static double table[ASIN_BINS];
  double deltaf, index;
  int index2;

  if (!have_lookup_table)
    {
      int i;
      for (i=0; i < ASIN_BINS; ++i)
	table[i] = asin( minval + (i * (maxval - minval))/(ASIN_BINS - 1));
      have_lookup_table = 1;
    }
  if ((x < minval) || (x > maxval)) return asin(x);
  else if (x == maxval) return table[ASIN_BINS - 1];
  else
    {
      index = (ASIN_BINS - 1) * (x - minval) / (maxval - minval);
      index2 = (int) index; 
      deltaf = table[index2 + 1] - table[index2];
      return table[index2] + (index - index2) * deltaf;
    }
}


#define ACOS_BINS 100
/* this is faster than acos on sun5 */
double table_acos(double x)
{
  static int have_lookup_table = 0;
  double minval = -1.0;
  double maxval = 1.0;
  static double table[ACOS_BINS];
  double deltaf, index;
  int index2;

  if (!have_lookup_table)
    {
      int i;
      for (i=0; i < ACOS_BINS; ++i)
        table[i] = acos( minval + (i * (maxval - minval))/(ACOS_BINS - 1));
      have_lookup_table = 1;
    }
  if ((x < minval) || (x > maxval)) return acos(x);
  else if (x == maxval) return table[ACOS_BINS - 1];
  else
    {
      index = (ACOS_BINS - 1) * (x - minval) / (maxval - minval);
      index2 = (int) index;
      deltaf = table[index2 + 1] - table[index2];
      return table[index2] + (index - index2) * deltaf;
    }
}


#define ATAN2_BINS 100
/* make two tables, depending on which is greater, x or y */
double table_atan2(double y, double x)
{
  static int have_lookup_table = 0;
  double minval = 0.0;
  double maxval = 1.0;
  static double table1[ATAN2_BINS];
  static double table2[ATAN2_BINS];
  double deltaf, index;
  int index2;
  double answer;
  double x2 = fabs(x);
  double y2 = fabs(y);
  double ratio;

  if (!have_lookup_table)
    {
      int i;
      for (i=0; i < ATAN2_BINS; ++i)
	{
	  /* y >= x */
	  table1[i] = atan2( 1.0, minval + (i * (maxval - minval))/(ATAN2_BINS - 1));
	  table2[i] = atan2( minval + (i * (maxval - minval))/(ATAN2_BINS - 1), 1.0);
	}
      have_lookup_table = 1;
    }

  if (y2 >= x2) /* use table 1 */
    {
      ratio = x2/y2;
      if ((ratio < minval) || (ratio > maxval)) answer =  atan2(y2,x2);
      else if (ratio == maxval) answer = table1[ATAN2_BINS - 1];
      else
	{
	  index = (ATAN2_BINS - 1) * (ratio - minval) / (maxval - minval);
	  index2 = (int) index; 
	  deltaf = table1[index2 + 1] - table1[index2];
	  answer = table1[index2] + (index - index2) * deltaf;
	}
    }
  else /* use table 2 */
    {
      ratio = y2/x2;
      if ((ratio < minval) || (ratio > maxval)) answer =  atan2(y2,x2);
      else if (ratio == maxval) answer = table2[ATAN2_BINS - 1];
      else
	{
	  index = (ATAN2_BINS - 1) * (ratio - minval) / (maxval - minval);
	  index2 = (int) index; 
	  deltaf = table2[index2 + 1] - table2[index2];
	  answer = table2[index2] + (index - index2) * deltaf;
	}
    }
  if (x < 0) answer = PI - answer;
  if (y < 0) answer = - answer;
  return answer;
}
