#include "timelice.hxx"
#include <ctime>

// $Id: timelice.cxx,v 1.3 2001/11/29 00:04:43 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/timelice.cxx,v 1.3 2001/11/29 00:04:43 bruns Exp $
// $Log: timelice.cxx,v $
// Revision 1.3  2001/11/29 00:04:43  bruns
// Removed ^M characters
// Added cvs tags
//

// Only allow operation of program on a range of dates
int timed_license()
{
  struct timeval tv;
  if (gettimeofday(&tv, NULL)) return 666;

  time_t ltarg = tv.tv_sec;
  struct tm * ltp;
  ltp = localtime(&ltarg);

  // cout << "Year = " << ltp->tm_year + 1900 << endl;
  // cout << "Month = " << ltp->tm_mon + 1 << endl;

  if ((ltp->tm_year + 1900) < START_YEAR) return -1;
  else if (((ltp->tm_year + 1900) == START_YEAR) && ((ltp->tm_mon + 1) < START_MONTH)) return -1;

  if ((ltp->tm_year + 1900) > END_YEAR) return 1;
  else if (((ltp->tm_year + 1900) == END_YEAR) && ((ltp->tm_mon + 1) > END_MONTH)) return 1;

  return 0;
}
