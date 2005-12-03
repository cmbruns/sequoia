#include "consens.hxx"

// $Id: consens.cxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/consens.cxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: consens.cxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

// Generate a consensus sequence from an alignment
Sequence consensus(const Alignment & a, 
		   const Comparison_matrix & cm, 
		   float cutoff)
{
  Sequence answer;
  double aa_score[256];
  if (a.n_seq() > 0)
    {
      answer = a[0];
      uint rnum;
      // for each residue in the sequence
      for (rnum=0; rnum < answer.length(); ++rnum)
	{
	  char rtype;
	  float are_gaps = 0;
	  float non_gaps = 0;
	  // for each amino acid type
	  for (rtype = 'A'; rtype <= 'Z'; ++rtype)
	    aa_score[rtype] = 0;
	  // for each sequence in alignment
	  uint snum;
	  for (snum = 0; snum < a.n_seq(); ++snum)
	    {
	      char ares = a[snum][rnum].letter();
	      double weight = a[snum].weight();
	      if (is_gap(ares)) are_gaps += weight;
	      else non_gaps +=weight;
	      for (rtype = 'A'; rtype <= 'Z'; ++rtype)
		aa_score[rtype] += weight * cm[rtype][ares];
	    }
	  // now make statistics on result
	  if (are_gaps > non_gaps) answer[rnum].letter() = '-';
	  else
	    {
	      char max = 'a';
	      double mean = 0;
	      int count = 0;
	      for (rtype = 'A'; rtype <= 'Z'; ++rtype)
		{
		  if (aa_score[rtype] > aa_score[max]) max = rtype;
		  mean += aa_score[rtype];
		  ++count;
		}
	      mean = mean/count;
	      double sigma = 0;
	      for (rtype = 'A'; rtype <= 'Z'; ++rtype)
		{
		  double df = aa_score[rtype] - mean;
		  sigma += df*df;
		}
	      sigma = sqrt(sigma / count);
	      if (sigma == 0) answer[rnum].letter() = max;
	      else if ((aa_score[max]/sigma) > cutoff) answer[rnum].letter() = max;
	      else answer[rnum].letter() = 'X';
	    }
	}
    }
  else answer = "";
  return answer;
}
