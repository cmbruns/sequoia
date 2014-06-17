#include "overlay.hxx"
#include "sort.h"
#include <vector>

// $Id: overlay.cxx,v 1.3 2002/02/19 18:40:16 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/overlay.cxx,v 1.3 2002/02/19 18:40:16 bruns Exp $
// $Log: overlay.cxx,v $
// Revision 1.3  2002/02/19 18:40:16  bruns
// Added throws of new exception ill_conditioned_rotation_error to overlay
// functions
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

// return a matrix specifying the optimal overlay of the first set of
// coordinates onto the second.  X and Y contain the coordinate sets
// to be compared.  w is an array of weights for the coordinates.
// Eventually perhaps this should default to all 1s. The
// answer matrix is the homogeneous transformation matrix, including a
// translation.  In 3D, the typical final row of this matrix will be
// {0, 0, 0, 1}.  This result is based upon the formulas of Wolfgang
// Kabsch (Acta Crystallographica A32: 922-923 (1976))

// 4-99 Attempt to incorporate case of Kabsch(1978) to prevent
// occaisional production of non proper rotations.

// *** NOTE! the cross products still make this routine 3-space specific ***
// Return zero matrix if not well-conditioned
RMat overlay(const Array<RVec> & X, const Array<RVec> & Y, const
	       Array<Real> & w) 
//   throw(ill_conditioned_rotation_error)
{
  if (X.dim() < 1) 
    {
      RMat answer(4,4,0);
      throw ill_conditioned_rotation_error();
      return answer; // assume 3 dimensions if no data
    }
  uint space = X[0].dim();
  RMat U(space + 1, space + 1, 0); // zero matrix
  uint i, j, k, n;

  // find center-of-mass of each coordinate set
  double sum = w[0];
  RVec t1 = X[0] * w[0];
  RVec t2 = Y[0] * w[0];
  for (n=1; n < X.dim(); ++n)
    {
      sum += w[n];
      t1 += X[n] * w[n];
      t2 += Y[n] * w[n];
    }
  t1 = t1 / sum;
  t2 = t2 / sum;

  // generate the symmetric positive definite matrix
  RMat R(space, space);
  for (i=0; i < space; ++i)
    for (j=0; j < space; ++j)
      {
	R[i][j] = 0;
	for (n=0; n < X.dim(); ++n)
	    R[i][j] += w[n] * (Y[n][i] - t2[i]) * (X[n][j] - t1[j]);
      }
  RMat spd = R.T() * R; // symmetric positive definite matrix

  // generate the eigenvectors and eigenvalues of the spd matrix
  RMat a1 = spd.sym_eigen(); // returns eigenvalues and eigenvectors

  // sort eigenvalues in increasing order
  // int order[space];
  std::vector<int> order(space);
  // Real mu1[space]; // sorted eigenvalues, increasing
  std::vector<Real> mu1(space);
  for (i = 0; i < space; ++i) 
    {
      order[i] = i;
      mu1[i] = a1[space][i];
    }
  quicksort(&mu1[0], space, &order[0]);

  // Real mu[space]; // sorted eigenvalues, decreasing
  std::vector<Real> mu(space);
  RMat a = a1;
  for (i = 0; i < space; ++i)
    {
      a[space - 1 - i] = a1[order[i]]; // sorted eigenvector
      mu[space - 1 - i] = a1[space][order[i]];
      a[space][space - i - 1] = mu[space - 1 - i]; // sorted eigenvalue
    }

  // If mu[space - 1] <= zero, there is no problem, since this
  // vector-value pair is never used.
  // BUT, if mu[space - 2] <= zero, the problem is ill-defined
  if (mu[space - 2] <= 0)
    {
      // Leave U as zero matrix
      throw ill_conditioned_rotation_error();
    }

  else
    {
      // Multiplier to correct for non proper rotations (Kabsch, 1978)
      // Doesn't actually do anything here!!!
      RVec sigma(space);
      for (i=0; i < space; ++i)
	sigma[i] = 1.0;
      
      // Cross product here makes this routine 3-space specific
      // generate vectors b, "effect of U on a"
      RMat b(space, space);
      
      a[2] = Vector3D(a[0]) % Vector3D(a[1]);
      b[0] = sigma[0] * (1/sqrt(mu[0])) * (R * a[0]);
      b[1] = sigma[1] * (1/sqrt(mu[1])) * (R * a[1]);
      b[2] = Vector3D(b[0]) % Vector3D(b[1]);
      
      // Generate the rotation matrix
      for (i=0; i < space; ++i)
	for (j=0; j < space; ++j)
	  {
	    U[i][j] = 0;
	    for (k=0; k < space; ++k)
	      U[i][j] += b[k][i] * a[k][j];
	    // make just rotation part for checking
	    R[i][j] = U[i][j];
	  }

      if (R.determinant() < 0)
	{
	  cout << "Negative determinant in R." << endl;
	  sigma[2] = -1.0;
	  RMat answer(space + 1, space + 1, 0);
	  return answer;
	  // remake rotation matrix ?
	}

      if ( (R.determinant() < 0.9)
	   || (R.determinant() > 1.1) )
	{
	  cerr << "Ill conditioned determinant for Kabsch rotation matrix!" 
	       << endl;
	  cerr << "Determinant = " << R.determinant() << endl;
	  
	  cerr << "eigenvalues = ";
	  for (i = 0; i < space; ++i)
	    cerr << mu[i] << ", ";
	  cerr << endl;
	  
	  for (i = 0; i < space; ++i)
	    {
	      cerr << "b[" << i << "] = (";
	      for (j = 0; j < space; ++j)
		cerr << b[i][j] << " ";
	      cerr << endl;
	    }
	  for (i = 0; i < space; ++i)
	    {
	      cerr << "a[" << i << "] = (";
	      for (j = 0; j < space; ++j)
		cerr << a[i][j] << " ";
	      cerr << endl;
	    }
	  
	  RMat answer(space + 1, space + 1, 0);
	  throw ill_conditioned_rotation_error();
	  return answer;
	  // exit(1);
	}
    }

  // add "fourth"(final) row and column
  for (i=0; i < space; ++i)
    {
      U[i][space] = 0;
      U[space][i] = 0;
    }
  U[space][space] = 1;

  // Restore translational components to make final transformation
  // matrix
  RMat T2 = trans_mat(t2);
  RMat T1 = trans_mat(-t1);
  return (T2 * U * T1);
}


// Matrix to transform residue2 onto residue1
RMat overlay(const SeqRes & r1, const SeqRes & r2)
//  throw(ill_conditioned_rotation_error)
{
  Array<RVec> v1, v2;
  Array<Real> w;

  RMat answer = eye(1) - eye(1); // zero matrix
  // r2 is failing - is this refAla?
  if ((r1.refatomp() == NULL) || (r2.refatomp() == NULL))
    {
      throw ill_conditioned_rotation_error(r1);
    }

  // Compare all the atoms
  const pdbatom * p1 = r1.refatom().first_in_res();
  while (p1 != NULL)
    {
      const pdbatom * p2 = r2.refatom().first_in_res();
      while (p2 != NULL)
	{
	  // Are they the same atom name?
	  if (!strcmp(p1->atom_name(), p2->atom_name()))
	    {
	      // The atoms are equivalent, in name anyway
	      w += 1.00;
	      v1 += p1->coord();
	      v2 += p2->coord();	      
	    }
	  p2 = p2->next_in_res();
	}
      p1 = p1->next_in_res();
    }

  // Not enough coordinates?
  if (v1.dim() < 2) 
    {
      throw ill_conditioned_rotation_error(r1);
    }
  answer = overlay(v2, v1, w);
  double determinant = answer.strip_trans().determinant();
  if (determinant < 0.9 || determinant > 1.1) {
    throw ill_conditioned_rotation_error(r1);
  }
  return answer;
}


RMat overlay(const pdbprotein & p1, 
		  const pdbprotein & p2, 
		  const Alignment & a)
{
  pdbprotein s1 = p1.ca();
  pdbprotein s2 = p2.ca();
  Array<RVec> v1, v2;
  Array<Real> w;
  uint i,j,r1,r2,j2;
  i = j = r1 = r2 = j2 = 0;
  // i is position in alignment
  // j is number of paired coordinates
  // r1 is position in structure
  while (r1 < s1.dim())
    {
      if ( (!is_gap(a[0][i].letter()))
	  && (!is_gap(a[1][i].letter())) )
	{
	  if ( (a[0][i].letter() == toupper(a[0][i].letter())) 
	      && (a[1][i].letter() == toupper(a[1][i].letter())) )
	    {
	      w += 1.00;
	      v1 += s1[r1].coord();
	      v2 += s2[r2].coord();
	      ++j2;
	    }
	  ++j;
	}
      if (!is_gap(a[0][i].letter())) ++r1;
      if (!is_gap(a[1][i].letter())) ++r2;
      ++i;
    }
  RMat answer = overlay(v2, v1, w);
  return answer;
}

// Reassign capital letters
Alignment struct_equiv(const pdbprotein & p1, 
		       const pdbprotein & p2,
		       const Alignment & a,
		       Real cutoff, uint runlength)
{
  Alignment answer = a.downcase();
  pdbprotein s1 = p1.ca();
  pdbprotein s2 = p2.ca();
  uint i,j,r1,r2,j2;
  i = j = r1 = r2 = j2 = 0;
  Real diff1;
  // i is position in alignment
  // j is number of paired coordinates
  // r1 is position in structure
  while (r1 < s1.dim())
    {
      if ( (!is_gap(answer[0][i].letter()))
	  && (!is_gap(answer[1][i].letter())) )
	{
	  diff1 = s1[r1].distance(s2[r2]);
	  if (diff1 < cutoff)
	    {
	      answer[0][i].letter() = toupper(answer[0][i].letter());
	      answer[1][i].letter() = toupper(answer[1][i].letter());
	      ++j2;
	    }
	  ++j;
	}
      if (!is_gap(answer[0][i].letter())) ++r1;
      if (!is_gap(answer[1][i].letter())) ++r2;
      ++i;
    }
  // remove runs of less than runlength
  uint run = 0;
  for (i=0; i < answer.length(); ++i)
    {
      if (isupper(answer[0][i].letter())) ++ run;
      else
	{
	  if ((run > 0) && (run < runlength))
	    {
	      for (j = 1; j <= run; ++j)
		{
		  answer[0][i-j].letter() = tolower(answer[0][i-j].letter());
		  answer[1][i-j].letter() = tolower(answer[1][i-j].letter());
		}
	    }
	  run = 0;
	}
    }
  // don't forget final run
  if ((run > 0) && (run < runlength))
    {
      for (j = 1; j <= run; ++j)
	{
	  answer[0][i-j].letter() = tolower(answer[0][i-j].letter());
	  answer[1][i-j].letter() = tolower(answer[1][i-j].letter());
	}
    }
  return answer;
}

// This version should do all pairwise comparisons
// Reassign capital letters
// When there are more than 3 structures, this method may be ambiguous!
Alignment struct_equiv(const Alignment & a,
		       Real dcutoff, Real acutoff, uint runlength)
{
  // first see if there are at least two coordinate sets
  uint i;
  uint ncoords = 0;
  for (i = 0; i < a.dim(); ++i) {
    cerr << i << "\n" << endl;
    if (a[i].have_coordinates()) ++ncoords;
  }
  if (ncoords < 2)
    {
      return a;
    }

  Alignment answer = a.downcase();
  Alignment subalign;
  Sequence sq1, sq2;

  // Compare all n*(n-1)/2 pairs of sequences
  uint ind1, ind2;
  for (ind1 = 0; ind1 < (a.dim() - 1); ++ind1)
    for (ind2 = ind1 + 1; ind2 < a.dim(); ++ind2)
      {
	if ( (a[ind1].have_coordinates()) && (a[ind2].have_coordinates()))
	  {
	    pdbprotein s1 = a[ind1].coordinates();
	    pdbprotein s2 = a[ind2].coordinates();
	    // make temporary sequences for pristine comparison
	    sq1 = a[ind1].downcase();
	    sq2 = a[ind2].downcase();
	    subalign = sq1;
	    subalign += sq2;

	    subalign = struct_equiv(s1, s2, subalign, dcutoff, runlength);
	    sq1 = subalign[0];
	    sq2 = subalign[1];
	    // update capitalization for actual answer
	    uint i;
	    for (i=0; i < answer.length(); ++i)
	      {
		if (isupper(sq1[i].letter()))
		  answer[ind1][i].letter() = sq1[i].letter();
		if (isupper(sq2[i].letter()))
		  answer[ind2][i].letter() = sq2[i].letter();
	      }
	  }
      }
  return answer;
}


RMat overlay(const pdboverlay & o1, 
	     const pdboverlay & o2, 
	     const Alignment & a)
{
  Array<RVec> v1, v2;
  Array<Real> w;

  pdbprotein s1, s2;
  Sequence sq1, sq2;
  uint ind1, ind2;
  // only compare BETWEEN overlays
  for (ind1 = 0; ind1 < o1.dim(); ++ind1)
    for (ind2 = o1.dim(); ind2 < o1.dim() + o2.dim(); ++ind2)
      {
	// choose two coordinate sets to compare
	s1 = o1[ind1].ca();
	s2 = o2[ind2 - o1.dim()].ca();
	// choose two sequences to compare
	sq1 = a[ind1];
	sq2 = a[ind2];

	uint i,j,r1,r2,j2;
	i = j = r1 = r2 = j2 = 0;
	// i is position in alignment
	// j is number of paired coordinates
	// r1 is position in structure
	while (r1 < s1.dim())
	  {
	    if ( (!is_gap(sq1[i].letter()))
		 && (!is_gap(sq2[i].letter())) )
	      {
		if ( (sq1[i].letter() == toupper(sq1[i].letter())) 
		     && (sq2[i].letter() == toupper(sq2[i].letter())) )
		  {
		    w += 1.00;
		    v1 += s1[r1].coord();
		    v2 += s2[r2].coord();
		    ++j2;
		  }
		++j;
	      }
	    if (!is_gap(sq1[i].letter())) ++r1;
	    if (!is_gap(sq2[i].letter())) ++r2;
	    ++i;
	  }
      }
  RMat answer = overlay(v2, v1, w);
  return answer;
}


// equivalences BETWEEN overlays
Alignment struct_equiv(const pdboverlay & o1, 
		       const pdboverlay & o2,
		       const Alignment & a,
		       Real cutoff, uint runlength)
{
  Alignment answer = a.downcase();
  Alignment subalign;
  Sequence sq1, sq2;
  pdbprotein s1, s2;

  uint ind1, ind2;
  // only compare BETWEEN overlays
  for (ind1 = 0; ind1 < o1.dim(); ++ind1)
    for (ind2 = o1.dim(); ind2 < o1.dim() + o2.dim(); ++ind2)
      {
	// select coordinates for comparison
	s1 = o1[ind1].ca();
	s2 = o2[ind2 - o1.dim()].ca();
	// make temporary sequences for pristine comparison
	sq1 = a[ind1].downcase();
	sq2 = a[ind2].downcase();
	subalign = sq1;
	subalign += sq2;

	subalign = struct_equiv(s1, s2, subalign, cutoff, runlength);
	sq1 = subalign[0];
	sq2 = subalign[1];
	// update capitalization for actual answer
	uint i;
	for (i=0; i < answer.length(); ++i)
	  {
	    if (isupper(sq1[i].letter()))
	      answer[ind1][i].letter() = sq1[i].letter();
	    if (isupper(sq2[i].letter()))
	      answer[ind2][i].letter() = sq2[i].letter();
	  }
      }
  return answer;
}
