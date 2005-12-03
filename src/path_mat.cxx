#include "path_mat.hxx"
#include "geoscore.hxx"
#include "ltable.hxx"

// $Id: path_mat.cxx,v 1.4 2002/06/21 00:44:06 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/path_mat.cxx,v 1.4 2002/06/21 00:44:06 bruns Exp $
// $Log: path_mat.cxx,v $
// Revision 1.4  2002/06/21 00:44:06  bruns
// Added condition to return 0 for mean and standard deviation, if the matrix
// has no cells
//
// Revision 1.3  2002/02/19 17:22:01  bruns
// Added check to avoid computing orientation when it is not well determined
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

#define debug 0

path_element::path_element()
{
  comp_score = 0;
  path_score = 0;
  best_step = 0;
}


residue_attribute::residue_attribute()
{
  //cout << "residue_attribute constructor\n";
  gpen = DEFAULT_GAP_PENALTY;
  score_offset = 0;
}


Score_matrix::Score_matrix() : Array2d<path_element>()
{
  // cout << "score_matrix constructor\n";
  row_attrib = NULL;
  column_attrib = NULL;
  epen() = DEFAULT_EXTENSION_PENALTY; // gap extension penalty
  gpen() = DEFAULT_GAP_PENALTY;
  sigma = 10; // rms score for gap penalty calculation
  random_seed() = 1;
  dcutoff() = DEFAULT_CUTOFF;
  acutoff() = DEFAULT_ACUTOFF;
  have_stable = false; // is comp_score up to date?
  have_table = false; // is comp_score up to date?
  have_path = false; // is path_score up to date?
  highroad = false;
  midroad = false;
  lowroad = false;
  randroad = true;
  useangle() = true;
  suboptimal() = DEFAULT_SUBOPTIMAL_FACTOR;
  runlength() = DEFAULT_RUNLENGTH;
  a1 = NULL;
  a2 = NULL;
}

Score_matrix::Score_matrix(uint r, uint c) : Array2d<path_element>(r,c) 
{
  // cout << "score matrix sized constructor\n";
  row_attrib = new residue_attribute[r];
  column_attrib = new residue_attribute[c];
  epen() = DEFAULT_EXTENSION_PENALTY; // gap extension penalty
  gpen() = DEFAULT_GAP_PENALTY;
  sigma = 10; // rms score for gap penalty calculation
  random_seed() = 1;
  dcutoff() = DEFAULT_CUTOFF;
  acutoff() = DEFAULT_ACUTOFF;
  have_table = false; // is comp_score up to date?
  have_stable = false; // is comp_score up to date?
  have_path = false; // is path_score up to date?
  highroad = false;
  lowroad = false;
  midroad = false;
  randroad = true;
  useangle() = true;
  suboptimal() = DEFAULT_SUBOPTIMAL_FACTOR;
  runlength() = DEFAULT_RUNLENGTH;
  a1 = NULL;
  a2 = NULL;
  uint i,j;
  for (i=0; i < (m() + 1); ++i)
    row_attrib[i].gpen = gpen();
  for (j=0; j < (n() + 1); ++j)
    column_attrib[j].gpen = gpen();
}

Score_matrix::Score_matrix
(const Alignment & a1, const Alignment & a2, const Comparison_matrix & cm)
{
  // cout << "score matrix fancy constructor\n";
  epen() = DEFAULT_EXTENSION_PENALTY; // gap extension penalty
  gpen() = DEFAULT_GAP_PENALTY;
  sigma = 10; // rms score for gap penalty calculation
  random_seed() = 1;
  dcutoff() = DEFAULT_CUTOFF;
  acutoff() = DEFAULT_ACUTOFF;
  have_path = false; // is path_score up to date?
  highroad = false;
  lowroad = false;
  randroad = true;
  useangle() = true;
  suboptimal() = DEFAULT_SUBOPTIMAL_FACTOR;
  runlength() = DEFAULT_RUNLENGTH;
  tabulate(a1, a2, cm);
}

Score_matrix::Score_matrix(const Score_matrix & sm)
{
  // cout << "score matrix copy constructor\n";
  *this = sm;
}

Score_matrix::~Score_matrix()
{
  // cout << "score matrix destructor\n";
  delete [] row_attrib;
  row_attrib = NULL;
  delete [] column_attrib;
  column_attrib = NULL;
}

void Score_matrix::update_gap_penalties() {
  uint i,j;
  for (i=0; i < (m() + 1); ++i)
    row_attrib[i].gpen = gpen();
  for (j=0; j < (n() + 1); ++j)
    column_attrib[j].gpen = gpen();
}

void Score_matrix::init_score_mat(uint r, uint c)
{
  // cout << "score matrix initializer\n";
  init_array(r, c);
  delete [] row_attrib;
  row_attrib = new residue_attribute[m() + 1];
  delete [] column_attrib;
  column_attrib = new residue_attribute[n() + 1];
  update_gap_penalties();
}

Score_matrix & Score_matrix::operator=(const Score_matrix & sm)
{
  // cout << "score matrix assignment operator\n";
  if (this == &sm) return *this;
  uint i, j;
  Score_matrix & s = *this;

  init_score_mat(sm.m(), sm.n());
  for (i=0; i < m(); ++i)
    {
      row_attrib[i] = sm.row_attrib[i];
      for (j=0; j < n(); ++j)
	s[i][j] = sm[i][j];
    }
  for (j=0; j < n(); ++j)
    column_attrib[j] = sm.column_attrib[j];

  epen() = sm.epen();
  gpen() = sm.gpen();
  sigma = sm.sigma;
  random_seed() = sm.random_seed();
  dcutoff() = sm.dcutoff();
  acutoff() = sm.acutoff();
  have_table = sm.have_table;
  have_stable = sm.have_stable;
  have_path = sm.have_path;
  highroad = sm.highroad;
  midroad = sm.midroad;
  lowroad = sm.lowroad;
  randroad = sm.randroad;
  useangle() = sm.useangle();
  suboptimal() = sm.suboptimal();
  runlength() = sm.runlength();
  a1 = sm.a1;
  a2 = sm.a2;

  return *this;
}

const MATRIX_TYPE & Score_matrix::operator()(uint i, uint j) const
{
  const Score_matrix & s = *this;
  return s[i][j].comp_score;
}

MATRIX_TYPE & Score_matrix::operator()(uint i, uint j)
{
  Score_matrix & s = *this;
  return s[i][j].comp_score;
}


// void Score_matrix::gpen(MATRIX_TYPE g)
// {
//   // cout << "gpen function\n";
//   uint i;
//   gpen() = g;
//   for (i=0; i < m(); ++i) row_attrib[i].gpen = gpen();
//   for (i=0; i < n(); ++i) column_attrib[i].gpen = gpen();
// }


MATRIX_TYPE Score_matrix::mean() const
{
  uint i,j;
  MATRIX_TYPE sum, answer;
  const Score_matrix & M = *this;
  
  sum = 0;
  for (i=0;i < m();++i) 
    for (j=0;j < n();++j)
      sum += M(i, j);
  if ((m() * n()) == 0) answer = 0;
  else answer = sum/ ( (MATRIX_TYPE) (m()*n()) );
  return answer;
}

MATRIX_TYPE Score_matrix::std_dev() const
{
  uint i,j;
  MATRIX_TYPE diff, answer;
  const Score_matrix & s = *this;
  MATRIX_TYPE sum = 0;
  MATRIX_TYPE mean = this->mean();

  for (i=0;i < m();++i) for (j=0;j < n();++j)
  {
    diff = mean - s(i, j);
    sum += (diff * diff);
  }

  if ((m() * n()) == 0) answer = 0;
  else answer = (MATRIX_TYPE) sqrt(sum / ( (MATRIX_TYPE) (m()*n()) ));
  return answer;
}

void Score_matrix::tabulate
(const Alignment & al1, const Alignment & al2, const Comparison_matrix & cm)
{
  uint r, c;
  if (al1.n_seq()) r = al1[0].length();
  else r = 0;
  if (al2.n_seq()) c = al2[0].length();
  else c = 0;
  a1 = & al1;
  a2 = & al2;

  Score_matrix & answer = *this;
  init_score_mat(r,c);

  uint i, j, s1, s2, t1, t2;
  double weight;
  Sequence sq1, sq2;
 
  /* loop through all pairings of sequences between the two */
  for (s1=0; s1 < al1.n_seq(); ++s1)
    {
      sq1 = al1[s1];
      for (s2=0; s2 < al2.n_seq(); ++s2)
	{
	  sq2 = al2[s2];
	  weight = al1[s1].weight() * al2[s2].weight();
	  /* loop through the residue positions in the matrix */
	  for (i=0; i < r; ++i) 
	    for (j=0; j < c; ++j)
	      {
		t1 = sq1[i].letter();
		t2 = sq2[j].letter();
		answer(i, j) += weight * cm[t1][t2];
	      }
	}
    }
 // offset table values so that the mean is zero
  MATRIX_TYPE offset = mean();  
  for (i=0; i < r; ++i)
    for (j=0; j < c; ++j)
      {
	answer(i, j) -= offset;
      }
  have_table = true;
  have_stable = false;
  have_path = false;
}

// FUNCTIONS FOR THE COMPARISON OF ATOMIC COORDINATES

void Score_matrix::stabulate
  (const pdbprotein & p1, 
   const pdbprotein & p2)
{
  uint r = p1.sequence().dim();
  uint c = p2.sequence().dim();
  Score_matrix & answer = *this;
  init_score_mat(r,c);
  double s;

  uint i, j;
  /* loop through the residue positions in the matrix */
  for (i=0; i < r; ++i)
    {
      for (j=0; j < c; ++j)
	{
	  if (p1.sequence()[i].have_ref() &&
	      p2.sequence()[j].have_ref())
	    {
	      const SeqRes & r1 = p1.sequence()[i];
	      const SeqRes & r2 = p2.sequence()[j];
	      s = res_3dscore(r1, r2);

	      answer(i, j) = s;
	    }
	  else answer(i,j) = 0; // gap or something
	}
    }
  have_stable = true;
  have_table = false;
  have_path = false;
}

void Score_matrix::stabulate
  (const pdboverlay & o1, 
   const pdboverlay & o2)
{
#define debug 0
  // This is a bit cheesy
  // It assumes alpha carbons only, and structure 1 is full length of all
  uint r = o1[0].sequence().dim();
  uint c = o2[0].sequence().dim();
  if (debug) 
    {
      uint i;
      cerr << "rows = " << r;
      for (i = 0; i < o1.dim(); ++i)
	cerr << " (" << i+1 << ": " << o1[i].sequence().dim() << ")";
      cerr << endl;
      cerr << "columns = " << c;
      for (i = 0; i < o2.dim(); ++i)
	cerr << " (" << i+1 << ": " << o2[i].sequence().dim() << ")";
      cerr << endl;
    }
  Score_matrix & answer = *this;
  init_score_mat(r,c);

  double s;

  uint i, j;
  uint ind1, ind2;
  /* loop through the residue positions in the matrix */
  for (i=0; i < r; ++i)
    {
      for (j=0; j < c; ++j)
	{
	  // initialize matrix position
	  s = 0;
	  // loop through structures
	  for (ind1 = 0; ind1 < o1.dim(); ++ind1)
	    for (ind2 = 0; ind2 < o2.dim(); ++ind2)
	      {
		const SeqRes & r1 = o1[ind1].sequence()[i];
		const SeqRes & r2 = o2[ind2].sequence()[j];
		s += res_3dscore(r1, r2);
	      }
	  answer(i, j) = s;
	  if (!(s >= 0))
	    {
	      cerr << " s < 0, s = " << s;
	      cerr << " i = " << i;
	      cerr << " j = " << j;
	      cerr << endl;
	    }
	  if (!(s <= 500))
	    {
	      cerr << " s < 0, s = " << s;
	      cerr << " i = " << i;
	      cerr << " j = " << j;
	      cerr << endl;
	    }
	  if (debug)
	    {
	      cerr << " ";
	      cerr.ios::precision(2);
	      cerr << s;
	    }
	}
    }
  if (debug) cerr << endl;
  have_stable = true;
  have_table = false;
  have_path = false;
}


// comparison score between two residues in PDB files
double Score_matrix::res_3dscore(const SeqRes & r1, const SeqRes & r2)
{
  double answer = 0;
  double Pnotequiv;
  double Pgap = 0.10; // expected probability of a gap, should be strictly
                      // less than or equal to Pnotequiv 
  static double logPgap = log(Pgap); // for scaling matrix values
  Pnotequiv = 0.10; // made up number for now: probability
                    // of non-equivalence, but being in the
                    // final alignment anyway 

  // Can we even align these things?
  if (!(is_gap(r1.letter())) && // not a gap
      (!is_gap(r2.letter())) && 
      (r1.refatomp() != NULL) && // and has an atom
      (r2.refatomp() != NULL))
    {
      // Use the DISTANCE as a criterion for equivalence
      double d = r1.refatom().distance(r2.refatom());
      double P1 = distance_score2(d, dcutoff(), 100.0);
      //double P1 = dcutoff()*dcutoff()/((d+dcutoff())*(d+dcutoff())); 

      // Use the ORIENTATION as a criterion for equivalence, too
      // (convert from radians to degrees)
      // Looks like this is really slow
      double P2;
      if (useangle()) 
	{
          // Make sure the orientation is well-determined before computing
          if ((r1.orientation().dim() == 3) && (r2.orientation().dim() == 3)) {
	    double o = ACOS(kappa_cos_angle(r1.orientation(), r2.orientation()));
  	    P2 = orientation_score2(o * RAD2DEG , acutoff());
	    if (!(P2 >= 0)) // NaN
	      {
	        cerr << "Internal error: orientation score NaN" << endl;
	        cerr << "P2 = " << P2 << endl;	
	        cerr << "o = " << o << endl;
	        cerr << "r1.orientation() = " << r1.orientation() << endl;
	        cerr << "r2.orientation() = " << r2.orientation() << endl;
	        cerr << "kappa_cos_angle(r1.orientation(), r2.orientation()) = " 
	  	     <<  kappa_cos_angle(r1.orientation(), r2.orientation()) 
		       << endl;
	        cerr << "sigma = " << sigma << endl;
	      }
            } else {
              // Orientation is not well determined?
              P2 = 0.2; // Just a guess, it should not be better than good orientations...
            }
        }
      else 
	P2 = 1.0;

      // New method, after epiphany Feb 1998
      // probability of equivalence, ranges from zero to one
      double Pequiv = P1 * P2;

      // true log-odds matrix
      answer = log(Pequiv + Pnotequiv) - logPgap;
      if (answer < 0) // should not happen
	  answer = 0;
      if (!(answer >= 0)) // NaN
	{
	  cerr << "Internal error: res3dscore NaN" << endl;
	  cerr << "answer = " << answer << endl;	
	  cerr << "logPgap = " << logPgap << endl;
	  cerr << "Pequiv = " << Pequiv << endl;
	  cerr << "P1 = " << P1 << endl;
	  cerr << "P2 = " << P2 << endl;
	  cerr << "Pnotequiv = " << Pnotequiv << endl;
	  cerr << "log(Pnotequiv + Pequiv) = " << log(Pequiv + Pnotequiv) << endl;
	}
    }
  else answer = 0; // Gap or something
  return answer;
}

// try to get it right this time
MATRIX_TYPE Score_matrix::randomize()
{
  Score_matrix & sm = *this;
  uint * order2 = new uint [n()];
  uint i, j, temp2, j0, vj, jm1;
  MATRIX_TYPE max;
  int k1, k2, best_step;

  // prevent core dumps for empty matrix
  if (!(m() * n())) return 0;
  /* in case gap penalties have changed... */
  update_gap_penalties();
  // choose random order for second sequence
  for (i=0;i < n();++i) order2[i]=i;
  for (i=0;i < n();++i)
    {
      j = rand() % n();
      temp2 = order2[j];
      order2[j] = order2[i];
      order2[i] = temp2;
    }
  j0 = order2[0];

  /* first initialize the first row and column to the
     values in the comparison table.  This assumes that
     'end gaps' are not penalized, but why would you ever
     want to penalize them? */
  for (i=0; i < m(); ++i) 
    {
      row_attrib[i].best_score = 0;
      sm[i][j0].path_score = sm(i, j0);
    }
  for (vj=0; vj < n(); ++vj) 
    {
      j = order2[vj];
      column_attrib[j].best_score = 0;
      sm[0][j].path_score = sm(0, j);
    }
  MATRIX_TYPE A, B, C, A1, B1;
  for (i=1; i < m(); ++i)
    { 
      for (vj=1; vj < n(); ++vj)
	{
	  // j is the physical memory column in the (original) DP matrix
	  // vj is the virtual column in the "reordered" DP matrix
	  j = order2[vj];
	  jm1 = order2[vj - 1];
	  C = sm[i-1][jm1].path_score; // path score with no gaps
	  k1 = i - 1 - column_attrib[jm1].best_score; // gap size in seq1
	  A1 = sm[i-k1-1][jm1].path_score - sigma * k1 * epen();
	  A = A1 - sigma * column_attrib[jm1].gpen; // path score with gap in seq1
	  k2 = vj - 1 - row_attrib[i-1].best_score; // gap size in seq2
	  B1 = sm[i-1][order2[vj-k2-1]].path_score - sigma * k2 * epen();
	  B = B1 - sigma * row_attrib[i-1].gpen; // path score with gap in seq2
	  max = MAX(MAX(A,B),C);

	  if (C >= max) best_step = 0;
	  else if (A >= max) best_step = k1;
	  else best_step = -k2;

	  sm[i][j].best_step = best_step;
	  if (!best_step) sm[i][j].path_score = sm(i,j) + C;
	  else if (best_step < 0) sm[i][j].path_score = sm(i,j) + B;
	  else sm[i][j].path_score = sm(i,j) + A;

	  // update column and row tables
	  if (A1 < C)
	    {
	      column_attrib[jm1].best_score = i - 1;
	    }
	  if (B1 < C)
	    {
	      row_attrib[i-1].best_score = vj - 1;
	    }
 	}
    }

  i = m() - 1; vj = n() - 1;
  j = order2[vj];
  max = sm[i][j].path_score;
  j = order2[n() - 1];
  for (i=0; i< (m()-1); ++i)
    {
      if (sm[i][j].path_score > max) {
	max = sm[i][j].path_score;
      }
    }  
  for (vj=0; vj<(n()-1); ++vj) {
    j = order2[vj];
    if (sm[m()-1][j].path_score > max)
      max = sm[m()-1][j].path_score;
  }
  return max;
}

Alignment Score_matrix::dyn_prog()
{
  MATRIX_TYPE max, max2, delta, delta2;
  // MATRIX_TYPE temp;
  uint i, j;
  int k, k1, k2;
  int best_step;

  Score_matrix & sm = *this;

  // prevent core dumps for empty matrix
  if (!(m() * n()))
    {
      Alignment nothing;
      return nothing;
    }

  /* in case gap penalties have changed... */
  update_gap_penalties();

  /* first initialize the first row and column to the
     values in the comparison table.  This assumes that
     'end gaps' are not penalized, but why would you ever
     want to penalize them? */
  for (i=0; i < m(); ++i) 
    {
      row_attrib[i].best_score = 0;
      row_attrib[i].best_redun = 1;
      sm[i][0].path_score = sm(i, 0);
      sm[i][0].n_best = 1;
    }
  for (j=0; j < n(); ++j) 
    {
      column_attrib[j].best_score = 0;
      column_attrib[j].best_redun = 1;
      sm[0][j].path_score = sm(0, j);
      sm[0][j].n_best = 1;
    }

  /* Now go down the matrix and tabulate the maximum path
     score possible up to that point */
  // start at upper left, go systematically through the matrix
  MATRIX_TYPE A, B, C, A1, B1;
  uint ways, ways1, ways2, ways3, choice;
  for (i=1; i < m(); ++i)
    { 
      for (j=1; j < n(); ++j)
	{
	  C = sm[i-1][j-1].path_score; // path score with no gaps
	  k1 = i - 1 - column_attrib[j-1].best_score; // gap size in seq1
	  A1 = sm[i-k1-1][j-1].path_score - sigma * k1 * epen();
	  A = A1 - sigma * column_attrib[j-1].gpen; // path score with gap in seq1
	  k2 = j - 1 - row_attrib[i-1].best_score; // gap size in seq2
	  B1 = sm[i-1][j-k2-1].path_score - sigma * k2 * epen();
	  B = B1 - sigma * row_attrib[i-1].gpen; // path score with gap in seq2
	  max = MAX(MAX(A,B),C);
          // using "suboptimal_factor" as an actual factor is bad, because
          // it has a big effect where the score is large (C-terminus)
          // relative to where it is small (N-terminus)
	  // delta = ABS(max * (1 - sm.suboptimal()));
          //
          // A saner method is to scale it to the actual score distribution
	  delta = ABS(sigma * sm.suboptimal()); // amount of slop permitted

	  if (A < (max - delta)) ways1 = 0;
	  else ways1 = column_attrib[j-1].best_redun;
	  if (B < (max - delta)) ways2 = 0;
	  else ways2 = row_attrib[i-1].best_redun;
	  if (C < (max - delta)) ways3 = 0;
	  else ways3 = sm[i-1][j-1].n_best;
          // how many permissible paths are there?
	  sm[i][j].n_best = ways1 + ways2 + ways3;

	  if (highroad)
	    {
	      if (ways1) best_step = k1;
	      else if (ways3) best_step = 0;
	      else best_step = -k2;
	    }
	  else if (lowroad)
	    {
	      if (ways2) best_step = -k2;
	      else if (ways3) best_step = 0;
	      else best_step = k1;
	    }
	  else if (midroad)
	    {
	      if (ways3) best_step = 0;
	      else if (ways1 > ways2) best_step = k1;
	      else best_step = -k2;
	    }
	  else // random
	    {
	      choice = rand() % (ways1 + ways2 + ways3);
	      if (ways1 > choice) best_step = k1;
	      else if ((ways1 + ways2) > choice) best_step = -k2;
	      else best_step = 0;
	    }
	  sm[i][j].best_step = best_step;
	  if (!best_step) sm[i][j].path_score = sm(i,j) + C;
	  else if (best_step < 0) sm[i][j].path_score = sm(i,j) + B;
	  else sm[i][j].path_score = sm(i,j) + A;

	  // update column and row tables
	  max2 = MAX(A1,C);
	  // delta2 = ABS(max2 * (1 - sm.suboptimal()));
	  delta2 = ABS(sigma * sm.suboptimal());
	  if (A1 < (C - delta2))
	    {
	      column_attrib[j-1].best_score = i-1;
	      column_attrib[j-1].best_redun = sm[i-1][j-1].n_best;
	    }
	  else if (A1 > (C + delta));
	  else
	    {
	      if (highroad) column_attrib[j-1].best_score = i-1;
	      else if (lowroad);
	      else if (midroad);
	      else // random
		{
		  ways = column_attrib[j-1].best_redun + sm[i-1][j-1].n_best;
		  if ((rand() % ways) < sm[i-1][j-1].n_best)
		    column_attrib[j-1].best_score = i-1;
		}
	      column_attrib[j-1].best_redun += sm[i-1][j-1].n_best;
	    }
	  max2 = MAX(B1,C);
	  // delta2 = ABS(max2 * (1 - sm.suboptimal()));
	  delta2 = ABS(sigma * sm.suboptimal());
	  if (B1 < (C - delta2))
	    {
	      row_attrib[i-1].best_score = j-1;
	      row_attrib[i-1].best_redun = sm[i-1][j-1].n_best;
	    }
	  else if (B1 > (C + delta));
	  else
	    {
	      if (highroad) row_attrib[i-1].best_score = j-1;
	      else if (lowroad);
	      else if (midroad);
	      else // random
		{
		  ways = row_attrib[i-1].best_redun + sm[i-1][j-1].n_best;
		  if ((rand() % ways) < sm[i-1][j-1].n_best)
		    row_attrib[i-1].best_score = j-1;
		}
	      row_attrib[i-1].best_redun += sm[i-1][j-1].n_best;
	    }
 	}
    }

  // next reconstruct the alignment from the tables we just made

  for (i=0; i <= m(); ++i) row_attrib[i].newgaps = 0;
  for (i=0; i <= n(); ++i) column_attrib[i].newgaps = 0;

  Alignment answer;
  max_score(i, j);
  answer.score_val = sm[i][j].path_score;
  row_attrib[m()].newgaps = n() - 1 - j;
  column_attrib[n()].newgaps = m() - 1 - i;

  /* now trace it back */
  while ( (i > 0) && (j > 0) )
    {
      k = sm[i][j].best_step;
      if (k == 0) {--i; --j;}
      else if (k > 0)           /* gap in seq 2 */
	{
	  column_attrib[j].newgaps = k;
	  --j; i -= (k+1);
	}
      else                      /* gap in seq 1 */
	{
	  row_attrib[i].newgaps = -k;
	  --i; j -= (1-k);
	}
    }
  if (j) row_attrib[0].newgaps = j;
  else if (i) column_attrib[0].newgaps = i;
  
  // finally use these gap arrays to produce the alignment

  uint s, sumgaps1, sumgaps2, newlength;

  sumgaps1 = 0; sumgaps2 = 0;
  for (i=0; i <= m(); ++i) sumgaps1 += row_attrib[i].newgaps;
  for (i=0; i <= n(); ++i) sumgaps2 += column_attrib[i].newgaps;

  newlength = sumgaps1 + m();
  if ( newlength != (sumgaps2 + n()) )
    {
      cerr << "inconsistent sequence lengths in alignment";
      exit(1);
    }

  const Alignment & al1 = *a1;
  const Alignment & al2 = *a2;
  // Mystring temp_seq = "";
  // Sequence sp;
  /* grind out new sequences */
  for (s=0; s < al1.n_seq(); ++s)
    {
      Sequence temp_seq = al1[s];
      temp_seq = ""; // this loses title information
      for (i=0; i < al1[s].length(); ++i)
	{
	  for (j=0; j < row_attrib[i].newgaps; ++j)
	    temp_seq += '-';
	  temp_seq += al1[s][i];
	}
      for (j=0; j < row_attrib[al1[s].length()].newgaps; ++j)
	temp_seq += '-';
      // sp = al1[s];
      // sp = temp_seq;
      answer += temp_seq;
    }
  for (s=0; s < al2.n_seq(); ++s)
    {
      Sequence temp_seq = al2[s];
      temp_seq = "";
      for (i=0; i < al2[s].length(); ++i)
	{
	  for (j=0; j < column_attrib[i].newgaps; ++j)
	    temp_seq += '-';
	  temp_seq += al2[s][i];
	}
      for (j=0; j < column_attrib[al2[s].length()].newgaps; ++j)
	temp_seq += '-';
      // sp = al2[s]; // to save title, etc
      // sp = temp_seq; // gapped sequence
      answer += temp_seq;
    }
  // count gaps
  answer.newgaps() = 0;
  for (i=1; i < m(); ++i) if (row_attrib[i].newgaps) ++answer.newgaps();
  for (i=1; i < n(); ++i) if (column_attrib[i].newgaps) ++ answer.newgaps();
  return answer;
}


/* return the best score and its matrix position from a path-score
   matrix */
MATRIX_TYPE Score_matrix::max_score(uint & im, uint & jm)
{
  MATRIX_TYPE max, delta;
  uint i, j, redundancy;
  const Score_matrix & sm = *this;

  // prevent core dumps for empty matrix
  if (!(m() * n())) return 0;
  im = m() - 1; jm = n() - 1;
  max = sm[m()-1][n()-1].path_score;
  redundancy = 1;
  // delta = ABS(max - max * suboptimal());
  delta = ABS(sigma * suboptimal());
  for (i=0; i< (m()-1); ++i)
    {
      if (sm[i][n()-1].path_score < (max - delta));
      else if (sm[i][n()-1].path_score > (max + delta))
	{
	  max = sm[i][n()-1].path_score;
	  // delta = ABS(max - max * suboptimal());
	  delta = ABS(sigma * suboptimal());
	  im = i; jm = n() - 1;
	  redundancy = 1;
	}
      else
	{
	  ++ redundancy;
	  max = MAX(max, sm[i][n()-1].path_score);
	  // delta = ABS(max - max * suboptimal());
	  delta = ABS(sigma * suboptimal());
	  if (highroad);
	  else if (lowroad) {im = i; jm = n() - 1;}
	  else if (midroad) {im = i; jm = n() - 1;}
	  else if (randroad)
	    if ((uint)rand()%1000 > (1000/redundancy)) 
	      {im = i; jm = n() - 1;}
	}
    }  
  for (j=0; j<(n()-1); ++j)
      if (sm[m()-1][j].path_score < (max - delta));
      else if (sm[m()-1][j].path_score > (max + delta))
	{
	  max = sm[m()-1][j].path_score;
	  // delta = ABS(max - max * suboptimal());
	  delta = ABS(sigma * suboptimal());
	  im = m() - 1; jm = j;
	  redundancy = 1;
	}
      else
	{
	  ++ redundancy;
	  max = MAX(max, sm[m() - 1][j].path_score);
	  // delta = ABS(max - max * suboptimal());
	  delta = ABS(sigma * suboptimal());
	  if (lowroad);
	  else if (highroad) {im = m() - 1; jm = j;}
	  else if (midroad) {im = m() - 1; jm = j;}
	  else if (randroad)
	    if ((uint) rand()%1000 > (1000/redundancy)) 
	      {im = m() - 1; jm = j;}
	}
  return sm[im][jm].path_score;
}


ostream & operator<<(ostream & os, const Score_matrix & s)
{
  uint i, j;
  for (i=0; i < s.m(); ++i)
    {
      for (j=0; j < s.n(); ++j)
	{
	  os.width(5);
	  os << s[i][j].comp_score;
	  os << "  ";
	}
      os << "\n";
    }
  return os;
}


ostream & Score_matrix::outsubmat(ostream & os, const Score_matrix & s, int x1, int y1, int x2, int y2)
{
  uint i, j;
  uint a, b, c, d;
  // convert relative addresses (negative) to absolute
  if (x1 >= 0) a = x1;
  else a = s.m() + x1;
  if (x2 >= 0) c = x2;
  else c = s.m() + x2;
  if (y1 >= 0) b = y1;
  else b = s.n() + y1;
  if (y2 >= 0) d = y2;
  else d = s.n() + y2;
  // First print out line with sequence letters of seq2
  // I should do some bounds checking eventually...
  for (i=a; i <= c; ++i)
    {
      for (j=b; j <= d; ++j)
	{
          // line 0 contains a blank or connection
	  os.width(5);
	  os << "   .   ";
	}
      for (j=b; j <= d; ++j)
	{
          // line 1 contains comp_score
          // Print out sequence letters of seq1
	  os.width(5);
	  os << s[i][j].comp_score;
	  os << "  ";
	}
      for (j=b; j <= d; ++j)
	{
          // line 2 contains path_score and possibly connection
	  os.width(5);
	  os << s[i][j].path_score;
	  os << "  ";
	}
      os << "\n";
    }
  return os;
}
