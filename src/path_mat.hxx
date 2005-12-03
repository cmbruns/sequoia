#ifndef PATH_MAT_HXX
#define PATH_MAT_HXX

// $Id: path_mat.hxx,v 1.2 2001/11/28 23:40:22 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/path_mat.hxx,v 1.2 2001/11/28 23:40:22 bruns Exp $
// $Log: path_mat.hxx,v $
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

#include "array2d.hxx"
#include "matrix.hxx"
#include "sequence.hxx"
#include "brookhav.hxx"

#define DEFAULT_GAP_PENALTY 10
#define DEFAULT_EXTENSION_PENALTY 0.01
#define DEFAULT_CUTOFF 4.5
#define DEFAULT_ACUTOFF 45
#define DEFAULT_SUBOPTIMAL_FACTOR 0.10
#define DEFAULT_RUNLENGTH 4

class Score_matrix;
class cli;

class path_element
{
friend class Score_matrix;
friend ostream & operator<<(ostream & os, const Score_matrix & s);
private:
  MATRIX_TYPE comp_score; // comparison score of res1 with res2
  MATRIX_TYPE path_score; // score of best path up to here
  int best_step; // signed int stores relative pointer to next residue
  uint n_best; // how many maximal paths use this element?
public:
  path_element();
};


// parameters which vary along rows and columns
class residue_attribute
{
friend class Score_matrix;
private:
  MATRIX_TYPE gpen;
  MATRIX_TYPE score_offset;
  uint best_score; // saves best score during dyn_prog algorithm
  uint best_redun; // save number of ties for best score
  uint newgaps;
  residue_attribute();
public:
};

class Score_matrix : public Array2d<path_element>
{
  friend class cli;
  friend class command_token;
  friend class sequoia_variable;
private:
  residue_attribute * row_attrib;
  residue_attribute * column_attrib;
  MATRIX_TYPE sigma; // rms score for gap penalty calculation
  bool have_table; // is comp_score up to date?
  bool have_stable;
  bool have_path; // is path_score up to date?
  bool have_seq1;
  bool have_seq2;
  bool have_align;
  bool have_strct1;
  bool have_strct2;
  bool have_overlay;
  bool highroad;
  bool lowroad;
  bool midroad;
  bool randroad;
  const Alignment * a1;
  const Alignment * a2;

  Score_matrix(uint r, uint c);
  void init_score_mat(uint r, uint c);
  void update_gap_penalties();

  // alignment behavior variables
  double internal_random_seed;
  double internal_useangle;  // use orientation in coordinate overlay?
  double internal_suboptimal; // cutoff for equivalent paths
  double internal_runlength; // minimum length of stretch of equivalent residues
  MATRIX_TYPE internal_epen; // gap extension penalty
  MATRIX_TYPE internal_gpen;
  double internal_acutoff; // angular distance cutoff, in degrees
  double internal_dcutoff; // c alpha distance cutoff

public:
  Score_matrix();
  Score_matrix(const Alignment & a1, const Alignment & a2, const Comparison_matrix & cm);
  Score_matrix(const Score_matrix & sm);
  ~Score_matrix();

  Score_matrix & operator=(const Score_matrix & s2);
  const MATRIX_TYPE & operator()(uint i, uint j) const; // indexing
  MATRIX_TYPE & operator()(uint i, uint j); // indexing

  // * alignment behavior variables *
  const MATRIX_TYPE & epen() const {return internal_epen;}
  MATRIX_TYPE & epen() {return internal_epen;}
  const MATRIX_TYPE & gpen() const {return internal_gpen;}
  MATRIX_TYPE & gpen() {return internal_gpen;}
  const double & acutoff() const {return internal_acutoff;}
  double & acutoff() {return internal_acutoff;}
  const double & dcutoff() const {return internal_dcutoff;}
  double & dcutoff() {return internal_dcutoff;}
  const double & random_seed() const {return internal_random_seed;}
  double & random_seed() {return internal_random_seed;}
  const double & runlength() const {return internal_runlength;}
  double & runlength() {return internal_runlength;}
  const double & suboptimal() const {return internal_suboptimal;}
  double & suboptimal() {return internal_suboptimal;}
  const double & useangle() const {return internal_useangle;}
  double & useangle() {return internal_useangle;}

  void seed(int s);
  MATRIX_TYPE mean() const;
  MATRIX_TYPE std_dev() const;
  void tabulate(const Alignment & a1, const Alignment & a2, const Comparison_matrix & cm);
  void stabulate(const pdbprotein & p1, const pdbprotein & p2);
  void stabulate(const pdboverlay & p1, const pdboverlay & p2);
  double res_3dscore(const SeqRes & r1, const SeqRes & r2);
  Alignment dyn_prog();
  MATRIX_TYPE Score_matrix::randomize();
  MATRIX_TYPE max_score(uint & im, uint & jm);
  ostream & outsubmat(ostream & os, const Score_matrix & s, int x1, int y1, int x2, int y2);

friend ostream & operator<<(ostream & os, const Score_matrix & s);
// for outputting sub-matrices
// fails because '<<' must take two arguments; so just rename it
// friend ostream & operator<<(ostream & os, const Score_matrix & s, int x1, int y1, int x2, int y2);
};

#endif
