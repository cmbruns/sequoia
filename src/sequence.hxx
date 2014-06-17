#ifndef _SEQUENCE_HXX_
#define _SEQUENCE_HXX_

// $Id: sequence.hxx,v 1.3 2002/06/28 17:17:46 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/sequence.hxx,v 1.3 2002/06/28 17:17:46 bruns Exp $
// $Log: sequence.hxx,v $
// Revision 1.3  2002/06/28 17:17:46  bruns
// Changed to avoid warnings with gcc 3.1
//   changed <iostream.h> type includes to <iostream>
//   removed redundant function parameter default values from .cxx files
//     (when they are already defined in the header)
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

#include <iostream>
#include "mystring.hxx"
#include "matrix.hxx"
#include "comp_mat.hxx"
#include "rmat.hxx"
// #include "brookhav.hxx"

class pdbatom;
class pdbprotein;
class pdboverlay;

bool is_gap(const char c);

class SeqRes // single residue of a sequence
{
  friend class Sequence;
private:
  char onelettercode;
  int residuenumber;
  pdbatom * refatom_ptr; // points to "reference" atom of structure
  RMat priv_orientation;
public:
  SeqRes(void) {onelettercode = '\0'; residuenumber = 0; refatom_ptr = NULL;}
  const int res_num() const {return residuenumber;}
  int & res_num() {return residuenumber;}
  const char letter() const {return onelettercode;}
  char & letter() {return onelettercode;}
  bool have_ref() const {return (refatom_ptr != NULL);}
  const pdbatom & refatom() const {return *refatom_ptr;}
  pdbatom & refatom() {return *refatom_ptr;}
  pdbatom *const & refatomp() const {return refatom_ptr;}
  pdbatom * & refatomp() {return refatom_ptr;}
  bool operator==(const SeqRes & s2) const;
  bool operator!=(const SeqRes & s2) const;
  bool have_coordinates() const {return have_ref();}
  RMat & orientation() {return priv_orientation;}
  const RMat & orientation() const {return priv_orientation;}
};

ostream & operator<<(ostream & os, const SeqRes & c);

class Sequence : public Array<SeqRes>
{
friend ostream & operator<<(ostream & os, const Sequence & s);
friend istream & operator>>(istream & is, Sequence & s);
friend class Alignment;
private:
  Mystring title_str;
  double weight_val;
  pdbprotein * coord_ptr_val; // points to atomic coordinates, if any
  void initialize_internals();
public:
  Sequence(void);
  Sequence(const Sequence &);
  // Sequence(const Sequence & s) : Array<SeqRes> (const Array<SeqRes> & s)
  //   {
  //     weight() = s.weight();
  //     title() = s.title();
  //   }
  Sequence(unsigned int n, uchar c = '-');
  Sequence(const char * c);
  ~Sequence();
  virtual Sequence & operator=(const Sequence &);
  Sequence & operator=(const char * cp);
  // operator const char *() const;
  //  Sequence & operator=(const Mystring &);
  Sequence & operator=(const pdbprotein & p);
  Sequence & operator+=(const char c);
  Sequence & operator+=(const SeqRes sr);
  bool operator==(const Sequence & s2) const;
  bool operator!=(const Sequence & s2) const;
  unsigned int length() const {return dim();}
  double identity(const Sequence & s) const;
  const double & weight() const {return weight_val;}
  double & weight() {return weight_val;}
  const Mystring title() const {return title_str;}
  Mystring & title() {return title_str;}

  pdbprotein * coord_ptr() const {return coord_ptr_val;}
  pdbprotein * & coord_ptr() {return coord_ptr_val;}
  pdbprotein & coordinates() {return *coord_ptr();}
  const pdbprotein & coordinates() const {return *coord_ptr();}

  Sequence remove_gaps() const;
  uint n_aligned(const Sequence & s2) const;
  uint n_cap_aligned(const Sequence & s2) const;
  Sequence upcase() const;
  Sequence downcase() const;
  Sequence & remove(const char c);
  int first() const // number of first residue
    {
      const Sequence & s = *this;
      if (length() < 1) return 0;
      else return s[0].res_num();
      // return first_ordinal;
    }
    int transfer_atoms(const pdbprotein & p);
  bool have_coordinates() const 
    {
      const Sequence & t = *this;
      uint i = 0; 
      for (i=0;i < dim(); ++i) 
	if (t[i].have_coordinates()) return true;
      return false;
    }
};

void seqcpy(char * s, const Sequence & seq);


class Alignment
{
friend ostream & operator<<(ostream & os, const Alignment & a);
friend istream & operator>>(istream & is, Alignment & a);
private:
  unsigned int nseq;
  Sequence * seq;
  MATRIX_TYPE mean;
  double z_score;
  int new_gaps;
public:
  MATRIX_TYPE score_val;
  Alignment(void);
  ~Alignment(void);
  Alignment(const Alignment & a);
  Alignment(const Sequence & s);
  Sequence & operator[](unsigned int m) {return seq[m];}
  const Sequence & operator[](unsigned int m) const {return seq[m];}
  Alignment & operator=(const Alignment & a);
  Alignment & operator=(const pdbprotein & p);
  Alignment & operator=(const pdboverlay & o);
  Alignment & operator+=(const Alignment & a);
  Alignment & operator+=(const Sequence & s);
  bool operator==(const Alignment & a) const;
  bool operator!=(const Alignment & a) const{return !(*this == a);}
  MATRIX_TYPE score(const Alignment & a2, const Comparison_matrix & cm) const;
  unsigned int n_seq() const {return nseq;}
  unsigned int dim() const {return n_seq();}
  ostream & summarize(ostream & os) const;
  ostream & pretty(ostream & os, uint line_len = 50) const;
  ostream & long_output(ostream & os) const;
  ostream & short_output(ostream & os) const;
  unsigned int length() const 
    {if (nseq) return seq[0].length(); else return 0;}
  void clear();
  int & newgaps() {return new_gaps;}
  void set_weights();
  double identity(uint posn) const;
  int nongaps(uint posn) const;
  Alignment remove_gaps() const;
  Real rms(const pdbprotein & p1, const pdbprotein & p2) const;
  Real rms(const pdboverlay & o1, const pdboverlay & o2) const;
  Alignment upcase() const;
  Alignment downcase() const;
  int transfer_atoms(const pdboverlay & o);
  int transfer_atoms(const pdboverlay & o1, const pdboverlay & o2);
  bool have_coordinates() const 
    {
      const Alignment & t = *this;
      uint i = 0; 
      for (i=0;i < dim(); ++i) 
	if (t[i].have_coordinates()) return true;
      return false;
    }
};

#endif
