// brookhav.hxx 1-26-96 Chris Bruns
// C++ routines for protein data bank files

#ifndef _BROOKHAV_HXX_
#define _BROOKHAV_HXX_

// $Id: brookhav.hxx,v 1.3 2002/06/28 17:21:54 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/brookhav.hxx,v 1.3 2002/06/28 17:21:54 bruns Exp $
// $Log: brookhav.hxx,v $
// Revision 1.3  2002/06/28 17:21:54  bruns
// Minor changes to avoid compiler warnings with gcc 3.1
//   changed <iostream.h> style includes to modern <iostream> style
//   removed redundant default function parameter values when already specified in
//    the header
//
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

using namespace std;

# include <iostream>
# include <fstream>
# include "mystring.hxx"
# include "rvec.hxx"
# include "rmat.hxx"
# include "sequence.hxx"
# include "overlay.hxx"

/*
00000000011111111112222222222333333333344444444445555555555666666666677777777778
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM      8  OE1 GLU     1      24.176   3.462  31.843  1.00 43.26      1PBP 131
 */

class Sequence;
class pdbprotein;
class pdboverlay;
class pdbatom;

typedef pdbatom * atomptr;

class pdbatom
{
friend class pdbprotein;
private:
  char recordtype[7]; // i.e. "ATOM  "
  int serialnumber;
  char atomname[5];
  char altloc[2]; // alternate location indicator
  char residuename[4];
  char chainID[2];
  int residuenumber;
  char insertioncode[2];
  Vector3D coordinate; // orthogonal coordinates
  float occupancy;
  float bvalue; // temperature factor
  char footnote[4];
  char segid[5];
  char element[3];
  char charge[3];
  // For easier traversal of individual residues
  pdbatom * priv_next_in_res; // next in residue
  pdbatom * priv_first_in_res; // first in residue
  pdbatom * priv_next_in_prot; // next in protein
public:
  pdbatom(); // constructor
  pdbatom(const char * line);
  pdbatom(const pdbatom & a2);

  pdbatom & operator=(const char * line);
  pdbatom & operator=(const pdbatom & a2);

  const char * res_name() const {return residuename;}
  const char * atom_name() const {return atomname;}
  const int & res_number() const {return residuenumber;}
  const int & res_num() const {return residuenumber;}
  char * res_name() {return residuename;}
  char * atom_name() {return atomname;}
  int & res_number() {return residuenumber;}
  int & res_num() {return residuenumber;}
  atomptr & next_in_res() {return priv_next_in_res;}
  atomptr & first_in_res() {return priv_first_in_res;}
  atomptr & next_in_prot() {return priv_next_in_prot;}
  const atomptr & next_in_res() const {return priv_next_in_res;}
  const atomptr & first_in_res() const {return priv_first_in_res;}
  const atomptr & next_in_prot() const {return priv_next_in_prot;}
  const double distance(const pdbatom & p2) const;
  const Vector3D coord() const;
  Vector3D & coord();
  pdbatom transform(const RMat & M) const;
  char one_letter_code() const;

friend ostream & operator<<(ostream & os, const pdbatom & a);
friend istream & operator>>(istream & is, pdbatom & a);
friend ostream & operator<<(ostream & os, const pdboverlay & a);
};

class pdbprotein
{
friend class pdboverlay;
protected:
  Sequence prot_seq; // To contain string of actual sequence
  Array<pdbatom> atoms;
  Mystring title_str;
public:
  void clear()
  {
    atoms.clear_array();
    // residues.clear_array();
    title() = "";
  }
  pdbatom & operator[](uint n) {return atoms[n];}
  const pdbatom & operator[](uint n) const {return atoms[n];}
  pdbprotein & operator+=(const pdbatom a);
  const uint dim() const {return atoms.dim();}
  const uint n_atoms() const {return dim();}
  pdbprotein ca(void) const;
  pdbprotein transform(const RMat & M) const;
  const Mystring & title() const {return title_str;}
  Mystring & title() {return title_str;}
  const pdbprotein residue(uint i);
  const Sequence & sequence() const {return prot_seq;}
  Sequence & sequence() {return prot_seq;}
  int update_refatoms();
  int update_residues(); // atomic reference pointers
  int update_orientations();
friend ostream & operator<<(ostream & os, const pdbprotein & p);
friend istream & operator>>(istream & is, pdbprotein & p);
friend istream & operator>>(istream & is, pdboverlay & p);
};

class pdboverlay : public Array<pdbprotein>
{
public:
  pdboverlay(int n) : Array<pdbprotein>(n) {}
  pdboverlay() : Array<pdbprotein>() {}
  const uint n_atoms() const 
  {
    const pdboverlay & p = *this;
    uint i;
    int answer = 0;
    if (dim() < 1) return 0;
    for (i=0; i<dim(); ++i)
      answer += p[i].n_atoms();
    return answer;
  }
  pdboverlay ca(void) const;
  pdboverlay transform(const RMat & M) const;
  int update_refatoms();
  int update_orientations();
};

int get_one_record(istream & is, pdbprotein & p);

#endif
