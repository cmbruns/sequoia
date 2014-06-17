#include "brookhav.hxx"
#include "parse.hxx"
#include <iostream>

using namespace std;

// $Id: brookhav.cxx,v 1.5 2002/02/19 18:38:15 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/brookhav.cxx,v 1.5 2002/02/19 18:38:15 bruns Exp $
// $Log: brookhav.cxx,v $
// Revision 1.5  2002/02/19 18:38:15  bruns
// Changed exception handling to catch ill conditioned orientations in the
// overlay function (overlay.cxx), to use C++ exeptions, and to not make
// so many warnings with a CA only overlay
//
// Revision 1.4  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.3  2001/12/08 22:12:30  cdputnam
// brookhav.cxx:  cdputnam 12/8/01  Added MSE (Selenomet), HID,HIE,HIP (His) to residue types
//
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

pdbatom::pdbatom(void)  // constructor
{
  strcpy(recordtype, "UNDEF ");
  serialnumber = 0;
  strcpy(atomname, "UNDF");
  atomname[4] = '\0';
  strcpy(altloc, " ");
  altloc[1] = '\0';
  strcpy(residuename, "UND");
  strcpy(chainID, " ");
  residuenumber = 0;
  strcpy(insertioncode, " ");
  insertioncode[1] = '\0';
  occupancy = 0;
  bvalue = 0; // temperature factor
  strcpy(footnote, "   ");
  strcpy(segid, "    ");
  strcpy(element, "  ");
  strcpy(charge, "  ");
}


pdbatom::pdbatom(const char * line)
{
  *this = line;
}

pdbatom::pdbatom(const pdbatom & a2)
{
  *this = a2;
}


pdbatom & pdbatom::operator=(const char * line)
{
  strncpy(recordtype, line, 6);
  recordtype[6] = '\0';
  serialnumber = atoi(line + 6);
  strncpy(atomname, line + 12, 4);
  atomname[4] = '\0';
  strncpy(altloc, line + 16, 1);
  if (strlen(altloc) != 1) strcpy(altloc, " ");
  altloc[1] = '\0';
  strncpy(residuename, line + 17, 3);
  residuename[3] = '\0';
  strncpy(chainID, line + 21, 1);
  chainID[1] = '\0';
  residuenumber = atoi(line + 22);
  strncpy(insertioncode, line + 26, 1);
  insertioncode[1] = '\0';
  coordinate.x() = atof(line + 30);
  coordinate.y() = atof(line + 38);
  coordinate.z() = atof(line + 46);
  occupancy = atof(line + 54);
  bvalue = atof(line + 60);
  if (strlen(line) > 69) strncpy(footnote, line + 67, 3);
  if (strlen(line) > 75) strncpy(segid, line + 72, 4);
  if (strlen(line) > 77) strncpy(element, line + 76, 2);
  if (strlen(line) > 79) strncpy(charge, line + 78, 2);

  if (strlen(footnote) != 3) strcpy(footnote, "   ");
  if (strlen(segid) != 4) strcpy(segid, "    ");
  if (strlen(element) != 2) strcpy(element, "  ");
  if (strlen(charge) != 2) strcpy(charge, "  ");
  footnote[3] = '\0';
  segid[4] = '\0';
  element[2] = '\0';
  charge[2] = '\0';

  return *this;
}


pdbatom & pdbatom::operator=(const pdbatom & a2)
{
  if (this == &a2) return *this;
  strcpy(recordtype, a2.recordtype);
  serialnumber = a2.serialnumber;
  strcpy(atomname, a2.atomname);
  atomname[4] = '\0';
  strcpy(altloc, a2.altloc);
  altloc[1] = '\0';
  strcpy(residuename, a2.residuename);
  strcpy(chainID, a2.chainID);
  residuenumber = a2.residuenumber;
  strcpy(insertioncode, a2.insertioncode);
  insertioncode[1] = '\0';
  coordinate = a2.coordinate;
  occupancy = a2.occupancy;
  bvalue = a2.bvalue;
  strcpy(footnote, a2.footnote);
  strcpy(segid, a2.segid);
  strcpy(element, a2.element);
  strcpy(charge, a2.charge);
  return *this;
}

const Vector3D pdbatom::coord() const
{
  return coordinate;
}

Vector3D & pdbatom::coord()
{
  return coordinate;
}

const double pdbatom::distance(const pdbatom & p2) const
{
  return coord().distance(p2.coord());
}

istream & operator>>(istream & is, pdbatom & a)
{
  char buffer[100];
  is.getline(buffer, 100);
  a = buffer;
  return is;
}


pdbatom pdbatom::transform(const RMat & M) const
{
  pdbatom answer = *this;
  answer.coord() = vec3_transform(M, answer.coord());
  return answer;
}

ostream & operator<<(ostream & os, const pdbatom & a)
{
  std::streamsize oldw;

  os << a.recordtype;
  oldw = os.width(5); os << a.serialnumber; os.width(oldw);
  os << " ";
  os << a.atomname;
  os << a.altloc;
  os << a.residuename;
  os << " ";
  os << a.chainID;
  os.width(4); os << a.residuenumber;
  os << a.insertioncode;
  os << "   ";
  ios::fmtflags old = os.setf(ios::floatfield, ios::fixed);
  os.width(8); os.precision(3); os << a.coordinate.x();
  os.width(8); os.precision(3); os << a.coordinate.y();
  os.width(8); os.precision(3); os << a.coordinate.z();
  os.width(6); os.precision(2); os << a.occupancy;
  os.width(6); os.precision(2); os << a.bvalue;
  os.setf(old, ios::fixed);
  os << " " << a.footnote;
  os << "  " << a.segid;
  os << a.element << a.charge;
  return os;
}


ostream & operator<<(ostream & os, const pdbprotein & p)
{
  int oldresnum = -100000; // start with something impossible
  pdbatom a;
  uint i;
  uint seqposn = 0; // placement within sequence data

  if (p.title().length() > 0)
    os << "TITLE     " << p.title() << "\n";
  // 1) Output any initial end gap cards
  while ((seqposn < p.sequence().dim()) && 
	 is_gap(p.sequence()[seqposn].letter()))
    {
      os << "REMARK GAP\n";
      ++seqposn;
    }
  // 2) Go though all atoms
  for (i=0; i < p.atoms.dim(); ++i)
    {
      a = p.atoms[i];
      os << a << "\n";
      oldresnum = a.res_number();
      // check the sequence data
      if (a.one_letter_code() != '\0')
	{
	  ++seqposn;
	}
      // write out the rest of this residue
      while ((i < p.atoms.dim() - 1) && 
	     (p.atoms[i+1].res_number() == oldresnum))
	{
	  ++i;
	  a = p.atoms[i];
	  os << a << "\n";
	}
      // Are we at a gap in the sequence?
      while ((seqposn < p.sequence().dim()) && 
	     is_gap(p.sequence()[seqposn].letter()))
	{
	  os << "REMARK GAP\n";
	  ++seqposn;
	}
    }
  // 3) Output any final end gap cards
  while ((seqposn < p.sequence().dim()) && 
	 is_gap(p.sequence()[seqposn].letter()))
    {
      os << "REMARK GAP\n";
      ++seqposn;
    }

  return os;
}


pdbprotein & pdbprotein::operator+=(const pdbatom a)
{
  pdbprotein & p = *this;
  atoms += a;
  if (a.one_letter_code() != '\0')
    {
      int oldresnum = -1000000; // impossible value
      if (p.n_atoms() > 1) oldresnum = p[p.n_atoms()-2].res_num();
      // Is this a new residue ?
      if (a.res_number() != oldresnum)
	{
	  SeqRes temp_res;
	  temp_res.letter() = a.one_letter_code();
	  temp_res.res_num() = a.res_num();
	  temp_res.refatomp() = NULL;
	  // Make certain the residue is "sequenceable"
	  if (temp_res.letter() != '\0')
	    {
	      p.sequence() += temp_res;
	    }
	}
    }
  return *this;
}


pdbprotein pdbprotein::ca(void) const
{
  pdbprotein answer;
  pdbprotein temp = *this;
  temp.update_refatoms();
  answer.sequence() = temp.sequence();
  uint i;
  for (i=0; i < temp.sequence().dim(); ++i)
    if (temp.sequence()[i].refatomp() != NULL)
      answer.atoms += temp.sequence()[i].refatom();
  answer.update_refatoms();
  return answer;
}

pdbprotein pdbprotein::transform(const RMat & M) const
{
  pdbprotein answer = *this;
  uint i;
  for (i=0; i < answer.dim(); ++i)
    {
      answer[i] = answer[i].transform(M);
    }  
  return answer;
}


// get one pdb record of a structure from a file
int get_one_record(istream & is, pdbprotein & p)
{
  if (!is.good()) return 1; // end-of-file or something

  char buffer[500];
  Mystring buf2;
  is.getline(buffer, 500); // get one line
  buf2 = buffer;

  // coordinate records get tacked onto our structure
  if ((buf2.at(0,4) == "ATOM") ||
      (buf2.at(0,6) == "HETATM")
      )
    {
      pdbatom currentatom = buffer;
      p += currentatom;
      return 0;
    }
  // Also examine my special gap remark
  // NULL pointer in residue list implies a gap;
  else if (buf2.at(0,10) == "REMARK GAP")
    {
      p.sequence() += '-';
      return 0;
    }
  // TITLE records
  else if (buf2.at(0,5) == "TITLE")
    {
      p.title() = advance_to_next_token(buf2.c_str());
      return 0;
    }
  // ENDMDL
  else if (buf2.at(0,6) == "ENDMDL")
    {
      return 2;
    }
  return 0;
}


istream & operator>>(istream & is, pdbprotein & answer)
{
  pdbprotein p;
  int result;
  while (!(result = get_one_record(is, p)));
  answer = p;
  answer.update_refatoms();
  return is;
}


istream & operator>>(istream & is, pdboverlay & answer)
{
  pdbprotein tempstrct;
  pdboverlay o;
  while (is.good())
    {
      is >> tempstrct;
      if (tempstrct.n_atoms() > 0) o += tempstrct;
    }
  answer = o;
  return is;
}

ostream & operator<<(ostream & os, const pdboverlay & p)
{
  uint i;
  // If there is only one structure, don't use MODEL stuff
  if (p.dim() < 2) os << p[0];
  else for (i=0; i < p.dim(); ++i)
    {
      os << "MODEL " << i + 1 << "\n";
      os << p[i];
      os << "ENDMDL\n";
    }
  return os;
}

pdboverlay pdboverlay::ca(void) const
{
  pdboverlay answer;
  const pdboverlay & o = *this;
  uint i;
  for (i=0; i < dim(); ++i)
    {
      answer += o[i].ca();
    }
  return answer;
}

pdboverlay pdboverlay::transform(const RMat & M) const
{
  pdboverlay answer = *this;
  uint i;
  for (i=0; i < answer.dim(); ++i)
    {
      answer[i] = answer[i].transform(M);
    }
  return answer;
}


int pdboverlay::update_refatoms()
{
  pdboverlay & o = * this;
  uint i;
  int result = 0;
  for (i = 0; i < dim(); ++i)
    result = o[i].update_refatoms();
  return result;
}

int pdboverlay::update_orientations()
{
  pdboverlay & o = * this;
  uint i;
  int result = 0;
  for (i = 0; i < dim(); ++i)
    result = o[i].update_orientations();
  return result;
}

// Put proper 'pointers to reference atoms' into the sequence data
// The expected behavior is to reference the first " CA " atom in the 
//  residue, if it exists; Otherwise ref the first atom of the residue.
int pdbprotein::update_refatoms()
{
  uint atompos = 0; // position in atom array
  uint seqpos = 0; // position in sequence array
  char seqr, atomr;

  // clear out reference array
  if (sequence().dim() < 1) return 0;
  for (seqpos = 0; seqpos < sequence().dim(); ++ seqpos)
    {
      sequence()[seqpos].refatomp() = NULL;
    }
  seqpos = 0; // rewind to beginning

  // check first residue
  if (atoms.dim() < 1) return 0;
  int oldresnum = atoms[atompos].res_num();

  // advance atoms to next intelligible residue
  while ((atompos < atoms.dim() - 1) && (atoms[atompos].one_letter_code() == '\0'))
    ++ atompos;
  if ((atoms[atompos].one_letter_code() == '\0') && (atompos >= atoms.dim()-1))
    return 0;

  // advance sequence to next non-gap
  while ((seqpos < sequence().dim() - 1) && is_gap(sequence()[seqpos].letter()))
    ++seqpos;
  if (is_gap(sequence()[seqpos].letter()) && seqpos >= sequence().dim()-1) 
    return 0;

  // sequence and atoms should now be in register
  seqr = toupper(sequence()[seqpos].letter());
  atomr = toupper(atoms[atompos].one_letter_code());
  if (seqr != atomr) 
    cerr << "WARNING: Sequence and coordinates do not agree. " 
	 << seqpos << seqr << " " << atomr << "\n";
  // Use first atom as refatom for now;
  sequence()[seqpos].refatomp() = & atoms[atompos];

  // Is this an alpha carbon?
  bool is_ca = (!strcmp(sequence()[seqpos].refatom().atomname, " CA "));
  bool was_ca = is_ca;

  oldresnum = atoms[atompos].res_num();
  // go through the rest of the atoms
  for (++atompos; atompos < atoms.dim(); ++atompos)
    {
      // is this a new residue?  If so, update pointer
      if (atoms[atompos].res_num() != oldresnum)
	{
	  // Does it have a sequence letter?, if not, save no refatom
	  if (atoms[atompos].one_letter_code() == '\0')
	    {
	      // advance atoms to next intelligible residue
	      while ((atompos < atoms.dim() - 1) && 
		     (atoms[atompos].one_letter_code() == '\0'))
		++ atompos;
	      if ((atoms[atompos].one_letter_code() == '\0') && (atompos >= atoms.dim()-1))
		return 0;
	      --atompos;
	    }
	  else
	    {
	      // advance sequence to next non-gap
	      ++seqpos;
	      while ((seqpos < sequence().dim() - 1) && is_gap(sequence()[seqpos].letter()))
		++seqpos;
	      if (is_gap(sequence()[seqpos].letter()) && seqpos == sequence().dim()-1)
		return 0;
	      
	      // sequence and atoms should now be in register
	      seqr = toupper(sequence()[seqpos].letter());
	      atomr = toupper(atoms[atompos].one_letter_code());
	      if (seqr != atomr) 
		cerr << "WARNING: Sequence and coordinates do not agree. " 
		  << seqpos << seqr << " " << atomr << "\n";
	      sequence()[seqpos].refatomp() = & atoms[atompos];
	      
	      // Is this an alpha carbon?
	      is_ca = (!strcmp(atoms[atompos].atom_name(), " CA "));
	      was_ca = is_ca;
	    }
	  oldresnum = atoms[atompos].res_num();
	}
      else // Not a new residue, see if reference atom might be better
	{
	  if (!was_ca) // if we already have CA, forget it
	    {
	      is_ca = (!strcmp(atoms[atompos].atom_name(), " CA "));
	      if (is_ca) // found a CA atom, update reference
		{
		  sequence()[seqpos].refatomp() = & atoms[atompos];
		  was_ca = true;
		}
	    }
	}
    }
  return 0;
}


// Make some pointers to make residue traversal easier
int pdbprotein::update_residues()
{
  uint i;
  pdbprotein & p = *this;
  int oldresnum = -1000000; // impossible value
  pdbatom * oldatom = NULL;

  for (i=0; i < n_atoms(); ++i)
    {
      pdbatom & a = p[i];

      // Put in default values
      a.next_in_res() = NULL;
      a.first_in_res() = NULL;
      a.next_in_prot() = NULL;
      // Is this a new residue?
      int resnum = a.res_num();
      if (resnum != oldresnum) // Yes, a new residue
	{
	  a.first_in_res() = &a; // This is the first atom of the
			       // residue
	  if (oldatom != NULL)
	    {
	      oldatom->next_in_prot() = &a;
	      oldatom->next_in_res() = NULL;
	    }
	}
      else // Same old residue
	{
	  if (oldatom != NULL)
	    {
	      a.first_in_res() = oldatom->first_in_res();
	      oldatom->next_in_prot() = &a;
	      oldatom->next_in_res() = &a;
	    }	
	}
      oldatom = &a;
      oldresnum = resnum;
    }

  return 0;
}


SeqRes ref_ala;
pdbprotein ref_prot;

// Make reference residue for orientation calculations
int init_ref_ala()
{
  ref_prot.clear();
  pdbatom a;
  strcpy(a.res_name(),"ALA");
  a.res_num() = 1;

  // No 'O' atom, because its local orientation varies
  // ATOM    626  N   ALA    98      29.675   3.008  18.066  1.00  9.09
  // ATOM    627  CA  ALA    98      28.932   2.345  19.125  1.00  8.94
  // ATOM    628  C   ALA    98      29.695   2.201  20.434  1.00 11.42
  // ATOM    630  CB  ALA    98      28.475   0.975  18.555  1.00  7.45

  strcpy(a.atom_name()," N  ");
  a.coord().x() = 29.68;
  a.coord().y() =  3.01;
  a.coord().z() = 18.07;
  ref_prot += a;
  
  strcpy(a.atom_name()," CA ");
  a.coord().x() = 28.93;
  a.coord().y() =  2.35;
  a.coord().z() = 19.13;
  ref_prot += a;
  
  strcpy(a.atom_name()," C  ");
  a.coord().x() = 29.70;
  a.coord().y() =  2.20;
  a.coord().z() = 20.43;
  ref_prot += a;
  
  strcpy(a.atom_name()," CB ");
  a.coord().x() = 28.48;
  a.coord().y() =  0.98;
  a.coord().z() = 18.56;
  ref_prot += a;
  
  ref_prot.update_residues();
  ref_prot.update_refatoms();

  ref_ala = ref_prot.sequence()[0];

  return 0;
}


// Calculate rotation matrices which place a reference alanine residue
// onto each residue of the structure; Uses only those atoms in common
// between the reference and each residue
int pdbprotein::update_orientations()
{
  pdbprotein & p = *this;
  // Make sure everything is kosher
  p.update_residues();
  p.update_refatoms();
  init_ref_ala();

  uint i;
  bool has_rotation_problem = false;
  for (i=0; i < p.sequence().dim(); ++i) {
    SeqRes & r = p.sequence()[i];
    // remove translational component, it is useless
    try {
      r.orientation() = overlay(r, ref_ala).strip_trans();
    } catch (ill_conditioned_rotation_error) {
      // Don't throw exceptions - some residues don't have orientations
      has_rotation_problem = true;
    }
  }
  if (has_rotation_problem)
    cerr << "*** One or more residues lack a well defined orientation (C-alpha only?) ***" << endl;
  return 0;
}


// *** NON MEMBER/FRIEND FUNCTIONS GO HERE ***

// Return NULL character if not a normal amino acid
// or something that should be included in a sequence
char pdbatom::one_letter_code() const
{
  char answer = '\0'; // default to NULL
  const pdbatom & p = *this;
  if (!strcmp(p.res_name(), "ALA")) answer = 'A';
  else if (!strcmp(p.res_name(), "CYS")) answer = 'C';
  else if (!strcmp(p.res_name(), "ASP")) answer = 'D';
  else if (!strcmp(p.res_name(), "GLU")) answer = 'E';
  else if (!strcmp(p.res_name(), "PHE")) answer = 'F';
  else if (!strcmp(p.res_name(), "GLY")) answer = 'G';
  else if (!strcmp(p.res_name(), "HIS")) answer = 'H';
  else if (!strcmp(p.res_name(), "HIP")) answer = 'H';
  else if (!strcmp(p.res_name(), "HIE")) answer = 'H';
  else if (!strcmp(p.res_name(), "HIP")) answer = 'H';
  else if (!strcmp(p.res_name(), "ILE")) answer = 'I';
  else if (!strcmp(p.res_name(), "LYS")) answer = 'K';
  else if (!strcmp(p.res_name(), "LEU")) answer = 'L';
  else if (!strcmp(p.res_name(), "MET")) answer = 'M';
  else if (!strcmp(p.res_name(), "MSE")) answer = 'M';
  else if (!strcmp(p.res_name(), "ASN")) answer = 'N';
  else if (!strcmp(p.res_name(), "PRO")) answer = 'P';
  else if (!strcmp(p.res_name(), "GLN")) answer = 'Q';
  else if (!strcmp(p.res_name(), "ARG")) answer = 'R';
  else if (!strcmp(p.res_name(), "SER")) answer = 'S';
  else if (!strcmp(p.res_name(), "THR")) answer = 'T';
  else if (!strcmp(p.res_name(), "VAL")) answer = 'V';
  else if (!strcmp(p.res_name(), "TRP")) answer = 'W';
  else if (!strcmp(p.res_name(), "TYR")) answer = 'Y';
  // ignore waters and other things
  // else answer = 'X';
  return answer;
}

// *** SEQUENCE/ALIGNMENT functions that involve coordinates ***

Sequence & Sequence::operator=(const pdbprotein & p)
{
  *this = p.sequence();
  title() = p.title();
  coord_ptr() = (pdbprotein *) &p; // discard const
  return *this;
}

Alignment & Alignment::operator=(const pdbprotein & p)
{
  Sequence s;
  s = p;
  *this = s;
  return *this;
}


// Calculate RMS deviation for CA atoms capitalized in alignment
Real Alignment::rms(const pdbprotein & p1, const pdbprotein & p2) const
{
  Real answer = 0;
  pdbprotein s1 = p1.ca();
  pdbprotein s2 = p2.ca();

  // Fudge error checking
  uint max = MIN(s1.n_atoms(),s2.n_atoms());
  uint max2 = MIN(seq[0].length(), seq[1].length());
  uint i, r1, r2;
  r1 = r2 = 0;
  uint count = 0;
  Real diff;
  for (i=0; i < max2; ++i)
    {
      if (!is_gap(seq[0][i].letter())) ++r1;
      if (!is_gap(seq[1][i].letter())) ++r2;
      if (r1 >= max) i = max2; // break loop
      else if (r2 >= max) i = max2; // break loop
      else if ((!is_gap(seq[0][i].letter())) 
	       && (!is_gap(seq[1][i].letter()))
	       && (seq[0][i].letter() == toupper(seq[0][i].letter()))
	       && (seq[1][i].letter() == toupper(seq[1][i].letter())))
	{
	  diff = s1[r1].coord().distance(s2[r2].coord());
	  answer += diff*diff;
	  ++ count;
	}
    }
  if (!count) return 0;
  else return sqrt(answer/double(count));
}


// Calculate RMS deviation for CA atoms capitalized in alignment
Real Alignment::rms(const pdboverlay & o1, const pdboverlay & o2) const
{
  Real answer = 0;
  pdbprotein s1, s2;
  Sequence sq1, sq2;
  uint ind1, ind2;
  const Alignment & a = *this;

  uint count = 0;
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

	// Fudge error checking
	uint max = MIN(s1.n_atoms(),s2.n_atoms());
	uint max2 = MIN(sq1.length(), sq2.length());
	uint i, r1, r2;
	r1 = r2 = 0;
	Real diff;
	for (i=0; i < max2; ++i)
	  {
	    if (!is_gap(sq1[i].letter())) ++r1;
	    if (!is_gap(sq2[i].letter())) ++r2;
	    if (r1 >= max) i = max2; // break loop
	    else if (r2 >= max) i = max2; // break loop
	    else if ((!is_gap(sq1[i].letter())) 
		     && (!is_gap(sq2[i].letter()))
		     && (sq1[i].letter() == toupper(sq1[i].letter()))
		     && (sq2[i].letter() == toupper(sq2[i].letter())))
	      {
		diff = s1[r1].coord().distance(s2[r2].coord());
		answer += diff*diff;
		++ count;
	      }
	  }
      }
  if (!count) return 0;
  else return sqrt(answer/double(count));
}


// Deposit atom pointer information directly into sequence data type
int Alignment::transfer_atoms(const pdboverlay & o1, 
			      const pdboverlay & o2)
{
  Alignment & a = *this;
  // only count smaller of align.dim() and overlay.dim()
  uint maxind = o1.dim() + o2.dim();
  if (dim() < o1.dim()) maxind = dim();

  uint i;
  // cycle through structures
  for (i = 0; i < maxind; ++i)
    {
      if (i < o1.dim()) a[i].transfer_atoms(o1[i]);
      else a[i].transfer_atoms(o2[i - o1.dim()]);
    }
  return 0;
}


// Deposit atom pointer information directly into sequence data type
int Alignment::transfer_atoms(const pdboverlay & o1)
{
  Alignment & a = *this;
  // only count smaller of align.dim() and overlay.dim()
  uint maxind = o1.dim();
  if (dim() < o1.dim()) maxind = dim();

  uint i;
  // cycle through structures
  for (i = 0; i < maxind; ++i)
    {
      a[i].transfer_atoms(o1[i]);
    }
  return 0;
}


// Deposit atom pointer information directly into sequence data type
int Sequence::transfer_atoms(const pdbprotein & p)
{
  Sequence & s = *this;
  uint i;
  uint str_res = 0;
  for (i = 0; i < s.dim(); ++i)
    {
      SeqRes & r1 = s[i];
      if (!is_gap(r1.letter()))
	{
	  const SeqRes & r2 = p.sequence()[str_res];
	  r1.refatomp() = r2.refatomp();
	  r1.orientation() = r2.orientation();
	  ++str_res;
	}
      else 
	{
	  r1.refatomp() = NULL;
	  r1.orientation() = eye(1) - eye(1); // zero matrix
	}
    }
  return 0;
}


Alignment & Alignment::operator=(const pdboverlay & o)
{
  Sequence s;
  if (o.dim() < 1) s = "";
  else s = o[0];
  *this = s;
  uint i;
  for (i=1; i<o.dim(); ++i)
    {
      s = o[i];
      *this += s;
    }
  return *this;
}

