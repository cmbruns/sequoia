// This file is part of the Sequoia package for macromolecular 
//  sequence/structure analysis
// Copyright (C) 2004  Christopher M. Bruns, Ph.D.
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// 
// See the accompanying file 'LICENSE' for details
// 
// To contact the author, write to cmbruns@comcast.net or bruns@scripps.edu
// In publications please cite: Bruns et al (1999), J.Mol.Biol. 288:427-439
// Please submit bug reports at http://bruns.homeip.net/bugzilla/
// 
// To obtain a non-GPL version of this program, see http://bruns.homeip.net/sequoia.html
// 
// $Id$
//
// $Header$
//
// $Log$
// Revision 1.1  2004/06/04 19:34:47  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.5  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.4  2002/09/13 22:41:16  bruns
// removed redundant parameter default
//
// Revision 1.3  2002/05/15 19:22:49  bruns
// Changed get_atom calls to return PDBAtom, not pointer, and to use
// modern exception handling.
// Added incomplete set_torsion() routine, which only works for leucine side
// chains
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "PDBResidue.h"
#include "Matrix3D.h"
#include "Matrix4D.h"
#include <cstdio>
#include "Hash3D.h"

double PDBResidue::store_accessible_surface_area(const Hash3D & hash3d) {
	PDBResidue & residue = *this;
	
	double answer = 0.0;
	
	for (unsigned int a = 0; a < residue.get_atom().size(); a++) {
		PDBAtom & atom = residue.atom()[a];
		// TODO - restrict to non hydrogen protein atoms
		double atom_area = atom.store_accessible_surface_area(hash3d);
		// cout << atom_area << endl;
		answer += atom_area;
	}
	
	p_accessible_surface_area = answer;
	return answer;
}

double PDBResidue::compute_accessible_surface_area(const Hash3D & hash3d) const {
	const PDBResidue & residue = *this;

	double answer = 0.0;
	
	for (unsigned int a = 0; a < residue.get_atom().size(); a++) {
		const PDBAtom & atom = residue.get_atom()[a];
		// TODO - restrict to non hydrogen protein atoms
		double atom_area = atom.compute_accessible_surface_area(hash3d);
		// cout << atom_area << endl;
		answer += atom_area;
	}
	
	return answer;
}

const vector<PDBAtom> & PDBResidue::get_atom() const {
  return private_atom;
}

vector<PDBAtom> & PDBResidue::atom() {
  return private_atom;
}

const PDBAtom PDBResidue::get_atom(string atom_name) const {
  const PDBAtom * answer = NULL;
  vector<PDBAtom>::const_iterator atom_ptr;
  for (atom_ptr = get_atom().begin(); 
       atom_ptr != get_atom().end();
       atom_ptr ++)
    if (atom_ptr->get_atom_name() == atom_name) {
      answer = &(*atom_ptr);
      break;
    }
  if (answer == NULL) throw no_such_atom_error(atom_name);
  return *answer;
}

PDBAtom & PDBResidue::atom(string atom_name) {
  PDBAtom * answer = NULL;
  vector<PDBAtom>::iterator atom_ptr;
  for (atom_ptr = atom().begin(); 
       atom_ptr != atom().end();
       atom_ptr ++)
    if (atom_ptr->get_atom_name() == atom_name) {
      answer = &(*atom_ptr);
      break;
    }
  if (answer == NULL) throw no_such_atom_error(atom_name);
  return *answer;
}

void PDBResidue::add_atom(const PDBAtom & atom) {
  if (atom.get_record_name() == "ATOM  ") private_is_hetatm_only = false;
  if (atom.get_alt_loc() != ' ') private_has_altloc = true;
  private_atom.push_back(atom);
}

bool PDBResidue::is_protein() const {
  // set_one_letter_code();
  return private_is_protein;
}

bool PDBResidue::is_hetatm_only() const {
  return private_is_hetatm_only;
}

bool PDBResidue::has_modres() const {
  return private_has_modres;
}

bool PDBResidue::is_terminus() const {
  return private_is_terminus;
}

bool PDBResidue::has_altloc() const {
  return private_has_altloc;
}

bool PDBResidue::is_redundant() const {
  return private_is_redundant;
}

void PDBResidue::check_redundancy(const PDBResidue * previous_residue) {
  private_is_redundant = false;

  // Residue with no atoms (TER only?)
  if (get_atom().size() < 1)
    private_is_redundant = true;
  
  if (previous_residue == NULL) return;

  // Microheterogeneity and alternate conformations
  if ((previous_residue->get_residue_number() == get_residue_number())) {
    // Insertion code is used to mark microheterogeneity in 1al4
    // && previous_residue->get_insertion_code() == get_insertion_code()) {

    // Both have non-blank alternate location fields? -> they are the same
    if (
	(previous_residue->has_altloc())
	&& (has_altloc())
	&& (previous_residue->get_insertion_code() == get_insertion_code())
	) {
      private_is_redundant = true;
    }

    // Both have close CA (protein) coordinates? -> they are the same
    try {
      const PDBAtom & ca1 = get_atom(" CA ");
      const PDBAtom & ca2 = previous_residue->get_atom(" CA ");
      // 1.0 Angstrom cutoff for "same" (I made this up)
      if (ca1.distance(ca2) < 1.0) 
	private_is_redundant = true;      
    } catch (no_such_atom_error) {}

    // Both have close C4* (nucleotide) coordinates? -> they are the same
    try {
      const PDBAtom & ca1 = get_atom(" C4*");
      const PDBAtom & ca2 = previous_residue->get_atom(" C4*");
      // 1.0 Angstrom cutoff for "same" (I made this up)
      if (ca1.distance(ca2) < 1.5) // pdb1fj.ent, chain F needs more than one Angstrom
	private_is_redundant = true;      
    } catch (no_such_atom_error) {}

    // Same residue number
    // Both have same standard_residue_name, and lack distant CA coordinates?
    if (get_standard_residue_name() == previous_residue->get_standard_residue_name()) {
      try {
	const PDBAtom & ca1 = get_atom(" CA ");
	const PDBAtom & ca2 = previous_residue->get_atom(" CA ");
	if (ca1.distance(ca2) < 1.0) 
	  private_is_redundant = true;
      } catch (no_such_atom_error) {
	private_is_redundant = true;
      }
    }

  }

  // Microheterogeneity and alternate conformations 2
  if ((previous_residue->get_residue_number() == get_residue_number())
      && (previous_residue->has_altloc())
      && (has_altloc()))
    private_is_redundant = true;

  // extra residue for C-terminus
  if ((get_atom().size() == 1) 
      && (get_atom()[0].get_atom_name() == " OXT")
      && (get_residue_name() == previous_residue->get_residue_name()))
    private_is_redundant = true;

}

string generate_unique_residue_key(string residue_name, 
				   int residue_number,
				   char insertion_code) {
  char key_buffer[50];
  sprintf(key_buffer, "%3s%4d%c", residue_name.c_str(), residue_number, insertion_code);
  string key = key_buffer;
  return key;
}

char three_letter_to_one_letter_other_code(string residue_name) {
  if (residue_name == "UNK") return 'X'; // explicitly uncertain
  
  // Ignore solvent, acetyl group, empty string
  if (residue_name == "H2O") return '\0';
  if (residue_name == "WAT") return '\0';
  if (residue_name == "HOH") return '\0';
  if (residue_name == "ACE") return '\0'; // acetyl group
  if (residue_name == "   ") return '\0';
  if (residue_name == "  ") return '\0';
  if (residue_name == " ") return '\0';
  if (residue_name == "") return '\0';
  
  return 'X'; // default
}

char three_letter_to_one_letter_nucleotide_code(string residue_name) {
  // Nucleotide residues
  if (residue_name == "  A") return 'A';
  if (residue_name == " +A") return 'A';
  if (residue_name == "  C") return 'C';
  if (residue_name == " +C") return 'C';
  if (residue_name == "  G") return 'G';
  if (residue_name == " +G") return 'G';
  if (residue_name == "  I") return 'I';
  if (residue_name == " +I") return 'I';
  if (residue_name == "  T") return 'T';
  if (residue_name == " +T") return 'T';
  if (residue_name == "  U") return 'U';
  if (residue_name == " +U") return 'U';

  if (residue_name == " A ") return 'A';
  if (residue_name == "+A ") return 'A';
  if (residue_name == " C ") return 'C';
  if (residue_name == "+C ") return 'C';
  if (residue_name == " G ") return 'G';
  if (residue_name == "+G ") return 'G';
  if (residue_name == " I ") return 'I';
  if (residue_name == "+I ") return 'I';
  if (residue_name == " T ") return 'T';
  if (residue_name == "+T ") return 'T';
  if (residue_name == " U ") return 'U';
  if (residue_name == "+U ") return 'U';

  if (residue_name == "A  ") return 'A';
  if (residue_name == "C  ") return 'C';
  if (residue_name == "G  ") return 'G';
  if (residue_name == "I  ") return 'I';
  if (residue_name == "I  ") return 'I';
  if (residue_name == "T  ") return 'T';
  if (residue_name == "U  ") return 'U';

  // Adenosine                                  A
  // Modified adenosine                        +A
  // Cytidine                                   C
  // Modified cytidine                         +C
  // Guanosine                                  G
  // Modified guanosine                        +G
  // Inosine                                    I
  // Modified inosine                          +I
  // Thymidine                                  T
  // Modified thymidine                        +T
  // Uridine                                    U
  // Modified uridine                          +U
  // Unknown                                  UNK
  if (residue_name == "UNK") return 'X'; // explicitly uncertain
  return 'X'; // default
}

char three_letter_to_one_letter_protein_code(string residue_name) {
  // Protein residues
  if (residue_name == "ALA") return 'A';
  if (residue_name == "DAL") return 'A'; // D-alanine
  if (residue_name == "ASX") return 'B'; // ASN/ASP ambiguous
  if (residue_name == "CYS") return 'C';
  if (residue_name == "DCY") return 'C'; // D-cysteine
  if (residue_name == "CSX") return 'C'; // S-oxy cysteine
  if (residue_name == "CSW") return 'C'; // cysteine S-dioxide
  if (residue_name == "ASP") return 'D';
  if (residue_name == "DAS") return 'D'; // D-aspartic acid
  if (residue_name == "DSP") return 'D'; // D-aspartic acid
  if (residue_name == "PAS") return 'D'; // phosphorylated aspartate
  if (residue_name == "GLU") return 'E';
  if (residue_name == "GGL") return 'E'; // gamma glutamic acid
  if (residue_name == "DGL") return 'E'; // D-glutamic acid
  if (residue_name == "PHE") return 'F';
  if (residue_name == "DPN") return 'F'; // D-phenylalanine
  if (residue_name == "GLY") return 'G';
  if (residue_name == "HIS") return 'H';
  if (residue_name == "DHI") return 'H'; // D-histidine
  if (residue_name == "ILE") return 'I';
  if (residue_name == "DIL") return 'I'; // D-isoleucine
  if (residue_name == "LYS") return 'K';
  if (residue_name == "DLY") return 'K'; // D-lysine
  if (residue_name == "LEU") return 'L';
  if (residue_name == "DLE") return 'L'; // D-leucine
  if (residue_name == "MET") return 'M';
  if (residue_name == "CXM") return 'M'; // N-carboxymethionine
  if (residue_name == "FME") return 'M'; // N-formylmethionine
  if (residue_name == "MF3") return 'M'; // trifluoromethionine
  if (residue_name == "MHO") return 'M'; // S-oxymethionine
  if (residue_name == "MSE") return 'M'; // selenomethionine
  if (residue_name == "MSO") return 'M'; // selenomethionine selenoxide
  if (residue_name == "OMT") return 'M'; // S-dioxymethionine
  if (residue_name == "SME") return 'M'; // methionine sulfoxide
  if (residue_name == "ASN") return 'N';
  if (residue_name == "PRO") return 'P';
  if (residue_name == "DPR") return 'P'; // D-proline
  if (residue_name == "HYP") return 'P'; // hydroxyproline
  if (residue_name == "GLN") return 'Q';
  if (residue_name == "PCA") return 'Q'; // pyroglutamic acid (cylic blocked N-terminal Gln)
  if (residue_name == "5HP") return 'Q'; // pyroglutamic acid (cylic blocked N-terminal Gln)
  if (residue_name == "DGN") return 'Q'; // D-glutamine
  if (residue_name == "ARG") return 'R';
  if (residue_name == "DAR") return 'R'; // D-arginine
  if (residue_name == "SER") return 'S';
  if (residue_name == "DSN") return 'S'; // D-serine
  if (residue_name == "SEP") return 'S'; // phosphoserine
  if (residue_name == "THR") return 'T';
  if (residue_name == "DTH") return 'T'; // D-threonine
  if (residue_name == "TPO") return 'T'; // phosphothreonine
  if (residue_name == "CSE") return 'U'; // selenocysteine
  if (residue_name == "SOC") return 'U'; // dioxy-selenocysteine
  if (residue_name == "VAL") return 'V';
  if (residue_name == "DVA") return 'V'; // D-valine
  if (residue_name == "TRP") return 'W';
  if (residue_name == "TRP") return 'W'; // D-tryptophan
  if (residue_name == "LTR") return 'W'; // L-tryptophan
  if (residue_name == "TYR") return 'Y';
  if (residue_name == "DTY") return 'Y'; // D-tyrosine
  if (residue_name == "PTR") return 'Y'; // phosphotyrosine
  if (residue_name == "GLX") return 'Z'; // GLN/GLX ambiguous  
  // NAME                    CODE           FORMULA                 MOL. WT.
  // -----------------------------------------------------------------------------
  // Alanine                 ALA            C3 H7 N1 O2             89.09
  // Arginine                ARG            C6 H14 N4 O2            174.20
  // Asparagine              ASN            C4 H8 N2 O3             132.12
  // Aspartic acid           ASP            C4 H7 N1 O4             133.10
  // ASP/ASN ambiguous       ASX            C4 H71/2 N11/2 O31/2    132.61
  // Cysteine                CYS            C3 H7 N1 O2 S1          121.15
  // Glutamine               GLN            C5 H10 N2 O3            146.15
  // Glutamic acid           GLU            C5 H9 N1 O4             147.13
  // GLU/GLN ambiguous       GLX            C5 H91/2 N11/2 O31/2    146.64
  // Glycine                 GLY            C2 H5 N1 O2             75.07
  // Histidine               HIS            C6 H9 N3 O2             155.16
  // Isoleucine              ILE            C6 H13 N1 O2            131.17
  // Leucine                 LEU            C6 H13 N1 O2            131.17
  // Lysine                  LYS            C6 H14 N2 O2            146.19
  // Methionine              MET            C5 H11 N1 O2 S1         149.21
  // Phenylalanine           PHE            C9 H11 N1 O2            165.19
  // Proline                 PRO            C5 H9 N1 O2             115.13
  // Serine                  SER            C3 H7 N1 O3             105.09
  // Threonine               THR            C4 H9 N1 O3             119.12
  // Tryptophan              TRP            C11 H12 N2 O2           204.23
  // Tyrosine                TYR            C9 H11 N1 O3            181.19
  // Valine                  VAL            C5 H11 N1 O2            117.15
  // Undetermined            UNK            C5 H6 N1 O3             128.16
  if (residue_name == "UNK") return 'X'; // explicitly uncertain
  return 'X'; // default
}

bool three_letter_is_nucleotide(string residue_name) {
  if (three_letter_to_one_letter_nucleotide_code(residue_name) == 'X')
    return false;
  else return true;
}

bool three_letter_is_protein(string residue_name) {
  if (three_letter_to_one_letter_protein_code(residue_name) == 'X')
    return false;
  else return true;
}

char three_letter_to_one_letter_code(string residue_name) {
  if (three_letter_is_protein(residue_name))
    return three_letter_to_one_letter_protein_code(residue_name);
  else if (three_letter_is_nucleotide(residue_name))
    return three_letter_to_one_letter_nucleotide_code(residue_name);
  else return three_letter_to_one_letter_other_code(residue_name);
}

void PDBResidue::set_one_letter_code() {
  // Don't consult MODRES unless we don't know the original name
  // (I believe I know better...)
  const string & res1 = private_residue_name;
  const string & res2 = private_standard_residue_name;
  const char olc1 = three_letter_to_one_letter_code(res1);
  const char olc2 = three_letter_to_one_letter_code(res2);

  if (three_letter_is_protein(res1)) {
    private_is_protein = true;
    private_one_letter_code = olc1;
  }
  else if (three_letter_is_protein(res2)) {
    private_is_protein = true;
    private_one_letter_code = olc2;
  }
  else if ((olc1 != '\0') && (olc1 != 'X')) {
    private_is_protein = false;
    private_one_letter_code = olc1;
  }
  else if ((olc2 != '\0') && (olc2 != 'X')) {
    private_is_protein = false;
    private_one_letter_code = olc2;
  }
  else {
    private_one_letter_code = olc1;
  }
}

const char & PDBResidue::one_letter_code() const {
  return private_one_letter_code;
}

int PDBResidue::set_modres(const Modres & modres) {
  // If we already had a useful MODRES, ignore new ones...
  if ((get_standard_residue_name() != get_residue_name()) // Have earlier MODRES
      && (three_letter_to_one_letter_code(get_standard_residue_name()) != 'X')
      && (three_letter_to_one_letter_code(get_standard_residue_name()) != '\0')) {
    return 1;
  }
  private_standard_residue_name = modres.get_standard_residue_name();
  private_has_modres = true;
  set_one_letter_code();
  return 0;
}

void PDBResidue::set_terminus() {private_is_terminus = true;}

// TODO - make this work for more than just leucine
// TODO - incorporate other angle names
void PDBResidue::set_torsion(const string & torsion_name, const double angle) {
  
  double original_angle = get_torsion(torsion_name);
  double rotation_angle = angle - original_angle;

  vector<PDBAtom *> atoms_to_move;
  Vector3D rotation_axis;
  Vector3D initial_translation;

  if ((torsion_name == "CHI1") || 
      (torsion_name == "Chi1") ||
      (torsion_name == "chi1") ||
      (torsion_name == "x1") ||
      (torsion_name == "X1")
      ) {
    // TODO Change atoms CG, CD1, CD2, and hydrogens

    // To get origin centered rotation, need to translate to atom on axis
    initial_translation = - get_atom(" CB ").get_coordinate();
    rotation_axis = get_atom(" CB ").get_coordinate() - 
      get_atom(" CA ").get_coordinate();    

    atoms_to_move.push_back(& atom(" CG "));
    atoms_to_move.push_back(& atom(" CD1"));
    atoms_to_move.push_back(& atom(" CD2"));
  }

  if ((torsion_name == "CHI2") || 
      (torsion_name == "Chi2") ||
      (torsion_name == "chi2") ||
      (torsion_name == "x2") ||
      (torsion_name == "X2")
      ) {
    // TODO Change atoms CG, CD1, CD2, and hydrogens

    // To get origin centered rotation, need to translate to atom on axis
    initial_translation = - get_atom(" CG ").get_coordinate();
    rotation_axis = get_atom(" CG ").get_coordinate() - 
      get_atom(" CB ").get_coordinate();    

    atoms_to_move.push_back(& atom(" CD1"));
    atoms_to_move.push_back(& atom(" CD2"));
  }

  // Accumulate full transformation
  Matrix4D T1; T1.set_translation(initial_translation);
  Matrix4D R; R.set_rotation(axis_angle(rotation_axis.unit(), -rotation_angle));
  Matrix4D T2; T2.set_translation(-initial_translation);
  Matrix4D torsion = T2 * R * T1;
    
  vector<PDBAtom *>::iterator atom_ptr;
  for (atom_ptr = atoms_to_move.begin(); 
       atom_ptr != atoms_to_move.end();
       ++ atom_ptr) {
    (*atom_ptr)->coordinate() = torsion * (*atom_ptr)->coordinate();
  }
}

PDBResidue::PDBResidue(char chain_id, string residue_name, int residue_number, char insertion_code) :
	p_accessible_surface_area(-1)
{
	private_chain_id = chain_id;
	private_residue_name = residue_name;
	private_standard_residue_name = residue_name;
	private_residue_number = residue_number;
	private_insertion_code = insertion_code;
	private_unique_key = generate_unique_residue_key(residue_name,
													 residue_number,
													 insertion_code);
	private_is_hetatm_only = true; // Assume non-canonical residue
	private_has_modres = false;
	private_is_terminus = false;
	private_has_altloc = false;
	private_is_protein = false;
	private_is_redundant = false;
	
	private_one_letter_code = '>'; // Something ridiculous...
	
	set_one_letter_code();
}

ostream & operator<<(ostream & os, const PDBResidue & residue) {
  vector<PDBAtom>::const_iterator atom_ptr;
  for (atom_ptr = residue.get_atom().begin(); atom_ptr != residue.get_atom().end(); ++ atom_ptr) {
    os << *atom_ptr << endl;
  }
  return os;
}
