/* Copyright (c) 2005 Christopher M. Bruns
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
 
// $Id$
//
// $Header$
//
// $Log$
// Revision 1.1  2004/06/04 19:34:46  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.7  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.6  2002/09/13 22:42:27  bruns
// Removed redundant parameter default
//
// Revision 1.5  2002/05/22 02:17:55  bruns
// Changed FastaSequence to BioSequence, to match new way
// Changed to_fasta() methods to use Residue(residue) constructors instead
// of Residue(char) constructors, and to carry residue_numbers along.
//
// Revision 1.4  2001/11/19 20:45:36  bruns
// Changed BaseSequence::add_residue calls to take character argument
//
// Revision 1.3  2001/11/16 01:07:59  bruns
// Fixed bug 5, so that BaseSequence stores an array of pointers to Residue
// objects allocated with new.
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "PDBChain.h"
#include "PDBAtom.h"
#include "BioSequence.h"
#include "ChemicalElement.h"
// #include "FastaResidue.h"

double PDBChain::store_accessible_surface_area() {
	PDBChain & chain = *this;
	
	double answer = 0.0;
	
	double probe_radius = OXYGEN.vdw_radius; // Angstroms
	double max_radius = (2.0 * probe_radius + 2.0 * CARBON.vdw_radius) / 2.0;
	
	// Hash all of the atoms in the structure
	Hash3D hash3d(max_radius);
	hash3d.add_pdb_chain(chain);
	
	// Now compute accessible surface for each residue
	for (unsigned int res = 0; res < chain.residue_count(); res++) {
		PDBResidue & residue = chain.residue(res);
		if (!residue.is_protein()) continue;
		double res_area = residue.store_accessible_surface_area(hash3d);
		answer += res_area;
	}

	p_accessible_surface_area = answer;
	return answer;
}

double PDBChain::compute_accessible_surface_area() const {
	const PDBChain & chain = *this;
	
	double answer = 0.0;
	
	double probe_radius = OXYGEN.vdw_radius; // Angstroms
	double max_radius = (2.0 * probe_radius + 2.0 * CARBON.vdw_radius) / 2.0;

	// Hash all of the atoms in the structure
	Hash3D hash3d(max_radius);
	hash3d.add_pdb_chain(chain);
	
	// Now compute accessible surface for each residue
	for (unsigned int res = 0; res < chain.residue_count(); res++) {
		const PDBResidue & residue = chain.get_residue(res);
		if (!residue.is_protein()) continue;
		double res_area = residue.compute_accessible_surface_area(hash3d);
		answer += res_area;
	}
	
	return answer;
}

PDBChain::PDBChain(char chain_id) {
  private_chain_id = chain_id;
}

void PDBChain::check_residues() {
  tweak_modres();

  vector<string>::const_iterator residue_id;
  residue_id = residue_order().begin();
  if (residue_id == residue_order().end()) return; // No residues to check
  PDBResidue * previous_residue = residue_by_key(*residue_id);
  // while we're at it, finalize one-letter code
  previous_residue->set_one_letter_code();
  previous_residue->check_redundancy(NULL);

  for (residue_id = residue_order().begin();
       residue_id != residue_order().end(); 
       residue_id ++) {
    PDBResidue * residue = residue_by_key(*residue_id);
    if (residue == NULL) continue;
    // while we're at it, finalize one-letter code
    residue->set_one_letter_code();
    residue->check_redundancy(&(*previous_residue));
    previous_residue = residue;
  }
}

// Use MODRES mapping on residues that are still unknown
// (without regard to residue number)

void PDBChain::tweak_modres() {
  vector<Modres>::const_iterator modres_ptr;
  for (modres_ptr = private_modres.begin();
       modres_ptr != private_modres.end();
       modres_ptr ++) {
  
    // Adjust residues that have no explicit MODRES, but match name
    // First examine coordinate residues
    map<string,PDBResidue>::iterator residue_ptr;
    for (residue_ptr = private_residue.begin(); 
	 residue_ptr != private_residue.end(); 
	 residue_ptr ++) {
      PDBResidue & residue = (*residue_ptr).second;
      // Already got MODRES?
      if (residue.get_residue_name() != residue.get_standard_residue_name()) continue;
      if (residue.get_residue_name() != modres_ptr->get_residue_name()) continue;
      int set_status = residue.set_modres(*modres_ptr);
      if (set_status == 0) private_used_modres.push_back(*modres_ptr);
    }
  }
}

void PDBChain::set_modres(const Modres & modres) {
  private_modres.push_back(modres);
  PDBResidue * residue = this->residue(modres.get_residue_name(),
				       modres.get_sequence_number(),
				       modres.get_insertion_code());
  if (residue == NULL) return;

  int set_status = residue->set_modres(modres);
  if (set_status == 0) private_used_modres.push_back(modres);
}

const char & PDBChain::get_chain_id() const {
  return private_chain_id;
}

bool PDBChain::is_protein(string residue_name) {
  int protein_residue_count = 0;
  int total_residue_count = 0;
  int nonx_residue_count = 0;

  // First examine coordinate residues
  vector<string>::const_iterator residue_id;
  for (residue_id = residue_order().begin(); 
       residue_id != residue_order().end(); 
       residue_id ++) {
    const PDBResidue & residue = get_residue_by_key(*residue_id);
    if (residue.one_letter_code() != '\0') total_residue_count ++;
    if ((residue.one_letter_code() != 'X')
	&& (residue.one_letter_code() != '\0')) nonx_residue_count ++;
    if (residue.is_protein()) protein_residue_count ++;
  }

  // Next examine SEQRES residues
  for (residue_id = get_seqres().get_sequence().begin(); 
       residue_id != get_seqres().get_sequence().end(); 
       residue_id ++) {
    if (three_letter_to_one_letter_code(*residue_id) == '\0') continue;
    total_residue_count ++;
    if (three_letter_to_one_letter_code(*residue_id) == 'X') continue;
    nonx_residue_count ++;
    if (three_letter_is_protein(*residue_id)) 
      protein_residue_count ++;
  }
  
  if ((protein_residue_count > 0)
      && ((((double)protein_residue_count) / ((double)nonx_residue_count)) > 0.95)
      && ((((double)protein_residue_count) / ((double)total_residue_count)) > 0.50)
      )
    return true;
  else return false;
}

void PDBChain::set_seqres(const Seqres & seqres) {
  private_seqres = seqres;
}

const vector<string> & PDBChain::residue_order() const {return private_residue_order;}

void PDBChain::add_atom(const PDBAtom & atom) {
  PDBResidue * res = this->residue(atom.get_residue_name(),
				   atom.get_residue_number(),
				   atom.get_insertion_code());

  // Generate a new residue if this one does not exist
  if (res == NULL) res = add_residue(atom.get_residue_name(),
				     atom.get_residue_number(),
				     atom.get_insertion_code());
  if (res == NULL) abort();
  res->add_atom(atom);
}

PDBResidue * PDBChain::add_residue(const string & res_name,
				   const int & res_num,
				   const char & i_code) {
  PDBResidue new_residue(get_chain_id(),
			 res_name, res_num, i_code);
  string residue_key = generate_unique_residue_key(res_name, res_num, i_code);
  private_residue[residue_key] = new_residue;
  private_residue_order.push_back(residue_key);
  return this->residue(res_name, res_num, i_code);
}

void PDBChain::add_terminus(const string & res_name,
			    const int & res_num,
			    const char & i_code) {
  PDBResidue * res = this->residue(res_name, res_num, i_code);
  if (res == NULL) res = add_residue(res_name, res_num, i_code);
  res->set_terminus();
}

PDBResidue * PDBChain::residue(const string & res_name,
			       const int & res_num,
			       const char & i_code) {
  string residue_id = generate_unique_residue_key(res_name, res_num, i_code);
  if (private_residue.find(residue_id) == private_residue.end()) {
    return NULL;
  }
  return &((*(private_residue.find(residue_id))).second);
}

// const map<string, PDBResidue> & PDBChain::residue() const {return private_residue;}

PDBResidue * PDBChain::residue_by_key(const string & key) {
  return &((*(private_residue.find(key))).second);
}

const PDBResidue & PDBChain::get_residue(const string & res_name,
					 const int & res_num,
					 const char & i_code) const {
  string residue_id = generate_unique_residue_key(res_name, res_num, i_code);
  return (*(private_residue.find(residue_id))).second;
}
const PDBResidue & PDBChain::get_residue_by_key(const string & key) const {
  return (*(private_residue.find(key))).second;
}

// Make a fasta entry, iff possible
BioSequence PDBChain::to_fasta(const string & id, const string & title) const {
  BioSequence coord_fasta = to_coord_fasta(id, title);
  BioSequence seqres_fasta = to_seqres_fasta(id, title);
  BioSequence best_fasta;

  // If there are no coordinates, then there is nothing to return
  if (coord_fasta.length() < 1) return coord_fasta;

  // FIXME - make upper-case for residues with coordinates
  // Use SEQRES for non-pathological cases
  if (seqres_fasta.is_empty() && !coord_fasta.is_empty()) {
    best_fasta = coord_fasta;
    best_fasta.set_title("match:coord_only\t" + best_fasta.get_title());
  }
  else if (coord_fasta.is_identical_to(seqres_fasta)) {
    best_fasta = seqres_fasta;
    best_fasta.set_title("match:exact\t" + best_fasta.get_title());
  }
  else if (coord_fasta.is_substring_of(seqres_fasta)) {
    best_fasta = seqres_fasta;
    best_fasta.set_title("match:substring\t" + best_fasta.get_title());
  }
  else if (coord_fasta.is_subsequence_of(seqres_fasta)) {
    best_fasta = seqres_fasta;
    best_fasta.set_title("match:subsequence\t" + best_fasta.get_title());
  }
  else if (coord_fasta.is_erroneous_subsequence_of(seqres_fasta, 2)) {
    best_fasta = seqres_fasta;
    best_fasta.set_title("match:1_or_2_mismatch\t" + best_fasta.get_title());
  }
  else {
    best_fasta = coord_fasta;
    best_fasta.set_title("match:mismatch\t" + best_fasta.get_title());
  }

  return best_fasta;
}

BioSequence PDBChain::to_coord_fasta() const {
	// Version to populate id and title automatically
	BioSequence answer = to_coord_fasta("This is not the ID!!!", "This is not the description!!!");
	// TODO - populate id and title
	answer.set_id(get_sequence_id());
	answer.set_title(get_title());
	return answer;
}

// Make a fasta entry using sequence deduced from the atomic coordinates 
// (as opposed to SEQRES records)
BioSequence PDBChain::to_coord_fasta(const string & id, const string & title) const {

  BioSequence fasta;

  fasta.set_id(id);
  fasta.set_title(title);
  
  // Loop over residues
  vector<string>::const_iterator residue_id;
  int n_good_residue = 0;
  for (residue_id = residue_order().begin(); 
       residue_id != residue_order().end(); 
       residue_id ++) {
    const PDBResidue & residue = get_residue_by_key(*residue_id);

    char olc = residue.one_letter_code();
    if (olc == '\0') continue; // skip null characters - this is my code for known non-sequence

    // skip HETATM without one letter code or SEQRES counterpart - probable non-sequence residue
    if (residue.is_hetatm_only() 
	&& (residue.one_letter_code() == 'X')
	&& (private_seqres.get_residue_count(residue.get_residue_name()) == 0) )
      continue;

    // Skip residues that know themselves to be redundant - see PDBResidue::check_redundancy()
    if (residue.is_redundant()) {
      // But still notice if they are terminal
      if (residue.is_terminus() && (n_good_residue > 0)) break;
      continue;
    }
    
    if (olc != 'X') n_good_residue ++;
    fasta.add_residue(residue);
    
    // Stop at "TER" records - but not before we got some residues
    if (residue.is_terminus() && (n_good_residue > 0)) break;
  }
  
  // Only return something if we really got some residues
  if (n_good_residue > 0)
    return fasta;
  else {
    BioSequence empty;
    return empty;
  }
}

// Make a fasta entry, iff possible
BioSequence PDBChain::to_seqres_fasta(const string & id, const string & title) const {
  BioSequence fasta;

  fasta.set_id(id);
  fasta.set_title(title);
  
  // Loop over seqres
  vector<string>::const_iterator residue_id;
  int line_pos = 0;
  int n_good_residue = 0;
  int n_letter_residue = 0;
  for (residue_id = private_seqres.get_sequence().begin(); 
       residue_id != private_seqres.get_sequence().end(); 
       residue_id ++) {
    
    char olc = three_letter_to_one_letter_code(*residue_id);

    // incorporate MODRES changes
    // FIXME - start with ones that were used for coordinates
    if ((olc == 'X') || (olc == '\0')) {
      vector<Modres>::const_iterator modres_ptr;
      for (modres_ptr = private_used_modres.begin();
	   modres_ptr != private_used_modres.end();
	   modres_ptr ++) {
	char olc2;
	if (modres_ptr->get_residue_name() == *residue_id) {
	  olc2 = three_letter_to_one_letter_code(modres_ptr->get_standard_residue_name());
	  if (olc2 != '\0') {
	    olc = olc2;
	    if (olc2 != 'X')
	      break; // Stop at first good one
	  }
	}
      }
    }
    if ((olc == 'X') || (olc == '\0')) {
      vector<Modres>::const_iterator modres_ptr;
      for (modres_ptr = private_modres.begin();
	   modres_ptr != private_modres.end();
	   modres_ptr ++) {
	char olc2;
	if (modres_ptr->get_residue_name() == *residue_id) {
	  olc2 = three_letter_to_one_letter_code(modres_ptr->get_standard_residue_name());
	  if (olc2 != '\0') {
	    olc = olc2;
	    if (olc2 != 'X')
	      break; // Stop at first good one
	  }
	}
      }
    }

    if (olc == '\0') continue; // skip null characters

    if (olc != 'X') n_good_residue ++;
    n_letter_residue ++;
    

	Residue * fasta_residue = new_protein_residue(olc);
	fasta_residue->set_residue_number(n_letter_residue);
    // FastaResidue fasta_residue(olc, n_letter_residue);
    fasta.add_residue(*fasta_residue);
	delete fasta_residue;
    
    line_pos ++;
  }
  
  // Only return something if we really got some residues
  // if (n_letter_residue > 0)
  if (n_good_residue > 0)
    return fasta;
  else {
    BioSequence empty;
    return empty;
  }
}

ostream & operator<<(ostream & os, const PDBChain & chain) {
  vector<string>::const_iterator residue_id;
  for (residue_id = chain.residue_order().begin(); 
       residue_id != chain.residue_order().end(); 
       ++ residue_id) {
    const PDBResidue & residue = chain.get_residue_by_key(*residue_id);
    os << residue;
  }
  return os;
}
