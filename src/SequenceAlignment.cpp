#include "SequenceAlignment.h"

template class map<const char, float>;

// ********* Conservidue methods *************

Conservidue::Conservidue(const Residue & residue) {
	initialize_members();
	weighted_sequence_count = 1.0;  // 
	residue_counts[residue.one_letter_code()] = 1.0;

	residues.push_back(&residue);
}

const Residue * Conservidue::sequence_residue(const BioSequence * seq_ptr) const {
	map<const BioSequence *, const Residue *>::const_iterator it;
	if ((it = sequence_residues.find(seq_ptr)) != sequence_residues.end()) {
		const Residue * seq_ptr = it->second;
		return seq_ptr;
	}
	else return NULL; // sequence not found
}

// ********* SequenceAlignment methods *************

SequenceAlignment::SequenceAlignment(const BioSequence & seq) :
	left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
	right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
	left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
	right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR) {

	conservidues.clear(); // delete previous contents
	begins.clear();
	ends.clear();

	macromolecule = seq.macromolecule; // 
	sequences.push_back(&seq); // Store pointer to sequence

	// Add "begin" conservidue
	Conservidue begin_conservidue;
	begin_conservidue.array_sequence_index = 0;
	begin_conservidue.is_initial = true;
	begin_conservidue.is_final = false;
	begin_conservidue.parent_alignment = this;
	
	conservidues.push_back(begin_conservidue); // Commit conservidue to array
	Conservidue * current_conservidue_pointer =  & conservidues.back(); // For storing sequential pointers
	begins.push_back(current_conservidue_pointer);
	
	// Loop over residues in the sequence
	Conservidue * previous_conservidue_pointer =  & conservidues.back();
	int i;
	for (i = 0; i < seq.length(); ++i) {
		
		Conservidue current_conservidue(seq[i]); // initialize from single residue

		// store sequence/residue relationship in conservidue
		current_conservidue.sequence_residues[&seq] = &(seq[i]);

		current_conservidue.array_sequence_index = i + 1;
		current_conservidue.is_initial = false; // in sequence, begin is defined above
		current_conservidue.is_final = false;
		current_conservidue.parent_alignment = this;
		// current_conservidue.residue_counts[seq[i].one_letter_code()] = 1.0;
			
		// Set up links for linear sequence
		ConserviduePredecessor back_link;
		back_link.predecessor_conservidue = previous_conservidue_pointer;
		back_link.transition_score = 0; // no penalty for following normal sequence		
		current_conservidue.predecessors.push_back(back_link); // Link to previous conservidue
		
		conservidues.push_back(current_conservidue); // Commit conservidue to array

		// Update for next round
		current_conservidue_pointer = & conservidues.back(); // Will this pointer remain valid?		
		previous_conservidue_pointer = current_conservidue_pointer;
	}

	// Perhaps final empty conservidue is a bad idea
	if (0) { // Add full final fake conservidue
			 // Add "final" conservidue
		Conservidue final_conservidue;
		final_conservidue.array_sequence_index = seq.length() + 1;
		final_conservidue.is_initial = false;
		final_conservidue.is_final = true;
		final_conservidue.parent_alignment = this;

		// Set up links for linear sequence
		ConserviduePredecessor back_link;
		back_link.predecessor_conservidue = previous_conservidue_pointer;
		back_link.transition_score = 0; // no penalty for following normal sequence		
		final_conservidue.predecessors.push_back(back_link); // Link to previous conservidue
	
		conservidues.push_back(final_conservidue); // Commit conservidue to array
		current_conservidue_pointer =  & conservidues.back(); // For storing sequential pointers
		ends.push_back(current_conservidue_pointer); // Note final conservidue of sequence
	}
	else { // Simply declare final actual conservidue final
		current_conservidue_pointer->is_final = true;
		ends.push_back(current_conservidue_pointer);
	}
}

// simple string output routine
ostream & operator<<(ostream & os, const SequenceAlignment & s) {
	// TODO - loop over line length, sequences, and conservidues
	vector<const BioSequence *>::const_iterator seq; // Sequences in alignment
	for (seq = s.sequences.begin(); seq != s.sequences.end(); ++seq) {
		const BioSequence * seq_ptr = *seq;
		vector<Conservidue>::const_iterator cons; // Conservidues in alignment
		for (cons = s.conservidues.begin(); 
			 cons != s.conservidues.end(); 
			 ++cons) {
			const Conservidue & conservidue = * cons;

			if (conservidue.is_initial) continue;
			if (conservidue.is_final) continue;

			const Residue * res_ptr = conservidue.sequence_residue(seq_ptr);
			os << res_ptr->one_letter_code();
		}
		os << endl;  // end of one sequence
	}

	return os;
}


