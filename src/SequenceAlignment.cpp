#include "SequenceAlignment.h"

template class map<const char, float>;

istream & operator>>(istream & is, SequenceAlignment & alignment) {
	alignment.load_fasta(is);
	return is;
}

// Load sequence alignment from fasta file
istream & SequenceAlignment::load_fasta(istream & is) {
	SequenceAlignment & alignment = *this;
	// Read a series of sequences, and paste them together into an alignment
	while (is.good()) {
		BioSequence sequence;
		is >> sequence;
		alignment.add_sequence_automatic(sequence);
	}
	return is;
}

// ********* Conservidue methods *************

Conservidue::Conservidue(const Residue & residue) :
	p_gap_parameter(protein_gap_model) 
{
	initialize_members();
	weighted_sequence_count = 1.0;  // 
	residue_counts[residue.one_letter_code()] = 1.0;

	const BioSequence * sequence = residue.sequence_pointer();
	int sequence_number = sequence->alignment_sequence_index;
	int residue_index = residue.sequence_residue_index;
	if (sequence_number == 0) {
		sequence_residues.clear();
		sequence_residues.push_back(residue_index);
	}
	else {
		sequence_residues.assign(sequence_number, -1);
		sequence_residues[sequence_number] = residue_index;
	}
}

const Residue * Conservidue::sequence_residue(unsigned int sequence_index) const {
	int residue_index = sequence_residues[sequence_index];
	if (residue_index < 0) return NULL; // sequence not found
	const BioSequence & seq = parent_alignment->sequences[sequence_index];
	const Residue * res_ptr = & seq[residue_index];
	return res_ptr;
}

// Create a combined conservidue from the alignment of two matching conservidues
Conservidue Conservidue::combine_conservidues(const Conservidue & conservidue2) const {
	const Conservidue & conservidue1 = *this;
	Conservidue answer = conservidue1;
	
	float weight1 = conservidue1.weighted_sequence_count;
	float weight2 = conservidue2.weighted_sequence_count;
	float total_weight = weight1 + weight2;
	weight1 = weight1 / total_weight;
	weight2 = weight2 / total_weight;

	answer.weighted_sequence_count += conservidue2.weighted_sequence_count;
	answer.gap_opening_penalty += conservidue2.gap_opening_penalty;
	answer.gap_closing_penalty += conservidue2.gap_closing_penalty;
	answer.p_gap_parameter = conservidue1.p_gap_parameter.combine(conservidue2.p_gap_parameter);
	answer.gap_deletion_penalty += conservidue2.gap_deletion_penalty;
	answer.gap_opening_penalty += conservidue2.gap_opening_penalty;
	answer.parent_alignment = NULL;
	
	// Fold in sequence to residues hash
	for (unsigned int i = 0; i < conservidue2.sequence_residues.size(); i++) {
		answer.sequence_residues.push_back(conservidue2.sequence_residues[i]);
	}

	// Fold in residue type counts hash
	map<char, float>::const_iterator res_count2;
	for (res_count2 = conservidue2.residue_counts.begin();
		 res_count2 != conservidue2.residue_counts.end();
		 res_count2 ++) {
		answer.residue_counts[res_count2->first] += res_count2->second;
	}
	
	return answer;
}

// ********* SequenceAlignment methods *************

SequenceAlignment::SequenceAlignment(const BioSequence & seq0) :
	left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
	right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
	left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
	right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR),
	pair_alignment_score (0) 
{
	sequences.clear();
	add_sequence(seq0);
	BioSequence & seq = sequences.back(); // Use actual stored sequence copy from now on
	int sequence_index = sequences.size() - 1;
		
	conservidues.clear(); // delete previous contents
	begins.clear();
	ends.clear();

	// Add "begin" conservidue - does not actually contain residues (but end does!)
	Conservidue begin_conservidue;
	add_conservidue(begin_conservidue);
	Conservidue * current_conservidue_pointer =  & conservidues.back(); // For storing sequential pointers
	
	// Loop over residues in the sequence
	for (unsigned int i = 0; i < seq.length(); ++i) {
		const Residue & residue = seq[i];
		Conservidue current_conservidue0(residue); // initialize from single residue
		add_conservidue(current_conservidue0);
		Conservidue & current_conservidue = conservidues.back();
		
		// store sequence/residue relationship in conservidue
		current_conservidue.sequence_residues[sequence_index] = i;

		// Set up links for linear sequence
		ConserviduePredecessor back_link;
		back_link.predecessor_conservidue = i;
		back_link.transition_score = 0; // no penalty for following normal sequence
		current_conservidue.predecessors.clear();
		current_conservidue.predecessors.push_back(back_link); // Link to previous conservidue
		
		// Update for next round
		current_conservidue_pointer = & current_conservidue; // Will this pointer remain valid?		
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
		back_link.predecessor_conservidue = seq.length();
		back_link.transition_score = 0; // no penalty for following normal sequence		
		final_conservidue.predecessors.push_back(back_link); // Link to previous conservidue
	
		conservidues.push_back(final_conservidue); // Commit conservidue to array
		current_conservidue_pointer =  & conservidues.back(); // For storing sequential pointers
		ends.push_back(final_conservidue.array_sequence_index); // Note final conservidue of sequence
	}
	else { // Simply declare final actual conservidue final
	}
}

// simple pretty alignment output routine
ostream & operator<<(ostream & os, const SequenceAlignment & s) {
	s.print_pretty(os);
	return os;
}

ostream & SequenceAlignment::print_pretty(ostream & os) const {
	const SequenceAlignment & s = *this;
	// Pretty alignment output

	unsigned int line_length = 60;
	unsigned int conservidue_number = 0;
	unsigned int sequence_number = 0;

	os << "Alignment score = " << s.pair_alignment_score << endl;
	
	// Print sequence summary
	for (sequence_number = 0; sequence_number < s.sequences.size(); ++sequence_number) {
		const BioSequence & seq = s.sequences[sequence_number];
		os << sequence_number + 1 << ":\t";
		os << seq.get_id() << "\t";
		os << seq.get_title() << endl;
	}
	os << endl;
	
	// Keep track of current residue number
	vector<int> sequence_residue_numbers(s.sequences.size(), 0);

	// loop over line length, sequences, and conservidues
	while (conservidue_number < s.conservidues.size()) {  // number at start of block
		unsigned int block_conservidue_number = conservidue_number;
		// final fake "sequence" indicates summary row for each block
		for (sequence_number = 0; sequence_number <= s.sequences.size(); ++sequence_number) {

			if (sequence_number == s.sequences.size())  // Summary row
			  os << "  " << " \t";
			else os << sequence_number + 1 << ":\t";
							
			const Residue * residue;
			block_conservidue_number = conservidue_number;
			for (unsigned int column = 0; column < line_length; column ++) {
				// Skip "begin" conservidues
				while ((block_conservidue_number < s.conservidues.size()) &&
					   (s.conservidues[block_conservidue_number].is_initial)) {
					block_conservidue_number ++;
				}
				if ((block_conservidue_number) >= s.conservidues.size()) break;
				const Conservidue & conservidue = s.conservidues[block_conservidue_number];

				if (sequence_number == s.sequences.size()) { // Summary row
					if ((conservidue.residue_counts.size() == 1) &&
						(conservidue.residue_counts.begin()->second > 1.0))
						os << '*'; // identical column
					else os << ' '; // non-identical column
				}
				else { // normal sequence row, not summary
					residue = conservidue.sequence_residue(sequence_number); // doesn't work
					if (residue == NULL)  os << '-';
					else {
						os << *residue;
						if (! residue->is_gap()) // Update sequence number for end of line
							sequence_residue_numbers[sequence_number] = residue->get_residue_number();
					}
				}

				block_conservidue_number ++;
			}
			
			// show residue number
			if (sequence_number == s.sequences.size()) {} // Summary row
			else os << "  " << sequence_residue_numbers[sequence_number];

			// End of one sequence line in alignment
			os << endl;
			
		}

		os << endl; // break between sections

		conservidue_number = block_conservidue_number;
	}
	
	return os;
}


