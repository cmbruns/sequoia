#include "SequenceAlignment.h"
#include "MutationMatrix.h"

template class map<const char, float>;

// ********* Conservidue methods *************

Conservidue::Conservidue(const Residue & residue) :
  parent_alignment(NULL),
  weighted_sequence_count(residue.sequence_pointer()->get_weight()),
  weighted_internal_gap_count(0),
  weighted_left_end_gap_count(0),
  weighted_right_end_gap_count(0),
  weighted_left_end_sequence_count(0),
  weighted_right_end_sequence_count(0),
  weighted_internal_sequence_count(weighted_sequence_count),
  gap_opening_penalty(OPENING_SCORE, 0),
  gap_closing_penalty(CLOSING_SCORE, 0),
  gap_deletion_penalty(DELETION_SCORE, 0),
  p_gap_parameter(protein_gap_model, weighted_sequence_count),
  array_sequence_index(0),
  is_initial(false),
  is_final(false)
{
	const BioSequence * sequence = residue.sequence_pointer();

	// initialize_members();
	// weighted_sequence_count = sequence->weight;
	// weighted_internal_gap_count = 0.0;
	// weighted_left_end_gap_count = 0.0;
	// weighted_right_end_gap_count = 0.0;
	
	residue_counts[residue.one_letter_code()] = weighted_sequence_count;
	gap_opening_penalty = residue.p_gap_opening_penalty * weighted_sequence_count;
	gap_closing_penalty = residue.p_gap_closing_penalty * weighted_sequence_count;
	gap_deletion_penalty = residue.p_gap_deletion_penalty * weighted_sequence_count;
	
	int sequence_number = sequence->alignment_sequence_index;
	int residue_index = residue.sequence_residue_index;
	if (sequence_number == 0) { // create new sequence_residues mapping
		sequence_residues.clear();
		sequence_residues.push_back(residue_index);
	}
	else {
		sequence_residues.assign(sequence_number, -1); // default residue number -1 for other sequences
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
	
	answer.weighted_sequence_count += conservidue2.weighted_sequence_count;
	answer.weighted_internal_gap_count += conservidue2.weighted_internal_gap_count;
	answer.weighted_left_end_gap_count += conservidue2.weighted_left_end_gap_count;
	answer.weighted_right_end_gap_count += conservidue2.weighted_right_end_gap_count;
	answer.weighted_left_end_sequence_count += conservidue2.weighted_left_end_sequence_count;
	answer.weighted_right_end_sequence_count += conservidue2.weighted_right_end_sequence_count;
	answer.weighted_internal_sequence_count += conservidue2.weighted_internal_sequence_count;

	answer.gap_opening_penalty += conservidue2.gap_opening_penalty;
	answer.gap_closing_penalty += conservidue2.gap_closing_penalty;
	answer.p_gap_parameter = conservidue1.p_gap_parameter.combine(conservidue2.p_gap_parameter);
	answer.gap_deletion_penalty += conservidue2.gap_deletion_penalty;
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

void Conservidue::initialize_members() {
	gap_opening_penalty.opening() = 0;
	gap_closing_penalty.closing() = 0;
	gap_deletion_penalty.deletion() = 0;
	parent_alignment = NULL;
	weighted_sequence_count = 0;
	weighted_internal_gap_count = 0;
	weighted_left_end_gap_count = 0;
	weighted_right_end_gap_count = 0;
	weighted_left_end_sequence_count = 0;
	weighted_right_end_sequence_count = 0;
	weighted_internal_sequence_count = 0;
	array_sequence_index = 0;
	is_initial = false;
	is_final = false;
}

AlignmentScore Conservidue::gap_extension_penalty(unsigned int gseg) const { // -log2 probability of this residue adding to a loop
	return p_gap_parameter.extension_penalty(gseg);
}
AlignmentScore Conservidue::gap_open_offset(unsigned int gseg) const { // -log2 probability of this residue adding to a loop
	return p_gap_parameter.open_offset(gseg);
}
AlignmentScore Conservidue::gap_close_offset(unsigned int gseg) const { // -log2 probability of this residue adding to a loop
	return p_gap_parameter.close_offset(gseg);
}

ostream & Conservidue::print_debug(ostream & os, unsigned int indent_size) const {
	string indent = "";
	for(unsigned int i=0;i<indent_size;i++)indent += " ";

	os << indent << "conservidue pointer = " << this << endl;
	os << indent << "parent_alignment = " << parent_alignment << endl;
	os << indent << "array_sequence_index = " << array_sequence_index << endl;

	os << indent << "weighted_sequence_count = " << weighted_sequence_count << endl;
	os << indent << "weighted_internal_gap_count = " << weighted_internal_gap_count << endl;
	os << indent << "weighted_left_end_gap_count = " << weighted_left_end_gap_count << endl;
	os << indent << "weighted_right_end_gap_count = " << weighted_right_end_gap_count << endl;
	os << indent << "weighted_left_end_sequence_count = " << weighted_left_end_sequence_count << endl;
	os << indent << "weighted_right_end_sequence_count = " << weighted_right_end_sequence_count << endl;
	os << indent << "weighted_internal_sequence_count = " << weighted_internal_sequence_count << endl;

	os << indent << "gap_opening_penalty = " << gap_opening_penalty << endl;
	os << indent << "gap_closing_penalty = " << gap_closing_penalty << endl;
	os << indent << "gap_deletion_penalty = " << gap_deletion_penalty << endl;
	os << indent << "is_initial = " << is_initial << endl;
	os << indent << "is_final = " << is_final << endl;

	// os << indent << " gap parameter weight = " << p_gap_parameter.sequence_weight << endl;

	os << indent << "Residue counts:" << endl;
	map<char,float>::const_iterator count;
	for (count = residue_counts.begin(); count != residue_counts.end(); count++) {
		os << indent << "  " << count->first << ": " << count->second << endl;
	}
	os << indent << "Sequence residues:" << endl;
	for (unsigned int i = 0; i < sequence_residues.size(); i++) {
		os << indent << "  sequence " << i << ", residue " << sequence_residues[i] << endl;
	}
		
	return os;
}
	
Conservidue::Conservidue() : 
parent_alignment(NULL),
weighted_sequence_count(0.0),
weighted_internal_gap_count(0),
weighted_left_end_gap_count(0),
weighted_right_end_gap_count(0),
weighted_left_end_sequence_count(0),
weighted_right_end_sequence_count(0),
weighted_internal_sequence_count(0),
gap_opening_penalty(OPENING_SCORE, 0),
gap_closing_penalty(CLOSING_SCORE, 0),
gap_deletion_penalty(DELETION_SCORE, 0),
p_gap_parameter(protein_gap_model, weighted_sequence_count),
array_sequence_index(0),
is_initial(false),
is_final(false)
{
	// initialize_members();
} // Simple constructor	

// ********* SequenceAlignment methods *************

SequenceAlignment::SequenceAlignment(const BioSequence & seq0) :
left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR),
pair_alignment_score () 
{
	sequences.clear();
	add_sequence(seq0);
	BioSequence & seq = sequences.back(); // Use actual stored sequence copy from now on
	int sequence_index = sequences.size() - 1;
	weighted_sequence_count = seq0.get_weight();
	
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
		Conservidue * previous_conservidue_pointer = & conservidues[conservidues.size() - 2];
	
		if (i == 0) {
			current_conservidue.weighted_left_end_sequence_count = weighted_sequence_count;
			current_conservidue.weighted_internal_sequence_count = 0;
			current_conservidue.weighted_right_end_sequence_count = 0;
			
			// Make sure begin conservidue has the penalties it needs for the initializaion
			// portion of the dynamic programming table
			Conservidue * begin_conservidue_pointer = & conservidues.front();
			begin_conservidue_pointer->gap_opening_penalty = current_conservidue.gap_opening_penalty;
			begin_conservidue_pointer->gap_deletion_penalty = current_conservidue.gap_deletion_penalty;
			begin_conservidue_pointer->gap_closing_penalty = current_conservidue.gap_closing_penalty;
			begin_conservidue_pointer->weighted_sequence_count = current_conservidue.weighted_sequence_count;  // Needed in initialization
			for (unsigned int gseg = 0; gseg < current_conservidue.gap_segment_count(); gseg ++) {
				begin_conservidue_pointer->gap_open_offset(gseg) = current_conservidue.gap_open_offset(gseg); // piecewise penalty
				begin_conservidue_pointer->gap_close_offset(gseg) = current_conservidue.gap_close_offset(gseg); // piecewise penalty
			}
		}
		else {
			current_conservidue.weighted_left_end_sequence_count = 0;
			current_conservidue.weighted_internal_sequence_count = weighted_sequence_count;
			current_conservidue.weighted_right_end_sequence_count = 0;
		}
		
		// store sequence/residue relationship in conservidue
		current_conservidue.sequence_residues[sequence_index] = i;
		
		// Set up links for linear sequence
		ConserviduePredecessor back_link;
		back_link.predecessor_conservidue = i;
		back_link.transition_score.transition() = 0; // no penalty for following normal sequence
		current_conservidue.predecessors.clear();
		current_conservidue.predecessors.push_back(back_link); // Link to previous conservidue
		
		// Update for next round
		current_conservidue_pointer = & current_conservidue; // Will this pointer remain valid?		
	}
	
	current_conservidue_pointer->weighted_right_end_sequence_count = weighted_sequence_count;
	current_conservidue_pointer->weighted_internal_sequence_count = 0;
	
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
		back_link.transition_score.transition() = 0; // no penalty for following normal sequence		
		final_conservidue.predecessors.push_back(back_link); // Link to previous conservidue
		
		conservidues.push_back(final_conservidue); // Commit conservidue to array
		current_conservidue_pointer =  & conservidues.back(); // For storing sequential pointers
		ends.push_back(final_conservidue.array_sequence_index); // Note final conservidue of sequence
	}
	else { // Simply declare final actual conservidue final
	}
	
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
				
				// Insert blank every ten characters
				if ((column > 0) && (column % 10 == 0)) os << ' ';
				
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

unsigned int SequenceAlignment::length() const {return conservidues.size();}
const Conservidue & SequenceAlignment::operator[](int i) const {return conservidues.at(i);}

void SequenceAlignment::set_gap_penalty(double penalty) {
	// Distribute evenly among open, close, and delete
	double each_penalty = penalty / 3.0;
	for (unsigned int i = 0; i < conservidues.size(); i++) {
		conservidues[i].gap_opening_penalty.opening() = each_penalty;
		conservidues[i].gap_closing_penalty.closing() = each_penalty;
		conservidues[i].gap_deletion_penalty.deletion() = each_penalty;
	}
}
void SequenceAlignment::set_extension_penalty(const GapModel & gap_model) {
	// Distribute evenly among open, close, and delete
	for (unsigned int i = 0; i < conservidues.size(); i++) {
		ResidueGapParameter parameter(gap_model, conservidues[i].weighted_sequence_count);
		conservidues[i].p_gap_parameter = parameter;
	}
}
	
float SequenceAlignment::report_accuracy(const SequenceAlignment & true_alignment) {
	const SequenceAlignment & test_alignment = *this;
	// start by assuming both are pairwise alignments
		
	// index is sequence1 residue, value is sequence2 residue
	// first read matches from trusted alignment
	vector<int> true_matches(true_alignment.sequences[1].length(), -1);
	int true_match_count = 0;
	for (unsigned int i = 0; i < true_alignment.conservidues.size(); i++) {
		const Conservidue & conservidue = true_alignment.conservidues[i];
		int sequence1_residue_index = conservidue.sequence_residues[0];
		int sequence2_residue_index = conservidue.sequence_residues[1];

		if (sequence1_residue_index < 0) continue;
		if (sequence2_residue_index < 0) continue;

		const Residue & residue1 = true_alignment.sequences[0][sequence1_residue_index];
		const Residue & residue2 = true_alignment.sequences[1][sequence2_residue_index];

		if (residue1.is_gap()) continue;
		if (residue2.is_gap()) continue;

		int residue1_number = residue1.get_residue_number();
		int residue2_number = residue2.get_residue_number();

		// cout << residue1.one_letter_code() << residue2.one_letter_code() << endl;
		// cout << residue1_number << ", " << residue2_number << endl;
						
		true_matches[residue1_number] = residue2_number;

		true_match_count ++;
	}
	// cout << true_match_count << " true equivalent residue pairs found" << endl;

	// next read matches from test alignment
	int true_positive_count = 0;
	int false_positive_count = 0;
	for (unsigned int i = 0; i < test_alignment.conservidues.size(); i++) {
		const Conservidue & conservidue = test_alignment.conservidues[i];
		int sequence1_residue_index = conservidue.sequence_residues[0];
		int sequence2_residue_index = conservidue.sequence_residues[1];
			
		if (sequence1_residue_index < 0) continue;
		if (sequence2_residue_index < 0) continue;
			
		const Residue & residue1 = test_alignment.sequences[0][sequence1_residue_index];
		const Residue & residue2 = test_alignment.sequences[1][sequence2_residue_index];
			
		if (residue1.is_gap()) continue;
		if (residue2.is_gap()) continue;
			
		// cout << residue1.one_letter_code() << residue2.one_letter_code() << endl;
		int residue1_number = residue1.get_residue_number();
		int residue2_number = residue2.get_residue_number();
			
		if (true_matches[residue1_number] == residue2_number) {
			true_positive_count ++;
			// cout << residue1.one_letter_code() << residue2.one_letter_code() << endl;
			// cout << residue1_number << ", " << residue2_number << endl;
		}
		else {
			false_positive_count ++;
		}
		
	}
	// cout << true_positive_count << " aligned pairs correctly predicted" << endl;
	// cout << false_positive_count << " aligned pairs incorrectly predicted" << endl;

	float Q_accuracy = 100.0 * (float)true_positive_count / (float)true_match_count;
	return Q_accuracy;
}

	void SequenceAlignment::add_sequence(const BioSequence & sequence) { // does not automatically update conservidues!!
		if (sequences.size() == 0)
			weighted_sequence_count = 0.0;
		sequences.push_back(sequence);
		int sequence_index = sequences.size() - 1;
		sequences.back().alignment_sequence_index = sequence_index;
		sequences.back().parent_alignment = this;
		weighted_sequence_count += sequence.get_weight();
	}
	void SequenceAlignment::add_sequence_automatic(const BioSequence & sequence) {
		// Unlike add_sequence, updates conservidues as well
		// Requires that the sequence length is the same as the current Alignment length
		//  (the sequence may include gaps for this)
		SequenceAlignment & alignment = *this;

		if (alignment.sequences.size() == 0) { // No sequences have been added yet
			SequenceAlignment new_alignment(sequence);
			alignment = new_alignment;
			// cout << new_alignment;
		}
		else {
			if (alignment.length() != (sequence.length() + 1))
				throw ALIGNMENT_LENGTH_MISMATCH_EXCEPTION;
			SequenceAlignment single_sequence_alignment(sequence);
			SequenceAlignment new_alignment(alignment);
			// 1) Add sequence structure
			new_alignment.sequences.clear();
			// From starting alignment
			for (unsigned int i = 0; i < alignment.sequences.size(); i ++)
				new_alignment.add_sequence(alignment.sequences[i]);
			// And from new sequence
			new_alignment.add_sequence(sequence);
			// 2) Add conservidues
			new_alignment.conservidues.clear();
			Conservidue begin_conservidue;
			new_alignment.add_conservidue(begin_conservidue);  // Everyone needs a single begin
			for (unsigned int i = 1; i < alignment.length(); i ++) {
				new_alignment.add_conservidue(alignment[i].combine_conservidues(single_sequence_alignment[i]));
			}
			alignment = new_alignment;
		}
	}
	
void SequenceAlignment::add_conservidue(const Conservidue & c, unsigned int sequence_index_offset) {
		conservidues.push_back(c);
		Conservidue & conservidue = conservidues.back();
		int conservidue_index = conservidues.size() - 1;
		conservidue.array_sequence_index = conservidue_index;
		conservidue.parent_alignment = this;
						
		// TODO - I need a smarter predecessor method for POA alignments
		if (conservidue_index == 0) {
			conservidue.is_initial = true;

			begins.clear();
			begins.push_back(conservidue_index);
			ends.clear();
			ends.push_back(conservidue_index);
		}
		else {
			conservidue.is_initial = false;
			conservidues[conservidue_index - 1].is_final = false;
		}
		conservidue.is_final = true;
		ends[0] = conservidue_index;

		// Predecessor link
		ConserviduePredecessor pred;
		pred.transition_score.transition() = 0;
		pred.predecessor_conservidue = conservidue_index - 1;
		vector<ConserviduePredecessor> pred_vector(1,pred);
		conservidue.predecessors = pred_vector;

		// Sanity check sequence_residue mapping
		conservidue.sequence_residues.assign(sequences.size(), -1);
		for (unsigned int sr = 0; sr < c.sequence_residues.size(); sr ++) {
			if ((sr + sequence_index_offset) < conservidue.sequence_residues.size()) {
				conservidue.sequence_residues[sr + sequence_index_offset] = c.sequence_residues[sr];
			} else {
				throw 99; // something is wrong with sequence_residues logic
			}
		}
	}

	ostream & SequenceAlignment::print_debug(ostream & os, unsigned int indent_size) const {
		string indent = "";
		for(unsigned int i=0;i<indent_size;i++)indent += " ";
		
		os << indent << "alignment pointer = " << this << endl;
		os << indent << "left_gap_factor = " << left_gap_factor << endl;
		os << indent << "right_gap_factor = " << right_gap_factor << endl;
		os << indent << "left_gap_extension_factor = " << left_gap_extension_factor << endl;
		os << indent << "right_gap_extension_factor = " << right_gap_extension_factor << endl;

		os << indent << "pair_alignment_score = " << pair_alignment_score << endl;
		
		os << indent << "Conservidues:" << endl;
		for (unsigned int i = 0; i < conservidues.size(); i++) {
			conservidues[i].print_debug(os, indent_size + 2);
			os << endl;
		}
		os << indent << "Sequences:" << endl;
		for (unsigned int i = 0; i < sequences.size(); i++) {
			sequences[i].print_debug(os, indent_size + 2);
		}
		
		return os;
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

AlignmentScore SequenceAlignment::sequence_pair_score(unsigned int seq_index1, unsigned int seq_index2) const { // TODO
	// Compute extension score at gap close time
	// TODO - adapt methods for different alignment granularities
	
	AlignmentScore answer;

	// Note final residue in each sequence
	const Residue * final_residue1 = NULL;
	const Residue * final_residue2 = NULL;
	for (unsigned int i = 0; i < conservidues.size(); i++) {
		const Residue * res1 = conservidues[i].sequence_residue(seq_index1);
		const Residue * res2 = conservidues[i].sequence_residue(seq_index2);
		if (res1 != NULL) final_residue1 = res1;
		if (res2 != NULL) final_residue2 = res2;
	}
	
	int gap1 = 0; // Keep track of continuty of gaps in first sequence
	int gap2 = 0;
	const Residue * latest_residue1 = NULL; // Non-gap position in sequence1
	const Residue * latest_residue2 = NULL; // Non-gap position in sequence2
	const Residue * res1 = NULL;
	const Residue * res2 = NULL;
	const BioSequence & seq1 = sequences[seq_index1];
	const BioSequence & seq2 = sequences[seq_index2];
	for (unsigned int i = 0; i < conservidues.size(); i++) {
		res1 = conservidues[i].sequence_residue(seq_index1);
		res2 = conservidues[i].sequence_residue(seq_index2);
		
		if ((res1 == NULL) && (res2 == NULL)) continue; // Both are gaps, no effect on score

		// MATCH
		if ((res1 != NULL) && (res2 != NULL)) {
			answer.match() += blosum62.get_score(res1->one_letter_code(), res2->one_letter_code());

			double end_gap_factor1 = 1.0;
			double end_gap_extension_factor1 = 1.0;
			if (latest_residue1 == NULL) {// closing left gap in sequence1
				end_gap_factor1 = seq1.left_gap_factor;
				end_gap_extension_factor1 = seq1.left_gap_extension_factor;
			}
			double end_gap_factor2 = 1.0;
			double end_gap_extension_factor2 = 1.0;
			if (latest_residue2 == NULL) { // closing left gap in sequence1
				end_gap_factor2 = seq2.left_gap_factor;
				end_gap_extension_factor2 = seq2.left_gap_extension_factor;
			}
			
			// Gap closing penalty
			if (gap1 > 0) {
				answer.closing() += res2->p_gap_closing_penalty * end_gap_factor1;
				answer.extension() += seq2.p_gap_model->total_penalty(gap1) * end_gap_extension_factor1;
			}
			if (gap2 > 0) {
				answer.closing() += res1->p_gap_closing_penalty * end_gap_factor2;
				answer.extension() += seq1.p_gap_model->total_penalty(gap2) * end_gap_extension_factor2;
			}
			
			gap1 = 0;
			gap2 = 0;
		}

		// TODO - add gap penalties
		// TODO - extension penalties
		// TODO - open, close, delete penalties

		// Insert
		else if (res1 == NULL) { // Gap in sequence 1

			double end_gap_factor = 1.0;
			if (latest_residue1 == NULL) end_gap_factor = seq1.left_gap_factor;
			
			if (latest_residue1 == final_residue1) // final residue end gap
				end_gap_factor = seq1.right_gap_factor;
			
			// gap opening penalty
			if (gap1 == 0) {
				answer.opening() += res2->p_gap_opening_penalty * end_gap_factor;
				if (latest_residue1 != NULL) // not Initial end gap
					answer.deletion() += latest_residue1->p_gap_deletion_penalty * end_gap_factor; // TODO final end factor
				else // Initial end gap
					answer.deletion() += DEFAULT_GAP_DELETION_PENALTY * end_gap_factor; // TODO final end factor
			}

			gap1 ++;
			gap2 = 0;
		}
		
		// Delete
		else if (res2 == NULL) { // Gap in sequence 2

			double end_gap_factor = 1.0;
			if (latest_residue2 == NULL) end_gap_factor = seq2.left_gap_factor;
			
			if (latest_residue2 == final_residue2) // final residue end gap
				end_gap_factor = seq2.right_gap_factor;

			// gap opening penalty
			if (gap2 == 0) {
				answer.opening() += res1->p_gap_opening_penalty * end_gap_factor;
				if (latest_residue2 != NULL) // not Initial end gap
					answer.deletion() += latest_residue2->p_gap_deletion_penalty * end_gap_factor;
				else // Initial end gap
					answer.deletion() += DEFAULT_GAP_DELETION_PENALTY * seq2.left_gap_factor;
			}
			
			gap2 ++;
			gap1 = 0;
		}

		if (res1 != NULL) latest_residue1 = res1;
		if (res2 != NULL) latest_residue2 = res2;
		
	}
	
	// closing end gaps
	if ((gap1 > 0) && (latest_residue2 != NULL)) { // Terminal gap in seq1
		answer.closing() += latest_residue2->p_gap_closing_penalty * seq1.right_gap_factor;
		answer.extension() += seq2.p_gap_model->total_penalty(gap1) * seq1.right_gap_extension_factor;
	}
	if ((gap2 > 0) && (latest_residue1 != NULL)) { // Terminal gap in seq1
		answer.closing() += latest_residue1->p_gap_closing_penalty * seq2.right_gap_factor;
		answer.extension() += seq1.p_gap_model->total_penalty(gap2) * seq2.right_gap_extension_factor;
	}
	
	// Incorporate sequence weights
	double weight = sequences[seq_index1].get_weight() * sequences[seq_index2].get_weight();
	answer *= weight;
	
	return answer;
}

AlignmentScore SequenceAlignment::sum_of_pairs_score() const {
	AlignmentScore answer;
	for (unsigned int i = 0; i < sequences.size(); i++) {
		for (unsigned int j = i + 1; j < sequences.size(); j++) {
			AlignmentScore pair_score = sequence_pair_score(i,j);
			answer += pair_score;
			cout << "Pair score (" << i + 1 << ", " << j + 1 << ") = " << pair_score << endl;
			pair_score.print_details(cout);
		}
	}
	return answer;
}

SequenceAlignment::SequenceAlignment() : 
	left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
	right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
	left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
	right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR),
	weighted_sequence_count(0), 
	pair_alignment_score()
{} // simple constructor
// need assignment operator and copy constructor so that conservidue.parent_alignment pointers are properly set

SequenceAlignment::SequenceAlignment(const SequenceAlignment & alignment2) {
		*this = alignment2;
	}
	
SequenceAlignment & SequenceAlignment::operator=(const SequenceAlignment & alignment2) {
	if (this == &alignment2) return *this;
	begins = alignment2.begins;
	ends = alignment2.ends;
	left_gap_factor = alignment2.left_gap_factor;
	right_gap_factor = alignment2.right_gap_factor;
	left_gap_extension_factor = alignment2.left_gap_extension_factor;
	right_gap_extension_factor = alignment2.right_gap_extension_factor;
	pair_alignment_score = alignment2.pair_alignment_score;
	weighted_sequence_count = alignment2.weighted_sequence_count;

	sequences.clear();
	for (unsigned int i = 0; i < alignment2.sequences.size(); i++) {
		add_sequence(alignment2.sequences[i]);
	}

	conservidues.clear();
	for (unsigned int i = 0; i < alignment2.conservidues.size(); i++) {
		add_conservidue(alignment2.conservidues[i]);
	}
		
	return *this;
}

istream & operator>>(istream & is, SequenceAlignment & alignment) {
	alignment.load_fasta(is);
	return is;
}

// simple pretty alignment output routine
ostream & operator<<(ostream & os, const SequenceAlignment & s) {
	s.print_pretty(os);
	return os;
}



