#include "Residue.h"

// static member string variables must be initialized in a .cpp file, not a header file

string AminoAcid::p_three_letter_code = "Xxx"; // Should this base class be used for "unknown?"
string GapResidue::p_three_letter_code = "---";

string       Alanine::p_three_letter_code = "Ala";
string Asparambiguous::p_three_letter_code = "Asx";
string      Cysteine::p_three_letter_code = "Cys";
string     Aspartate::p_three_letter_code = "Asp";
string     Glutamate::p_three_letter_code = "Glu";
string Phenylalanine::p_three_letter_code = "Phe";
string       Glycine::p_three_letter_code = "Gly";
string     Histidine::p_three_letter_code = "His";
string    Isoleucine::p_three_letter_code = "Ile";
string        Lysine::p_three_letter_code = "Lys";
string       Leucine::p_three_letter_code = "Leu";
string	  Methionine::p_three_letter_code = "Met";
string    Asparagine::p_three_letter_code = "Asn";
string       Proline::p_three_letter_code = "Pro";
string     Glutamine::p_three_letter_code = "Gln";
string      Arginine::p_three_letter_code = "Arg";
string        Serine::p_three_letter_code = "Ser";
string     Threonine::p_three_letter_code = "Thr";
string Selenocysteine::p_three_letter_code = "Sec";
string        Valine::p_three_letter_code = "Val";
string    Tryptophan::p_three_letter_code = "Trp";
string      Tyrosine::p_three_letter_code = "Tyr";
string Glutambiguous::p_three_letter_code = "Glx";

// Functions 

AminoAcid * new_protein_residue(const char one_letter_code) {
	switch (one_letter_code) {
		case 'A': return new Alanine;
		case 'B': return new Asparambiguous;
		case 'C': return new Cysteine;
		case 'D': return new Aspartate;
		case 'E': return new Glutamate;
		case 'F': return new Phenylalanine;
		case 'G': return new Glycine;
		case 'H': return new Histidine;
		case 'I': return new Isoleucine;
		case 'K': return new Lysine;
		case 'L': return new Leucine;
		case 'M': return new Methionine;
		case 'N': return new Asparagine;
		case 'P': return new Proline;
		case 'Q': return new Glutamine;
		case 'R': return new Arginine;
		case 'S': return new Serine;
		case 'T': return new Threonine;
		case 'U': return new Selenocysteine;
		case 'V': return new Valine;
		case 'W': return new Tryptophan;
		case 'Y': return new Tyrosine;
		case 'Z': return new Glutambiguous;
		default: return new AminoAcid; // unknown
	}
}


