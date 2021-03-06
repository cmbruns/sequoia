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

// 

// $Id$
// $Header$
// $Log$
// Revision 1.6  2004/06/14 16:50:53  cmbruns
// Renamed Score.h to AlignmentScore.h
//
// Revision 1.5  2004/06/04 19:13:45  cmbruns
// Updated GPL header
//
// Moved initialization of one_letter_code values from header, here, to cpp file, so it will compile on baxter (Linux)
//
// Made residue_number accessors virtual, so that the structure residue numbers will work correctly
//

#ifndef __RESIDUE_H__
#define __RESIDUE_H__

#include <string>
#include <iostream>
#include "AlignmentScore.h"

#define DEFAULT_GAP_OPENING_PENALTY -0.34
#define DEFAULT_GAP_CLOSING_PENALTY -0.35
#define DEFAULT_GAP_DELETION_PENALTY -0.36

using namespace std;

// Residue.h
// Classes describing platonic amino acids or nucleotides

enum Macromolecule {
	PROTEIN_MOLECULE,
	NUCLEIC_MOLECULE,
	SPLICEDRNA_MOLECULE,
	GENOMIC_MOLECULE,
	PROTEINNUCLEIC_MOLECULE // Mixed hybrid alignment
};

class BioSequence;

///////// Base class ////////////
class Residue {
	friend class SequenceAlignment;
	friend class Conservidue;
protected:
	int p_residue_number;  // Only one true instance variable (non static)
	const static char p_one_letter_code; // should not be accessed, every subclass should override this...
	Macromolecule p_macromolecule;
	const BioSequence * p_sequence_pointer;
	
	// Alignment specific parameters
	AlignmentScore p_gap_deletion_penalty; // -log2 probability of gap after this Conservidue
	AlignmentScore p_gap_closing_penalty; // -log2 probability of loop before this Conservidue
	AlignmentScore p_gap_opening_penalty; // -log2 probability of loop beginning with this Conservidue
public:
	unsigned int sequence_residue_index; // Actual index in parent sequence

	// Need accessor functions for variables so "virtual" functionality will work
    virtual const char & one_letter_code() const {return p_one_letter_code;}

	const BioSequence * sequence_pointer() const {return p_sequence_pointer;}
	virtual void set_sequence_pointer(BioSequence * p) {p_sequence_pointer = p;}
	int & residue_number() {return p_residue_number;}
	virtual int get_residue_number() const {return p_residue_number;}
	Macromolecule & macromolecule() {return p_macromolecule;}
		
	virtual bool is_gap() const {return false;} // only GapResidue can be a gap
	virtual Residue * new_clone() const = 0; // so subclasses can send the correct pointer type

	ostream & print_debug(ostream & os, unsigned int indent_size = 0) const {
		string indent = "";
		for(unsigned int i=0;i<indent_size;i++)indent += " ";
		
		os << indent << "residue pointer = " << this << endl;
		os << indent << "residue number = " << p_residue_number << endl;
		os << indent << "one letter code = " << one_letter_code() << endl;
			
		return os;
	}
	void set_residue_number(int n) {p_residue_number = n;}

	Residue() 
		: 
		p_residue_number(-1),
		p_sequence_pointer(NULL), 
		p_gap_deletion_penalty(DELETION_SCORE, DEFAULT_GAP_DELETION_PENALTY),
		p_gap_closing_penalty(CLOSING_SCORE, DEFAULT_GAP_CLOSING_PENALTY),
		p_gap_opening_penalty(OPENING_SCORE, DEFAULT_GAP_OPENING_PENALTY),
		sequence_residue_index(0)		
	{}
};

////////// Gap instance class //////////////
class GapResidue : public Residue {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code;
public:
	GapResidue() {
		p_three_letter_code = "---";
	}
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	
	GapResidue * new_clone() const {return new GapResidue(*this);}
	bool is_gap() const {return true;} // override default 'false' value
};

////////// Intermediate classes for protein and nucleotides //////////
class Nucleotide : public Residue {	
protected:
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	virtual const char & one_letter_code() const {return p_one_letter_code;}

	Nucleotide() {
		macromolecule() = NUCLEIC_MOLECULE;
	}
	Nucleotide * new_clone() const {return new Nucleotide(*this);}
	bool is_gap() const {return false;}
};

// Use this class for "unknown amino acid"
class AminoAcid : public Residue {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	virtual const char & one_letter_code() const {return p_one_letter_code;}
	virtual const string & three_letter_code() const {return p_three_letter_code;}

	AminoAcid() {
		macromolecule() = PROTEIN_MOLECULE;
	}
	AminoAcid * new_clone() const {return new AminoAcid(*this);}
	bool is_gap() const {return false;}
};

// construct an object of the correct class
Residue * new_protein_residue(const char one_letter_code);

//////////// Instance classes for DNA /////////////////////////
// string members must be defined in .cpp file //

class Adenylate : public Nucleotide {	
protected:
	const static char p_one_letter_code; // default to ambiguity character
public:
	Adenylate * new_clone() const {return new Adenylate(*this);}
	const static string name;
};

class Cytidylate : public Nucleotide {	
protected:
	const static char p_one_letter_code; // default to ambiguity character
public:
	Cytidylate * new_clone() const {return new Cytidylate(*this);}
	const static string name;
};

class Guanylate : public Nucleotide {	
protected:
	const static char p_one_letter_code; // default to ambiguity character
public:
	Guanylate * new_clone() const {return new Guanylate(*this);}
	const static string name;
};

class Thymidylate : public Nucleotide {	
protected:
	const static char p_one_letter_code; // default to ambiguity character
public:
	Thymidylate * new_clone() const {return new Thymidylate(*this);}
	const static string name;
};

//////////// Instance classes for protein /////////////////////////
// string members must be defined in .cpp file //

class Alanine : public AminoAcid {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}	
	Alanine * new_clone() const {return new Alanine(*this);}
	const static string name;
};

class Asparambiguous : public AminoAcid {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Asparambiguous * new_clone() const {return new Asparambiguous(*this);}
	const static string name;
};

class Cysteine : public AminoAcid {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Cysteine * new_clone() const {return new Cysteine(*this);}
	const static string name;
};

class Aspartate : public AminoAcid {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Aspartate * new_clone() const {return new Aspartate(*this);}
	const static string name;
};

class Glutamate : public AminoAcid {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Glutamate * new_clone() const {return new Glutamate(*this);}
	const static string name;
};

class Phenylalanine : public AminoAcid {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Phenylalanine * new_clone() const {return new Phenylalanine(*this);}
	const static string name;
};

class Glycine : public AminoAcid {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Glycine * new_clone() const {return new Glycine(*this);}
	const static string name;
};

class Histidine : public AminoAcid {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Histidine * new_clone() const {return new Histidine(*this);}
	const static string name;
};

class Isoleucine : public AminoAcid {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Isoleucine * new_clone() const {return new Isoleucine(*this);}
	const static string name;
};

class Lysine : public AminoAcid {	
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Lysine * new_clone() const {return new Lysine(*this);}
	const static string name;
};

class Leucine : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Leucine * new_clone() const {return new Leucine(*this);}
	const static string name;
};

class Methionine : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Methionine * new_clone() const {return new Methionine(*this);}
	const static string name;
};

class Asparagine : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	
public:
	Asparagine * new_clone() const {return new Asparagine(*this);}
	const static string name;
};

class Proline : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Proline * new_clone() const {return new Proline(*this);}
	const static string name;
};

class Glutamine : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Glutamine * new_clone() const {return new Glutamine(*this);}
	const static string name;
};

class Arginine : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Arginine * new_clone() const {return new Arginine(*this);}
	const static string name;
};

class Serine : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Serine * new_clone() const {return new Serine(*this);}
	const static string name;
};

class Threonine : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Threonine * new_clone() const {return new Threonine(*this);}
	const static string name;
};

class Selenocysteine : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Selenocysteine * new_clone() const {return new Selenocysteine(*this);}
	const static string name;
};

class Valine : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Valine * new_clone() const {return new Valine(*this);}
	const static string name;
};

class Tryptophan : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Tryptophan * new_clone() const {return new Tryptophan(*this);}
	const static string name;
};

class Tyrosine : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Tyrosine * new_clone() const {return new Tyrosine(*this);}
	const static string name;
};

class Glutambiguous : public AminoAcid {
protected:
	static string p_three_letter_code;	
	const static char p_one_letter_code; // default to ambiguity character
public:
	// Need accessor functions for variables so "virtual" functionality will work
	const char & one_letter_code() const {return p_one_letter_code;}
	const string & three_letter_code() const {return p_three_letter_code;}
	Glutambiguous * new_clone() const {return new Glutambiguous(*this);}
	const static string name;
};

ostream & operator<<(ostream & os, const Residue & residue);

#endif
