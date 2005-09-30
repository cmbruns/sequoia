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
// $Header$
// $Log$
// Revision 1.4  2004/06/04 19:11:41  cmbruns
// Updated GPL header
//
// Set one_letter_code values here instead of header, so it will compile on baxter (Linux)
//

#include "Residue.h"

// static member string variables must be initialized in a .cpp file, not a header file

const char Residue::p_one_letter_code = '?'; // Should this base class be used for "unknown?"

string AminoAcid::p_three_letter_code = "Xxx"; // Should this base class be used for "unknown?"
const char AminoAcid::p_one_letter_code = 'X'; // Should this base class be used for "unknown?"

string GapResidue::p_three_letter_code = "---";
const char GapResidue::p_one_letter_code = '-';

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

const char       Alanine::p_one_letter_code = 'A';
const char Asparambiguous::p_one_letter_code = 'B';
const char      Cysteine::p_one_letter_code = 'C';
const char     Aspartate::p_one_letter_code = 'D';
const char     Glutamate::p_one_letter_code = 'E';
const char Phenylalanine::p_one_letter_code = 'F';
const char       Glycine::p_one_letter_code = 'G';
const char     Histidine::p_one_letter_code = 'H';
const char    Isoleucine::p_one_letter_code = 'I';
const char        Lysine::p_one_letter_code = 'K';
const char       Leucine::p_one_letter_code = 'L';
const char	  Methionine::p_one_letter_code = 'M';
const char    Asparagine::p_one_letter_code = 'N';
const char       Proline::p_one_letter_code = 'P';
const char     Glutamine::p_one_letter_code = 'Q';
const char      Arginine::p_one_letter_code = 'R';
const char        Serine::p_one_letter_code = 'S';
const char     Threonine::p_one_letter_code = 'T';
const char Selenocysteine::p_one_letter_code = 'U';
const char        Valine::p_one_letter_code = 'V';
const char    Tryptophan::p_one_letter_code = 'W';
const char      Tyrosine::p_one_letter_code = 'Y';
const char Glutambiguous::p_one_letter_code = 'Z';

// Functions 

Residue * new_protein_residue(const char one_letter_code) {
	char olc = toupper(one_letter_code);
	switch (olc) {
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
		case '-': return new GapResidue;
		default: return new AminoAcid; // unknown
	}
}

ostream & operator<<(ostream & os, const Residue & residue) {
	os << residue.one_letter_code();
	return os;
}


