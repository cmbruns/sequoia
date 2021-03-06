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
// Revision 1.3  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "Modres.h"
#include <cstdio>

Modres::Modres(string line) {
  private_record_name           = line.substr(0,6);
  private_id_code               = line.substr(7,4);
  private_residue_name          = line.substr(12,3);
  private_chain_id              = line.substr(16,1)[0];
  private_sequence_number       = atoi(line.substr(18,4).c_str());
  private_insertion_code        = line.substr(22,1)[0];
  private_standard_residue_name = line.substr(24,3);
  private_comment               = line.substr(29,40);
}
string Modres::get_record_name() const           {return private_record_name;}
string Modres::get_id_code() const               {return private_id_code;}
string Modres::get_residue_name() const          {return private_residue_name;}
char   Modres::get_chain_id() const              {return private_chain_id;}
int    Modres::get_sequence_number() const       {return private_sequence_number;}
char   Modres::get_insertion_code() const        {return private_insertion_code;}
string Modres::get_standard_residue_name() const {return private_standard_residue_name;}
string Modres::get_comment() const               {return private_comment;}

// COLUMNS        DATA TYPE       FIELD          DEFINITION
// --------------------------------------------------------------------------------
//  1 -  6        Record name     "MODRES"
// 
//  8 - 11        IDcode          idCode         ID code of this entry.
// 
// 13 - 15        Residue name    resName        Residue name used in this entry.
// 
// 17             Character       chainID        Chain identifier.
// 
// 19 - 22        Integer         seqNum         Sequence number.
// 
// 23             AChar           iCode          Insertion code.
// 
// 25 - 27        Residue name    stdRes         Standard residue name.
// 
// 30 - 70        String          comment        Description of the residue
//                                               modification.
