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
