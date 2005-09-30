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
// Hash3D.h
// Utility class to help efficiently answer questions like:
//  "find all pairs of atoms in this structure that are within
//   6 Angstroms of each other"
// Geometric hashing can make this process O(n) complexity instead
// of O(n^2) complexity

// $Id$
// $Header$
// $Log$
// Revision 1.1  2004/06/04 19:34:46  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.3  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.2  2002/05/23 18:04:44  bruns
// Added CVS tags
//

#ifndef _HASH_3D_H
#define _HASH_3D_H

#include <map>
#include <vector>
#include "Vector3D.h"

class PDBAtom;
class PDBChain;
class PDBResidue;

// Index to define one cubelet in a 3D hash
class Index3D {
 public:
  int first;
  int second;
  int third;
};

// Comparison function for Index3D
struct ltindex
{
  bool operator()(const Index3D i1, const Index3D i2) const
  {
    if (i1.first < i2.first) return true;
    if (i1.first > i2.first) return false;
    if (i1.second < i2.second) return true;
    if (i1.second > i2.second) return false;
    if (i1.third < i2.third) return true;
    if (i1.third > i2.third) return false;
    return false; // Keys are equal
  }
};

// The object contained in a cubelet
// Each type of thing we want to stuff into the hash should make
//  a descendant of this class
class BaseHash3DContent {
	friend class Hash3D;
protected:
	virtual BaseHash3DContent * new_clone() const = 0; // Needed for memory management
 public:
  virtual const void * get_object() const = 0;
  virtual Vector3D get_position() const = 0;
};

class PDBAtomHash3DContent : public BaseHash3DContent {
	friend class Hash3D;
protected:
	const PDBAtom * p_atom;
	PDBAtomHash3DContent * new_clone() const {return new PDBAtomHash3DContent(*this);}
public:
	const PDBAtom * get_object() const {return p_atom;}
	Vector3D get_position() const;
	PDBAtomHash3DContent(const PDBAtom & atom) {
		p_atom = & atom;
	}
};

// One volume element in a 3D hash
class Cubelet {
private:
	Index3D p_index;
public:
	const Index3D & get_index() const {return p_index;}
	void set_index(const Index3D & i) {p_index = i;}
	vector<const BaseHash3DContent *> contents;
	~Cubelet() {
		for (unsigned int c = 0; c < contents.size(); c++) 
			delete contents[c];
	}
};

class Hash3D {
 private:
  map<Index3D, Cubelet, ltindex> p_cubelet;
  double p_edge_length;
 public:
  const Index3D get_index(const Vector3D & v) const {
    int x, y, z;

    x = (int)(v.get_x() / p_edge_length);
    y = (int)(v.get_y() / p_edge_length);
    z = (int)(v.get_z() / p_edge_length);

    if (v.get_x() < 0) x -= 1;
    if (v.get_y() < 0) y -= 1;
    if (v.get_z() < 0) z -= 1;

    Index3D i;
    i.first = x;
    i.second = y;
    i.third = z;

    return i;
  }
  void add_object(const BaseHash3DContent & c) {
    Index3D index = get_index(c.get_position());
    cubelet(index).contents.push_back(c.new_clone());
  }

  // atom must be persistent!!!!
  void add_pdb_atom(const PDBAtom & atom) {
	  PDBAtomHash3DContent content(atom);
	  add_object(content);
  }
  void add_pdb_residue(const PDBResidue & residue);
  void add_pdb_chain(const PDBChain & chain);
	
  const Cubelet * get_cubelet(const Index3D & index) const {
	  if (p_cubelet.find(index) != p_cubelet.end())
		  return &(p_cubelet.find(index)->second);
	  return NULL; // Not found
  }
  const Cubelet * get_cubelet(const Vector3D & v) const {
   return get_cubelet(get_index(v));
  }

  Cubelet & cubelet(const Index3D & index) {
	  p_cubelet[index].set_index(index); // In case it's new
	  return p_cubelet[index];
  }
  Cubelet & cubelet(const Vector3D & v) {
	  return cubelet(get_index(v));
  }
  Hash3D(double edge_length) {p_edge_length = edge_length;}
};

#endif
