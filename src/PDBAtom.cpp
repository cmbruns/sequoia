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
// Revision 1.5  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.4  2002/05/15 19:27:41  bruns
// Added get_dihedral_angle (non-member) subroutine
//
// Revision 1.3  2001/12/14 19:21:55  bruns
// Changes to compile with gcc 3.0.2
//  - "using namespace std;"  to use standard classes
//  - "ios::fmtflags" instead of "fmtflags"
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//

#include "PDBAtom.h"
#include <cstdio> // atoi(), atof()
#include <set>
#include "ChemicalElement.h"
#include "Hash3D.h"

// return a vector containing 20 points on a sphere around an atom
vector<Vector3D> sphere_points(float radius, Vector3D center) {
	// 20 points on a unit sphere for crude surface area computation
	static vector<Vector3D> unit_sphere_points;
	if (unit_sphere_points.size() < 1) {
		// Initialize unit sphere
		// These are the vertices of a dodecahedron
		unit_sphere_points.push_back(Vector3D(0.35682, 0.00000, -0.93417));
		unit_sphere_points.push_back(Vector3D(0.35682, 0.00000, 0.93417));
		unit_sphere_points.push_back(Vector3D(-0.35682, 0.00000, -0.93417));
		unit_sphere_points.push_back(Vector3D(-0.35682, 0.00000, 0.93417));
		unit_sphere_points.push_back(Vector3D(0.00000, -0.93417, 0.35682));
		unit_sphere_points.push_back(Vector3D(0.00000, 0.93417, 0.35682));
		unit_sphere_points.push_back(Vector3D(0.00000, -0.93417, -0.35682));
		unit_sphere_points.push_back(Vector3D(0.00000, 0.93417, -0.35682));
		unit_sphere_points.push_back(Vector3D(-0.93417, 0.35682, 0.00000));
		unit_sphere_points.push_back(Vector3D(0.93417, 0.35682, 0.00000));
		unit_sphere_points.push_back(Vector3D(-0.93417, -0.35682, 0.00000));
		unit_sphere_points.push_back(Vector3D(0.93417, -0.35682, 0.00000));
		unit_sphere_points.push_back(Vector3D(0.57735, 0.57735, 0.57735));
		unit_sphere_points.push_back(Vector3D(0.57735, 0.57735, -0.57735));
		unit_sphere_points.push_back(Vector3D(0.57735, -0.57735, 0.57735));
		unit_sphere_points.push_back(Vector3D(0.57735, -0.57735, -0.57735));
		unit_sphere_points.push_back(Vector3D(-0.57735, 0.57735, 0.57735));
		unit_sphere_points.push_back(Vector3D(-0.57735, 0.57735, -0.57735));
		unit_sphere_points.push_back(Vector3D(-0.57735, -0.57735, 0.57735));
		unit_sphere_points.push_back(Vector3D(-0.57735, -0.57735, -0.57735));
	}

	vector<Vector3D> answer;
	for (unsigned int v = 0; v < unit_sphere_points.size(); v++) {
		Vector3D new_point = unit_sphere_points[v] * radius + center;
		answer.push_back(new_point);
	}
	return answer;
}

double PDBAtom::store_accessible_surface_area(const Hash3D & hash3d) {
	p_accessible_surface_area = compute_accessible_surface_area(hash3d);
	return p_accessible_surface_area;
}

// How wide is the cube of cubelets that must be searched?
#define CUBELET_OFFSET_COUNT 5

double PDBAtom::compute_accessible_surface_area(const Hash3D & hash3d) const {
	const PDBAtom & atom = *this;
	
	// Reject non-protein atoms
	if (atom.get_element() == NULL) return 0.0; // Unknown element, no surface area
	if (atom.get_element()->vdw_radius == 0.0) return 0.0; // Not part of solvent model
	
	// Set up points on accessible surface shell
	double probe_radius = 1.40;
	double atom_radius = atom.get_element()->vdw_radius;
	double shell_radius = atom_radius + probe_radius;
	const Vector3D atom_center = atom.get_coordinate();
	vector<Vector3D> shell_points = sphere_points(shell_radius, atom_center);

	// Now find other atoms that might occlude each point
	double min_distance_squared = 0.5;  // ignore self comparisons, distance to atom
	
	set<int> accessible_points; // start full
	for (unsigned int p = 0; p < shell_points.size(); p++)
		accessible_points.insert(p);
	
	// loop over center cubelet and neighbors
	int x_offset, y_offset, z_offset;
	Index3D center_index = hash3d.get_index(atom_center);
	// It would be most efficient to start with central cubelet
	int offset_order[CUBELET_OFFSET_COUNT] = {0, -1, 1, -2, 2};
	for (x_offset = 0; x_offset < CUBELET_OFFSET_COUNT; x_offset ++) {
		if (accessible_points.empty()) break;
		for (y_offset = 0; y_offset < CUBELET_OFFSET_COUNT; y_offset ++) {
			if (accessible_points.empty()) break;
			for (z_offset = 0; z_offset < CUBELET_OFFSET_COUNT; z_offset ++) {
				if (accessible_points.empty()) break;
				Index3D cubelet_index = center_index;
				cubelet_index.first += offset_order[x_offset];
				cubelet_index.second += offset_order[y_offset];
				cubelet_index.third += offset_order[z_offset];
				const Cubelet * cubelet = hash3d.get_cubelet(cubelet_index);
				if (cubelet == NULL) continue;
				// loop over points in cubelet
				for (unsigned int c = 0; c < cubelet->contents.size(); c++) {
					const PDBAtom * subject_atom = (const PDBAtom *)(cubelet->contents[c]->get_object());
					// TODO - only consider protein atoms
					if (subject_atom->get_record_name() != "ATOM  ") continue;
					if (subject_atom->get_element() == NULL) continue;
					if (subject_atom->get_element()->vdw_radius <= 0) continue;
					
					const Vector3D subject_center = subject_atom->get_coordinate();
					// ignore atoms that are the same as the center atom
					if (subject_center.distance_squared(atom_center) < min_distance_squared) continue;
					
					double subject_max_distance = subject_atom->get_element()->vdw_radius + probe_radius;
					double subject_max_distance_squared = subject_max_distance * subject_max_distance;
					
					// compare to shell point
					for (unsigned int point = 0; point < shell_points.size(); point ++) {
						// Don't do more work on already discarded points
						if (accessible_points.find(point) == accessible_points.end()) continue;
						
						const Vector3D & surf_point = shell_points[point];
						
						if (subject_center.distance_squared(surf_point) < subject_max_distance_squared) {							
							// handle case of occluded atom
							accessible_points.erase(point);
							if (accessible_points.empty()) break;
						} // end if close atom
					} // end for shell points
				} // end for cubelet atoms
			} // end for cubelet z
		} // end for cubelet y
	} // end for cubelet x
	
	double PI = 3.14159;
	// Surface area associated with one point on the shell
	// (don't) Compute the area on the vdW surface, vs. the outer "accessible surface"
	double point_area = 4 * PI * (atom_radius + probe_radius) * (atom_radius + probe_radius) / shell_points.size();	
	double answer = accessible_points.size() * point_area;
	
	return answer;
}

void PDBAtom::set_element() {
	if ((get_record_name() != "ATOM  ") && (get_record_name() != "HETATM")) {
		p_element = NULL;
		return;
	}
	string esymbol = get_element_symbol();
	if ((esymbol.length() == 2) && 
		(esymbol.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ") != string::npos)) {
		p_element = find_element(get_element_symbol());
		return;
	}
	if (get_atom_name().length() > 0) {
		p_element = find_element(get_atom_name());
		return;
	}
}

PDBAtom::PDBAtom(const string & s) : p_element(NULL), p_accessible_surface_area(-1) {
	if (s.length() < 6) {set_element();return;}
  private_record_name        = s.substr(0,  6);
  if (s.length() < 11) {set_element();return;}
  private_serial_number      = atoi( s.substr(6,  5).c_str() );
  if (s.length() < 16) {set_element();return;}
  private_atom_name          = s.substr(12, 4);
  if (s.length() < 17) {set_element();return;}
  private_alt_loc            = s.substr(16, 1)[0];
  if (s.length() < 20) {set_element();return;}
  private_residue_name       = s.substr(17, 3);
  if (s.length() < 22) {set_element();return;}
  private_chain_id           = s.substr(21, 1)[0];
  if (s.length() < 26) {set_element();return;}
  private_residue_number     = atoi( s.substr(22, 4).c_str() );
  if (s.length() < 27) {set_element();return;}
  private_insertion_code     = s.substr(26, 1)[0];
  if (s.length() < 38) {set_element();return;}
  private_coordinate.set_x() = atof( s.substr(30, 8).c_str() );
  if (s.length() < 46) {set_element();return;}
  private_coordinate.set_y() = atof( s.substr(38, 8).c_str() );
  if (s.length() < 54) {set_element();return;}
  private_coordinate.set_z() = atof( s.substr(46, 8).c_str() );
  if (s.length() < 60) {set_element();return;}
  private_occupancy          = atof( s.substr(54, 6).c_str() );
  if (s.length() < 66) {set_element();return;}
  private_temperature_factor = atof( s.substr(60, 6).c_str() );
  if (s.length() < 76) {set_element();return;}
  private_segment_identifier = s.substr(72, 4);
  if (s.length() < 78) {set_element();return;}
  private_element_symbol     = s.substr(76, 2);
  if (s.length() < 80) {set_element();return;}
  private_charge             = s.substr(78, 2);

  set_element();
  return;
}

string   PDBAtom::get_record_name() const {return private_record_name;}
int      PDBAtom::get_serial_number() const {return private_serial_number;}
string   PDBAtom::get_atom_name() const {return private_atom_name;}
char     PDBAtom::get_alt_loc() const {return private_alt_loc;}
string   PDBAtom::get_residue_name() const {return private_residue_name;}
char     PDBAtom::get_chain_id() const {return private_chain_id;}
int      PDBAtom::get_residue_number() const {return private_residue_number;}
char     PDBAtom::get_insertion_code() const {return private_insertion_code;}
Vector3D PDBAtom::get_coordinate() const {return private_coordinate;}
float    PDBAtom::get_occupancy() const {return private_occupancy;}
float    PDBAtom::get_temperature_factor() const {return private_temperature_factor;}
string   PDBAtom::get_segment_identifier() const {return private_segment_identifier;}
string   PDBAtom::get_element_symbol() const {return private_element_symbol;}
string   PDBAtom::get_charge() const {return private_charge;}

double   PDBAtom::distance(const PDBAtom & atom2) const {
  return get_coordinate().distance(atom2.get_coordinate());
}

ostream & operator<<(ostream & os, const PDBAtom & atom) {
  // Remember initial format flags settings
  ios::fmtflags old_format = os.flags(); // new iostream way
  // long old_format = os.flags(); // old iostream way

  os << atom.get_record_name();
  os.width(5); os << atom.get_serial_number();
  os << " ";
  os << atom.get_atom_name();

  os << atom.get_alt_loc();
  os << atom.get_residue_name();
  os << " ";
  os << atom.get_chain_id();
  os.width(4); os << atom.get_residue_number();
  os << atom.get_insertion_code();
  os << "   ";
  os.setf(ios::fixed); // Makes precision be what I want?
  os.precision(3); os.width(8); os << atom.get_coordinate().get_x();
  os.precision(3); os.width(8); os << atom.get_coordinate().get_y();
  os.precision(3); os.width(8); os << atom.get_coordinate().get_z();
  os.precision(2); os.width(6); os << atom.get_occupancy();
  os.precision(2); os.width(6); os << atom.get_temperature_factor();
  os << "      ";
  os << atom.get_segment_identifier();
  os << atom.get_element_symbol();
  os << atom.get_charge();

  // Restore initial settings
  os.flags(old_format);

  return os;
}


// COLUMNS        DATA TYPE       FIELD         DEFINITION
// ---------------------------------------------------------------------------------
//  1 -  6        Record name     "ATOM  "
// 
//  7 - 11        Integer         serial        Atom serial number.
// 
// 13 - 16        Atom            name          Atom name.
// 
// 17             Character       altLoc        Alternate location indicator.
// 
// 18 - 20        Residue name    resName       Residue name.
// 
// 22             Character       chainID       Chain identifier.
// 
// 23 - 26        Integer         resSeq        Residue sequence number.
// 
// 27             AChar           iCode         Code for insertion of residues.
// 
// 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
//                                              Angstroms.
// 
// 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
//                                              Angstroms.
// 
// 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
//                                              Angstroms.
// 
// 55 - 60        Real(6.2)       occupancy     Occupancy.
// 
// 61 - 66        Real(6.2)       tempFactor    Temperature factor.
// 
// 73 - 76        LString(4)      segID         Segment identifier, left-justified.
// 
// 77 - 78        LString(2)      element       Element symbol, right-justified.
// 
// 79 - 80        LString(2)      charge        Charge on the atom.

double get_dihedral_angle(const PDBAtom & atom1,
		   const PDBAtom & atom2,
		   const PDBAtom & atom3,
		   const PDBAtom & atom4) {

  Vector3D center = atom3.get_coordinate() - atom2.get_coordinate();
  Vector3D end1 = atom1.get_coordinate() - atom2.get_coordinate();
  Vector3D end2 = atom4.get_coordinate() - atom3.get_coordinate();

  end2 = center.cross(end2).unit();
  end1 = center.cross(end1).unit();

  double cos_angle = end1.dot(end2);
  if (cos_angle > 1.0) cos_angle = 1.0;
  if (cos_angle < -1.0) cos_angle = -1.0;
  
  double angle = acos(cos_angle);

  double sign = end1.cross(end2).dot(center);
  if (sign < 0) {angle *= -1.0;}

  return angle;
}
