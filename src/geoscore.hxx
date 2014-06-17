#ifndef _GEOSCORE_HXX_
#define _GEOSCORE_HXX_

// $Id: geoscore.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/geoscore.hxx,v 1.2 2001/11/28 23:08:17 bruns Exp $
// $Log: geoscore.hxx,v $
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

#include "rmat.hxx"

double distance_score1(double distance, double sigma);
double distance_score2(double distance, double sigma, double sphere_rad);
double orientation_score1(double angle, double sigma);
double orientation_score2(double angle, double sigma);
double kappa_cos_angle(const RMat & o1, const RMat & o2);

#endif
