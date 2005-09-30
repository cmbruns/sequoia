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
// Revision 1.2  2004/06/14 16:37:59  cmbruns
// Added minus operator
// renamed AlignmentScore to DebugAlignmentScore
// sketched mechanism for reverting AlignmentScore to double once testing is complete
//
// Revision 1.1  2004/06/07 18:47:02  cmbruns
// Renamed Score object to AlignmentScore
//
// Revision 1.5  2004/06/04 19:15:07  cmbruns
// Updated GPL header
//
// Moved conservidue pair score function from here to be a member function of conservidue
//

#include "AlignmentScore.h"

DebugAlignmentScore operator*(const double r, const DebugAlignmentScore & score) {
	DebugAlignmentScore answer = score;
	answer.match() *= r;
	answer.extension() *= r;
	answer.opening() *= r;
	answer.closing() *= r;
	answer.deletion() *= r;
	answer.transition() *= r;
	return answer;
}

