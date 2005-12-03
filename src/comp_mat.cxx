#include "comp_mat.hxx"
// #define DEBUG1
#undef DEBUG1

// $Id: comp_mat.cxx,v 1.3 2001/12/12 18:26:27 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/comp_mat.cxx,v 1.3 2001/12/12 18:26:27 bruns Exp $
// $Log: comp_mat.cxx,v $
// Revision 1.3  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

Comparison_matrix::Comparison_matrix(void) : Matrix<MATRIX_TYPE>(256,256,(MATRIX_TYPE)0)
{
}

Comparison_matrix::Comparison_matrix(int flavor) : Matrix<MATRIX_TYPE>(256,256,(MATRIX_TYPE)0)
{
  int i, j, i2, j2;
//  Matrix<MATRIX_TYPE> temp(256, 256, 0);
//  this->element = temp;
  MATRIX_TYPE blo62trans[27][27] =
    {/* A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z  -         */
      { 8, 0, 4, 2, 3, 2, 4, 2, 3, 0, 3, 3, 3, 2, 0, 3, 3, 3, 5, 4, 0, 4, 1, 0, 2, 0, 0 }, /* A */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* B */
      { 4, 0,13, 1, 0, 2, 1, 1, 3, 0, 1, 3, 3, 1, 0, 1, 1, 1, 3, 3, 0, 3, 2, 0, 2, 0, 0 }, /* C */
      { 2, 0, 1,10, 6, 1, 3, 3, 1, 0, 3, 0, 1, 5, 0, 3, 4, 2, 4, 3, 0, 1, 0, 0, 2, 0, 0 }, /* D */
      { 3, 0, 0, 6, 9, 1, 2, 4, 1, 0, 5, 1, 2, 4, 0, 3, 6, 4, 4, 3, 0, 2, 1, 0, 2, 0, 0 }, /* E */
      { 2, 0, 2, 1, 1,10, 1, 3, 4, 0, 1, 4, 4, 1, 0, 0, 1, 1, 2, 2, 0, 3, 5, 0, 7, 0, 0 }, /* F */
      { 4, 0, 1, 3, 2, 1,10, 2, 0, 0, 2, 0, 1, 4, 0, 2, 2, 2, 4, 2, 0, 1, 2, 0, 1, 0, 0 }, /* G */
      { 2, 0, 1, 3, 4, 3, 2,12, 1, 0, 3, 1, 2, 5, 0, 2, 4, 4, 3, 2, 0, 1, 2, 0, 6, 0, 0 }, /* H */
      { 3, 0, 3, 1, 1, 4, 0, 1, 8, 0, 1, 6, 5, 1, 0, 1, 1, 1, 2, 3, 0, 7, 1, 0, 3, 0, 0 }, /* I */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* J */
      { 3, 0, 1, 3, 5, 1, 2, 3, 1, 0, 9, 2, 3, 4, 0, 3, 5, 6, 4, 3, 0, 2, 1, 0, 2, 0, 0 }, /* K */
      { 3, 0, 3, 0, 1, 4, 0, 1, 6, 0, 2, 8, 6, 1, 0, 1, 2, 6, 2, 3, 0, 5, 2, 0, 3, 0, 0 }, /* L */
      { 3, 0, 3, 1, 2, 4, 1, 2, 5, 0, 3, 6, 9, 2, 0, 2, 4, 3, 3, 3, 0, 5, 3, 0, 3, 0, 0 }, /* M */
      { 2, 0, 1, 5, 4, 1, 4, 5, 1, 0, 4, 1, 2,10, 0, 2, 4, 4, 5, 4, 0, 1, 0, 0, 2, 0, 0 }, /* N */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* O */
      { 3, 0, 1, 3, 3, 0, 2, 2, 1, 0, 3, 1, 2, 2, 0,11, 3, 2, 3, 3, 0, 2, 0, 0, 1, 0, 0 }, /* P */
      { 3, 0, 1, 4, 6, 1, 2, 4, 1, 0, 5, 2, 4, 4, 0, 3, 9, 5, 4, 3, 0, 2, 2, 0, 3, 0, 0 }, /* Q */
      { 3, 0, 1, 2, 4, 1, 2, 4, 1, 0, 6, 6, 3, 4, 0, 2, 5, 9, 3, 3, 0, 1, 1, 0, 2, 0, 0 }, /* R */
      { 5, 0, 3, 4, 4, 2, 4, 3, 2, 0, 4, 2, 3, 5, 0, 3, 4, 3, 8, 5, 0, 2, 1, 0, 2, 0, 0 }, /* S */
      { 4, 0, 3, 3, 3, 2, 2, 2, 3, 0, 3, 3, 3, 4, 0, 3, 3, 3, 5, 9, 0, 4, 2, 0, 2, 0, 0 }, /* T */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* U */
      { 4, 0, 3, 1, 2, 3, 1, 1, 7, 0, 2, 5, 5, 1, 0, 2, 2, 1, 2, 4, 0, 8, 1, 0, 3, 0, 0 }, /* V */
      { 1, 0, 2, 0, 1, 5, 2, 2, 1, 0, 1, 2, 3, 0, 0, 0, 2, 1, 1, 2, 0, 1,15, 0, 6, 0, 0 }, /* W */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* X */
      { 2, 0, 2, 2, 2, 7, 1, 6, 3, 0, 2, 3, 3, 2, 0, 1, 3, 2, 2, 2, 0, 3, 6, 0,11, 0, 0 }, /* Y */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* Z */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  /* - */
    };
  MATRIX_TYPE blo62[27][27] =
    {/* A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z  -         */
      { 4, 0, 0,-2,-1,-2, 0,-2,-1, 0,-1,-1,-1,-2, 0,-1,-1,-1, 1, 0, 0, 0,-3, 0,-2, 0, 0 }, /* A */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* B */
      { 0, 0, 9,-3,-4,-2,-3,-3,-1, 0,-3,-1,-1,-3, 0,-3,-3,-3,-1,-1, 0,-1,-2, 0,-2, 0, 0 }, /* C */
      {-2, 0,-3, 6, 2,-3,-1,-1,-3, 0,-1,-4,-3, 1, 0,-1, 0,-2, 0,-1, 0,-3,-4, 0,-2, 0, 0 }, /* D */
      {-1, 0,-4, 2, 5,-3,-2, 0,-3, 0, 1,-3,-2, 0, 0,-1, 2, 0, 0,-1, 0,-2,-3, 0,-2, 0, 0 }, /* E */
      {-2, 0,-2,-3,-3, 6,-3,-1, 0, 0,-3, 0, 0,-3, 0,-4,-3,-3,-2,-2, 0,-1, 1, 0, 3, 0, 0 }, /* F */
      { 0, 0,-3,-1,-2,-3, 6,-2,-4, 0,-2,-4,-3, 0, 0,-2,-2,-2, 0,-2, 0,-3,-2, 0,-3, 0, 0 }, /* G */
      {-2, 0,-3,-1, 0,-1,-2, 8,-3, 0,-1,-3,-2, 1, 0,-2, 0, 0,-1,-2, 0,-3,-2, 0, 2, 0, 0 }, /* H */
      {-1, 0,-1,-3,-3, 0,-4,-3, 4, 0,-3, 2, 1,-3, 0,-3,-3,-3,-2,-1, 0, 3,-3, 0,-1, 0, 0 }, /* I */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* J */
      {-1, 0,-3,-1, 1,-3,-2,-1,-3, 0, 5,-2,-1, 0, 0,-1, 1, 2, 0,-1, 0,-2,-3, 0,-2, 0, 0 }, /* K */
      {-1, 0,-1,-4,-3, 0,-4,-3, 2, 0,-2, 4, 2,-3, 0,-3,-2, 2,-2,-1, 0, 1,-2, 0,-1, 0, 0 }, /* L */
      {-1, 0,-1,-3,-2, 0,-3,-2, 1, 0,-1, 2, 5,-2, 0,-2, 0,-1,-1,-1, 0, 1,-1, 0,-1, 0, 0 }, /* M */
      {-2, 0,-3, 1, 0,-3, 0, 1,-3, 0, 0,-3,-2, 6, 0,-2, 0, 0, 1, 0, 0,-3,-4, 0,-2, 0, 0 }, /* N */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* O */
      {-1, 0,-3,-1,-1,-4,-2,-2,-3, 0,-1,-3,-2,-2, 0, 7,-1,-2,-1,-1, 0,-2,-4, 0,-3, 0, 0 }, /* P */
      {-1, 0,-3, 0, 2,-3,-2, 0,-3, 0, 1,-2, 0, 0, 0,-1, 5, 1, 0,-1, 0,-2,-2, 0,-1, 0, 0 }, /* Q */
      {-1, 0,-3,-2, 0,-3,-2, 0,-3, 0, 2, 2,-1, 0, 0,-2, 1, 5,-1,-1, 0,-3,-3, 0,-2, 0, 0 }, /* R */
      { 1, 0,-1, 0, 0,-2, 0,-1,-2, 0, 0,-2,-1, 1, 0,-1, 0,-1, 4, 1, 0,-2,-3, 0,-2, 0, 0 }, /* S */
      { 0, 0,-1,-1,-1,-2,-2,-2,-1, 0,-1,-1,-1, 0, 0,-1,-1,-1, 1, 5, 0, 0,-2, 0,-2, 0, 0 }, /* T */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* U */
      { 0, 0,-1,-3,-2,-1,-3,-3, 3, 0,-2, 1, 1,-3, 0,-2,-2,-3,-2, 0, 0, 4,-3, 0,-1, 0, 0 }, /* V */
      {-3, 0,-2,-4,-3, 1,-2,-2,-3, 0,-3,-2,-1,-4, 0,-4,-2,-3,-3,-2, 0,-3,11, 0, 2, 0, 0 }, /* W */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* X */
      {-2, 0,-2,-2,-2, 3,-3, 2,-1, 0,-2,-1,-1,-2, 0,-3,-1,-2,-2,-2, 0,-1, 2, 0, 7, 0, 0 }, /* Y */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* Z */
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  /* - */
    };
  switch (flavor)
    {
    case BLOSUM62:
      i2 = 0;
      for (i = 'A'; i <= 'Z'; ++i)
	{
	  j2 = 0;
	  for (j='A'; j <= 'Z'; ++j)
	    {
	      this->operator[](i)[j] = blo62[i2][j2];
	      ++j2;
	    }
	  ++i2;
	}
      i2 = 0;
      for (i = 'a'; i <= 'z'; ++i)
	{
	  j2 = 0;
	  for (j='a'; j <= 'z'; ++j)
	    {
	      this->operator[](i)[j] = blo62[i2][j2];
	      ++j2;
	    }
	  ++i2;
	}
      i2 = 0;
      for (i = 'a'; i <= 'z'; ++i)
	{
	  j2 = 0;
	  for (j='A'; j <= 'Z'; ++j)
	    {
	      this->operator[](i)[j] = this->operator[](j)[i] = blo62[i2][j2];
	      ++j2;
	    }
	  ++i2;
	}
      this->used = "CSTPAGNDEQHRKMILVFYWcstpagndeqhrkmilvfyw-";
      this->title = "BLOSUM 62 matrix (PNAS 89:10915)";
    }
}

Comparison_matrix::~Comparison_matrix(void)
{
}

// const Matrix_Row<MATRIX_TYPE> & Comparison_matrix::operator[](unsigned int m) const 
// {return element[m];}

// Matrix_Row<MATRIX_TYPE> & Comparison_matrix::operator[](unsigned int m) 
// {return element[m];}

ostream & operator<<(ostream & os, const Comparison_matrix & cm)
{
  int i, j, i2, count;
  int fieldwidth = 3;
  Mystring phlerb;
  
  os << cm.title << "\n";
  count = cm.used.length();
  os << "  ";
  for (j=0; j < count; ++j)
    {
      phlerb = cm.used[j]; // to make fieldwidth work
      os.width(fieldwidth);
      os << phlerb;
    }
  os << "\n";
  for (i=0; i < count; ++i)
    {
      i2 = cm.used[i];
      os << (char) i2 << " ";
      for (j=0; j < count; ++j)
	{
	  os.width(fieldwidth);
	  os << cm[i2][cm.used[j]];
	}
    os << "\n";
    }
  return os;
}

istream & operator>>(istream & is, Comparison_matrix & cm)
{
  char c2;
  unsigned int c;
  int  i;
  int count = 0; /* number of residues used in the table */
  char buffer[2000]; /* storage of lines from file */
  Mystring buf_string;
  Array<Mystring> tokens;

  /* read in first line */
  is.getline(buffer, 2000);
  cm.title = buffer;

  /* now read in characters that define the rows and columns */
  is.getline(buffer, 2000);
  buf_string = buffer;
  buf_string = buf_string.remove(' '); // remove spaces
  cm.used = buf_string;
  count = cm.used.length();
#ifdef DEBUG1
  cerr << "cm.used = " << cm.used << "***\n";
#endif
  if (is.eof()) cerr << "end of file before data\n";

 /* now read in the actual data */
 /* the values are in free format and must be separated by spaces */
 while (is.good())
   {
     is >> c2;
     c = c2; // to avoid compiler warnings
     is.getline(buffer, 2000);
     buf_string = buffer;
     tokens = buf_string.split();
     count = tokens.dim();
#ifdef DEBUG1
cerr << "tokens[0] = " << tokens[0] << "\n";
cerr << "tokens[1] = " << tokens[1] << "\n";
#endif
     // for (i=1; i < count; ++i) // tokens[0] seems to be "" (?!)
     for (i=0; i < count; ++i)
       {
        const char * cp = tokens[i].c_str();
	strcpy(buffer, cp);
	// cm[cm.used[i-1]][c] = atof(buffer);
	// cm[c][cm.used[i-1]] = atof(buffer);
        // want to ignore final row of letters, if there
        if (isdigit(buffer[0]))
          {
            cm[cm.used[i]][c] = atof(buffer);
            cm[c][cm.used[i]] = atof(buffer);
          }
      }
   }
 return is;
}
