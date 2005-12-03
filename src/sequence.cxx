#include "sequence.hxx"

// $Id: sequence.cxx,v 1.4 2002/06/21 00:45:00 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/sequence.cxx,v 1.4 2002/06/21 00:45:00 bruns Exp $
// $Log: sequence.cxx,v $
// Revision 1.4  2002/06/21 00:45:00  bruns
// Made Alignment score_val, and z_score initialize to 0, in case nothing
// is aligned
//
// Revision 1.3  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

#define debug 0

void Sequence::initialize_internals(void)
{
  title() = "Sequence with no title";
  weight() = 1.0;
  coord_ptr() = NULL;
  *this = "";
}

Sequence::Sequence(void)
{
  initialize_internals();
}

Sequence::~Sequence()
{
}

Sequence::Sequence(const char * cp)
{
  initialize_internals(); // use void constructor first
  *this = cp;
}

Sequence::Sequence(unsigned int n, uchar c)
{
  initialize_internals();
  unsigned int i;
  for (i=0; i < n; ++i) 
    *this += c; // fill array with string of char
}

Sequence::Sequence(const Sequence & s) : Array<SeqRes>(s) // use array copy constructor first
{
  *this = s;
}

Sequence & Sequence::operator=(const Sequence & s)
{
  if (this == &s) return *this;
  this->Array<SeqRes>::operator=(s);
  weight() = s.weight();
  title() = s.title();
  coord_ptr() = s.coord_ptr();
  return *this;
}

Sequence & Sequence::operator=(const char * cp)
{
  SeqRes sr;
  Sequence & seq = *this;
  uint i;
  uint slen = strlen(cp);
  init_array(0);
  for (i=0; i < slen; ++i)
    {
      sr.letter() = cp[i];
      sr.res_num() = i+1; // assign residue numbers
      seq += sr;
    }
  return *this;
}

Sequence & Sequence::operator+=(const char c)
{
  Sequence & s = *this;
  SeqRes sr;
  sr.letter() = c;
  // Guess at residue number
  if (length() < 1) 
    {
      if (is_gap(c)) sr.res_num() = 0;
      else sr.res_num() = 1;
    }
  // make gap have same number as previous actual residue
  else if (is_gap(c)) sr.res_num() = s[length()-1].res_num();
  else sr.res_num() = s[length()-1].res_num() + 1;
  s += sr;
  return *this;
}


Sequence & Sequence::operator+=(const SeqRes sr)
{
  this->Array<SeqRes>::operator+=(sr);
  return *this;
}

bool SeqRes::operator==(const SeqRes & s2) const
{
  if (letter() == s2.letter()) return true;
  else return false;
}

bool SeqRes::operator!=(const SeqRes & s2) const
{
  return (!(*this == s2));
}

bool Sequence::operator==(const Sequence & s2) const
{
  if (length() != s2.length()) return false;
  else if (length() < 1) return true;
  const Sequence & s = *this;
  uint i;
  for (i=0; i < length(); ++i)
    {
      if (s[i] != s2[i]) return false;
    }
  return true;
}

bool Sequence::operator!=(const Sequence & s2) const
{
  return !(*this == s2);
}

double Sequence::identity(const Sequence & s) const
{
  uint len, i;
  double id, count;
  const Sequence & t = * this;

  id = 0;
  count = 0;
  len = length();
  if (s.length() < len) len = s.length();
  for (i=0; i<len; ++i)
    {
      if (t[i].letter() == '-');
      else if (s[i].letter() == '-');
      else if (t[i].letter() == '.');
      else if (s[i].letter() == '.');
      else
	{
	  if (toupper(t[i].letter()) == toupper(s[i].letter()))
	    ++ id;
	  ++ count;
	}
    }
  if (count == 0) return 0;
  else return (id/count);
}

Sequence Sequence::remove_gaps() const
{
  Sequence answer = *this;
  answer.remove('-').remove('.').remove(' ');
  return answer;
}

uint Sequence::n_aligned(const Sequence & s2) const
{
  uint answer = 0;
  const Sequence & t = *this;
  uint max = MIN(length(), s2.length());
  uint i;
  for (i=0; i < max; ++i)
    if ((!is_gap(t[i].letter())) && (!is_gap(s2[i].letter()))) ++ answer;
  return answer;
}

uint Sequence::n_cap_aligned(const Sequence & s2) const
{
  uint answer = 0;
  const Sequence & t = *this;
  uint max = MIN(length(), s2.length());
  uint i;
  for (i=0; i < max; ++i)
    if (
	(!is_gap(t[i].letter())) 
	&& (!is_gap(s2[i].letter()))
	&& (t[i].letter() == toupper(t[i].letter()))
	&& (s2[i].letter() == toupper(s2[i].letter())))
      ++ answer;
  return answer;
}

Sequence Sequence::downcase() const
{
  Sequence answer = *this;
  uint i;
  for (i=0; i<length(); ++i)
    answer[i].letter() = tolower(answer[i].letter());
  return answer;
}

Sequence Sequence::upcase() const
{
  Sequence answer = *this;
  uint i;
  for (i=0; i < length(); ++i)
    answer[i].letter() = toupper(answer[i].letter());
  return answer;
}

Sequence & Sequence::remove(const char c)
{
  Sequence temp = *this;
  temp.init_array(0);
  Sequence & t = *this;
  uint i;
  for (i=0; i < length(); ++i)
    if (t[i].letter() != c) temp += t[i];
  *this = temp;
  return *this;
}

void seqcpy(char * s, const Sequence & seq)
{
  uint i;
  for (i=0; i < seq.length(); ++i)
    s[i] = seq[i].letter();
  s[seq.length()] = '\0';
}

ostream & operator<<(ostream & os, const Sequence & s)
{
  os << ">" << s.title() << "\n\n";
  // Should we output a numbering cue?
  // First make sure there are residues to count
  if (s.length() > 0)
    // Next, see if first character has an unusual number
    if ((is_gap(s[0].letter()) && s.first() != 0) ||
	((!is_gap(s[0].letter())) && (s.first() != 1)))
      {
	// put the numbering card
	os << "#" << s.first();
      }
  uint i;
  for (i=0; i<s.length(); ++i)
    os << s[i].letter();
  os << "*";
  return os;
}

istream & operator>>(istream & is, Sequence & s)
{
  uchar c;
  char buffer[2001];
  SeqRes sr;

  if (debug)
    {
      cout << "Reading in a sequence\n";
    }
  s = "";
  // 1. find a '>' character
  while ( (c != '>') && is.good() )
    is >> c;
  if (!is.good())
      s.title() = "";
  else
    // 2. read title until carriage return
    {
      is.getline(buffer, 2000);
      s.title() = buffer;
      if (debug) cout << "Title = " << s.title() << "\n";
      // Now read in the sequence data
      if (!is.good());
      else
	{
	// first look for #XXX identifier
	  is >> c;
	  if (is.good() && (c == '#'))
	    {
	      Mystring thing = "";
	      // check for minus sign
	      is >> c;
	      if ((c == '-') && (is.good())) {thing +=c; is >> c;}
	      while (isdigit(c) && is.good()) {thing += c; is >> c;}
	      sr.res_num() = atoi(thing.c_str());
	    }
	  else 
	    {
	      if (is_gap(c)) sr.res_num() = 0;
	      else sr.res_num() = 1;
	    }

	// 3. read sequence until '*' or end of file or '>'
	  while ( (c != '*') && is.good() && (c != '>'))
	    {
	      sr.letter() = c;
	      if (s.length() < 1) 
		{
		  if (is_gap(c)) sr.res_num() = 0;
		  else sr.res_num() = 1;
		}
	      // make gap have same number as previous actual residue
	      else if (is_gap(c)) sr.res_num() = s[s.length()-1].res_num();
	      else sr.res_num() = s[s.length()-1].res_num() + 1;
	      s += sr;
	      is >> c; // (skips whitespace)
	      if (debug) cout << c;
	    }
	  if (c == '>') // return it to the stream, belongs to next sequence
	    is.putback(c);
	  if (debug)
	    cout << "sequence = (" << s.dim() << ")" << s << "\n";
	}
    }
  return is;
}

bool is_gap(const char c)
{
  if (c == '-') return true;
  else if (c == ' ') return true;
  else if (c == '.') return true;
  else return false;
}

Alignment::Alignment(void)
{
  nseq = 0;
  seq = NULL;
  new_gaps = 0;
  score_val = 0; // to avoid undefined value for null alignments
  z_score = 0;
  mean = 0;
}

Alignment::~Alignment(void)
{
  nseq = 0;
  delete [] seq;
  seq = NULL;
}

Alignment::Alignment(const Alignment & a)
{
  unsigned int i;
  nseq = a.nseq;
  new_gaps = a.new_gaps;
  seq = new Sequence[nseq];
  for (i=0; i < nseq; ++i)
    seq[i] = a.seq[i];
  z_score = a.z_score;
  score_val = a.score_val;
  mean = a.mean;
}

Alignment::Alignment(const Sequence & s)
{
  nseq = 1;
  seq = new Sequence[1];
  seq[0] = s;
  new_gaps = 0;
  score_val = 0; // to avoid undefined value for null alignments
  z_score = 0;
  mean = 0;
}

Alignment & Alignment::operator=(const Alignment & a)
{
  if (this == &a) return *this;
  unsigned int i;
  new_gaps = a.new_gaps;
  z_score = a.z_score;
  score_val = a.score_val;
  mean = a.mean;
  nseq = a.nseq;
  delete [] seq;
  seq = new Sequence[nseq];
  for (i=0; i < nseq; ++i)
    seq[i] = a.seq[i];

  return *this;
}

Alignment & Alignment::operator+=(const Alignment & a)
{
  uint i;
  for (i=0; i < a.n_seq(); ++i)
    *this += a[i];
  return *this;
}

Alignment & Alignment::operator+=(const Sequence & s)
{
  Sequence * sp = new Sequence[n_seq() + 1];
  unsigned int i;
  for (i=0; i < n_seq(); ++i)
    sp[i] = seq[i];
  sp[n_seq()] = s;
  delete [] seq;
  seq = sp;
  ++nseq;
  return *this;
}

bool Alignment::operator==(const Alignment & a) const
{
  if (n_seq() != a.n_seq()) return false;
  const Alignment & t = *this;
  uint i;
  for (i=0; i < n_seq(); ++i)
    if (t[i] != a[i]) return false;
  return true;
}

ostream & Alignment::summarize(ostream & os) const
{
  uint i;
  const Alignment & t = *this;

  for (i=0; i < n_seq(); ++i)
    {
      os.width(3);
      os << i + 1 << ") ";
      os << t[i].title();
      if (t[i].weight() != 1)
        {
          os << " (w=";
          os.precision(2);
          os << t[i].weight() << ")";
        }
      os << endl;
    }
  return os;
}

ostream & Alignment::pretty(ostream & os, uint line_len) const
{
  const  Alignment & t = *this;
  summarize(os);
  os << "\n";

  uint i, j, count;

  count = 0;
  uint len;
  if (n_seq()) len = t[0].length();
  else len = 0;
  char c;
  while(count < len)
    {
      // Top summary line: dots and a useless number
      os << "      ";
      for (j=0; j < line_len; ++j)
	{
	  if (! ((j + count + 1) % 10) ) os << ".";
	  else os << " ";
	}
      os << count + line_len << "\n";
      // Block of sequences
      for (i=0; i < n_seq(); ++i)
	{
	  os.width(4);
	  os << i + 1 << ") ";
	  for (j=0; j < line_len; ++j)
	    {
	      if (t[i].length() > (count + j))
		{
		  c = t[i][count+j].letter();
		  os << c;
		}
	    }
	  if (t[i].length() > (count + line_len))
            /* too large by one ? */
	    /* os << "  " << t[i][count + line_len].res_num() << "\n"; */
	    os << "  " << t[i][count + line_len - 1].res_num() << "\n";
	  else os << "  " << t[i][t[i].length() - 1].res_num() << "\n";
	}
      // print out summary line ("*" for identities)
      os << "      ";
      for (j=0; j < line_len; ++j)
	{
	  if ((count+j) >= length()) os << " "; // past end of sequence
          else if (t.nongaps(count+j) < 2) os << " "; // just 1 sequence
	  else if (t.identity(count+j) < 1) os << " ";
	  // else if (t.identity(count+j) < 1) os << t.nongaps(count+j);
          // else if (t.nongaps(count+j) < 2) os << "1"; // just 1 sequence
	  else os << "*";
	}
      os << "\n";
      os << "\n";
      count += line_len;
    }
  return os;
}

void Alignment::clear()
{
  nseq = 0;
  delete [] seq;
  seq = NULL;
}

void Alignment::set_weights()
{
  uint i,j;
  double proximity;
  double min, scale;
  Alignment & a = *this;
  
  min = 1000000;
  for (i=0; i < n_seq(); ++i)
    {
      a[i].weight() = 0; /* initialize weight to zero */
      
      /* Down-weight sequence based upon proximity to other sequences */
      for (j=0; j < n_seq(); ++j)
	{
	  /* Other indices of proximity between two sequences might 
	     be chosen.  Any such index should vary between zero and 
	     some well defined upper limit.  Here I choose the very 
	     simple proximity index of (sequence_identity - 5%).  
	     Hopefully this will never go below zero, and cannot 
	     exceed 95%) */
	  proximity = a[i].identity(a[j]) - 0.04;
	  a[i].weight() += proximity;
	}

      /* I said DOWN-weight, you fool */
      a[i].weight() = 1 / a[i].weight();

      if (min > a[i].weight()) min = a[i].weight();
    }
  scale = 1/min; /* make SMALLEST weight 1 */
  for (i=0; i < n_seq(); ++i) a[i].weight() *= scale;
}


// overall identity at a given position
double Alignment::identity(uint posn) const
{
  const Alignment & a = *this;
  double answer = 0.0;
  if (n_seq() < 2) return 1;
  uint i, j;
  for (i=0; i < n_seq()-1; ++i)
   {
    if (a[i][posn].letter() == '-'); // ignore gaps
    else for (j=i+1; j < n_seq(); ++j)
     {
      if (a[j][posn].letter() == '-'); // ignore gaps
      else if (toupper(a[i][posn].letter()) == toupper(a[j][posn].letter())) answer += 1;
     }
   }
  // only count identities with no gaps
  double nseq2 = a.nongaps(posn);
  answer /= (nseq2 * (nseq2 - 1) * 0.5);
  // answer /= (n_seq() * (n_seq() - 1) * 0.5);
  return answer;
}


// how many actual sequences at a given position
// (for use with determining where to place stars in pretty)
int Alignment::nongaps(uint posn) const
{
  const Alignment & a = *this;
  int answer = 0;
  uint i;
  for (i=0; i < n_seq(); ++i)
    if (a[i][posn].letter() != '-') answer += 1;
  return answer;
}


Alignment Alignment::remove_gaps() const
{
  const Alignment & t = *this;
  Alignment answer = *this;
  if (!n_seq());
  else if (n_seq() == 1) answer[0] = answer[0].remove_gaps();
  else
    {
      bool isgap;
      uint posn, seq;
      for (seq=0; seq < n_seq(); ++ seq) answer[seq] = "";
      for (posn=0; posn < length(); ++ posn)
	{
	  isgap = true;
	  for (seq=0; seq < n_seq(); ++seq) 
	    if (!is_gap(t[seq][posn].letter())) isgap = false;
	  if (!isgap) 
	    for (seq=0; seq < n_seq(); ++seq)
	      answer[seq] += t[seq][posn];
	}
    }
  return answer;
}

ostream & Alignment::short_output(ostream & os) const
{
  unsigned int i;
  for (i=0; i < nseq; ++i)
    os << seq[i] << "\n";
  return os;
}

ostream & operator<<(ostream & os, const Alignment & a)
{
  a.short_output(os);
  return os;
}

istream & operator>>(istream & is, Alignment & a)
{
  Sequence s;
  while (is.good())
    {
      is >> s;
      if (s.length()) a += s;
    }
  return is;
}

Alignment Alignment::upcase() const
{
  Alignment answer = *this;
  uint i;
  for (i=0; i < n_seq(); ++i)
    answer[i] = answer[i].upcase();
  return answer;
}

Alignment Alignment::downcase() const
{
  Alignment answer = *this;
  uint i;
  for (i=0; i < n_seq(); ++i)
    answer[i] = answer[i].downcase();
  return answer;
}

