#include "mystring.hxx"

// $Id: mystring.cxx,v 1.3 2001/12/12 18:26:27 bruns Exp $
// $Header: /usr/data/cvs/sequoia1/mystring.cxx,v 1.3 2001/12/12 18:26:27 bruns Exp $
// $Log: mystring.cxx,v $
// Revision 1.3  2001/12/12 18:26:27  bruns
// Many changes to get compilation with gcc 3.0.2
// Removed all templated friend functions (no friends needed)
// Added .c_str() method to Mystring, and removed auto conversion to const char *
//
// Revision 1.2  2001/11/28 23:40:22  bruns
// Added cvs header tags
// Removed ^M characters
//

// Start with at least 8 bytes
Mystring::Mystring() : Array<char> (1, '\0')
{}

Mystring::Mystring(const char * s) : Array<char>(s, strlen(s) + 1)
{
}

Mystring::Mystring(const char c) : Array<char> (2)
{
  Mystring & t = *this;
  t[0] = c;
  t[1] = '\0';
}

Mystring::Mystring(const Mystring & s) : Array<char> (s.length() + 1)
{
  *this = s;
}

Mystring::~Mystring()
{
}

Mystring & Mystring::operator=(const Mystring & ms)
{
  if (this == &ms) return *this;
  Array<char> & ac = *this;
  ac = (Array<char>) ms;
  return *this;
}

Mystring & Mystring::operator=(const char * s)
{
  Mystring temp(s);
  *this = temp;
  return *this;
}

Abstract_array<char> & Mystring::operator+=(const char & c)
{
  Array<char> & ac = *this;
  ac.Array<char>::operator+=('\0');
  ac[ac.dim() - 2] = c;
  return *this;
}

bool Mystring::operator==(const char * s2) const
{
  if (!strcmp(this->c_str(), s2)) return true;
  else return false;
}

bool Mystring::operator!=(const char * s2) const
{
  if (strcmp(this->c_str(), s2)) return true;
  else return false;
}

bool Mystring::operator==(const Mystring & s2) const
{
  if (!strcmp(this->c_str(), s2.c_str())) return true;
  else return false;
}

Mystring & Mystring::remove(const char c)
{
  Mystring temp = "";
  Mystring & t = *this;
  uint i;
  for (i=0; i < length(); ++i)
    if (t[i] != c) temp += t[i];
  *this = temp;
  return *this;
}

Array<Mystring> Mystring::split(char delim) const
{
  const Mystring & t = *this;
  Array<Mystring> answer (0);
  uint k;
  Mystring nextword;
  for (k=0; k < length(); ++k)
    {
      nextword = "";
      while ((t[k] != delim) && (t[k] != '\0'))
	{
	  nextword += t[k];
	  ++ k;
	}
      if (nextword.length()) answer += nextword;
    }
  return answer;
}

Mystring Mystring::at(uint start, uint len) const
{
  const Mystring & t = *this;
  Mystring answer = "";
  uint k;
  for (k=start; k < start + len; ++k)
    answer += t[k];
  return answer;
}

Mystring Mystring::after(const uint & start) const
{
  return this->at(start + 1, length() - start - 1);
}

Mystring Mystring::after(const char & c) const
{
  const Mystring & t = *this;
  uint i = 0;
  while ((i < (length()-1)) && (t[i] != c)) ++i;
  return this->at(i + 1, length() - i - 1);
}

Mystring Mystring::upcase() const
{
  Mystring answer = *this;
  uint i;
  for (i=0; i < length(); ++ i)
    answer[i] = toupper((int) answer[i]);
  return answer;
}

Mystring Mystring::downcase() const
{
  Mystring answer = *this;
  uint i;
  for (i=0; i < length(); ++ i)
    answer[i] = tolower((int) answer[i]);
  return answer;
}

Mystring upcase(const Mystring & s)
{
  return s.upcase();
}

ostream & operator<<(ostream & os, const Mystring & s)
{
  os << s.c_str();
  return os;
}

istream & operator>>(istream & is, Mystring & s)
{
  s = "";
  char c;
  is.get(c);
  while (is.good() && c != '\n')
    {
      s += c;
      is.get(c);
    }
  return is;
}
