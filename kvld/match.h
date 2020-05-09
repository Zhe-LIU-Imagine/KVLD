// Copyright (c) 2010 Pascal Monasse
// Copyright (c) 2011 Pierre Moulon: encapsulation
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

#ifndef MATCH_H
#define MATCH_H

#include <vector>
#include <fstream>
#include <sstream>

/// Save matching position between two points.
struct Match
{
    Match() {}
    Match(float ix1, float iy1, float ix2, float iy2)
    : x1(ix1), y1(iy1), x2(ix2), y2(iy2) {}
    float x1, y1, x2, y2;

    /**
    * Load the corresponding matches from file.
    * \param nameFile   The file where matches were saved.
    * \param vec_match  The loaded corresponding points.
    * \return bool      True if everything was ok, otherwise false.
    */
    static bool loadMatch(const char* nameFile, std::vector<Match>& vec_match)
    {
      vec_match.clear();
      std::ifstream f(nameFile);
      while( f.good() ) {
        std::string str;
        std::getline(f, str);
        if( f.good() ) {
          std::istringstream s(str);
          Match m;
          s >> m;
          if(!s.fail() )
            vec_match.push_back(m);
        }
      }
      return !vec_match.empty();
    }
    
    /**
    * Save the corresponding matches to file.
    * \param nameFile   The file where matches will be saved.
    * \param vec_match  The matches that we want to save.
    * \return bool True if everything was ok, otherwise false.
    */
    static bool saveMatch(const char* nameFile, const std::vector<Match>& vec_match)
    {
      std::ofstream f(nameFile);
      if( f.is_open() ) {
        std::vector<Match>::const_iterator it = vec_match.begin();
        for(; it != vec_match.end(); ++it)
          f << *it;
      }
      return f.is_open();
    }

    /// Lexicographical ordering of matches. Used to remove duplicates.
    friend bool operator<(const Match& m1, const Match& m2)
    {
      if(m1.x1 < m2.x1) return true;
      if(m1.x1 > m2.x1) return false;

      if(m1.y1 < m2.y1) return true;
      if(m1.y1 > m2.y1) return false;

      if(m1.x2 < m2.x2) return true;
      if(m1.x2 > m2.x2) return false;

      return (m1.y2 < m2.y2);
    }

    /// Comparison Operator
    friend bool operator==(const Match& m1, const Match& m2)
    {
      return (m1.x1==m2.x1 && m1.y1==m2.y1 &&
              m1.x2==m2.x2 && m1.y2==m2.y2);
    }

    friend std::ostream& operator<<(std::ostream& os, const Match & m)
    {
      return os << m.x1 << " " << m.y1 << " "
        << m.x2 << " " << m.y2 << std::endl;
    }

    friend std::istream& operator>>(std::istream & in, Match & m)
    {
      return in >> m.x1 >> m.y1 >> m.x2 >> m.y2;
    }
};

#endif // MATCH_H
