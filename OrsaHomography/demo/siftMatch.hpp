//Copyright (C) 2011 Lionel Moisan, Pascal Monasse, Pierre Moulon
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef LIBS_SIFT_MATCH_H
#define LIBS_SIFT_MATCH_H

#include "extras/sift/demo_lib_sift.h"

/// SIFT matches
static void SIFT(const Image<unsigned char> &im1,
                 const Image<unsigned char> &im2,
                 std::vector<Match>& vec_matchings,
                 float fMatchRatio=0.6f) {
    //Convert images to float
    Image<float> If1, If2;
    libs::convertImage(im1, &If1);
    libs::convertImage(im2, &If2);

    siftPar param;
    default_sift_parameters(param);
    param.MatchRatio = fMatchRatio;
    param.DoubleImSize=0;

    keypointslist keyp1, keyp2;
    compute_sift_keypoints(If1.data(), keyp1, If1.Width(), If1.Height(), param);
    std::cout<< "sift:: 1st image: " << keyp1.size() << " keypoints"<<std::endl;
    compute_sift_keypoints(If2.data(), keyp2, If2.Width(), If2.Height(), param);
    std::cout<< "sift:: 2nd image: " << keyp2.size() << " keypoints"<<std::endl;

    // Find putatives matches
    compute_sift_matches(keyp1, keyp2, vec_matchings, param);
    std::cout << "sift:: matches: " << vec_matchings.size() <<std::endl;
}

/// Remove multiple "same position" matches
static void rm_duplicates(std::vector<Match>& m) {
  std::sort(m.begin(), m.end());
  std::vector<Match>::iterator end = std::unique(m.begin(), m.end());
  if(end != m.end()) {
    std::cout << "Remove " << std::distance(end, m.end())
      << "/" << m.size() << " duplicate matches, "
      << " keeping " << std::distance(m.begin(), end) <<std::endl;
    m.erase(end, m.end());
  }
}

#endif
