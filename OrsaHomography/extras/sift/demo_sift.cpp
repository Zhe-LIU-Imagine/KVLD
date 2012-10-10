// WARNING: 
// This file implements an algorithm possibly linked to the patent
//
// David Lowe  "Method and apparatus for identifying scale invariant 
// features in an image and use of same for locating an object in an 
// image",  U.S. Patent 6,711,293.
//
// This file is made available for the exclusive aim of serving as
// scientific tool to verify of the soundness and
// completeness of the algorithm description. Compilation,
// execution and redistribution of this file may violate exclusive
// patents rights in certain countries.
// The situation being different for every country and changing
// over time, it is your responsibility to determine which patent
// rights restrictions apply to you before you compile, use,
// modify, or redistribute this file. A patent lawyer is qualified
// to make this determination.
// If and only if they don't conflict with any patent terms, you
// can benefit from the following license terms attached to this
// file.
//
// This program is provided for scientific and educational only:
// you can use and/or modify it for these purposes, but you are
// not allowed to redistribute this work or derivative works in
// source or executable form. A license must be obtained from the
// patent right holders for any other use.

#include "demo_lib_sift.h"
#include "library.h"

#include "libImage/image.hpp"
#include "libImage/image_io.hpp"
#include "libImage/image_converter.hpp"
#include "libImage/sample.hpp"
#include "libImage/image_drawing.hpp"
#include "libImage/image_concat.hpp"

#include <iostream>
#include <fstream>

using namespace libs;

int main(int argc, char **argv)
{	
    if(argc != 4 && argc != 5 && argc != 6) {
        std::cerr << "Usage: " << argv[0] << " imgIn imgIn2 fileMatchesOut [imgMatchesHorizontal] [imgMatchesVertical]"
                  << std::endl;
        return 1;
    }

	//////////////////////////////////////////////// Input

  Image<RGBColor> image1, image2;
  if(! ReadImage(argv[1], &image1))
    return 1;
  if(! ReadImage(argv[2], &image2))
    return 1;

  Image<unsigned char> image1Gray, image2Gray;
  libs::convertImage(image1, &image1Gray);
  libs::convertImage(image2, &image2Gray);

	int w1 = image1Gray.Width(), h1 = image1Gray.Height();
  int w2 = image2Gray.Width(), h2 = image2Gray.Height();

  //Convert images to float
  Image<float> If1;
  convertImage(image1Gray , &If1);
  Image<float> If2;
  convertImage(image2Gray , &If2);
  
	///////////////////////////////////////// Applying Sift
	siftPar siftparameters;
	default_sift_parameters(siftparameters);
    siftparameters.DoubleImSize=0;

	keypointslist keyp1, keyp2;
	compute_sift_keypoints(If1.data(), keyp1, w1, h1, siftparameters);
    std::cout<< "sift:: 1st image: " << keyp1.size() << " keypoints"<<std::endl;
	compute_sift_keypoints(If2.data(), keyp2, w2, h2, siftparameters);
    std::cout<< "sift:: 2nd image: " << keyp2.size() << " keypoints"<<std::endl;

	matchingslist matchings;
	compute_sift_matches(keyp1,keyp2,matchings,siftparameters);
    std::cout << "sift:: matches: " << matchings.size() <<std::endl;

	//////////////////////////////////////////////////////////////// Save file with matches
  Match::saveMatch(argv[3], matchings);

	//////////////////////////////////////////////// Output image containing line matches
  if( argc >= 5 )
  {

    Image<unsigned char> concat;
    ConcatH( image1Gray, image2Gray, concat );
    //-- Draw link between features :
    matchingslist::iterator ptr = matchings.begin();
    for(; ptr != matchings.end(); ++ptr)
    {
      libs::DrawLine( (int) ptr->x1, (int) ptr->y1, w1 + (int) ptr->x2, (int) ptr->y2, (unsigned char)(255), &concat);
    }
    std::ostringstream os;
    os << argv[4];
    libs::WriteImage(os.str().c_str(), concat);
	}

  if( argc >= 6 )
  {

    Image<unsigned char> concat;
    ConcatV( image1Gray, image2Gray, concat );
    //-- Draw link between features :
    matchingslist::iterator ptr = matchings.begin();
    for(; ptr != matchings.end(); ++ptr)
    {
      libs::DrawLine( (int) ptr->x1, (int) ptr->y1, (int) ptr->x2, h1 + (int) ptr->y2, (unsigned char)(255), &concat);
    }
    std::ostringstream os;
    os << argv[5];
    libs::WriteImage(os.str().c_str(), concat);
  }

  return 0;
}
