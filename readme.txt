                      KVLD Library
               Zhe Liu and Pierre Moulon

ABOUT
  The KVLD open source library implements the Kth virtual line descriptor matching method introduced in the paper: 

Z,Liu and R,Marlet. Virtual Line Descriptor and Semi-Local Matching Method for Reliable Feature Correspondence in BMVC 2012

  KVLD is distributed under the BSD license (see the COPYING file).

INSTALLING
  This implementation is on C++ and depends on openCV library, whose installation guild is available online. You will need Cmake 2.6 or later to compile the program.
  However the kvld library is independent of openCV, and if you develop your own applications using kvld, please 
  *include folders as OrsaHomography and kvld 
  *modify the CMakeLists.txt file to remove openCV library (see CMakeLists.txt)
  *include kvld in your code by adding #include "kvld/kvld.h"

FOLDERS:
  Kvld: containing all KVLD algorithm, some structures depend on
   OrsaHomography library, so please include both of them to make KVLD running.
  OrsaHomography: containing ORSA algorithm implemented by Pierre Moulon. It also offers basic structures for KVLD algorithm.
  demo_image: some illustrating pairs of images for different applications.
  demo_output: results of demos are sent here, including 
    (? means the index)
    * IMG_?_Detectors1: detectors in the first image
    * IMG_?_Detectors2: detectors in the second image
    * IMG_?_initial_matches: input matches of KVLD(pairs of indices)
    * IMG_?_kvld_matches: output matches of KVLD(pairs of indices)
    * IMG_?_kvld_filtered: visual output matches of KVLD
    * IMG_?_initial: visual input matches of KVLD
    * (optional)IMG_?_matrix: fundametal or homography matrix

APPLICATIONS
  The code contains three main applications (KVLD_Deformable, KVLD_Calibrate and KVLD_Interface). They come with demonstration using images in demo_image folder. All the outputs are sent to the demo_output folder.

KVLD_Deformable: 
  To establish correspondence of points in a deformable object, KVLD can efficiently retrieve correct matches among a large number of false matches.

KVLD_Calibration: 
  For image calibration, KVLD generates filtered matches, then ORSA (Ransac variant) estimate homography or fundamental matrix using those (less contaminated) matches.

KVLD_Interface:
  This allow you to feed KVLD with your own detectors and initial matches with the original images and let KVLD generates filtered matches. You need to provide the following three files.
 * IMG_?_Detectors1: detectors in the first image
 * IMG_?_Detectors2: detectors in the second image
 * IMG_?_initial_matches: input matches of KVLD(pairs of indices)

For more information, please contact zhe.liu@enpc.fr



