/** @Converting functions between structures of openCV and KVLD suitable structures
 ** @author Zhe Liu
 **/

/*
Copyright (C) 2011-12 Zhe Liu and Pierre Moulon.
All rights reserved.

This file is part of the KVLD library and is made available under
the terms of the BSD license (see the COPYING file).
*/
#ifndef CONVERT_H
#define CONVERT_H
#include <vector>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"


//#include "demo/libImage/image.hpp"

#include "kvld/algorithm.h"

typedef std::pair<size_t,size_t> Pair;

int Convert_image(const cv::Mat& In, Image<float> & imag);//convert only gray scale image of opencv

int Convert_detectors(const  std::vector<cv::KeyPoint>& feat1,std::vector<keypoint>& F1);//convert openCV detectors to KVLD suitable detectors
int Convert_detectors(const  std::vector<keypoint>& F1,std::vector<cv::KeyPoint>& feat1);

int Convert_matches(const std::vector<cv::DMatch>& matches, std::vector<Pair>& matchesPair);

int read_detectors(const std::string& filename ,  std::vector<keypoint>& feat);//reading openCV style detectors


int read_matches(const std::string& filename , std::vector<cv::DMatch>& matches);//reading openCV style matches

#endif //CONVERT_H