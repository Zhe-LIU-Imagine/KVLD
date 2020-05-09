/** @Image calibration application
 ** @Estimate fundamental or homography matrix
 ** @author Zhe Liu
 **/

/*
Copyright (C) 2011-12 Zhe Liu and Pierre Moulon.
All rights reserved.

This file is part of the KVLD library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#include <algorithm>
#include <memory>

#include "kvld/kvld.h"
#include "convert.h"

#include "opencv2/opencv.hpp"
#include "opencv2/features2d.hpp"

#include "opencv2/calib3d.hpp"

#include "opencv2/core.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"

#include "opencv2/xfeatures2d.hpp"
#include "opencv2/xfeatures2d/nonfree.hpp"

using namespace cv::xfeatures2d;

const float sift_matching_criterion=0.98;


int main(int argc,char*argv[]) {
	std::cout<<"Warming: K-VLD may suffer performance degradation under Linux OS!"<<std::endl
		<<"Please first check existing result in the output folder as a reference!"<<std::endl;

	//================= load images ======================//
	cv::Mat image1, image2;
	int imageID=2;// index of a pair of images you want to use in folder demo_images
	std::string index;
	std::stringstream f;
	f<<imageID;
	f>>index;
	std::string input=std::string(SOURCE_DIR)+"/demo_image/IMG_";
	image1= cv::imread(input+index+".jpg", cv::IMREAD_GRAYSCALE);
	image2= cv::imread(input+index+"bis.jpg", cv::IMREAD_GRAYSCALE);

    cv::Mat image1color, image2color, concat;//for visualization
	image1color= cv::imread(input+index+".jpg", cv::IMREAD_COLOR);
	image2color= cv::imread(input+index+"bis.jpg", cv::IMREAD_COLOR);

	//=============== compute SIFT points =================//
	std::cout<<"Extracting SIFT features"<<std::endl;
	std::vector<cv::KeyPoint> feat1,feat2;

	//cv::SiftFeatureDetector* detectortype=new  cv::SiftFeatureDetector() ;
	//cv::PyramidAdaptedFeatureDetector detector2(detectortype,3);// 3 levels of image scale
	cv::Ptr<SIFT> detector = SIFT::create(0.5);//default setting is ok, 5 levels generate too much features
	cv::Ptr<SiftDescriptorExtractor> extractor = SiftDescriptorExtractor::create();
	cv::Mat descriptors1,descriptors2;

	detector->detect(image1,feat1);
	extractor->compute(image1,feat1,descriptors1);
	std::cout<< "sift:: 1st image: " << feat1.size() << " keypoints"<<std::endl;

	detector->detect(image2,feat2);
	extractor->compute(image2,feat2,descriptors2);
	std::cout<< "sift:: 2nd image: " << feat2.size() << " keypoints"<<std::endl;

	//=============== compute matches using brute force matching ====================//
	std::vector<cv::DMatch> matches;
	bool bSymmetricMatches = false;//caution, activate this with knn matching will cause errors.
	cv::BFMatcher matcher(cv::NORM_L2, bSymmetricMatches);
	if (bSymmetricMatches){
		matcher.match(descriptors1,descriptors2,matches);
	}
	else
	{
		std::vector<std::vector<cv::DMatch>> knnmatches;
		matcher.knnMatch(descriptors1,descriptors2,knnmatches,2);
		for (std::vector<std::vector<cv::DMatch>>::const_iterator it=knnmatches.begin();it!=knnmatches.end();it++){
			if (it->at(0).distance<sift_matching_criterion*it->at(1).distance) 
				matches.push_back((*it)[0]);
		}
	}
	//=============== convert openCV sturctures to KVLD recognized elements
	Image<float> If1, If2;
	Convert_image(image1, If1);
	Convert_image(image2, If2);

	std::vector<keypoint> F1, F2;
	Convert_detectors(feat1,F1);//we only need detectors for KVLD
	Convert_detectors(feat2,F2);//we only need detectors for KVLD
	std::vector<Pair> matchesPair;
	Convert_matches(matches,matchesPair);

	//===============================  KVLD method ==================================//
	std::cout<<"K-VLD starts with "<<matches.size()<<" matches"<<std::endl;

	std::vector<Pair> matchesFiltered;
	std::vector<double> vec_score;

	//In order to illustrate the gvld(or vld)-consistant neighbors, the following two parameters has been externalized as inputs of the function KVLD.
	//libNumerics::matrix<float> E = libNumerics::matrix<float>::ones(matches.size(),matches.size())*(-1);
	std::vector<std::vector<float>> E(matches.size(), std::vector<float>(matches.size(), -1));
	// gvld-consistency matrix, intitialized to -1,  >0 consistency value, -1=unknow, -2=false  

	std::vector<bool> valide(matches.size(), true);// indices of match in the initial matches, if true at the end of KVLD, a match is kept.

	size_t it_num=0;
	KvldParameters kvldparameters;//initial parameters of KVLD

	while (it_num < 5 && kvldparameters.inlierRate>KVLD(If1, If2,F1,F2, matchesPair, matchesFiltered, vec_score,E,valide,kvldparameters)) {
		kvldparameters.inlierRate/=2;
		kvldparameters.rang_ratio=sqrt(2.0f);
		std::cout<<"low inlier rate, re-select matches with new rate="<<kvldparameters.inlierRate<<std::endl;
		if (matchesFiltered.size()==0) kvldparameters.K=2;
		it_num++;
	}
	std::cout<<"K-VLD filter ends with "<<matchesFiltered.size()<<" selected matches"<<std::endl;
	//====================fundamental matrix verification================//
	cv::Mat points1, points2;
	for (auto & p1p2 : matchesFiltered) {
		int id1 = p1p2.first;
		int id2 = p1p2.second;
		cv::Mat M1 = (cv::Mat_<float>(1, 2) << feat1[id1].pt.x, feat1[id1].pt.y);
		cv::Mat M2 = (cv::Mat_<float>(1, 2) << feat2[id2].pt.x, feat2[id2].pt.y);
		points1.push_back(M1);
		points2.push_back(M2);

	}
	if (points1.empty()) {
		std::cout << "kep points empty " << std::endl;
		return false;
	}
	// relation estimation 
	std::vector<uchar> mask;
	cv::Mat matrix = cv::findFundamentalMat(points1, points2, cv::FM_RANSAC, 3.0, 0.98, mask);

	//================= write files to output folder ==================//
	std::cout<<"Writing results to the output folder..."<<std::endl;
	std::string output=std::string(SOURCE_DIR)+"/demo_output/IMG_"+index+"_";
	writeResult(output,F1, F2, matchesPair, matchesFiltered, vec_score);

	std::ofstream matrix_str((output+"matrix.txt").c_str());
	matrix_str<<matrix;
	matrix_str.close();

	//================= Visualize matching result ====================//
	cv::vconcat(image1color, image2color,concat);

	for (auto & p1p2 : matches) {
		cv::KeyPoint start = feat1[p1p2.queryIdx];
		cv::KeyPoint end = feat2[p1p2.trainIdx];
		cv::line(concat, start.pt, end.pt + cv::Point2f(0, image1.rows), cv::Scalar(255, 0, 0));	
	}

	cv::imwrite(output+"initial.png",concat);

	//========== KVLD result =============//
	cv::vconcat(image1color, image2color,concat);

	//draw gvld-consistant neighbors (not exhostive), may include outliers rejected by ORSA
	for (int it1=0; it1<matchesPair.size()-1;it1++){
		for (int it2=it1+1; it2<matchesPair.size();it2++){
			if (valide[it1] && valide[it2] && E[it1][it2]>=0 ){

				cv::KeyPoint l1 = feat1[matchesPair[it1].first];
				cv::KeyPoint l2 = feat1[matchesPair[it2].first];

				cv::KeyPoint r1 = feat2[matchesPair[it1].second];
				cv::KeyPoint r2 = feat2[matchesPair[it2].second];

				cv::line(concat,l1.pt, l2.pt,cv::Scalar(255,0,255),2);
				cv::line(concat,r1.pt+cv::Point2f(0,image1.rows), r2.pt+cv::Point2f(0,image1.rows),cv::Scalar(255,0,255),2);

			}
		}
	}
	for (int i = 0; i < mask.size(); i++) {
		cv::KeyPoint start = feat1[matchesFiltered[i].first];
		cv::KeyPoint end = feat2[matchesFiltered[i].second];
		if (mask[i] > 0) {
			cv::line(concat, start.pt, end.pt + cv::Point2f(0, image1.rows), cv::Scalar(0, 255, 0));
		}
		else {
			cv::line(concat, start.pt, end.pt + cv::Point2f(0, image1.rows), cv::Scalar(255, 0, 0));
		}
	}

	cv::imwrite(output+"kvld_filtered.png",concat);
	std::cout<<"Please check the output folder for results."<<std::endl;
	return 0;
}
