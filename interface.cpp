/** @KVLD interface with external detectors and matches
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
#include <vector>

#include "opencv2/opencv.hpp"
#include "opencv2/features2d.hpp"

#include "opencv2/calib3d.hpp"

#include "opencv2/core.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
//#include "opencv2/nonfree/features2d.hpp" 


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
	std::string output=std::string(SOURCE_DIR)+"/demo_output/IMG_"+index+"_";

	image1= cv::imread(input+index+".jpg", cv::IMREAD_GRAYSCALE);
	image2= cv::imread(input+index+"bis.jpg", cv::IMREAD_GRAYSCALE);

	cv::Mat image1color, image2color, concat;//for visualization
	image1color= cv::imread(input+index+".jpg", cv::IMREAD_COLOR);
	image2color= cv::imread(input+index+"bis.jpg", cv::IMREAD_COLOR);

	//=============== Read SIFT points =================//
	std::cout<<"Loading SIFT features"<<std::endl;
	std::vector<keypoint> F1, F2;
	read_detectors(output+"Detectors1.txt",F1);
	read_detectors(output+"Detectors2.txt",F2);

	std::vector<cv::KeyPoint> feat1,feat2;
	Convert_detectors(F1,feat1);//we only need detectors for KVLD
	Convert_detectors(F2,feat2);//we only need detectors for KVLD

	std::cout<< "sift:: 1st image: " << F1.size() << " keypoints"<<std::endl;
	std::cout<< "sift:: 2nd image: " << F2.size() << " keypoints"<<std::endl;
	//=============== compute matches using brute force matching ====================//
	std::vector<cv::DMatch> matches;
	read_matches(output+"initial_matches.txt",matches);
	//=============== convert openCV sturctures to KVLD recognized elements
	Image<float> If1, If2;
	Convert_image(image1, If1);
	Convert_image(image2, If2);

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

	while (it_num < 5 && kvldparameters.inlierRate>KVLD(If1, If2,F1,F2, matchesPair, matchesFiltered, vec_score,E,valide,kvldparameters)) 
	{
		kvldparameters.inlierRate/=2;
		kvldparameters.rang_ratio=sqrt(2.0f);
		std::cout<<"low inlier rate, re-select matches with new rate="<<kvldparameters.inlierRate<<std::endl;
		if (matchesFiltered.size()==0) kvldparameters.K=2;
		it_num++;
	}
	std::cout<<"K-VLD filter ends with "<<matchesFiltered.size()<<" selected matches"<<std::endl;

	//================= write files to output folder ==================//
	std::cout<<"Please check the output folder for results"<<std::endl;
	writeResult(output,F1, F2, matchesPair, matchesFiltered, vec_score);
//================= Visualize matching result ====================//
	//  for (int it1=0; it1<matchesPair.size()-1;it1++){
	//  for (int it2=it1+1; it2<matchesPair.size();it2++){
	//    if (valide[it1] && valide[it2] && E(it1,it2)>=0){

	//      cv::KeyPoint l1 = feat1[matchesPair[it1].first];
	//      cv::KeyPoint l2 = feat1[matchesPair[it2].first];

	//      cv::KeyPoint r1 = feat2[matchesPair[it1].second];
	//      cv::KeyPoint r2 = feat2[matchesPair[it2].second];

	//      //cv::line(image1color,l1.pt, l2.pt,cv::Scalar(255,0,255),2);
	//      cv::circle(image1color,l1.pt,5,cv::Scalar(0,255,0));
	//    }
	//  }
	//}
	//cv::imwrite(output+"output.png",image1color);



	cv::vconcat(image1color, image2color,concat);
	for(  std::vector<cv::DMatch>::const_iterator ptr = matches.begin(); ptr != matches.end(); ++ptr)
	{
		cv::KeyPoint start = feat1[ptr->queryIdx]; 
		cv::KeyPoint end = feat2[ptr->trainIdx];
		cv::line( concat,start.pt, end.pt+cv::Point2f(0,image1.rows),cv::Scalar(0, 255,0 ));

	}
	cv::imwrite(output+"initial.png",concat);

	//========== KVLD result =============//
	cv::vconcat(image1color, image2color,concat);

	//draw gvld-consistant neighbors (not exhostive)
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
	for( std::vector<Pair >::const_iterator ptr = matchesFiltered.begin(); ptr != matchesFiltered.end(); ++ptr)
	{
		size_t i = ptr->first;
		size_t j = ptr->second;
		cv::KeyPoint start = feat1[i]; 
		cv::KeyPoint end = feat2[j];

		cv::line( concat,start.pt, end.pt+cv::Point2f(0,image1.rows),cv::Scalar(0, 255, 0 ));

	}
	cv::imwrite(output+"kvld_filtered.png",concat);

	return 0;
}
