/** @basic structures implementation
 ** @author Zhe Liu
 **/

/*
Copyright (C) 2011-12 Zhe Liu and Pierre Moulon.
All rights reserved.

This file is part of the KVLD library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#include "algorithm.h"



//template<typename T>
IntegralImages::IntegralImages(const Image<float>& I){
		map.Resize(I.Width()+1,I.Height()+1);
		map.fill(0);
		for (int y=0;y<I.Height();y++)
			for (int x=0;x<I.Width();x++){
				map(y+1,x+1)=double(I(y,x))+map(y,x+1)+map(y+1,x)-map(y,x);
			}
	}

float getRange(const Image<float>& I,int a,const float p){
  float range=sqrt(float(3*I.Height()*I.Width())/(p*a*PI));
  std::cout<<"range ="<<range<<std::endl;
  return range;
}

///================================= ORSA functions ======================================//
static void display_stats(const std::vector<Match>& vec_matchings,
  const std::vector<size_t>& vec_inliers,
   libNumerics::matrix<double>& F, bool homography) {
     double l2=0, linf=0;

     if (homography){
       std::vector<size_t>::const_iterator it=vec_inliers.begin();
       for(; it!=vec_inliers.end(); ++it) {
         const Match& m=vec_matchings[*it];
         Matrix x(3,1);
         x(0)=m.x1;
         x(1)=m.y1;
         x(2)=1.0;
         x = F*x;
         x /= x(2);
         double e = (m.x2-x(0))*(m.x2-x(0)) + (m.y2-x(1))*(m.y2-x(1));
         l2 += e;
         if(linf < e)
           linf = e;
       }
     }else{
       std::vector<size_t>::const_iterator it=vec_inliers.begin();

       for(; it!=vec_inliers.end(); ++it) {
         const Match& m=vec_matchings[*it];
         double a = F(0,0) * m.x1 + F(0,1) * m.y1 + F(0,2);
         double b = F(1,0) * m.x1 + F(1,1) * m.y1 + F(1,2);
         double c = F(2,0) * m.x1 + F(2,1) * m.y1 + F(2,2);
         double d = a*m.x2 + b*m.y2 + c;
         // double e =  (d*d) / (a*a + b*b);
         double e =  (d*d) / (a*a + b*b);
         l2 += e;
         if(linf < e)
           linf = e;
       }
     }
     std::cout << "Average/max error: "
       << sqrt(l2/vec_inliers.size()) << "/"
       << sqrt(linf) <<std::endl;
}

bool ORSA(const std::vector<Match>& vec_matchings, int w1,int h1, int w2,int h2,
          double& precision,
          libNumerics::matrix<double>& H, std::vector<size_t>& vec_inliers,bool homo)
{
  const size_t n = vec_matchings.size();
  if(n < 5)
  {
      std::cerr << "Error: ORSA needs 5 matches or more to proceed" <<std::endl;
      return false;
  }
  libNumerics::matrix<double> xA(2,n), xB(2,n);

  for (size_t i=0; i < n; ++i)
  {
    xA(0,i) = vec_matchings[i].x1;
    xA(1,i) = vec_matchings[i].y1;
    xB(0,i) = vec_matchings[i].x2;
    xB(1,i) = vec_matchings[i].y2;
  }

  std::auto_ptr< orsa::OrsaModel > modelEstimator;
  if(homo){
	  modelEstimator = std::auto_ptr< orsa::HomographyModel  >(new orsa::HomographyModel(xA, w1, h1, xB, w2, h2, true));
  }else{
	  //Fundamental
	  modelEstimator = std::auto_ptr< orsa::FundamentalModel >(new orsa::FundamentalModel(xA, w1, h1, xB, w2, h2, true));
  }
  if(modelEstimator->orsa(vec_inliers, 2000, &precision, &H, false)>0.0)
	  return false;

  //std::cout << "Before refinement: ";
	//display_stats(vec_matchings, vec_inliers, H,homo);
 
  if( modelEstimator->ComputeModel(vec_inliers,&H) ) // Re-estimate with all inliers
  {
	  std::cout << "After  refinement: ";
    display_stats(vec_matchings, vec_inliers, H,homo);
  } else
	  std::cerr << "Warning: error in refinement, result is suspect" <<std::endl;
  return true;
}

FCrit Find_Model(const Image<float>& I1,const Image<float>& I2,
	const std::vector<keypoint>& F1,const std::vector<keypoint>& F2,
			const std::vector<Pair>& matches,double& precision, bool homograghy){
				//precision=0 default threshold search, else put given precision ex. 2.0
		std::cout<<"==========Runing Orsa==========="<<std::endl;

	std::vector<Match> vec_matches;
  for (std::vector<Pair>::const_iterator it=matches.begin();it!=matches.end();it++){
    vec_matches.push_back(Match(F1[it->first].x,F1[it->first].y,F2[it->second].x,F2[it->second].y));
  }
	rm_duplicates(vec_matches);

	libNumerics::matrix<double> H(3,3);
	std::vector<size_t> vec_inliers;


	bool bRes= ORSA(vec_matches, I1.Width(), I1.Height(), I2.Width(), I2.Height(),precision, H, vec_inliers,homograghy);

	std::cout << std::endl << "Orsa estimation :\n"
		<< "putative: " << vec_matches.size() <<std::endl
		<< "validated:" << vec_inliers.size()<<std::endl;

	if (homograghy)
		H/=H(2,2);
	return FCrit (H,precision,homograghy);
}

template <typename T>
static void rm_duplicates(T& m) {
  std::sort(m.begin(), m.end());
  typename T::iterator end = std::unique(m.begin(), m.end());
  if(end != m.end()) {
    std::cout << "Remove " << std::distance(end, m.end())
      << "/" << m.size() << " duplicate matches, "
      << " keeping " << std::distance(m.begin(), end) <<std::endl;
    m.erase(end, m.end());
  }
}

//=============================IO interface, convertion of object types======================//

std::ofstream& writeDetector(std::ofstream& out, const keypoint& feature){
  out<<feature.x<<" "<<feature.y<<" "<<feature.scale<<" "<<feature.angle<<std::endl;
  /*for(int i=0;i<128;i++)  
    out<<feature.vec[i]<<" ";
  out<<std::endl;*/
return out;
}

std::ifstream& readDetector(std::ifstream& in,keypoint& point){
  in>>point.x>>point.y>>point.scale>>point.angle;
  //for(int i=0;i<128;i++)  {
  //  in>>point.vec[i];
  //}
return in;
}
