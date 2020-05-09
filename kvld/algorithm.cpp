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

float getRange(const Image<float>& I,int a,const float p, const float ratio){
  float range=ratio*sqrt(float(3*I.Height()*I.Width())/(p*a*PI));
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
