#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <CImg.h>
#include <assert.h>
#include <cfloat>

using namespace cimg_library;
using namespace std;


const int win_size = 1;
const int classes = 100;
const double alpha = 1e0; 

double sqr(double a) { return a*a; }
double msg_cost(int i, int j, int c, int d, CImg<double>& total_cost){
    double cost = 0;
    if(d == 0){
        cost = total_cost(i+1, j, c) + total_cost(i, j-1, c) + total_cost(i, j+1, c);
    } else if(d == 1) {
        cost = total_cost(i-1, j, c) + total_cost(i, j-1, c) + total_cost(i, j+1, c);        
    } else if(d == 2) {
        cost = total_cost(i-1, j, c) + total_cost(i+1, j, c) + total_cost(i, j+1, c);        
    } else if(d == 3) {
        cost = total_cost(i-1, j, c) + total_cost(i+1, j, c) + total_cost(i, j-1, c);        
    }
    return cost;
}


double msg_cost(int i, int j, int c, int d, CImg<double>& total_cost, bool vertical){
    double cost = 0;
    
    /* horizontal spanning tree */
    if(d == 0 && !vertical) {
        cost = total_cost(i+1, j, c);
    }
    if(d == 1 && !vertical) {
        cost = total_cost(i-1, j, c);
    }

    /* vertical spanning tree */
    if(d == 0 && vertical) {
        cost = total_cost(i, j+1, c);
    }
    if(d == 1 && vertical) {
        cost = total_cost(i, j-1, c);
    }

    return cost;
}

CImg<double> mrf_stereo(const CImg<double> &img1, const CImg<double> &img2){
  //CImg<double> result(img1.width(), img1.height());
  CImg<double> result(img1.height(), img1.width());
  CImg<double> total_cost(img1.width(), img1.height(), classes);
  //CImg<double> total_cost(img1.height(), img1.width(), classes);
  CImg<double> d_cost(img1.width(), img1.height(), classes);
  //CImg<double> d_cost(img1.height(), img1.width(), classes);
  int iterations = 1;

  //cout << "MRF calculating d_cost..." << endl;
  for(int i=0; i<img1.width(); i++){
      for(int j=0; j<img1.height(); j++){
          for(int c=0; c<classes;c++){
              float cost = 0;
              for(int p = max(i-win_size, 0); p <= min(i+win_size, img1.height()-1); p++){
                for(int q = max(j-win_size, 0); q <= min(j+win_size, img1.width()-1); q++){
                  cost += sqr(img1(min(q+c, img1.width()-1), p) - img2(q, p));
                }
              }
              d_cost(i, j, c) = cost;
          }
      }
  }

  //cout << "loopy BP MRF iterating..." << endl;
  while(iterations--){
      for(int i=1; i<img1.width()-1; i++){
          for(int j=1; j<img1.height()-1; j++){
              for(int c=0; c<classes;c++){
                double tmp_cost = 0; 
                double min_cost = 1e7;
                for(int c1=0; c1<classes;c1++){
                    for(int d=0; d<4; d++){
                    tmp_cost = d_cost(i, j, c1);
                    tmp_cost += alpha*(c1 != c ? 1 : 0);
                        tmp_cost += msg_cost(i, j, c1, d, total_cost);
                    if(tmp_cost < min_cost){
                        min_cost = tmp_cost;
                    }
                    }
                }
                total_cost(i, j, c) = min_cost;
              }
          }
      }
  }
  //cout << "loop BP MRF predicting..."<< endl;

  for(int i=1; i<img1.width()-1; i++){
      for(int j=1; j<img1.height()-1; j++){
          double tmp_cost = 0;
          double min_cost = 1e7;
          int idx = -1;
          for(int c=0; c<classes;c++){
              tmp_cost = d_cost(i, j, c);
              tmp_cost += total_cost(i-1, j, c) + total_cost(i+1, j, c) + total_cost(i, j+1, c) + total_cost(i, j-1, c);
              if(tmp_cost < min_cost){
                  min_cost = tmp_cost;
                  idx = c;
              }
          }
          result(j, i) = idx;
          //result(i, j) = idx;
      }
  }
  return result;
}


CImg<double> TRWS_stereo(const CImg<double> &img1, const CImg<double> &img2){
  CImg<double> result(img1.height(), img1.width());
  //CImg<double> result(img1.width(), img1.height());
  CImg<double> total_cost(img1.width(), img1.height(), classes);
  CImg<double> d_cost(img1.width(), img1.height(), classes);
  int iterations = 1;

  //cout << "MRF calculating d_cost..." << endl;
  for(int i=0; i<img1.width(); i++){
      for(int j=0; j<img1.height(); j++){
          for(int c=0; c<classes;c++){
              float cost = 0;
              for(int p = max(i-win_size, 0); p <= min(i+win_size, img1.height()-1); p++){
                for(int q = max(j-win_size, 0); q <= min(j+win_size, img1.width()-1); q++){
                  cost += sqr(img1(min(q+c, img1.width()-1), p) - img2(q, p));
                }
              }
              d_cost(i, j, c) = cost;
          }
      }
  }

  //cout << "TRWS on MRF iterating..." << endl;
  while(iterations--){
      for(int i=1; i<img1.width()-1; i++){
          for(int j=1; j<img1.height()-1; j++){
              for(int c=0; c<classes;c++){
                double tmp_cost = 0; 
                double min_cost = 1e10;
                double min_cost1 = 1e10;
                for(int c1=0; c1<classes;c1++){
                    // vertical spanning tree
                    tmp_cost = d_cost(i, j, c1);
                    tmp_cost += alpha*(c1 != c ? 1 : 0);
                    for(int d=0; d<2; d++){
                        tmp_cost += msg_cost(i, j, c1, d, total_cost, true);
                    }
                    if(tmp_cost < min_cost){
                        min_cost = tmp_cost;
                    }
                    
                    // horizontal spannng tree
                    tmp_cost = d_cost(i, j, c1);
                    tmp_cost += alpha*(c1 != c ? 1 : 0);
                    for(int d=0; d<2; d++){
                        tmp_cost += msg_cost(i, j, c1, d, total_cost, false);
                    }
                    if(tmp_cost < min_cost1){
                        min_cost1 = tmp_cost;
                    }
                }
                total_cost(i, j, c) = 0.5 * (min_cost1 + min_cost);
              }
          }
      }
  }

  //cout << "TRWS on MRF predicting..."<< endl;
  for(int i=1; i<img1.width()-1; i++){
      for(int j=1; j<img1.height()-1; j++){
          double tmp_cost = 0;
          double min_cost = 1e10;
          int idx = -1;
          for(int c=0; c<classes;c++){
              tmp_cost = d_cost(i, j, c);
              tmp_cost += total_cost(i-1, j, c) + total_cost(i+1, j, c) + total_cost(i, j+1, c) + total_cost(i, j-1, c);
              if(tmp_cost < min_cost){
                  min_cost = tmp_cost;
                  idx = c;
              }
          }
          //result(i, j) = idx;
          result(j, i) = idx;
      }
  }
  return result;
}


int main(int argc, char *argv[])
{
  if(argc != 4 && argc != 3)
    {
      cerr << "usage: " << argv[0] << " image_file1 image_file2 [gt_file]" << endl;
      return 1;
    }
  string input_filename1 = argv[1], input_filename2 = argv[2];
  string gt_filename;
  if(argc == 4)
    gt_filename = argv[3];

  CImg<double> image1(input_filename1.c_str());
  CImg<double> image2(input_filename2.c_str());
  CImg<double> gt;

  if(gt_filename != "")
  {
    gt = CImg<double>(gt_filename.c_str());

    for(int i=0; i<gt.height(); i++)
      for(int j=0; j<gt.width(); j++)
        gt(j,i) = gt(j,i) / 3.0;
  }

  CImg<double> mrf_disp = mrf_stereo(image1, image2);
  mrf_disp.get_normalize(0,255).save((input_filename1 + "-disp_mrf.png").c_str());
  CImg<double> mrf_disp1 = TRWS_stereo(image1, image2);
  mrf_disp1.get_normalize(0,255).save((input_filename1 + "-disp_mrf1.png").c_str());

  if(gt_filename != "")
    {
      cout << "loopy BP, MRF stereo technique mean error = " << (mrf_disp-gt).sqr().sum()/gt.height()/gt.width() << endl;
      cout << "TRWS, MRF stereo technique mean error = " << (mrf_disp1-gt).sqr().sum()/gt.height()/gt.width() << endl;

    }

  return 0;
}


