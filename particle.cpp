/*
 *   CPP header template
 *   Developed by Dr. Prasanta Pal
 */

#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <thread>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include "particle.hpp"
int main(int argc,char** argv){
double val=M_PI;
double value =5*M_PI + 0.1;

std::cout << std::fmod(value , val) << std::endl;

  Particle p;
p.set_system_len(2.0*M_PI);
p.set_len(0.1);

  auto x = 6.0*M_PI + 0.28;
 x = -0.0500;

std::cout << p.does_belong_to(x) << std::endl;

return 0;
  std::cout << p.mod_len(x) << std::endl;


  std::string data_file = "random_walk.txt";
  p.set_filename(data_file);
  int n = 100;
  p.time_evolve(n);

  return 0;
}
