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
#include <thread>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include "particle.hpp"
int main(int argc,char** argv){
  Particle p;
  std::string data_file = "random_walk.txt";
  p.set_filename(data_file);
  int n = 1000000;
p.time_evolve(n);

	return 0;
}
