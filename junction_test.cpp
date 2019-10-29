/*
 *   CPP header template
 *   Developed by Dr. Prasanta Pal
*/

#include <cstdio>
#include <iomanip>
#include <cstdlib>
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
#include "junction.hpp"
int main(int argc,char** argv){
 	std::cout << "hello world " << std::endl;
        Junction j[10];
        Junction j1;
        auto x = 1.0;
        auto len = 0.2;

        j1.set_x(x);
        j1.set_len(len);
        auto pos = x - 0.1000000000001;
        std::cout << j1.does_belong_to( pos) << std::endl;
        auto position = 1;
        j1.set_occupation(position);
        j1.set_occupation(position-1);
        std::cout << j1.is_occupied() << std::endl;
        j1.reset_occupation(position);
        std::cout << j1.is_occupied() << std::endl;

        j1.reset_occupation(position-1);
        std::cout << j1.is_occupied() << std::endl;



	return 0;
}
