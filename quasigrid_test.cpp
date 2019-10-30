#include "quasigrid.hpp"
int main(int argc, char** argv){
  std::cerr << "testing quasi grid" << std::endl;
  QuasiGrid q;
 // int n = 100;
 // for(auto i:std::vector<bool>(n,false)){
 //   std::cout << q.get_random() << std::endl;
 // }
 auto n= 1;
 auto j=3;
 q.set_num_particles(n);
 q.set_num_intersections(j);

 std::cout << q.get_num_particles() << std::endl;
q.set_particle_selection_distribution();
 std::cout << "critical density:" << q.calculate_critical_density() << std::endl;
for(auto i=0;i<n;i++)
 std::cout << q.get_random_particle() << std::endl;

q.initialize_system();



return 0;
}
