#include "quasigrid.hpp"
int main(int argc, char** argv){
  std::cerr << "testing quasi grid" << std::endl;
  QuasiGrid q;
 // int n = 100;
 // for(auto i:std::vector<bool>(n,false)){
 //   std::cout << q.get_random() << std::endl;
 // }
 auto n= 3;
 auto j=3;
 q.set_num_particles(n);
 q.set_num_intersections(j);
q.initialize_system();

 q.set_system_len(2.0*M_PI);
 auto len = 0.1;
 q.set_len(len);

 std::cout << q.get_num_particles() << std::endl;
q.set_particle_selection_distribution();
 std::cout << "critical density:" << q.calculate_critical_density_c() << std::endl;
for(auto i=0;i<n;i++)
 std::cout << q.get_random_particle() << std::endl;

q.initialize_system();
int k =1;
std::cout << n << " "<< " neighbor of  " << k << " " << q.get_neighbor_particle(k,QuasiGrid::Direction::Forward) << std::endl;
q.move_particle(n-1,0.01);
return 0;
}
