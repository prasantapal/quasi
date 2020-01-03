#include "quasigrid.hpp"
int main(int argc, char** argv){
  /*
  Junction J;
  Particle p;
  Particle p2;
  std::cerr << "occupation status:"<< J.is_occupied() << std::endl;
  std::cerr << "occupying..." << std::endl;
  J.occupy(&p);
  std::cerr << "occupation status:"<< J.is_occupied() << std::endl;
  p2.print();
  J.occupy(&p2,false);
  std::cerr << "occupation status:"<< J.is_occupied() << std::endl;

  J.print_occupation();

  std::cerr << "unoccupying..." << std::endl;

  J.unoccupy();

  J.print_occupation();

  std::cerr << "occupation status:"<< J.is_occupied() << std::endl;

  J.unoccupy();

  std::cerr << "occupation status:"<< J.is_occupied() << std::endl;

  J.unoccupy();
  std::cerr << "occupation status:"<< J.is_occupied() << std::endl;
  */
  QuasiGrid q;
  // int n = 100;
  // for(auto i:std::vector<bool>(n,false)){
  //   std::cout << q.get_random() << std::endl;
  // }
  auto n = 3;
  auto j = 3;
  q.set_num_particles(n);
  q.set_num_intersections(j);
  q.initialize_system();

  q.print_intersection_labels();
  q.set_system_len(2.0*M_PI);
  auto len = 0.1;
  q.set_len(len);

return 0;
  std::cout << q.get_num_particles() << std::endl;
  q.set_particle_selection_distribution();
  std::cout << "critical density:" << q.calculate_critical_density_c() << std::endl;
  for(auto i=0;i<n;i++)
    std::cout << q.get_random_particle() << std::endl;

  q.initialize_system();

  std::cout << "system initialized:" << std::endl;

  int k =1;
  std::cout << n << " "<< " neighbor of  " << k << " " << q.get_neighbor_particle(k,QuasiGrid::Direction::Forward) << std::endl;
  q.move_particle(n-1,0.01);
  std::cout << "printing junction labels" << std::endl;

  return 0;
}
