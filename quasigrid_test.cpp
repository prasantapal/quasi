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
  auto n = 3;
  auto j = 2;
  QuasiGrid q;
  q.set_num_particles(n);
  q.set_num_intersections(j);
  auto system_length = 2.0*M_PI;
  q.set_system_len(system_length);
  q.set_arm_length();
  q.initialize_system();
  q.print_intersection_labels();
  q.print_junction_labels();
  q.set_junction_coordinates();
  auto inters = q.get_intersection(1);
  inters.print_intersection_coordinates();
  q.print_junction_coordinates();
  auto len = 0.1;
  q.set_len(len);
  ///all keyword means there is a uniformity
q.set_all_intersection_length();
q.print_all_intersection_length();


  std::cout << q.get_num_particles() << std::endl;
  q.set_particle_selection_distribution();
  auto  x=0.7854;
  auto val = q.does_belong_to_junction(x);
  if(val == nullptr){
    std::cout << x << " does not belongs to junction " << std::endl;
  }else {
    std::cout << x << " belongs to junction " << val->get_label() << " " << val->get_x() << std::endl;
  }
  return 0;
  std::cout << "critical density:" << q.calculate_critical_density_c() << std::endl;
  for(auto i=0;i<n;i++)
    std::cout << q.get_random_particle() << std::endl;
  std::cout << std::endl;
  int k =1;
  std::cout << n << " "<< " neighbor of  " << k << " " << q.get_neighbor_particle(k,QuasiGrid::Direction::Forward) << std::endl;
  q.move_particle(n-1,0.01);
  std::cout << "printing junction labels" << std::endl;
  return 0;
}
