#include "quasigrid.hpp"
#define TEST_REFERENCE_WRAPPER
#undef TEST_REFERENCE_WRAPPER
int main(int argc, char** argv){

#ifdef TEST_REFERENCE_WRAPPER ///Tes section
//double f = std::numeric_limits<double>::quiet_NaN();
//f = NAN;
//if(std::isnan(f))
//std::cout << "f:" << f << std::endl;
//
//f = 100.0;
//
//std::cout << "f:" << f << std::endl;
//exit(0);
int s = std::nan(0);
std::cout << "s:" << s << std::endl;



  std::string hello = "Hello, ";
  std::string world = "everyone!";
  typedef std::vector<std::reference_wrapper<std::string>> vec_t;
  vec_t vec = {hello, world};
  std::string& hello_clone = hello;
  vec.push_back(hello_clone);
  vec.at(1).get() = "world!";
  std::cout << hello << world << std::endl;
  std::copy(vec.begin(),vec.end(),std::ostream_iterator<std::string>(std::cout," "));
#endif

  //     Junction J;
  //     Particle p;
  //     Particle p2;
  //     std::cerr << "occupation status:"<< J.is_occupied() << std::endl;
  //     std::cerr << "occupying..." << std::endl;
  //     J.occupy(&p);
  //     std::cerr << "occupation status:"<< J.is_occupied() << std::endl;
  //     p2.print();
  //     J.occupy(&p2,false);
  //     std::cerr << "occupation status:"<< J.is_occupied() << std::endl;
  //     J.print_occupation();
  //     std::cerr << "unoccupying..." << std::endl;
  //     J.unoccupy();
  //     J.print_occupation();
  //     std::cerr << "occupation status:"<< J.is_occupied() << std::endl;
  //     J.unoccupy();
  //     std::cerr << "occupation status:"<< J.is_occupied() << std::endl;
  //     J.unoccupy();
  //     std::cerr << "occupation status:"<< J.is_occupied() << std::endl;

  auto n = 3;
  auto j = 3;
  QuasiGrid q;
  q.set_num_particles(n);
  q.set_num_intersections(j);
  auto system_length = 2.0*M_PI;
  q.set_system_len(system_length);
  q.set_arm_length();
  int num_particles_per_middle_arm_at_max_packing = 1;
  q.set_num_particles_per_middle_arm_at_max_packing(num_particles_per_middle_arm_at_max_packing);
  std::cout << q.get_num_particles_per_middle_arm_at_max_packing() << std::endl;
  q.calculate_and_set_max_allowed_particles_in_end_lobes_at_max_packing();
  std::cout << q.get_max_allowed_particles_in_end_lobes_at_max_packing() << std::endl;
  q.calculate_and_set_min_num_particles_at_kinetic_arrest();
  std::cout << "get_min_no_of_particles_at_kinetic_arrest:" << q.get_min_num_particles_at_kinetic_arrest() << std::endl;
  q.calculate_and_set_num_particles_at_closed_packing();
  std::cout << q.get_num_particles_at_closed_packing()  << std::endl;
  unsigned int m = 2;
  q.set_num_particles_above_min_no_particles_at_kinetic_arrest(m);
  std::cout << "no_of_particles_above_min_no_particles_at_kinetic_arrest:" << q.get_num_particles_above_min_no_particles_at_kinetic_arrest() << std::endl;

  q.calculate_and_set_num_particles_from_filling_mode();
  q.calculate_and_set_num_particles_possible_in_system();
  std::cout << "get_num_particles_possible_in_system:" << q.get_num_particles_possible_in_system () << std::endl;

  std::cout << "get_num_particles_at_closed_packing:" << q.get_num_particles_at_closed_packing()  << std::endl;


  std::cout << "get_min_num_particles_at_kinetic_arrest:" << q.get_min_num_particles_at_kinetic_arrest() << std::endl;

  q.calculate_and_set_max_num_particles_at_kinetic_arrest();
  std::cout << "get_max_num_particles_at_kinetic_arrest:" << q.get_max_num_particles_at_kinetic_arrest() << std::endl;


  double density_delta_log_cale = 2.0;
  q.set_density_delta_log_scale(density_delta_log_cale);

  std::cout << "num particles " << q.get_num_particles() << std::endl;

  q.calculate_and_set_arm_void_length();

  std::cout << "get_arm_void_length:" << q.get_arm_void_length() << std::endl;
  std::cout << "get_arm_length:" << q.get_arm_length() << std::endl;
  std::cout << "get_particle_length:" << q.get_particle_len() << std::endl;
q.calculate_num_particles_possible_in_system();
  q.calculate_and_set_density_kinetic_arrest();

  q.calculate_and_set_density_from_density_delta_and_kinetic_arrest_density();
  std::cout << "density kinetic arrest:" << q.get_density_kinetic_arrest() << std::endl;
  std::cout << "density :" << q.get_density() << std::endl;
  q.calculate_and_set_particle_len();

  std::cout << "particle_len :" << q.get_particle_len() << std::endl;


q.print_system();


  exit(0);

  q.initialize_system();
  q.print_intersection_labels();
  q.print_junction_labels();
  q.set_junction_coordinates();
  auto& intersections= q.get_intersection(1);
  intersections.print_intersection_coordinates();
  q.print_junction_coordinates();
  auto len = 0.1;
  auto eps = -0.00000001;
  q.set_particle_len(len);
  ///all keyword means there is a uniformity
  q.set_all_intersection_length();
  q.print_all_intersection_length();
  auto intersection_ptr = intersections.get_intersection_ptr();
  auto& jn = q.get_junction(1);
  std::cout << "jn location:" << jn.get_x() << std::endl;

  std::cout << q.get_num_particles() << std::endl;
  q.set_particle_selection_distribution();
  auto  x = jn.get_x() - len - eps;
  auto val = q.does_belong_to_junction(x);
  if(val == nullptr){
    std::cout << x << " does not belongs to junction " << std::endl;
  }else {
    std::cout << x << " belongs to junction " << val->get_label() << " " << val->get_x() << std::endl;
  }
  return 0;
  std::cout << "critical density:" << q.calculate_density_close_packing() << std::endl;
  for(auto i=0;i<n;i++)
    std::cout << q.get_random_particle() << std::endl;
  std::cout << std::endl;
  int k =1;
  std::cout << n << " "<< " neighbor of  " << k << " " << q.get_neighbor_particle(k,QuasiGrid::Direction::Forward) << std::endl;
  q.move_particle(n-1,0.01);
  std::cout << "printing junction labels" << std::endl;
  return 0;
}
