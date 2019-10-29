#include "particle.hpp"
#include "intersection.hpp"
#include <list>
#include <regex>
#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
class QuasiGrid{
  public:
    QuasiGrid();
    ~QuasiGrid();
    inline double get_random()   { return  normal_dist_(particle_motion_generator_); }
    inline double get_random_particle()   { return  particle_selection_dist_(particle_selection_generator_); }
    void set_particle_selection_distribution();
    double calculate_critical_density() const;
    void set_num_particles( const int n);
    int get_num_particles() const ;
    void set_num_intersections( const int n);
    int get_num_intersections() const ;
    void initialize_system();
  private:
    void quasigrid_helper();
    static std::string class_name_;
    std::random_device particle_motion_{};
    std::mt19937 particle_motion_generator_{particle_motion_()};
    std::normal_distribution<> normal_dist_{0,1};
    std::random_device particle_selection_{};
    std::mt19937 particle_selection_generator_{particle_selection_()};
    std::uniform_int_distribution<> particle_selection_dist_{1, 6};
    double density_;
    double critical_density_;
    int num_intersections_;
    int num_particles_;
    double particle_length_;
    std::list<Particle> particles_;
    std::list<Intersection> intersections_;
};
std::string QuasiGrid::class_name_ = {""};
QuasiGrid::QuasiGrid(){
  quasigrid_helper();
  std::cerr << class_name_ << " ctor" << std::endl;
}
QuasiGrid::~QuasiGrid(){
  std::cerr << class_name_<< " dtor" << std::endl;
}
void QuasiGrid::quasigrid_helper(){
  class_name_ = static_cast<std::string>(typeid(this).name());
  std::cerr << class_name_ << " helper" << std::endl;
  /*
  ////Fill this in to get the automatied class name
  std::regex re(".*[[:digit:]].*");
  std::smatch m;
  class_name_="hello 2 world 4";
  auto match_status=std::regex_match (class_name_,re);
  std::regex_search ( class_name_, m, re );
  for (unsigned i=0; i<m.size(); ++i) {
  std::cout << "match " << i << " (" << m[i] << ") ";
  std::cout << "at position " << m.position(i) << std::endl;
  }
  */
}
double QuasiGrid::calculate_critical_density() const{
  double density = {0};
  switch (num_intersections_){
    case 1:
      density = num_intersections_%2 == 0?static_cast<double>(num_particles_/(num_particles_ + 2.0)):static_cast<double>(num_particles_/(num_particles_ + 1.0));
      break;
    case 2:
      break;
    default:
      break;
  }
  return density;
}
void QuasiGrid::set_num_particles( const int n) {
  num_particles_ = n;
}
int QuasiGrid::get_num_particles() const {
  return num_particles_;
}
void QuasiGrid::set_num_intersections( const int n){
  num_intersections_ = n;
}
int QuasiGrid::get_num_intersections() const {
  return num_intersections_;
}
void QuasiGrid::set_particle_selection_distribution(){
  std::cout << "setting particle selection dist" << std::endl;
  particle_selection_dist_ = std::uniform_int_distribution<> (0, num_particles_);
}
void QuasiGrid::initialize_system() {
  particles_ = std::list<Particle>(num_particles_);
  intersections_ = std::list<Intersection>(num_intersections_);
}
