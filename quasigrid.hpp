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
///Junction is the unit of intersection
class QuasiGrid{
  public:
    enum Direction {Forward,Backward};
    QuasiGrid();
    ~QuasiGrid();
    void distribute_junction_labels();
    void reassign_intersection_labels();
    int get_neighbor_particle(const int& n,const int direction);

int mod_num_particles(const int& n);

    inline double get_random()   { return  normal_dist_(particle_motion_generator_); }
    inline double get_random_particle()   { return  particle_selection_dist_(particle_selection_generator_); }
    void set_particle_selection_distribution();
    double calculate_critical_density() const;
    void set_num_particles( const int n);
    int get_num_particles() const ;
    void set_num_intersections( const int n);
    int get_num_intersections() const ;
    void initialize_system();
    void populate_blocked_intersection(const int n);
    void depopulate_blocked_intersection(const int n);
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
    std::vector<Particle> particles_;
    std::vector<Intersection> intersections_;
    std::list<int> blocked_intersections_; ///This contains the labels for blocked intersections
};
std::string QuasiGrid::class_name_ = {""};
QuasiGrid::QuasiGrid(){
  quasigrid_helper();
  std::cerr << class_name_ << " ctor" << std::endl;
  initialize_system();
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

/**
 * distribute junction labels to each intersections
 */
void QuasiGrid::distribute_junction_labels(){
  std::cout << "distributing labels:" << std::endl;
  for(auto it = intersections_.begin(); it != intersections_.end();++it) {
    std::cout << "counter:" << it->get_counter() << " label " << it->get_label() << " ";
    it->allocate_labels(0,it->get_label());
    it->allocate_labels(1,2*intersections_.size() -1 -it->get_label());

  }
  std::cout << std::endl;
  std::cout << "printing junction labels:" << std::endl;
  for(auto it:intersections_) {
    std::cout << "intersection: " << it.get_label() << std::endl;
    it.print_junction_labels();
  }

}
void QuasiGrid::initialize_system() {
  particles_ = std::vector<Particle>(num_particles_);
  intersections_ = std::vector<Intersection>(num_intersections_);
  distribute_junction_labels();
}
void QuasiGrid::populate_blocked_intersection(const int n){
  blocked_intersections_.push_back(n);
}
void QuasiGrid::depopulate_blocked_intersection(const int n){
  blocked_intersections_.remove(n);
}
void QuasiGrid::reassign_intersection_labels() {//Use it optionally if the intersection labels need reassignments

}

int QuasiGrid::mod_num_particles(const int& n){
if(n>=0){
  return n%num_particles_;
}
return num_particles_ - ((-n)%num_particles_);

}
int QuasiGrid::get_neighbor_particle(const int& n,const int direction){///zero is backward, 1 is forward
if(direction == Direction::Forward) {
  return mod_num_particles(n+1);
}

  return mod_num_particles(n-1);
}
