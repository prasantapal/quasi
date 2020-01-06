#ifndef QUASI_GRID_HPP
#define QUASI_GRID_HPP
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
using  JunctionRef = std::reference_wrapper<Junction> ;
// Outlier(Outlier const& other) = default;
//    Outlier& operator=(Outlier const& other) = default;
//    Outlier(Outlier&& other) = default;
//    Outlier& operator=(Outlier&& other) = default;
//
///Junction is the unit of intersection
class QuasiGrid{
  public:
    enum Direction {Forward,Backward};
    QuasiGrid();
    QuasiGrid(QuasiGrid&) = delete;
    QuasiGrid& operator=(QuasiGrid const& ) = delete;
    QuasiGrid& operator=(QuasiGrid&& ) = delete;
    ~QuasiGrid();
    bool is_blocked(const int&) const;///Check whether or not a junction is blocked
    void set_arm_length();
    void set_all_intersection_length() const;
    void print_all_intersection_length();
    void print_intersection_labels();
    void print_junction_labels();
    void print_junction_coordinates();
    void set_len(const double& len);
    void set_junction_coordinates() const;
    void set_system_len(const double& len);
    int get_next_junction_index(const double& x);
    const Junction* does_belong_to_junction(const double&) const;
    inline double mod_len(const double x) {
      if(x>=0)
        return std::fmod(x,system_len_);
      return (system_len_ - std::fmod(std::fabs(x),system_len_));
    }
    void move_particle(const int& n,const double& dx);
    void distribute_junction_labels();
    void reassign_intersection_labels();
    int get_neighbor_particle(const int& n,const int direction);
    int mod_num_particles(const int& n);
    inline double get_random()   { return  normal_dist_(particle_motion_generator_); }
    inline double get_random_particle()   { return  particle_selection_dist_(particle_selection_generator_); }
    void set_particle_selection_distribution();
    double calculate_critical_density_c() const;
    void set_num_particles( const int n);
    int get_num_particles() const ;
    void set_num_intersections( const int n);
    int get_num_intersections() const ;
    void initialize_system();
    void populate_blocked_intersection(const int n);
    void depopulate_blocked_intersection(const int n);
    Intersection& get_intersection(const int index);
    const Junction& get_junction(const int index) const;
    static double get_arm_length();
    void populate_junction_conjugates();
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
    double critical_density_c_;
    double critical_density_g_;
    int num_intersections_;
    int num_particles_;
    double len_;
    double system_len_;;
    std::vector<Particle> particles_;
    std::vector<Intersection> intersections_;
    std::list<int> blocked_intersections_; ///This contains the labels for blocked intersections
    ///Make a vector tuple based on junction index with the contents being
    std::vector< Junction const*> junctions_;
    std::vector<JunctionRef> junction_refs_;
    std::multimap<int,int> junction_conjugates_;
    static double arm_length_;
    static int default_arm_offset_;
};
void QuasiGrid::populate_junction_conjugates(){
  for(const auto& intersection:intersections_){
    auto& intersection_ref = intersection.get_intersection_ref();
    auto& conjugate_map= intersection.get_junction_conjugate();
    for(const auto& pair:conjugate_map){
      std::cout << pair.first << " " << pair.second << std::endl;
      junction_conjugates_.emplace(pair);
    }
  }
}
double QuasiGrid::arm_length_ = {0.0};
int QuasiGrid::default_arm_offset_ = {2};
std::string QuasiGrid::class_name_ = {""};
QuasiGrid::QuasiGrid(){
  quasigrid_helper();
  std::cerr << class_name_ << " ctor" << std::endl;
}
QuasiGrid::~QuasiGrid(){
  std::cerr << class_name_<< " dtor" << std::endl;
}
void QuasiGrid::quasigrid_helper(){
  ///TEST
  //
  arm_length_ = 0.1;
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
double QuasiGrid::calculate_critical_density_c() const{
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
  for(auto& it:intersections_){
    it.allocate_labels(0,it.get_label());
    it.allocate_labels(1,2*intersections_.size() -1 -it.get_label());
  }
  print_intersection_labels();
}
void QuasiGrid::print_intersection_labels(){
  std::cerr << "printing intersection labels:" << std::endl;
  for(auto& it:intersections_){
    std::cout << it.get_label() << std::endl;
  }
  std::cout << std::endl;
}
void QuasiGrid::print_junction_labels(){
  std::cerr << "printing junction labels:" << std::endl;
  for(auto& it:junctions_){
    std::cout << it->get_label() << " ";
  }
  std::cout << std::endl;
}
void QuasiGrid::print_junction_coordinates(){
  std::cerr << "printing intersection coordinates:" << std::endl;
  for(auto& it:junctions_){
    std::cout << it->get_x() << " ";
  }
  std::cout << std::endl;
}
void QuasiGrid::initialize_system() {
  std::cerr << "initializing system " <<  std::endl;
  particles_ = std::vector<Particle>(num_particles_);
  intersections_ = std::move(std::vector<Intersection>(num_intersections_));
  distribute_junction_labels();
  for(auto& it:intersections_){
    auto& intersection_ref = it.get_intersection_ref();
    for(auto& its:intersection_ref) {
      junctions_.push_back(const_cast<Junction*>(&its));
    }
    it.set_junction_conjugate();
    it.print_conjugate();
  }
  populate_junction_conjugates();
  std::cout << std::endl;
  std::sort(junctions_.begin(),junctions_.end(),[](const auto& a,const auto&b){ return a->get_label()<b->get_label();});
  std::cerr << "printing junction labels" << std::endl;
  for(auto& jun:junctions_){
    std::cout << jun->get_label() << " ";
  }
  std::cout << std::endl;
  std::cout << "checking is blocked" << std::endl;
  for(auto& it:junctions_){
    this->is_blocked(it->get_label());
  }
  exit(0);
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
///Find the blocked junction closer than the next particle
int QuasiGrid::get_next_junction_index(const double& x){
  ///return the next junction with respect to a given point
  auto  last_junction_Index = 1;
}
///Check if there is a blocked junction before or after the next particle
void QuasiGrid::move_particle(const int& n,const double& dx){
  auto x = particles_[n].get_x();
  if(dx>0){
    auto neighbor_particle  = mod_num_particles(n+1);
    if(!particles_[neighbor_particle].does_belong_to(x + dx)){///Particle constraint is cleared now the junction constraint
      bool is_blocked_by_junction = {false};
      for(auto it:intersections_){
        auto& intersection_ref = it.get_intersection_ref();
        for(auto& its:intersection_ref){
          std::cout << "blocking status " << its.get_is_blocked() << std::endl;
        }
      }
    }else {///Well do nothing
    }
  }else {
    auto neighbor_particle  = mod_num_particles(n-1);
    if(!particles_[neighbor_particle].does_belong_to(x - dx)){///Particle constraint is cleared now the junction constraint
    }else {///Well do nothing
    }
  }
}
void QuasiGrid::set_len(const double& len){
  len_ = len;
}
void QuasiGrid::set_system_len(const double& len){
  system_len_ = len;
}
const Junction* QuasiGrid::does_belong_to_junction(const double& x) const {
  for(auto& it:junctions_){
    if(it->does_belong_to(x))
      return it;
  }
  return nullptr;
}
void QuasiGrid::print_all_intersection_length(){
  std::cerr << "printing all intersection length:" << std::endl;
  for(auto& it:intersections_){
    auto& ref = *it.get_intersection_ptr();
    std::cout << "intersection " << it.get_label() << ": ";
    for(auto& its:ref){
      std::cout << its.get_len() << " ";
    }
    std::cout << std::endl;
  }
}
void QuasiGrid::set_all_intersection_length()const {
  std::cerr << "setting all junction length" << std::endl;
  for(const auto& it:junctions_){
    (const_cast<Junction*>(it))->set_len(this->len_);
  }
}
void QuasiGrid::set_junction_coordinates() const{
  std::cerr << "setting junction coordinates:" << std::endl;
  std::vector<double> junction_coordinates(Intersection::get_junction_count());
  for(  auto const& it:junctions_){
    auto label = it->get_label();
    std::cout << "label:"<< label << " ";
    double x = {0.0};
    if(label < intersections_.size()) {
      x = (label + 1)*arm_length_;
    }else {
      x = (default_arm_offset_ + label )*arm_length_;
    }
    (const_cast<Junction*>(it))->set_x(x);
  }
  std::cout << std::endl;
  std::cout << std::endl;
}
const Junction& QuasiGrid::get_junction(const int n) const{
  return *junctions_.at(n-1);
}
Intersection& QuasiGrid::get_intersection(const int n){
  if(n<=intersections_.size() && n>=1)
    return intersections_.at(n-1);
  else{
    std::cerr << "index " << n << " out of range " << std::endl;
    throw "index out of range";
  }
}
void QuasiGrid::set_arm_length() {
  arm_length_ = system_len_ / (2.0*(num_intersections_ + 1));
  std::cerr << "arm_length:" << arm_length_ << std::endl;
}
double QuasiGrid::get_arm_length(){
  return arm_length_;
}
typedef std::multimap<int, int>::iterator MMAPIterator;
bool QuasiGrid::is_blocked(const int& index) const{
  std::cout << "index:" << index  << std::endl;
  auto matches  = junction_conjugates_.equal_range(index);
  // Iterate over the range
  for (auto it = matches.first; it != matches.second; it++){
    std::cout << it->second << std::endl;
    if((junctions_.at(it->second))->is_occupied()){
      return true;
    }
  }
  return false;
}
#endif
