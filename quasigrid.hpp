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
using USH = unsigned short;
using UIN = unsigned int;
using  JunctionRef = std::reference_wrapper<Junction> ;///This is a proposal to wrap a junction ref in a container
  template<class TupType, size_t... I>
void print_tuple(const TupType& _tup, std::index_sequence<I...>)
{
  std::cout << "(";
  (..., (std::cout << (I == 0? "" : ", ") << std::get<I>(_tup)));
  std::cout << ")\n";
}
  template<class... T>
void print_tuple (const std::tuple<T...>& _tup)
{
  print_tuple(_tup, std::make_index_sequence<sizeof...(T)>());
}
// Outlier(Outlier const& other) = default;
//    Outlier& operator=(Outlier const& other) = default;
//    Outlier(Outlier&& other) = default;
//    Outlier& operator=(Outlier&& other) = default;
//
///Junction is the unit of intersection
class QuasiGrid{
  public:
    enum Direction {Forward,Backward};
    enum StateLabels {EndLobe,JunctionLocation,MiddleLobe};
    QuasiGrid();
    QuasiGrid(QuasiGrid&) = delete;
    QuasiGrid& operator=(QuasiGrid const& ) = delete;
    QuasiGrid& operator=(QuasiGrid&& ) = delete;
    ~QuasiGrid();
    double calculate_mid_lobe_location(const int index) const;
    void set_system_length_scales();
    void set_system_void_spaces();
    std::map<std::string,double> get_system_length_scales() const;
    std::map<std::string,double> get_system_void_spaces() const;

    double get_mid_lobe_location(const int&);
    std::tuple<USH,USH> which_region_particle_belong_to(const Particle& p);
    int does_belong_to_end_lobe(const Particle& p);
    int does_belong_to_junction(const Particle& p);
    int does_belong_to_middle_lobe(const Particle& p) ;
    USH calculate_arm_index(const Particle& p);
    void populate_end_lobe_boundaries() ;
    void print_particles() const;
    void print_system() const;
    bool is_blocked(const Junction&) const;///Check whether or not a junction is blocked
    //bool is_blocked(const int&) const;///Check whether or not a junction is blocked
    void set_arm_length();
    double get_arm_length() const;
    void set_all_intersection_length() const;
    void print_all_intersection_length();
    void occupy_junction(const int&,Particle * const& particle,const bool&  is_forward_direction=true) const;
    void print_intersection_labels();
    void print_junction_labels();
    void print_junction_coordinates();
    void set_particle_len(const double& len);
    double get_particle_len() const;
    double calculate_particle_len()const;
    void calculate_and_set_particle_len();
    void set_junction_coordinates() const;
    void set_system_len(const double& len);
    double get_system_len() const;
    int get_next_junction_index(const double& x);
    const Junction* does_belong_to_junction(const double&) const;
    inline double mod_len(const double x) { if(x>=0) return std::fmod(x,system_len_); return (system_len_ - std::fmod(std::fabs(x),system_len_)); }
    void move_particle(const int& n,const double& dx);
    void distribute_junction_labels();
    void reassign_intersection_labels();///\brief assign the intersection labels starting from zero go max intersections
    int get_neighbor_particle(const int& n,const int direction);
    int mod_num_particles(const int& n);///\brief get particle as a module of total number of particles
    inline double get_random()   { return  normal_dist_(particle_motion_generator_); }
    inline double get_random_particle()   { return  particle_selection_dist_(particle_selection_generator_); }
    void set_particle_selection_distribution();
    double calculate_density_close_packing() const;
    double calculate_density_kinetic_arrest();
    void calculate_and_set_density_kinetic_arrest();
    double calculate_density_from_density_delta_and_kinetic_arrest_density() const ;
    void calculate_and_set_density_from_density_delta_and_kinetic_arrest_density();
    /************* SYSTEM DEFINITIONS  ***************/
    void set_num_particles( const int n);
    int get_num_particles() const ;
    UIN calculate_num_particles_from_filling_mode();
    UIN calculate_and_set_num_particles_from_filling_mode();
    void set_num_intersections( const int n);
    UIN get_num_intersections() const ;
    void set_num_junctions( const int n);
    UIN get_num_junctions()const;
    /************* INITIALIZATION *******************/
    void initialize_system();
    /************** DYNAMICS ***********************/
    void populate_blocked_intersection(const int n);
    void print_end_lobe_boundaries() const;
    void depopulate_blocked_intersection(const int n);
    Intersection& get_intersection(const int index);
    Junction& get_junction(const int index) ;
    void set_num_arms(const UIN&);
    UIN calculate_num_arms() const;
    void calculate_and_set_num_arms();
    UIN get_num_arms() const;
    UIN get_num_arms_half() const;
    void populate_junction_conjugates();
    void set_min_num_particles_at_kinetic_arrest(const int&);
    int calculate_min_num_particles_at_kinetic_arrest() const ;
    void calculate_and_set_min_num_particles_at_kinetic_arrest()  ;
    int get_min_num_particles_at_kinetic_arrest() const ;
    void set_num_particles_per_middle_arm_at_max_packing(const int& k); ///This is a design choice related to the topology of a particular instance of quasi
    int get_num_particles_per_middle_arm_at_max_packing() const; ///This is a design choice related to the topology of a particular instance of quasi
    void set_max_allowed_particles_in_end_lobes_at_max_packing( const int& val);
    int get_max_allowed_particles_in_end_lobes_at_max_packing() const;
    int calculate_max_allowed_particles_in_end_lobes_at_max_packing() const;
    void calculate_and_set_max_allowed_particles_in_end_lobes_at_max_packing();
    void set_num_particles_at_closed_packing(const int& n);
    unsigned int get_num_particles_at_closed_packing() const ;
    UIN calculate_num_particles_at_closed_packing() const;
    void calculate_and_set_num_particles_at_closed_packing();
    void set_max_num_particles_at_kinetic_arrest( const UIN);
    UIN get_max_num_particles_at_kinetic_arrest()const;
    UIN calculate_max_num_particles_at_kinetic_arrest() const;
    void calculate_and_set_max_num_particles_at_kinetic_arrest();
    void set_num_particles_above_min_no_particles_at_kinetic_arrest(const int& n);
    UIN get_num_particles_above_min_no_particles_at_kinetic_arrest() const;
    void set_num_particles_possible_in_system( const int&);
    UIN get_num_particles_possible_in_system() const;
    UIN calculate_num_particles_possible_in_system() const;
    void calculate_and_set_num_particles_possible_in_system();
    /************* STATE CALCULATIONS ******************/
    void AssignStateLablesToParticles() ;
    void set_density(const int& phi);
    double get_density() const;
    double calculate_density() const;
    void calculate_and_set_density();
    void set_arm_void_length(const double&);
    double get_arm_void_length() const;
    double calculate_arm_void_length() const;
    void calculate_and_set_arm_void_length();
    void set_density_delta_log_scale(const double& density);
    double get_density_delta_log_scale()const;
    void set_density_close_packing(const double&);
    double get_density_close_packing() const;
    void set_density_kinetic_arrest(const double& density);
    double get_density_kinetic_arrest() const;
    double calculate_density_kinetic_arrest()const;
    void set_trajectory_output_filename(const std::string& filename);
    std::string get_trajectory_output_filename() const;
    void set_trajectory_input_filename(const std::string& filename);
    std::string get_trajectory_input_filename() const;
    UIN calculate_system_state_size() const;
    void calculate_and_set_system_state_size();
    void set_system_state_size(const UIN&);
    UIN get_system_state_size() const;
    std::vector<UIN> calculate_system_state();
    UIN calculate_lower_end_lobe_occupation() const;
    UIN calculate_upper_end_lobe_occupation() const;
    std::vector<Particle>&  get_particles();
    Particle&  get_particle(const USH& index);
  private:
    void quasigrid_helper();
    static std::string class_name_;
    /*************** Dynamics and Distribution...later to be moved to a class ********/
    std::random_device particle_motion_{};
    std::mt19937 particle_motion_generator_{particle_motion_()};
    std::normal_distribution<> normal_dist_{0,1};
    std::random_device particle_selection_{};
    std::mt19937 particle_selection_generator_{particle_selection_()};
    std::uniform_int_distribution<> particle_selection_dist_{1, 6};
    double density_;
    double density_close_packing_;
    double density_kinetic_arrest_;
    double density_delta_log_scale_;
    double density_delta_;
    int num_intersections_;
    int num_intersections_half_;
    int num_junctions_;
    int num_junctions_half_;
    int num_particles_;
    int num_particles_possible_in_system_at_closed_packing_;
    int num_particle_size_voids_at_kinetic_arrest_;
    double particle_len_;
    double particle_len_half_;
    double system_len_;
    double system_len_half_;
    double len_at_kinetic_arrest_;
    std::vector<Particle> particles_;
    std::vector<Intersection> intersections_;
    std::list<int> blocked_intersections_; ///This contains the labels for blocked intersections
    ///Make a vector tuple based on junction index with the contents being
    std::vector<Junction *> junctions_;
    //    std::vector<JunctionRef> junction_refs_;
    std::multimap<int,int> junction_conjugates_;
    static double arm_length_;
    static double arm_void_length_;
    static int default_arm_offset_;
    static USH num_particles_per_middle_arm_at_max_packing_;///\brief K parameter in the paper
    static UIN num_particles_at_closed_packing_;///\brief $N^{cp}$ parameter in the paper
    static UIN num_particles_possible_in_system_;///\brief $N^{cp}$ parameter in the paper
    static UIN max_no_of_particles_at_kinetic_arrest_;///\brief this is the max number of particles when the system can still move
    static UIN min_num_particles_at_kinetic_arrest_;///\brief min number of particles when the sytem can move and also KA happens
    static UIN max_num_particles_at_kinetic_arrest_;///\brief min number of particles when the sytem can move and also KA happens
    static UIN no_of_particles_above_min_no_particles_at_kinetic_arrest_;
    static UIN range_of_particles_between_min_to_max_at_kinetic_arrest_;
    static USH max_allowed_particles_in_end_lobes_at_max_packing_;
    static UIN min_particle_size_voids_needed_for_kinetic_arrest_;
    static double kinetic_arrest_length_scale_;
    std::string trajectory_output_filename_ = {""};
    std::string trajectory_input_filename_ = {""};
    std::string input_file_path_ = {""};///\brief input folder path
    std::string trajectory_input_file_complete_path_ = {""};
    std::string output_file_path_ = {""};///\brief output folder path
    std::string trajectory_output_file_complete_path_ = {""};
    std::vector<UIN> system_state_;
    std::string system_state_string_;
    UIN system_state_size_;
    static std::string system_state_delimiter_;
    static std::vector<std::tuple<double,double>> lower_end_lobe_boundaries_;
    static std::vector<std::tuple<double,double>> upper_end_lobe_boundaries_;
    static UIN  num_arms_;
    static UIN  num_arms_half_;
    static std::map<std::string,double> system_length_scales_; ///Different length scales
    static std::map<std::string,double> system_void_spaces_; ///Different void spaces
};

std::map<std::string,double> QuasiGrid::system_void_spaces_; ///Different length scales
std::map<std::string,double> QuasiGrid::system_length_scales_;
UIN  QuasiGrid::num_arms_ = {0};
UIN  QuasiGrid::num_arms_half_ = {0};
std::vector<std::tuple<double,double>> QuasiGrid::lower_end_lobe_boundaries_;
std::vector<std::tuple<double,double>> QuasiGrid::upper_end_lobe_boundaries_;
std::string QuasiGrid::system_state_delimiter_ = {"-"};
double QuasiGrid::kinetic_arrest_length_scale_ = {NAN};
UIN QuasiGrid::min_particle_size_voids_needed_for_kinetic_arrest_ = {2};
UIN QuasiGrid::max_num_particles_at_kinetic_arrest_ = {0};
UIN QuasiGrid::num_particles_possible_in_system_ = {0};///\brief this is the max number of particles when the system can still move
UIN QuasiGrid::min_num_particles_at_kinetic_arrest_ = {0};
UIN QuasiGrid::max_no_of_particles_at_kinetic_arrest_ = {0};
UIN QuasiGrid::no_of_particles_above_min_no_particles_at_kinetic_arrest_ = {0};
UIN QuasiGrid::range_of_particles_between_min_to_max_at_kinetic_arrest_ = {0};
UIN QuasiGrid::num_particles_at_closed_packing_ = {0};///\brief $N^{cp}$ parameter in the paper
USH QuasiGrid::num_particles_per_middle_arm_at_max_packing_ = {0};
USH QuasiGrid::max_allowed_particles_in_end_lobes_at_max_packing_ ={0};
double QuasiGrid::arm_void_length_ = NAN;
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
void QuasiGrid::set_num_particles_at_closed_packing(const int& n){
  num_particles_at_closed_packing_ = n;
}
unsigned int QuasiGrid::get_num_particles_at_closed_packing() const {
  return num_particles_at_closed_packing_;
}
UIN QuasiGrid::calculate_num_particles_at_closed_packing() const{
  //  return  2*(num_intersections_ + 1)*(num_particles_per_middle_arm_at_max_packing_ + 1) - num_intersections_;
  return  2*(num_particles_per_middle_arm_at_max_packing_ + 1) + num_intersections_*(2*num_particles_per_middle_arm_at_max_packing_ + 1) ;
}
void QuasiGrid::calculate_and_set_num_particles_at_closed_packing() {
  this->set_num_particles_at_closed_packing(std::move(this->calculate_num_particles_at_closed_packing()));
}
double QuasiGrid::calculate_density_close_packing() const{
  double density = {0};
  switch (num_intersections_){
    case 1:
      density = num_intersections_%2 == 0?static_cast<double>(num_particles_)/(num_particles_ + 2.0):static_cast<double>(num_particles_/(num_particles_ + 1.0));
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
  num_intersections_half_ = num_intersections_/2;
  num_junctions_half_  = num_intersections_;

}
UIN QuasiGrid::get_num_intersections() const {
  return num_intersections_;
}
void QuasiGrid::set_num_junctions( const int n){
  num_junctions_ = n;
  num_junctions_half_ = num_junctions_/2;
}
UIN QuasiGrid::get_num_junctions() const {
  return num_junctions_;
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
    this->is_blocked(*junctions_.at(it->get_label()));
  }
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
  short int sign_of_dx = dx>0?1:-1;
bool is_forward_move = sign_of_dx ==1?true:false;
  auto x = particles_.at(n).get_x();
  auto x_next = x + sign_of_dx*dx;

    auto neighbor_particle  = mod_num_particles(n + sign_of_dx);
    if(!particles_.at(neighbor_particle).does_belong_to(x_next)){///Particle constraint is cleared now the junction constraint
      bool is_blocked_by_junction = {false};
      bool should_break_intersection_loop = {false};
      for(auto it:intersections_){
        auto& intersection_ref = it.get_intersection_ref();

        for(const Junction& its:intersection_ref){///its is a junction
          bool did_belong_to_junction = its.does_belong_to(x);
          bool does_belong_to_junction = its.does_belong_to(x_next);
          if(did_belong_to_junction) {
            if(does_belong_to_junction) {//regular update

            } else {//exit from junction
              //its.unoccupy(true);
#warning ENSURE NEW OBJECTS ARE NOT BEING CREATED HERE
              std::remove_const<Junction>::type its_non_const(its);
              its_non_const.unoccupy(is_forward_move);
            }
          }else {
            if(does_belong_to_junction){ //new belongning to junction
              std::remove_const<Junction>::type its_non_const(its);
              its_non_const.occupy(&particles_.at(n),is_forward_move);

            } else {///neighter belonged to junction nor it does now!
              ///regular jump
            }
          }
        }
        if(should_break_intersection_loop)
          break;

      }
    }else {///Well do nothing
    }

}

void QuasiGrid::set_particle_len(const double& len){
  particle_len_ = len;
}
double QuasiGrid::get_particle_len()const{
  return  particle_len_;
}
void QuasiGrid::set_system_len(const double& len){
  system_len_ = len;
  std::cerr << "setting system length " << system_len_ << std::endl;
  system_len_half_ = system_len_/2.0;
  std::cerr << "setting system_len_half_ " << system_len_half_ << std::endl;
}
double QuasiGrid::get_system_len() const {
  return  system_len_;
}
/**
 *@returns ptr to junction else nullptr
 */
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
    (const_cast<Junction*>(it))->set_len(this->particle_len_);
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
/**
 *
 */
Junction& QuasiGrid::get_junction(const int n) {
  return *junctions_.at(n-1);
}
Intersection& QuasiGrid::get_intersection(const int n){
  if(n<=intersections_.size() && n>=1)
    return intersections_.at(n-1);
  else{
    std::cerr << "intersection index " << n << " out of range " << std::endl;
    throw "index out of range";
  }
}
void QuasiGrid::set_arm_length() {
  arm_length_ = system_len_ / (2.0*(num_intersections_ + 1));
  std::cerr << "arm_length:" << arm_length_ << std::endl;
}
double QuasiGrid::get_arm_length() const{
  return arm_length_;
}
bool QuasiGrid::is_blocked(const Junction& j) const{
  std::cout << "index:" << index  << std::endl;
  auto matches  = junction_conjugates_.equal_range(j.get_label());
  // Iterate over the range
  for (auto it = matches.first; it != matches.second; it++){
    //std::cout << it->second << std::endl;
    if((junctions_.at(it->second))->is_occupied()){
      return true;
    }
  }
  return false;
}
void QuasiGrid::occupy_junction(const int& junction_index, Particle * const& particle,const bool& is_forward_direction)const{
  const_cast<Junction *>(this->junctions_.at(junction_index))->occupy(particle,is_forward_direction);
}
void QuasiGrid::set_num_particles_per_middle_arm_at_max_packing(const int& k){
  num_particles_per_middle_arm_at_max_packing_ = k;
}
int QuasiGrid::get_min_num_particles_at_kinetic_arrest() const {
  return this->min_num_particles_at_kinetic_arrest_;
}
void QuasiGrid::set_min_num_particles_at_kinetic_arrest(const int& n)  {
  this->min_num_particles_at_kinetic_arrest_ = n;
}
int QuasiGrid::calculate_min_num_particles_at_kinetic_arrest() const{ ///This is the minimum number of particles where kinetic arrest will happen
  if(num_intersections_ > 1 ){
    auto val =  2*(max_allowed_particles_in_end_lobes_at_max_packing_ -1 -1) + num_intersections_ * num_particles_per_middle_arm_at_max_packing_;
    return val;
  }else {
    std::cerr << "sorry this can be calculated only for num_intersections > 1" << std::endl;
    throw "num intersections > 1 ";
  }
  return 0;
}
void QuasiGrid::calculate_and_set_min_num_particles_at_kinetic_arrest()  {
  this->set_min_num_particles_at_kinetic_arrest(std::move(this->calculate_min_num_particles_at_kinetic_arrest()));
}
void QuasiGrid::set_max_allowed_particles_in_end_lobes_at_max_packing( const int& val) {
  this->max_allowed_particles_in_end_lobes_at_max_packing_ = val;
}
int QuasiGrid::calculate_max_allowed_particles_in_end_lobes_at_max_packing() const{
  return 2*num_particles_per_middle_arm_at_max_packing_ + 1;
}
int QuasiGrid::get_max_allowed_particles_in_end_lobes_at_max_packing() const{
  return this->max_allowed_particles_in_end_lobes_at_max_packing_;
}
void QuasiGrid::calculate_and_set_max_allowed_particles_in_end_lobes_at_max_packing(){
  auto n(std::move(calculate_max_allowed_particles_in_end_lobes_at_max_packing()));
  set_max_allowed_particles_in_end_lobes_at_max_packing(n);
}
int QuasiGrid::get_num_particles_per_middle_arm_at_max_packing() const{
  return this->num_particles_per_middle_arm_at_max_packing_;
}
void QuasiGrid::set_num_particles_above_min_no_particles_at_kinetic_arrest(const int& n) {
  this->no_of_particles_above_min_no_particles_at_kinetic_arrest_ = n;
}
UIN QuasiGrid::get_num_particles_above_min_no_particles_at_kinetic_arrest() const{
  return no_of_particles_above_min_no_particles_at_kinetic_arrest_;
}
UIN QuasiGrid::calculate_num_particles_from_filling_mode() {
  return      (this->min_num_particles_at_kinetic_arrest_ + this->no_of_particles_above_min_no_particles_at_kinetic_arrest_);
}
///K param defines a filling mode
UIN QuasiGrid::calculate_and_set_num_particles_from_filling_mode(){
  this->set_num_particles(std::move(this->calculate_num_particles_from_filling_mode()));
}
void QuasiGrid::set_density(const int& phi){
  this->density_ = phi;
}
double QuasiGrid::get_density() const{
  return this->density_;
}
double QuasiGrid::calculate_density_from_density_delta_and_kinetic_arrest_density() const {
  return density_kinetic_arrest_ - std::pow(10.0,-density_delta_log_scale_);
}
void QuasiGrid::calculate_and_set_density_from_density_delta_and_kinetic_arrest_density()  {
  this->density_ = std::move(this->calculate_density_from_density_delta_and_kinetic_arrest_density());
}
double QuasiGrid::calculate_density() const {
  return static_cast<double>(num_particles_*particle_len_)/system_len_;
}
void QuasiGrid::calculate_and_set_density() {
  this->density_ = std::move(this->calculate_density());
}
UIN QuasiGrid::get_num_particles_possible_in_system() const{
  return this->num_particles_possible_in_system_;
}
//void QuasiGrid::set_num_particles_at_closed_packing(const int& n){
//  this->num_particles_possible_in_system_at_closed_packing_ = n;
//}
void QuasiGrid::calculate_and_set_num_particles_possible_in_system(){
  this->num_particles_possible_in_system_  = std::move(this->calculate_num_particles_possible_in_system() );
  //this->max_no_of_particles_at_kinetic_arrest_ = std::move(this->calculate_num_particles_possible_in_system_at_closed_packing() );
}
UIN QuasiGrid::calculate_num_particles_possible_in_system() const{
  return 2*(2*num_particles_per_middle_arm_at_max_packing_ + 1 ) + 2*num_intersections_ + 2*(num_intersections_-1)*num_particles_per_middle_arm_at_max_packing_;
}
void QuasiGrid::print_system() const{
  std::cout << "J:" << num_intersections_ << " K:" << num_particles_per_middle_arm_at_max_packing_  << " N:" << num_particles_ << std::endl;
  std::cout << "N_min:" << min_num_particles_at_kinetic_arrest_ << " N_max:" << max_no_of_particles_at_kinetic_arrest_  << std::endl;
  std::cout << "N_above_min:" <<  no_of_particles_above_min_no_particles_at_kinetic_arrest_ << " N_voids:" << num_particle_size_voids_at_kinetic_arrest_  << std::endl;
}
void QuasiGrid::set_max_num_particles_at_kinetic_arrest( const UIN n){
  max_num_particles_at_kinetic_arrest_ = n;
}
UIN QuasiGrid::get_max_num_particles_at_kinetic_arrest()const{
  return max_num_particles_at_kinetic_arrest_;
}
UIN QuasiGrid::calculate_max_num_particles_at_kinetic_arrest() const{
  return num_particles_at_closed_packing_ - this->min_particle_size_voids_needed_for_kinetic_arrest_;
}
void QuasiGrid::calculate_and_set_max_num_particles_at_kinetic_arrest(){
  this->max_num_particles_at_kinetic_arrest_ = std::move(this->calculate_max_num_particles_at_kinetic_arrest());
}
void QuasiGrid::set_arm_void_length(const double& len){
  arm_void_length_ = len;
}
double QuasiGrid::get_arm_void_length() const{
  return arm_void_length_ ;
}
double QuasiGrid::calculate_arm_void_length() const{
  return arm_length_ - 2.0*particle_len_;
}
void QuasiGrid::calculate_and_set_arm_void_length(){
  this->arm_void_length_ = std::move(this->calculate_arm_void_length());
}
void QuasiGrid::set_density_close_packing(const double& density){
  density_close_packing_ = density;
}
double QuasiGrid::get_density_close_packing() const{
  return this->density_close_packing_;
}
void QuasiGrid::set_density_kinetic_arrest(const double& density){
  density_kinetic_arrest_ = density;
}
double QuasiGrid::get_density_kinetic_arrest() const{
  return this->density_kinetic_arrest_;
}
double QuasiGrid::calculate_density_kinetic_arrest(){//phi_g
  return static_cast<double>(num_particles_)/num_particles_possible_in_system_;
}
void QuasiGrid::calculate_and_set_density_kinetic_arrest(){//phi_g
  this->set_density_kinetic_arrest(std::move(this->calculate_density_kinetic_arrest()));
}
void QuasiGrid::set_density_delta_log_scale(const double& density){
  this->density_delta_log_scale_ = density;
}
double QuasiGrid::get_density_delta_log_scale()const{
  return  this->density_delta_log_scale_;
}
double QuasiGrid::calculate_particle_len()const{
  return  (density_*system_len_)/num_particles_;
}
void QuasiGrid::calculate_and_set_particle_len(){
  this->particle_len_ = std::move(this->calculate_particle_len());
}
void QuasiGrid::set_trajectory_output_filename(const std::string& filename){
  this->trajectory_output_filename_ = filename;
}
std::string QuasiGrid::get_trajectory_output_filename() const{
  return this->trajectory_output_filename_;
}
void QuasiGrid::set_trajectory_input_filename(const std::string& filename){
  this->trajectory_input_filename_ = filename;
}
std::string QuasiGrid::get_trajectory_input_filename() const{
  return this->trajectory_input_filename_;
}
UIN QuasiGrid::calculate_system_state_size() const{
  std::cout << "num_junctions_:" << num_intersections_ << " " << 2 + num_intersections_ + 2*num_intersections_ << std::endl;
  return 2 + num_intersections_ + 2*(num_intersections_ -1);
}
void QuasiGrid::set_system_state_size( const UIN& size){
  this->system_state_size_ = size;
}
UIN QuasiGrid::get_system_state_size() const{
  return this->system_state_size_;
}
void QuasiGrid::calculate_and_set_system_state_size(){
  this->system_state_size_ = std::move(this->calculate_system_state_size());
}
UIN QuasiGrid::calculate_lower_end_lobe_occupation() const{
}
UIN QuasiGrid::calculate_upper_end_lobe_occupation() const{
}

void QuasiGrid::AssignStateLablesToParticles() {
#ifdef  VERBOSE
  std::cerr << "in function " << __func__ << std::endl;
#endif
  for(auto& it:particles_){
    std::cerr << "particle label:"<< it.get_label() << std::endl;
  }
}

std::vector<UIN> QuasiGrid::calculate_system_state() {
  std::cerr << "calculating system states for " << num_particles_ << std::endl;
  decltype(system_state_) system_state(this->system_state_size_,0);


  for( auto& particle:particles_){
    auto region = which_region_particle_belong_to(particle);
    auto region_label = std::get<0>(region);
    auto region_index = std::get<1>(region);
    short vector_location = {0};

    switch(region_label){
      case StateLabels::EndLobe:
        std::cout << "end lobe:" << region_index << " " <<  system_state.at(vector_location)<< std::endl;
        vector_location = (region_index == 0)?region_index:(system_state_size_-1);
        system_state.at(vector_location)++;

        break;
      case StateLabels::JunctionLocation:
        std::cout << "junction location:" << region_index << std::endl;
        if(region_index < num_junctions_half_){
          vector_location = 1 + 3*region_index;
          system_state.at(vector_location)++;

        }else {
          std::cout << "else condition " << region_index << std::endl;
          //    vector_location =
          region_index -= num_junctions_half_;
          region_index = num_junctions_half_ - region_index -1;
          std::cout << "else condition " << region_index << std::endl;
          vector_location = 1 + 3*region_index;
          if(system_state.at(vector_location) == 0)
            system_state.at(vector_location) = 2;
          system_state.at(vector_location)++;

        }
        break;
      case StateLabels::MiddleLobe:
        std::cout << "mid lobe location:" << region_index << " num arms " << num_arms_ << " " << num_arms_half_  << std::endl;
        if(region_index < num_arms_half_ ){
          vector_location = 2 + 3*(region_index-1);
          system_state.at(vector_location)++;

        }else {
          std::cout << "else condition " << region_index << std::endl;
          //    vector_location =
          region_index -= num_arms_half_;
          region_index = num_arms_half_ - region_index -1;
          std::cout << "mid lobe else condition " << region_index << std::endl;
          vector_location = 3 + 3*(region_index-1);
          system_state.at(vector_location)++;

        }

        break;
      default:
        std::cout << "unknown region label" << std::endl;
        throw "unknown region label";
        break;
    }
  }
  std::cout << "system state " << std::endl;
  for(const auto& it:system_state) {
    std::cout << it;
  }

  std::cout << std::endl;
  ///First lower end lobe this.StateLabels.EndLobe with value 0
  //junctions
  //middle lobes
  ///last upper end lobe this.StateLabels.EndLobe with value 1
  return system_state;
}


void QuasiGrid::print_particles() const{
  std::cout << "label" << " " << " x " << " uncorrected x" << " state label " << " state label id" << std::endl;
  for(const auto& it:particles_) {
    std::cout << it.get_label() << " " << it.get_x() << " " << it.get_X() << " " << it.get_state_label() << " " << it.get_state_label_id() << std::endl;
  }
}

void QuasiGrid::print_end_lobe_boundaries() const{
  std::cout << "lower end lobe boundaries" << std::endl;
  for(const auto& it:lower_end_lobe_boundaries_){
    print_tuple(it);
  }
  std::cout << "upper end lobe boundaries" << std::endl;
  for(const auto& it:upper_end_lobe_boundaries_){
    print_tuple(it);
  }
}
/**
 *
 */
std::tuple<USH,USH> QuasiGrid::which_region_particle_belong_to(const Particle& p){
  ///END LOBE CALCULATION
  int location;
  location = does_belong_to_end_lobe(p);
  std::tuple<USH,USH> belongs_to;
  if( location >= 0){
    belongs_to = decltype(belongs_to)(StateLabels::EndLobe,location);
  }else {

    location = does_belong_to_junction(p);
    if(location >= 0){

      belongs_to = decltype(belongs_to)(StateLabels::JunctionLocation,location);
    }else {
      location = calculate_arm_index(p);
      belongs_to = decltype(belongs_to)(StateLabels::MiddleLobe,location);
    }
  }

  return belongs_to;
}
int QuasiGrid::does_belong_to_junction(const Particle& p) {
  auto junction_ref = does_belong_to_junction(p.get_x());
  if(junction_ref != nullptr){
    return junction_ref->get_label();
  }
  return -1;
}
int QuasiGrid::does_belong_to_middle_lobe(const Particle& p) {
}
int QuasiGrid::does_belong_to_end_lobe(const Particle& p) {
  auto x = p.get_x();
  for(const auto& it:lower_end_lobe_boundaries_){
    if(it == *lower_end_lobe_boundaries_.begin()){
      if(x>= std::get<0>(it) && x<std::get<1>(it) ){
        return 0;
      }
    }else {
      if(x> std::get<0>(it) && x<=std::get<1>(it) ){
        return 0;
      }
    }
  }
  std::cout << "upper end lobe boundaries" << std::endl;
  for(const auto& it:upper_end_lobe_boundaries_){
    if(x> std::get<0>(it) && x<std::get<1>(it) ){
      return 1;
    }
  }
  return -1;
}
void QuasiGrid::populate_end_lobe_boundaries() {
  lower_end_lobe_boundaries_.push_back(std::tuple<double,double>(0,arm_length_ - particle_len_half_));
  lower_end_lobe_boundaries_.push_back(std::tuple<double,double>(system_len_-(arm_length_ - particle_len_half_),system_len_));
  auto& last_intersection = intersections_.at(num_intersections_-1).get_intersection_ref();
  std::cout << last_intersection.at(0).get_x() << std::endl;
  std::cout << last_intersection.at(1).get_x() << std::endl;
  upper_end_lobe_boundaries_.push_back(std::tuple<double,double>(last_intersection.at(0).get_x() + particle_len_half_,last_intersection.at(1).get_x() - particle_len_half_));
  //lower_end_lobe_boundaries_
  //
  //upper_end_lobe_boundaries_;
  //
}
/**
 *
 *dependencie s
 */
USH QuasiGrid::calculate_arm_index(const Particle& p){
  return static_cast<USH>(p.get_x()/arm_length_);
}
Particle&  QuasiGrid::get_particle(const USH& index) {
  if(index>0 && index <= particles_.size())
    return particles_.at(index-1);
  else {
    std::cerr << "sorry index is out of bound" << std::endl;
    throw "index out of bound";
  }

}
std::vector<Particle>&  QuasiGrid::get_particles() {
  return particles_;
}
void QuasiGrid::set_num_arms(const UIN& num_arms){
  num_arms_ = num_arms;
  num_arms_half_ = num_arms/2;
}
UIN QuasiGrid::get_num_arms() const {
  return this->num_arms_;
}
UIN QuasiGrid::get_num_arms_half() const{
  return this->num_arms_half_;
}
UIN QuasiGrid::calculate_num_arms() const{
  return 2*(num_intersections_ + 1);
}
void QuasiGrid::calculate_and_set_num_arms(){
  this->set_num_arms(std::move(this->calculate_num_arms()));
}


double QuasiGrid::calculate_mid_lobe_location(const int index) const{
  std::cout << "arm_length:" << arm_length_ << std::endl;
  return arm_length_*(static_cast<double>(index) - 0.5);

}
//total system void
//total system void minus filled junctions

void QuasiGrid::set_system_void_spaces() {
}
std::map<std::string,double>  QuasiGrid::get_system_void_spaces() const{

  return this->system_void_spaces_;
}

void QuasiGrid::set_system_length_scales() {
  //particle length
  //squeezed middle lobe length scale
  //total system void minus filled junctions scale
}

std::map<std::string,double> QuasiGrid::get_system_length_scales() const{
  return this->system_length_scales_;
}
#endif
