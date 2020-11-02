#ifndef QUASI_GRID_HPP
#define QUASI_GRID_HPP

#include <string_view>
#include "particle.hpp"
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
#include <cxxabi.h>
#include <json/json.h>
#include "intersection.hpp"
using USH = unsigned short;
using UIN = unsigned int;
using  JunctionRef = std::reference_wrapper<Junction> ;
///This is a proposal to wrap a junction ref in a container

/// @brief print_tuple :
//
///
/// @tparam TupType
/// @tparam I
/// @param _tup
/// @param std::index_sequence
template<class TupType, size_t... I>
void print_tuple(const TupType& _tup, std::index_sequence<I...>) {
  std::cout << "(";
  (..., (std::cout << (I == 0? "" : ", ") << std::get<I>(_tup)));
  std::cout << ")\n";
}


/**
 *
 **/
  template<class... T>
void print_tuple (const std::tuple<T...>& _tup)
{
  print_tuple(_tup, std::make_index_sequence<sizeof...(T)>());
}
///Junction is the unit of intersection
// -------------------------------
/// @Synopsis  <QuasiGrid Base Class>
// ---------------------------------
class QuasiGrid{
  public:
    enum Direction {
      Forward, /// \brief differentiates between forward and backward movement
      Backward
    };
    enum StateLabels {
      EndLobeLocation,
      JunctionLocation,
      MiddleLobeLocation
    };
    enum SimulationMode {///\brief defines  the mode of simulation
      FIXED_STEP_SIZE,
      ADAPTIVE_STEP_SIZE
    };
    enum ParticleFillingStrategy{
      direct_from_input,
      above_below_K_g,
      above_below_K_c
    };
    enum ParticleFillingMode{
      below_N_g,
      above_N_g,
      below_N_c,
      above_N_c
    };

    enum DensityMode{
      below_phi_g,
      above_phi_g,
      below_phi_c,
      above_phi_c
    };

    enum DensityDefinition{
      direct_from_input_phi,///Taken as a direct input
      above_below_phi_g_direct,
      above_below_phi_g_logarithmic, /// e.g. phi_g - 10^-10
      above_below_phi_c_logarithmic
    };

    ///ctor and the rule of 5
    QuasiGrid();///\brief default ctor
    QuasiGrid(const QuasiGrid&) = delete; ///\brief copy ctor deleted to avoid extraneous copies
    QuasiGrid& operator=(QuasiGrid const& ) = delete; ///\brief same as above
    QuasiGrid& operator=(QuasiGrid&& ) = delete; ///\brief same as above
    void init();
    void un_init();


    std::string calculate_demangled_classname();
    /// @brief create filename at ease
    /// @tparam Args
    /// @param 
    /// @returns void
    template <typename... Args>
      static decltype(auto) filename_creator(Args&&... args);
    /// \brief debug print service
    template <typename... Args>
      static void debug_print_service(Args&&... args);
    /// \brief error print service
    template <typename... Args>
      static  void error_print_service(Args&&... args);
    /// \brief log service
    template <typename... Args, typename T>
      static void log_service(Args&&... args, T& file_stream);
    /**
     *
     *      STUFF about Quasi system    
     *
     */

    UIN get_num_particle_size_voids_at_kinetic_arrest() const;
    template<typename T>
      void set_num_particle_size_voids_at_kinetic_arrest(const T) ;
    template<typename T>
      void set_particle_len_at_complete_kinetic_arrest(const T);
    double  get_particle_len_at_complete_kinetic_arrest() const;
    double  calculate_particle_len_at_complete_kinetic_arrest() const;
    void  calculate_and_set_particle_len_at_complete_kinetic_arrest() const;
    void calculate_and_set_num_particle_size_voids_at_kinetic_arrest() ;
    UIN calculate_num_particle_size_voids_at_kinetic_arrest() const;
    template<typename T>
      void print_map(const T t);
    virtual ~QuasiGrid();
    void time_update();
    double calculate_mid_lobe_location(const int index) const;///\brief calculates the location of the midlobe
    void set_system_length_scales();///\brief set the different length scales of the system
    void set_system_void_spaces();///\brief set the different void spaces of the system
    std::map<std::string,double> get_system_length_scales() const;
    std::map<std::string,double> get_system_void_spaces() const;
    double get_mid_lobe_location(const int&);
    std::tuple<USH,USH> which_region_particle_belong_to(const Particle& p);
    int does_belong_to_end_lobe(const Particle& p);
    int does_belong_to_junction(const Particle& p);
    const Junction* does_belong_to_junction_location(const double&) const;
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
    inline double mod_len(const double x) { if(x>=0) return std::fmod(x,system_len_); return (system_len_ - std::fmod(std::fabs(x),system_len_)); }
    /**
     *
     **/
    void move_particle(const int& n,const double& dx);
    void distribute_junction_labels();
    void reassign_intersection_labels();///\brief assign the intersection labels starting from zero go max intersections
    int get_neighbor_particle(const int& n,const int direction);
    int mod_num_particles(const int& n);///\brief get particle as a module of total number of particles
    inline double get_random()   { return  normal_dist_(particle_motion_generator_); }
    /**
     *
     **/
    inline double get_random_particle()   { return  particle_selection_dist_(particle_selection_generator_); }
    /**
     *
     **/
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

    /**
     * Helper Functions
     */

    template<typename T>
      Json::Value  read_json(const T& filaneme );
    /// @brief parse_config_json
    void parse_config_json();
    template<typename T>
      static bool is_empty_string(const T& input) ;
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
    void increment_simulation_time_elapsed(const  double incr = 1.0);
    void increment_simulation_steps_elapsed(const UIN incr = 1.0);


    template<typename T>
      decltype(auto) set_log_file_name(const T& log_file_name);
    decltype(auto) get_log_file_name() const;
    template<typename T>
      decltype(auto) set_data_file_name(const T& data_file_name);
    decltype(auto) get_data_file_name() const;

    template<typename T>
      decltype(auto) set_output_file_name(const T& output_file_name);
    decltype(auto) get_output_file_name() const;

    template<typename T>
      void set_output_folder_root(const T&);
    decltype(auto) get_output_folder_root() const;


    template<typename T>
      decltype(auto) set_runtime_config_file_name(const T& runtime_config_file_name);
    decltype(auto) get_runtime_config_file_name() const ;

    template<typename T>
      decltype(auto) set_config_file_name(const T& config_file_name) ;
    decltype(auto) get_config_file_name() const;

  private:
    void particle_length_dependencies();
    void quasigrid_helper();
    static std::string class_name_;
    /*************** Dynamics and Distribution...later to be moved to a class ********/
    std::random_device particle_motion_{};
    std::mt19937 particle_motion_generator_{particle_motion_()};
    std::normal_distribution<> normal_dist_{0,1};
    std::random_device particle_selection_{};
    std::mt19937 particle_selection_generator_{particle_selection_()};
    std::uniform_int_distribution<> particle_selection_dist_{1, 6};
/**
 * system params
 *
 */
    int num_particles_; ///\brief N parameter
    int num_intersections_; ///\brief J parameter
    static USH num_particles_per_middle_arm_at_max_packing_;///\brief K parameter in the paper
    double density_;
    double density_close_packing_; ///\brief density at closed packing
    double density_kinetic_arrest_; ///\brief Density at the kinetic arrest

    double density_delta_log_scale_;
    double density_delta_;
    int num_intersections_half_;
    int num_junctions_;
    int num_junctions_half_;
    int num_particles_possible_in_system_at_closed_packing_;
    int num_particle_size_voids_at_kinetic_arrest_;
    /**
     * Length scales
     */
    double system_len_;
    double particle_len_;
    double particle_len_at_critical_packing_;
    double particle_len_half_;
    double system_len_half_;
    double particle_len_at_kinetic_arrest_;
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
    static UIN num_particles_at_closed_packing_;///\brief $N^{g}$ parameter in the paper
    static UIN num_particles_possible_in_system_;///\brief $N^{total}$ parameter in the paper
    static UIN min_num_particles_at_kinetic_arrest_;///\brief min number of particles when the sytem can move and also KA happens
    static UIN no_of_particles_above_min_no_particles_at_kinetic_arrest_;
    static UIN max_num_particles_at_kinetic_arrest_;///\brief min number of particles when the sytem can move and also KA happens
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
    /**
     * System state
     */
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
    static short int simulation_mode_;
    static double simulation_time_elapsed_;
    static long long simulation_steps_elapsed_;
    static std::string length_scale_tag_particle_length_;
    static std::string length_scale_tag_arm_length_;
    static std::string length_scale_tag_end_lobe_length_;
    static std::string length_scale_tag_middle_lobe_length_;
    static std::string void_space_tag_total_system_wide_void_;
    static std::string void_space_tag_total_junction_corrected_system_wide_void_;
    static std::string void_space_tag_squeezed_end_lobe_void_;
    static std::string void_space_tag_squeezed_middle_lobe_void_;
    static std::string void_space_tag_average_squeezed_end_lobe_void_;
    static std::string void_space_tag_average_squeezed_middle_lobe_void_;
    static double particle_len_at_complete_kinetic_arrest_;
    static std::string log_file_name_;
    static std::string data_file_name_;
    static std::string output_file_name_;
    static std::string runtime_config_file_name_;
    static std::string config_file_name_;

    static std::string home_directory_;
    static std::string_view constexpr log_file_name_default_ = {"quasi_log.txt"};
    static std::string_view constexpr data_file_name_default_ = {"quasi_data.txt"};
    static std::string_view constexpr output_file_name_default_  = {"quasi_data_output.txt"};
    static std::string_view constexpr runtime_config_file_name_default_ = {"quasi_runtime_config_file_"};
    static std::string_view constexpr config_file_name_default_ = {"config.json"};
    static Json::Reader json_reader_;
    static bool is_config_json_set_;
    static   Json::Value config_json_root_;
    static std::string output_folder_root_;
    static std::string_view constexpr output_folder_root_default_ = {"output"};
};
std::string QuasiGrid::output_folder_root_ = {std::string(output_folder_root_default_)};
bool QuasiGrid::is_config_json_set_ = {false};

Json::Value QuasiGrid::config_json_root_;
/// KEY
std::string QuasiGrid::log_file_name_ = {"quasi_log.txt"};
std::string QuasiGrid::data_file_name_ = {"quasi_data.txt"};
std::string QuasiGrid::output_file_name_ = {"quasi_data_output.txt"};
std::string QuasiGrid::runtime_config_file_name_ = {"quasi_runtime_config_file_"};
std::string QuasiGrid::config_file_name_ = {"config.json"};


std::string QuasiGrid::home_directory_ = {""};

double QuasiGrid::particle_len_at_complete_kinetic_arrest_ = {0.0};
std::string QuasiGrid::length_scale_tag_particle_length_ = {"particle_length"};
std::string QuasiGrid::length_scale_tag_arm_length_= {"arm_length"};
std::string QuasiGrid::length_scale_tag_end_lobe_length_= {"end_lobe_length"};
std::string QuasiGrid::length_scale_tag_middle_lobe_length_ = {"middle_lobe_length"};
std::string QuasiGrid::void_space_tag_total_system_wide_void_ = {"system_wide_void"};
std::string QuasiGrid::void_space_tag_total_junction_corrected_system_wide_void_ = {"junction_corrected_system_wide_void"};
std::string QuasiGrid::void_space_tag_squeezed_end_lobe_void_ = {"squeezed_end_lobe_void"};
std::string QuasiGrid::void_space_tag_squeezed_middle_lobe_void_ = {"squeezed_middle_lobe_void"};
std::string QuasiGrid::void_space_tag_average_squeezed_end_lobe_void_ = {"average_squeezed_end_lobe_void"};
std::string QuasiGrid::void_space_tag_average_squeezed_middle_lobe_void_ = {"average_squeezed_middle_lobe_void"};
long long QuasiGrid::simulation_steps_elapsed_ = {0};
double QuasiGrid::simulation_time_elapsed_ = {0.0};
short int QuasiGrid::simulation_mode_ = {QuasiGrid::SimulationMode::FIXED_STEP_SIZE};
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
//UIN QuasiGrid::max_no_of_particles_at_kinetic_arrest_ = {0};
UIN QuasiGrid::no_of_particles_above_min_no_particles_at_kinetic_arrest_ = {0};
UIN QuasiGrid::range_of_particles_between_min_to_max_at_kinetic_arrest_ = {0};
UIN QuasiGrid::num_particles_at_closed_packing_ = {0};///\brief $N^{cp}$ parameter in the paper
USH QuasiGrid::num_particles_per_middle_arm_at_max_packing_ = {0};
USH QuasiGrid::max_allowed_particles_in_end_lobes_at_max_packing_ ={0};
double QuasiGrid::arm_void_length_ = NAN;
/***
 *
 */
void QuasiGrid::increment_simulation_time_elapsed(const double incr ){
  simulation_time_elapsed_ += std::fabs(incr);
}
/**
 *
 **/
void QuasiGrid::increment_simulation_steps_elapsed(const UIN incr){
  this->simulation_steps_elapsed_ += incr;
}
/**
 *
 **/
void QuasiGrid::populate_junction_conjugates(){
  for(const auto& intersection:intersections_){
    auto& intersection_ref = intersection.get_intersection_ref();
    auto& conjugate_map= intersection.get_junction_conjugate();
    for(const auto& pair:conjugate_map){
      std::cout << pair.first << " " << pair.second << std::endl;
      junction_conjugates_.emplace(pair);
    }
    /**
     *
     **/
  }
  /**
   *
   **/
}
void QuasiGrid::init(){
  if(is_config_json_set_){
    config_json_root_ = read_json(config_file_name_);
  }else {
    std::cout << "sorry config json is not set...quit()" << std::endl;
    exit(0);
  }
  parse_config_json();

}
/**
 *
 **/
double QuasiGrid::arm_length_ = {0.0};
int QuasiGrid::default_arm_offset_ = {2};
std::string QuasiGrid::class_name_ = {""};
QuasiGrid::QuasiGrid(){
  quasigrid_helper();
  std::cerr << class_name_ << " ctor" << std::endl;

}
/**
 *
 **/
QuasiGrid::~QuasiGrid(){
  std::cerr << class_name_<< " dtor" << std::endl;
}
/**
 *
 **/
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
/**
 *
 **/
void QuasiGrid::set_num_particles_at_closed_packing(const int& n){
  num_particles_at_closed_packing_ = n;
}
/**
 *
 **/
unsigned int QuasiGrid::get_num_particles_at_closed_packing() const {
  return num_particles_at_closed_packing_;
}
/**
 *
 **/
UIN QuasiGrid::calculate_num_particles_at_closed_packing() const{
  //  return  2*(num_intersections_ + 1)*(num_particles_per_middle_arm_at_max_packing_ + 1) - num_intersections_;
  return  2*(num_particles_per_middle_arm_at_max_packing_ + 1) + num_intersections_*(2*num_particles_per_middle_arm_at_max_packing_ + 1) ;
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_num_particles_at_closed_packing() {
  this->set_num_particles_at_closed_packing(std::move(this->calculate_num_particles_at_closed_packing()));
}
/**
 *
 **/
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
  /**
   *
   **/
  return density;
}
/**
 *
 **/
void QuasiGrid::set_num_particles( const int n) {
  num_particles_ = n;
}
/**
 *
 **/
int QuasiGrid::get_num_particles() const {
  return num_particles_;
}
/**
 *
 **/
void QuasiGrid::set_num_intersections( const int n){
  num_intersections_ = n;
  num_intersections_half_ = num_intersections_/2;
  num_junctions_half_  = num_intersections_;
}
/**
 *
 **/
UIN QuasiGrid::get_num_intersections() const {
  return num_intersections_;
}
/**
 *
 **/
void QuasiGrid::set_num_junctions( const int n){
  num_junctions_ = n;
  num_junctions_half_ = num_junctions_/2;
}
/**
 *
 **/
UIN QuasiGrid::get_num_junctions() const {
  return num_junctions_;
}
/**
 *
 **/
void QuasiGrid::set_particle_selection_distribution(){
  std::cout << "setting particle selection dist" << std::endl;
  particle_selection_dist_ = std::uniform_int_distribution<> (0, num_particles_);
}
/**
 *
 **/
/**
 * distribute junction labels to each intersections
 */
void QuasiGrid::distribute_junction_labels(){
  std::cout << "distributing labels:" << std::endl;
  for(auto& it:intersections_){
    it.allocate_labels(0,it.get_label());
    it.allocate_labels(1,2*intersections_.size() -1 -it.get_label());
  }
  /**
   *
   **/
  print_intersection_labels();
}
/**
 *
 **/
void QuasiGrid::print_intersection_labels(){
  std::cerr << "printing intersection labels:" << std::endl;
  for(auto& it:intersections_){
    std::cout << it.get_label() << std::endl;
  }
  /**
   *
   **/
  std::cout << std::endl;
}
/**
 *
 **/
void QuasiGrid::print_junction_labels(){
  std::cerr << "printing junction labels:" << std::endl;
  for(auto& it:junctions_){
    std::cout << it->get_label() << " ";
  }
  /**
   *
   **/
  std::cout << std::endl;
}
/**
 *
 **/
void QuasiGrid::print_junction_coordinates(){
  std::cerr << "printing intersection coordinates:" << std::endl;
  for(auto& it:junctions_){
    std::cout << it->get_x() << " ";
  }
  /**
   *
   **/
  std::cout << std::endl;
}
/**
 *
 **/
void QuasiGrid::initialize_system() {///Initialize the system
  std::cerr << "initializing system " <<  std::endl;
  particles_ = std::vector<Particle>(num_particles_);
  intersections_ = std::move(std::vector<Intersection>(num_intersections_));
  distribute_junction_labels();
  for(auto& it:intersections_){
    auto& intersection_ref = it.get_intersection_ref();
    for(auto& its:intersection_ref) {
      junctions_.push_back(const_cast<Junction*>(&its));
    }
    /**
     *
     **/
    it.set_junction_conjugate();
    it.print_conjugate();
  }
  /**
   *
   **/
  populate_junction_conjugates();///For each junction populate the conjugate junction
  std::cout << std::endl;
  std::sort(junctions_.begin(),junctions_.end(),[](const auto& a,const auto&b){ return a->get_label()<b->get_label();});
  std::cerr << "printing junction labels" << std::endl;
  for(auto& jun:junctions_){
    std::cout << jun->get_label() << " ";
  }
  /**
   *
   **/
  std::cout << std::endl;
  std::cout << "checking is blocked" << std::endl;
  for(auto& it:junctions_){
    this->is_blocked(*junctions_.at(it->get_label()));
  }
  calculate_num_particles_possible_in_system();
  calculate_and_set_max_num_particles_at_kinetic_arrest();
  calculate_and_set_num_particle_size_voids_at_kinetic_arrest();
  calculate_and_set_particle_len_at_complete_kinetic_arrest();
}
/**
 *
 **/
void QuasiGrid::populate_blocked_intersection(const int n){
  blocked_intersections_.push_back(n);
}
/**
 *
 **/
void QuasiGrid::depopulate_blocked_intersection(const int n){
  blocked_intersections_.remove(n);
}
/**
 *
 **/
void QuasiGrid::reassign_intersection_labels() {//Use it optionally if the intersection labels need reassignments
}
/**
 *
 **/
int QuasiGrid::mod_num_particles(const int& n){
  if(n>=0){
    return n%num_particles_;
  }
  /**
   *
   **/
  return num_particles_ - ((-n)%num_particles_);
}
/**
 *
 **/
int QuasiGrid::get_neighbor_particle(const int& n,const int direction){///zero is backward, 1 is forward
  if(direction == Direction::Forward) {
    return mod_num_particles(n+1);
  }
  /**
   *
   **/
  return mod_num_particles(n-1);
}
/**
 *
 **/
///Find the blocked junction closer than the next particle
int QuasiGrid::get_next_junction_index(const double& x){
  ///return the next junction with respect to a given point
  auto  last_junction_Index = 1;
}
/**
 *
 **/
///Check if there is a blocked junction before or after the next particle
void QuasiGrid::move_particle(const int& n,const double& dx){
  short int sign_of_dx = dx>0?1:-1;
  bool is_forward_move = sign_of_dx ==1?true:false;
  auto x = particles_.at(n).get_x();
  auto x_next = x + sign_of_dx*dx;
  auto neighbor_particle  = mod_num_particles(n + sign_of_dx);
  if(!particles_.at(neighbor_particle).does_belong_to(x_next)){///Particle constraint is cleared now the junction constraint
    bool is_blocked_by_junction = {false};
    bool does_move_involve_junction = {false};
    for(auto it:intersections_){
      auto& intersection_ref = it.get_intersection_ref();
      for(const Junction& its:intersection_ref){///its is a junction
        bool did_belong_to_junction = its.does_belong_to(x);
        bool does_belong_to_junction = its.does_belong_to(x_next);
        if(did_belong_to_junction || does_belong_to_junction)
          does_move_involve_junction = true;
        if(does_move_involve_junction){
          if(did_belong_to_junction) {
            if(does_belong_to_junction) {//regular update
            } else {//exit from junction
#warning ENSURE NEW OBJECTS ARE NOT BEING CREATED HERE
              std::remove_const<Junction>::type its_non_const(its);
              its_non_const.unoccupy(is_forward_move);
            }
            /**
             *
             **/
          }else {
            if(does_belong_to_junction){ //new belongning to junction
              std::remove_const<Junction>::type junction_non_const(its);
              junction_non_const.occupy(&particles_.at(n),is_forward_move);
            } else {///neighter belonged to junction nor it does now!
              //simple update with regular jump
            }
            /**
             *
             **/
          }
          /**
           *
           **/
        }else {
          ///simple jump here
        }
        /**
         *
         **/
      }///junction loop per intersection
    }///intersection loop
  }else {//nothing to do as blocked by neighbors
  }
  /**
   *
   **/
  time_update();
}
/**
 *
 **/
void QuasiGrid::time_update(){
}
/**
 *
 **/
void QuasiGrid::set_particle_len(const double& len){
  particle_len_ = len;
}
/**
 *
 **/
double QuasiGrid::get_particle_len()const{
  return  particle_len_;
}
/**
 *
 **/
void QuasiGrid::set_system_len(const double& len){
  system_len_ = len;
  std::cerr << "setting system length " << system_len_ << std::endl;
  system_len_half_ = system_len_/2.0;
  std::cerr << "setting system_len_half_ " << system_len_half_ << std::endl;
}
/**
 *
 **/
double QuasiGrid::get_system_len() const {
  return  system_len_;
}
/**
 *
 **/
/**
 *@returns ptr to junction else nullptr
 */
const Junction* QuasiGrid::does_belong_to_junction_location(const double& x) const {
  for(auto& it:junctions_){
    if(it->does_belong_to(x))
      return it;
  }
  /**
   *
   **/
  return nullptr;
}
/**
 *
 **/
void QuasiGrid::print_all_intersection_length(){
  std::cerr << "printing all intersection length:" << std::endl;
  for(auto& it:intersections_){
    auto& ref = *it.get_intersection_ptr();
    std::cout << "intersection " << it.get_label() << ": ";
    for(auto& its:ref){
      std::cout << its.get_len() << " ";
    }
    /**
     *
     **/
    std::cout << std::endl;
  }
  /**
   *
   **/
}
/**
 *
 **/
void QuasiGrid::set_all_intersection_length()const {
  std::cerr << "setting all junction length" << std::endl;
  for(const auto& it:junctions_){
    (const_cast<Junction*>(it))->set_len(this->particle_len_);
  }
  /**
   *
   **/
}
/**
 *
 **/
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
    /**
     *
     **/
    (const_cast<Junction*>(it))->set_x(x);
  }
  /**
   *
   **/
  std::cout << std::endl;
  std::cout << std::endl;
}
/**
 *
 **/
/**
 *
 */
Junction& QuasiGrid::get_junction(const int n) {
  return *junctions_.at(n-1);
}
/**
 *
 **/
Intersection& QuasiGrid::get_intersection(const int n){
  if(n<=intersections_.size() && n>=1)
    return intersections_.at(n-1);
  else{
    std::cerr << "intersection index " << n << " out of range " << std::endl;
    throw "index out of range";
  }
  /**
   *
   **/
}
/**
 *
 **/
void QuasiGrid::set_arm_length() {
  arm_length_ = system_len_ / (2.0*(num_intersections_ + 1));
  std::cerr << "arm_length:" << arm_length_ << std::endl;
}
/**
 *
 **/
double QuasiGrid::get_arm_length() const{
  return arm_length_;
}
/**
 *
 **/
bool QuasiGrid::is_blocked(const Junction& j) const{
  std::cout << "index:" << index  << std::endl;
  auto matches  = junction_conjugates_.equal_range(j.get_label());
  // Iterate over the range
  for (auto it = matches.first; it != matches.second; it++){
    //std::cout << it->second << std::endl;
    if((junctions_.at(it->second))->is_occupied()){
      return true;
    }
    /**
     *
     **/
  }
  /**
   *
   **/
  return false;
}
/**
 *
 **/
void QuasiGrid::occupy_junction(const int& junction_index, Particle * const& particle,const bool& is_forward_direction)const{
  const_cast<Junction *>(this->junctions_.at(junction_index))->occupy(particle,is_forward_direction);
}
/**
 *
 **/
void QuasiGrid::set_num_particles_per_middle_arm_at_max_packing(const int& k){
  num_particles_per_middle_arm_at_max_packing_ = k;
}
/**
 *
 **/
int QuasiGrid::get_min_num_particles_at_kinetic_arrest() const {
  return this->min_num_particles_at_kinetic_arrest_;
}
/**
 *
 **/
void QuasiGrid::set_min_num_particles_at_kinetic_arrest(const int& n)  {
  this->min_num_particles_at_kinetic_arrest_ = n;
}
/**
 *
 **/
int QuasiGrid::calculate_min_num_particles_at_kinetic_arrest() const{ ///This is the minimum number of particles where kinetic arrest will happen
  if(num_intersections_ > 1 ){
    auto val =  2*(max_allowed_particles_in_end_lobes_at_max_packing_ -1 -1) + num_intersections_ * num_particles_per_middle_arm_at_max_packing_;
    return val;
  }else {
    std::cerr << "sorry this can be calculated only for num_intersections > 1" << std::endl;
    throw "num intersections > 1 ";
  }
  /**
   *
   **/
  return 0;
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_min_num_particles_at_kinetic_arrest()  {
  this->set_min_num_particles_at_kinetic_arrest(std::move(this->calculate_min_num_particles_at_kinetic_arrest()));
}
/**
 *
 **/
void QuasiGrid::set_max_allowed_particles_in_end_lobes_at_max_packing( const int& val) {
  this->max_allowed_particles_in_end_lobes_at_max_packing_ = val;
}
/**
 *
 **/
int QuasiGrid::calculate_max_allowed_particles_in_end_lobes_at_max_packing() const{
  return 2*num_particles_per_middle_arm_at_max_packing_ + 1;
}
/**
 *
 **/
int QuasiGrid::get_max_allowed_particles_in_end_lobes_at_max_packing() const{
  return this->max_allowed_particles_in_end_lobes_at_max_packing_;
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_max_allowed_particles_in_end_lobes_at_max_packing(){
  auto n(std::move(calculate_max_allowed_particles_in_end_lobes_at_max_packing()));
  set_max_allowed_particles_in_end_lobes_at_max_packing(n);
}
/**
 *
 **/
int QuasiGrid::get_num_particles_per_middle_arm_at_max_packing() const{
  return this->num_particles_per_middle_arm_at_max_packing_;
}
/**
 *
 **/
void QuasiGrid::set_num_particles_above_min_no_particles_at_kinetic_arrest(const int& n) {
  this->no_of_particles_above_min_no_particles_at_kinetic_arrest_ = n;
}
/**
 *
 **/
UIN QuasiGrid::get_num_particles_above_min_no_particles_at_kinetic_arrest() const{
  return no_of_particles_above_min_no_particles_at_kinetic_arrest_;
}
/**
 *
 **/
UIN QuasiGrid::calculate_num_particles_from_filling_mode() {
  return      (this->min_num_particles_at_kinetic_arrest_ + this->no_of_particles_above_min_no_particles_at_kinetic_arrest_);
}
/**
 *
 **/
///K param defines a filling mode
UIN QuasiGrid::calculate_and_set_num_particles_from_filling_mode(){
  this->set_num_particles(std::move(this->calculate_num_particles_from_filling_mode()));
}
/**
 *
 **/
void QuasiGrid::set_density(const int& phi){
  this->density_ = phi;
}
/**
 *
 **/
double QuasiGrid::get_density() const{
  return this->density_;
}
/**
 *
 **/
double QuasiGrid::calculate_density_from_density_delta_and_kinetic_arrest_density() const {
  return density_kinetic_arrest_ - std::pow(10.0,-density_delta_log_scale_);
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_density_from_density_delta_and_kinetic_arrest_density()  {
  this->density_ = std::move(this->calculate_density_from_density_delta_and_kinetic_arrest_density());
}
/**
 *
 **/
double QuasiGrid::calculate_density() const {
  return static_cast<double>(num_particles_*particle_len_)/system_len_;
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_density() {
  this->density_ = std::move(this->calculate_density());
}
/**
 *
 **/
UIN QuasiGrid::get_num_particles_possible_in_system() const{
  return this->num_particles_possible_in_system_;
}
/**
 *
 **/
//void QuasiGrid::set_num_particles_at_closed_packing(const int& n){
//  this->num_particles_possible_in_system_at_closed_packing_ = n;
//}
/**
 *
 **/
void QuasiGrid::calculate_and_set_num_particles_possible_in_system(){
  std::cout << "calculate_num_particles_possible_in_system:" << num_particles_possible_in_system_ << std::endl;
  this->num_particles_possible_in_system_  = std::move(this->calculate_num_particles_possible_in_system() );
  //this->max_no_of_particles_at_kinetic_arrest_ = std::move(this->calculate_num_particles_possible_in_system_at_closed_packing() );
}
/**
 *
 **/
UIN QuasiGrid::calculate_num_particles_possible_in_system() const{
  return 2*(2*num_particles_per_middle_arm_at_max_packing_ + 1 ) + 2*num_intersections_ + 2*(num_intersections_-1)*num_particles_per_middle_arm_at_max_packing_;
}
/**
 *
 **/
void QuasiGrid::print_system() const{
  std::cout << "J:" << num_intersections_ << " K:" << num_particles_per_middle_arm_at_max_packing_  << " N:" << num_particles_ << std::endl;
  std::cout << "N_min:" << min_num_particles_at_kinetic_arrest_ << " N_max:" << max_num_particles_at_kinetic_arrest_  << std::endl;
  std::cout << "N_above_min:" <<  no_of_particles_above_min_no_particles_at_kinetic_arrest_ << " N_voids:" << num_particle_size_voids_at_kinetic_arrest_  << std::endl;
}
/**
 *
 **/
void QuasiGrid::set_max_num_particles_at_kinetic_arrest( const UIN n){
  max_num_particles_at_kinetic_arrest_ = n;
}
/**
 *
 **/
UIN QuasiGrid::get_max_num_particles_at_kinetic_arrest()const{
  return max_num_particles_at_kinetic_arrest_;
}
/**
 *
 **/
UIN QuasiGrid::calculate_max_num_particles_at_kinetic_arrest() const{
  return num_particles_at_closed_packing_ - this->min_particle_size_voids_needed_for_kinetic_arrest_;
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_max_num_particles_at_kinetic_arrest(){
  this->max_num_particles_at_kinetic_arrest_ = std::move(this->calculate_max_num_particles_at_kinetic_arrest());
  std::cout << "max_num_particles_at_kinetic_arrest:" << max_num_particles_at_kinetic_arrest_ << std::endl;
}
/**
 *
 **/
void QuasiGrid::set_arm_void_length(const double& len){
  arm_void_length_ = len;
}
/**
 *
 **/
double QuasiGrid::get_arm_void_length() const{
  return arm_void_length_ ;
}
/**
 *
 **/
double QuasiGrid::calculate_arm_void_length() const{
  return arm_length_ - 2.0*particle_len_;
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_arm_void_length(){
  this->arm_void_length_ = std::move(this->calculate_arm_void_length());
}
/**
 *
 **/
void QuasiGrid::set_density_close_packing(const double& density){
  density_close_packing_ = density;
}
/**
 *
 **/
double QuasiGrid::get_density_close_packing() const{
  return this->density_close_packing_;
}
/**
 *
 **/
void QuasiGrid::set_density_kinetic_arrest(const double& density){
  density_kinetic_arrest_ = density;
}
/**
 *
 **/
double QuasiGrid::get_density_kinetic_arrest() const{
  return this->density_kinetic_arrest_;
}
/**
 *
 **/
double QuasiGrid::calculate_density_kinetic_arrest(){//phi_g
  return static_cast<double>(num_particles_)/num_particles_possible_in_system_;
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_density_kinetic_arrest(){//phi_g
  this->set_density_kinetic_arrest(std::move(this->calculate_density_kinetic_arrest()));
}
/**
 *
 **/
void QuasiGrid::set_density_delta_log_scale(const double& density){
  this->density_delta_log_scale_ = density;
}
/**
 *
 **/
double QuasiGrid::get_density_delta_log_scale()const{
  return  this->density_delta_log_scale_;
}
/**
 *
 **/
double QuasiGrid::calculate_particle_len()const{
  return  (density_*system_len_)/num_particles_;
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_particle_len(){
  this->particle_len_ = std::move(this->calculate_particle_len());
  particle_length_dependencies();
}
/**
 *
 **/
void QuasiGrid::set_trajectory_output_filename(const std::string& filename){
  this->trajectory_output_filename_ = filename;
}
/**
 *
 **/
std::string QuasiGrid::get_trajectory_output_filename() const{
  return this->trajectory_output_filename_;
}
/**
 *
 **/
void QuasiGrid::set_trajectory_input_filename(const std::string& filename){
  this->trajectory_input_filename_ = filename;
}
/**
 *
 **/
std::string QuasiGrid::get_trajectory_input_filename() const{
  return this->trajectory_input_filename_;
}
/**
 *
 **/
UIN QuasiGrid::calculate_system_state_size() const{
  std::cout << "num_junctions_:" << num_intersections_ << " " << 2 + num_intersections_ + 2*num_intersections_ << std::endl;
  return 2 + num_intersections_ + 2*(num_intersections_ -1);
}
/**
 *
 **/
void QuasiGrid::set_system_state_size( const UIN& size){
  this->system_state_size_ = size;
}
/**
 *
 **/
UIN QuasiGrid::get_system_state_size() const{
  return this->system_state_size_;
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_system_state_size(){
  this->system_state_size_ = std::move(this->calculate_system_state_size());
}
/**
 *
 **/
UIN QuasiGrid::calculate_lower_end_lobe_occupation() const{
}
/**
 *
 **/
UIN QuasiGrid::calculate_upper_end_lobe_occupation() const{
}
/**
 *
 **/
void QuasiGrid::AssignStateLablesToParticles() {
#ifdef  VERBOSE
  std::cerr << "in function " << __func__ << std::endl;
#endif
  for(auto& it:particles_){
    std::cerr << "particle label:"<< it.get_label() << std::endl;
  }
  /**
   *
   **/
}
/**
 *
 **/
std::vector<UIN> QuasiGrid::calculate_system_state() {
  std::cerr << "calculating system states for " << num_particles_ << std::endl;
  decltype(system_state_) system_state(this->system_state_size_,0);
  for( auto& particle:particles_){
    auto region = which_region_particle_belong_to(particle);
    auto region_label = std::get<0>(region);
    auto region_index = std::get<1>(region);
    short vector_location = {0};
    switch(region_label){
      case StateLabels::EndLobeLocation:
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
        /**
         *
         **/
        break;
      case StateLabels::MiddleLobeLocation:
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
        /**
         *
         **/
        break;
      default:
        std::cout << "unknown region label" << std::endl;
        throw "unknown region label";
        break;
    }
    /**
     *
     **/
  }
  /**
   *
   **/
  std::cout << "system state " << std::endl;
  for(const auto& it:system_state) {
    std::cout << it;
  }
  /**
   *
   **/
  std::cout << std::endl;
  ///First lower end lobe this.StateLabels.EndLobe with value 0
  //junctions
  //middle lobes
  ///last upper end lobe this.StateLabels.EndLobe with value 1
  return system_state;
}
/**
 *
 **/
void QuasiGrid::print_particles() const{
  std::cout << "label" << " " << " x " << " uncorrected x" << " state label " << " state label id" << std::endl;
  for(const auto& it:particles_) {
    std::cout << it.get_label() << " " << it.get_x() << " " << it.get_X() << " " << it.get_state_label() << " " << it.get_state_label_id() << std::endl;
  }
  /**
   *
   **/
}
/**
 *
 **/
void QuasiGrid::print_end_lobe_boundaries() const{
  std::cout << "lower end lobe boundaries" << std::endl;
  for(const auto& it:lower_end_lobe_boundaries_){
    print_tuple(it);
  }
  /**
   *
   **/
  std::cout << "upper end lobe boundaries" << std::endl;
  for(const auto& it:upper_end_lobe_boundaries_){
    print_tuple(it);
  }
  /**
   *
   **/
}
/**
 *
 **/
/**
 *
 */
std::tuple<USH,USH> QuasiGrid::which_region_particle_belong_to(const Particle& p){
  ///END LOBE CALCULATION
  int location;
  location = does_belong_to_end_lobe(p);
  std::tuple<USH,USH> belongs_to;
  if( location >= 0){
    belongs_to = decltype(belongs_to)(StateLabels::EndLobeLocation,location);
  }else {
    location = does_belong_to_junction(p);
    if(location >= 0){
      belongs_to = decltype(belongs_to)(StateLabels::JunctionLocation,location);
    }else {
      location = calculate_arm_index(p);
      belongs_to = decltype(belongs_to)(StateLabels::MiddleLobeLocation,location);
    }
    /**
     *
     **/
  }
  /**
   *
   **/
  return belongs_to;
}
/**
 * @returns {Junction index or -1 if it does not}
 **/
int QuasiGrid::does_belong_to_junction(const Particle& p) {
  //auto junction_ref = does_belong_to_junction(p.get_x());
  const double& x = p.get_x();
  auto junction_ref = does_belong_to_junction_location(x);
  if(junction_ref != nullptr){
    return junction_ref->get_label();
  }
  /**
   *
   **/
  /// pre defined garbage value. potentially can be defined through an enum value
  return -1;
}
/**
 *
 **/
int QuasiGrid::does_belong_to_middle_lobe(const Particle& p) {
}
/**
 *
 **/
int QuasiGrid::does_belong_to_end_lobe(const Particle& p) {
  auto x = p.get_x();
  for(const auto& it:lower_end_lobe_boundaries_){
    if(it == *lower_end_lobe_boundaries_.begin()){
      if(x>= std::get<0>(it) && x<std::get<1>(it) ){
        return 0;
      }
      /**
       *
       **/
    }else {
      if(x> std::get<0>(it) && x<=std::get<1>(it) ){
        return 0;
      }
      /**
       *
       **/
    }
    /**
     *
     **/
  }
  /**
   *
   **/
  std::cout << "upper end lobe boundaries" << std::endl;
  for(const auto& it:upper_end_lobe_boundaries_){
    if(x> std::get<0>(it) && x<std::get<1>(it) ){
      return 1;
    }
    /**
     *
     **/
  }
  /**
   *
   **/
  return -1;
}
/**
 *
 **/
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
 **/
/**
 *
 *dependencie s
 */
USH QuasiGrid::calculate_arm_index(const Particle& p){
  return static_cast<USH>(p.get_x()/arm_length_);
}
/**
 *
 **/
Particle&  QuasiGrid::get_particle(const USH& index) {
  if(index>0 && index <= particles_.size())
    return particles_.at(index-1);
  else {
    std::cerr << "sorry index is out of bound" << std::endl;
    throw "index out of bound";
  }
  /**
   *
   **/
}
/**
 *
 **/
std::vector<Particle>&  QuasiGrid::get_particles() {
  return particles_;
}
/**
 *
 **/
void QuasiGrid::set_num_arms(const UIN& num_arms){
  num_arms_ = num_arms;
  num_arms_half_ = num_arms/2;
}
/**
 *
 **/
UIN QuasiGrid::get_num_arms() const {
  return this->num_arms_;
}
/**
 *
 **/
UIN QuasiGrid::get_num_arms_half() const{
  return this->num_arms_half_;
}
/**
 *
 **/
UIN QuasiGrid::calculate_num_arms() const{
  return 2*(num_intersections_ + 1);
}
/**
 *
 **/
void QuasiGrid::calculate_and_set_num_arms(){
  this->set_num_arms(std::move(this->calculate_num_arms()));
}
/**
 * get the location of the arm based on 1 based arm index
 **/
double QuasiGrid::calculate_mid_lobe_location(const int index) const{
  std::cout << "arm_length:" << arm_length_ << std::endl;
  return arm_length_*(static_cast<double>(index) - 0.5);
}
/**
 *
 **/
//total system void
//total system void minus filled junctions
///Should be called before system length scale set function
void QuasiGrid::set_system_void_spaces() {
  system_void_spaces_[void_space_tag_total_system_wide_void_] = system_len_ -  num_particles_*particle_len_;
  system_void_spaces_[void_space_tag_total_junction_corrected_system_wide_void_] = system_void_spaces_[void_space_tag_total_system_wide_void_] - num_intersections_*particle_len_;
  system_void_spaces_[void_space_tag_squeezed_end_lobe_void_] = system_length_scales_[length_scale_tag_end_lobe_length_ ] - (2*num_particles_per_middle_arm_at_max_packing_ + 1)*particle_len_ ;
  system_void_spaces_[void_space_tag_squeezed_middle_lobe_void_] = system_length_scales_[length_scale_tag_middle_lobe_length_ ] - (num_particles_per_middle_arm_at_max_packing_*particle_len_);
  system_void_spaces_[void_space_tag_average_squeezed_end_lobe_void_] = system_void_spaces_[void_space_tag_squeezed_end_lobe_void_]/((2.0*num_particles_per_middle_arm_at_max_packing_ +1) + 1.0);
  system_void_spaces_[void_space_tag_average_squeezed_middle_lobe_void_] = system_void_spaces_[void_space_tag_average_squeezed_middle_lobe_void_]/(num_particles_per_middle_arm_at_max_packing_ + 1.0);
#ifdef DEBUG
  print_map(system_void_spaces_);
#endif
}
/**
 *
 **/
std::map<std::string,double>  QuasiGrid::get_system_void_spaces() const{
  return this->system_void_spaces_;
}
/**
 *
 **/
template<typename T>
void QuasiGrid::print_map(const T t){
  std::cout << "printing map:" << std::endl;
  for(const auto& it:t) {
    std::cout << it.first << " " << it.second << std::endl;
  }
}
//Should be set after particle length has been set
void QuasiGrid::set_system_length_scales() {
  //particle length
  //squeezed middle lobe length scale
  //total system void minus filled junctions scale
  system_length_scales_[length_scale_tag_particle_length_] = particle_len_;
  system_length_scales_[length_scale_tag_arm_length_] = arm_length_;
  system_length_scales_[length_scale_tag_end_lobe_length_ ] = 2.0*arm_length_;
  system_length_scales_[length_scale_tag_middle_lobe_length_ ] = arm_length_;
#ifdef DEBUG
  print_map(system_length_scales_);
#endif
}
/**
 *
 **/
std::map<std::string,double> QuasiGrid::get_system_length_scales() const{
  return this->system_length_scales_;
}
void QuasiGrid::particle_length_dependencies(){
  this->set_system_length_scales();
  this->set_system_void_spaces();
}
template<typename T>
void QuasiGrid::set_num_particle_size_voids_at_kinetic_arrest( const T n){
  num_particle_size_voids_at_kinetic_arrest_ = n;
}
UIN QuasiGrid::get_num_particle_size_voids_at_kinetic_arrest() const{
  return num_particle_size_voids_at_kinetic_arrest_;
}
void QuasiGrid::calculate_and_set_num_particle_size_voids_at_kinetic_arrest() {
  num_particle_size_voids_at_kinetic_arrest_ =   std::move(calculate_num_particle_size_voids_at_kinetic_arrest());
  std::cout << "num_particle_size_voids_at_kinetic_arrest:" << num_particle_size_voids_at_kinetic_arrest_ << std::endl;
  exit(0);
}
/***
 *
 */

UIN QuasiGrid::calculate_num_particle_size_voids_at_kinetic_arrest() const{
  std::cout << num_particles_possible_in_system_ << "  " <<  num_intersections_ << " " << num_particles_ << std::endl;
  return  num_particles_possible_in_system_  - num_intersections_ - num_particles_;
}
/***
 *
 */

template<typename T>
void QuasiGrid::set_particle_len_at_complete_kinetic_arrest(const T l){
  particle_len_at_complete_kinetic_arrest_ = l;
}
/***
 *
 */

double  QuasiGrid::get_particle_len_at_complete_kinetic_arrest() const{
  return particle_len_at_complete_kinetic_arrest_;
}
/***
 *
 */

double  QuasiGrid::calculate_particle_len_at_complete_kinetic_arrest() const{
  return system_len_/num_particles_possible_in_system_;
}
/***
 *
 */
void    QuasiGrid::calculate_and_set_particle_len_at_complete_kinetic_arrest() const{
  particle_len_at_complete_kinetic_arrest_  = std::move(calculate_particle_len_at_complete_kinetic_arrest());
}
/**
 *
 **/
std::string QuasiGrid::calculate_demangled_classname() {
  int status;
  std::string classname = std::string(abi::__cxa_demangle(typeid(*this).name(),0,0,&status));
  if(!status){
#ifdef DEBUG
    std::cerr << "demangled status " << status << std::endl;
#endif
  }else {
  }

  return classname;
}


template<typename T>
bool QuasiGrid::is_empty_string(const T& input) {
  if(std::all_of(input.begin(), input.end(), [](const auto& c){ return (c == ' ') || c == '\t' || c == '\n';}))
    return true;
}


template<typename T>
decltype(auto) QuasiGrid::set_log_file_name(const T& log_file_name) {
  log_file_name_ = log_file_name;
}
decltype(auto) QuasiGrid::get_log_file_name() const{

  return log_file_name_;
}
template<typename T>
decltype(auto) QuasiGrid::set_data_file_name(const T& data_file_name){

  data_file_name_ = data_file_name;
}
decltype(auto) QuasiGrid::get_data_file_name() const{
  return data_file_name_;
}

template<typename T>
decltype(auto) QuasiGrid::set_output_file_name(const T& output_file_name){
  output_file_name_ = output_file_name;

}

decltype(auto) QuasiGrid::get_output_file_name() const{

  return output_file_name_;
}
template<typename T>
decltype(auto) QuasiGrid::set_runtime_config_file_name(const T& runtime_config_file_name){
  runtime_config_file_name_ = runtime_config_file_name;
}

decltype(auto) QuasiGrid::get_runtime_config_file_name() const {
  return runtime_config_file_name_;

}

template<typename T>
decltype(auto) QuasiGrid::set_config_file_name(const T& config_file_name) {
  config_file_name_ = config_file_name;
  is_config_json_set_ = true;
}

decltype(auto) QuasiGrid::get_config_file_name() const {
  return config_file_name_;

}


/// @brief parse_config_json 
void QuasiGrid::parse_config_json(){
  if(!is_config_json_set_){
    //error_print_service("please set config.json file...exit");
    exit(0);
  }
  Json::Value json = config_json_root_;
  std::string key = "";
  // Get the value of the member of root named 'encoding',
  // and return 'UTF-8' if there is no such member.
  std::string json_encoding = "UTF-8";
  std::string encoding = json.get("encoding", json_encoding.c_str() ).asString();
  std::string project_name = json.get("project_name", json_encoding.c_str() ).asString();
  std::cout << "project_name:" << project_name << std::endl;
  key="output_folder_root";
  std::string output_folder_root = json.get(key,json_encoding.c_str()).asString();
  set_output_folder_root(output_folder_root);
  //  key = "is_multi_threaded";
  //  auto is_multi_threaded = json.get(key, json_encoding.c_str() ).asBool();
  //  set_is_multi_threaded(is_multi_threaded);
  //  debug_print_service("set_is_multi_threaded:",is_multi_threaded_);
  //  key = "is_recursive_processing";
  //  auto is_recursive_processing = json.get(key, json_encoding.c_str() ).asBool();
  //  set_is_recursive_processing(is_recursive_processing);
  //  debug_print_service("is_recursive_processing:",is_recursive_processing_);
  //  key = "image_kernel_radius";
  //  auto image_kernel_radius = json.get(key, json_encoding.c_str() ).asInt();
  //  set_image_kernel_radius(image_kernel_radius);
  //  debug_print_service("image_kernel_radius:",image_kernel_radius_);
  //  std::cout << "setting image kernel radius:" <<  image_kernel_radius_ << std::endl;
  //  key = "data_block_window_size_increment_step";
  //  auto data_block_window_size_increment_step = json.get(key, json_encoding.c_str() ).asInt();
  //  set_data_block_window_size_increment_step(data_block_window_size_increment_step);
  //  debug_print_service("is_recursive_processing:",is_recursive_processing_);
  //  key = "recursion_depth_max";
  //  auto recursion_depth_max = json.get(key, json_encoding.c_str() ).asInt();
  //  set_recursion_depth_max(recursion_depth_max);
  //  debug_print_service("recursion_depth_max:",recursion_depth_max);
  //  std::cout << "recursion_depth_max:" << recursion_depth_max_ << std::endl;
  //  key = "data_processing_domain_type";
  //  auto data_processing_domain_type = json.get(key, json_encoding.c_str() ).asInt();
  //  set_data_processing_domain_type(data_processing_domain_type);
  //  debug_print_service("data_processing_domain_type:",data_processing_domain_type);
  //  std::cout << "data_processing_domain_type:" << data_processing_domain_type << std::endl;
  //  Json::Value  data_processing = json.get("data_processing", json_encoding.c_str());
  //  ///\brief set outlier threshold factor. the factor by which 3*MAD would be multiplied to get the final threshold value
  //  key = "outlier_threshold_factor";
  //  auto outlier_threshold_factor = data_processing[key].asFloat();
  //  set_outlier_threshold_factor(outlier_threshold_factor);
  //  debug_print_service("outlier_threshold_factor:",outlier_threshold_factor);
  //  key = "raw_data_file_start_column_input";
  //  auto raw_data_file_start_column_input = data_processing[key].asInt();
  //  set_raw_data_file_start_column_input(raw_data_file_start_column_input);
  //  key = "raw_data_file_end_column_input";
  //  auto data_file_end_column_input = data_processing[key].asInt();
  //  set_data_file_end_column_input(data_file_end_column_input);
  //  debug_print_service("data_file_end_column_input:", raw_data_file_end_column_input_," raw_data_file_start_column_input:",raw_data_file_start_column_input_);
  //  key = "correction_type";
  //  auto correction_type = data_processing[key].asString();
  //  set_correction_type(correction_type);
  //  USH val = get_correction_type();
  //  key = "data_sampling_rate";
  //  auto data_sampling_rate = data_processing[key].asFloat();
  //  set_data_sampling_rate(data_sampling_rate);
  //  key = "time_unit";
  //  auto time_unit = data_processing[key].asString();
  //  ///\brief ignore the case sensitivity from user input
  //  transform(time_unit.begin(), time_unit.end(), time_unit.begin(), ::tolower);
  //  /*
  //   *
  //   */
  //  auto time_unit_enum = TimeUnitDescriptor[time_unit]; /// \breif set the time unit in terms of descriptor
  //  debug_print_service("time_unit_enum:",time_unit_enum);
  //  set_time_unit(time_unit_enum);
  //  /**
  //   * window and overlap keys
  //   */
  //  key = "spatial_correction_zonal_coverage_depth";
  //  auto spatial_correction_zonal_coverage_depth = data_processing[key].asInt(); /// \brief how much of the spatial correction to cover
  //  set_spatial_correction_zonal_coverage_depth(spatial_correction_zonal_coverage_depth);
  //  debug_print_service("get_spatial_correction_zonal_coverage_depth:",get_spatial_correction_zonal_coverage_depth());
  //  /**
  //   * Data blocks and overlaps
  //   * *****************************************************************
  //   */
  //  key = "data_block_window_size_common_dual";
  //  auto data_block_window_size_common_dual = data_processing[key].asInt();
  //  set_data_block_window_common_dual(data_block_window_size_common_dual);
  //  debug_print_service("get_data_block_window_size_common_dual:",get_data_block_window_size_common_dual());
  //  /**
  //   *
  //   */
  //  key = "data_block_window_size_common_dual_overlap_percent";
  //  auto data_block_window_size_common_dual_overlap_percent = data_processing[key].asFloat();
  //  set_data_block_window_size_common_dual_overlap_percent(data_block_window_size_common_dual_overlap_percent);
  //  debug_print_service("get_data_block_window_overlap_percent:",get_data_block_window_size_common_dual_overlap_percent());
  //  /**
  //   *
  //   */
  //  key = "data_block_window_temporal";
  //  auto data_block_window_temporal = data_processing[key].asInt();
  //  set_data_block_window_temporal(data_block_window_temporal);
  //  debug_print_service("get_data_block_window_temporal:",get_data_block_window_temporal());
  //  /**
  //   *
  //   */
  //  key = "data_block_window_temporal_common_overlap_percent";
  //  auto data_block_window_temporal_common_overlap_percent = data_processing[key].asFloat();
  //  set_data_block_window_temporal_common_overlap_percent(data_block_window_temporal_common_overlap_percent);
  //  debug_print_service("get_data_block_window_temporal_common_overlap_percent:",get_data_block_window_temporal_common_overlap_percent());
  //  /**
  //   *
  //   */
  //  key = "data_block_window_spatial";
  //  auto data_block_window_spatial = data_processing[key].asInt();
  //  set_data_block_window_spatial(data_block_window_spatial);
  //  debug_print_service ( "get_data_block_window_spatial():", get_data_block_window_spatial() );
  //  /**
  //   *
  //   */
  //  key = "data_block_window_spatial_common_overlap_percent";
  //  auto data_block_window_spatial_overlap_percent = data_processing[key].asFloat();
  //  set_data_block_window_spatial_common_overlap_percent(data_block_window_spatial_overlap_percent);
  //  debug_print_service("get_data_block_window_spatial_overlap_percent:",get_data_block_window_spatial_common_overlap_percent());
  //  /**
  //   *
  //   *
  //   */
  //  key = "is_temporal_polyfit_applied";
  //  auto is_temporal_polyfit_applied = data_processing[key].asFloat();
  //  set_is_temporal_polyfit_applied(is_temporal_polyfit_applied);
  //  debug_print_service("get_is_temporal_polyfit_applied:",get_is_temporal_polyfit_applied());
  //  key = "temporal_polyfit_order";
  //  auto temporal_polyfit_order = data_processing[key].asFloat();
  //  set_temporal_polyfit_order(temporal_polyfit_order);
  //  debug_print_service( "get_temporal_polyfit_order:", get_temporal_polyfit_order());
  //  key = "is_half_gaussian_smoothing_applied";
  //  auto is_half_gaussian_smoothing_applied = data_processing[key].asFloat();
  //  set_is_half_gaussian_smoothing_applied(is_half_gaussian_smoothing_applied);
  //  debug_print_service("get_is_half_gaussian_smoothing_applied:",get_is_half_gaussian_smoothing_applied());
  exit(0);

}

/**
 * \brief read the config.json file
 */
template<typename T>
Json::Value   QuasiGrid::read_json(const T& filename ){
  Json::Value root;
  Json::Reader reader;
  std::ifstream file(filename.c_str());
  if(!reader.parse(file, root, true)){
    //for some reason if it fails to parse
    //        error_print_service("Failed to parse configuration for file ",filename,reader.getFormattedErrorMessages());
  }
  return root;
}
template<typename T>
void QuasiGrid::set_output_folder_root(const T& output_folder_root) {
  output_folder_root_ = output_folder_root;
}

/// @brief
///
/// @param {void}
decltype(auto) QuasiGrid::get_output_folder_root() const{
  return output_folder_root_;
}


/**
 *
 *
 */
template <typename... Args>
void QuasiGrid::error_print_service(Args&&... args) {
  static std::mutex error_print_mutex;
#ifdef ERROR_OUTPUT
  const std::lock_guard<std::mutex> lock(error_print_mutex);
  ((std::cerr << std::forward<Args>(args) << std::endl), ...);
#endif
}
/**
 *A custom filename creator based on arguments
 */
template <typename... Args>
decltype(auto) QuasiGrid::filename_creator(Args&&... args){
  std::stringstream ss;
  ((ss << std::forward<Args>(args)), ...);
  return ss.str();
}
/**
 * \brief  printing service for  dubuging purpose
 */
  template <typename... Args>
void QuasiGrid::debug_print_service(Args&&... args)
{
  static std::mutex debug_print_mutex;
#ifdef IS_DEBUG
  const std::lock_guard<std::mutex> lock(debug_print_mutex);
  ((std::cout << std::forward<Args>(args) << " "), ...);
#ifdef IS_DEBUG_WITH_ENDL
  std::cout << std::endl;
#endif
#endif
}
/// @brief A log printing service
/// @tparam Args
/// @tparam T
/// @param args
/// @param file_stream
template <typename... Args, typename T>
void QuasiGrid::log_service(Args&&... args, T& file_stream) {
  std::mutex log_mutex;
  const std::lock_guard<std::mutex> lock(log_mutex);
  ((file_stream << std::forward<Args>(args) << " "), ...);
  file_stream << std::endl;
}

#endif
