#ifndef JUNCTION_HPP
#define JUNCTION_HPP
///  @brief Solves non-linear equation with Newton method.
///
///  @tparam T   Any float-point type such as float, double or long double
///  @param fun  Non-linear function f(x)
///  @param dfun Derivative of non-linear function df(x) = d/dx f(x)
///  @param x0   Initial guess
///  @param eps  Tolerance for stopping criteria.
///  @return     Equation result as a float point type T.
///
///  @details
///  Solves non-linear equation using Newton method. This function needs two
///  functions, the function to be solved @p fun and its derivate @p dfun
///
///  @note     The function f(x) must be continues and differentiable.
///  @warning  Throws NonCoverge exception when the root is not found.
///
///  @see NewtonSolver
///  @see https://en.wikipedia.org/wiki/Newton%27s_method
///
///  Example:
///  @code
///    // Solve f(x) = x^2 - 25.0 , df(x) = 2x around x0 = 10.0
///    auto fun = [](double x){ return x * x -  25.0 };
///    auto dfun = [](double x){ return 2 * x; }
///
///    double root = GenericNewtonsolver(fun, dfun, 10.0, 0.001);
///    std::cout << "Root = " << root << std::endl;
///  @endcode
///
/**
 * @brief A Junction  class to implement quasi one dimensional models
 *
 */
#include <list>
#include <iostream>
#include <string>
#include <vector>
#include "particle.hpp"
class Junction{
  public:
    Junction();
//    Junction(Junction&) = delete;
    Junction(const double& x);
    Junction( double&& x);
    ~Junction();
    static void set_length_fixed_size(const double&);
    bool operator <(const Junction & junction) const
    {
      return x_ < junction.get_x();
    }
    bool get_is_blocked() const;
    int get_counter() const;
    void set_label(const int label);
    int get_label() const;
    inline bool does_belong_to(const double& x) { return (x>= lower_boundary_ && x <= upper_boundary_); }
    inline Junction* does_belong_to_junction(const double& x) { return (x>= lower_boundary_ && x <= upper_boundary_)?this:nullptr; } ///return the pointer to the junction if it belongs else return nullptr
    void set_x(const double  x);
    void set_len(const double  len);
    double get_x() const;
    unsigned get_id() const;
    double get_len() const;
    void occupy(Particle*,bool is_forward = true);///occupy a junction
    void unoccupy(bool is_forward = true);///occupy a junction
    void print_occupation() const;
    bool is_occupied() const;
    void set_occupation(const short pos);
    void reset_occupation(const short pos);
    void set_is_blocked(const bool is_blocked);
    static int get_junction_count() ;
  private:
    void set_len_dependency();
    void ctor_helper();
    double x_;
    unsigned id_;
    double len_;
    double len_half_;
    double lower_boundary_;
    double upper_boundary_;
    int state_;
    std::vector<bool> occupation_list_status_ = {false,false}; ///Store the list of ids occupying the junction
    std::list<Particle*> occupations_;
    bool is_blocked_= {false};
    bool is_active_ = {true}; //whether or not the junction is is_active_ive
    static unsigned  counter_;
    static int max_occupation_;
    static double len_fixed_size_;///This is the length when the junction is of fixed size
};
double Junction::len_fixed_size_ = {0.0};

///**
// *return whether or not the junction is occupied
// */
//bool Junction::is_occupied() const {
//bool is_occupied = {false};
//for(const auto it:occupation_list_status_) {
//is_occupied |= it;
//}
//return is_occupied;
//}
/**
 *return whether or not the junction is occupied
 */
bool Junction::is_occupied() const {
  return (occupations_.size() > 0?true:false);
}
int Junction::max_occupation_ = {2}; ///The maximumum number of particles that can occupy a junction
//This number is decided by the underlying geometry but in the simplest case it is 2
unsigned Junction::counter_ = {0}; ///Initialize the counter with zero
void Junction::ctor_helper(){
  id_ = counter_++;
  std::cerr << "setting id:" << id_ << std::endl;
  ///Initialize the occupation status list....this is more efficient than vector push back or pop methods
  occupation_list_status_ = std::vector<bool>(max_occupation_,false);
}
unsigned Junction::get_id() const {
  return id_;
}
Junction::Junction(const double& x):x_(x){
  ctor_helper();
  std::cerr << __func__ << " ctor" << std::endl;
}
Junction::Junction( double&& x):x_(std::move(x)){
  ctor_helper();
  std::cerr << __func__ << " ctor" << std::endl;
}
Junction::Junction():x_(0){
  ctor_helper();
  std::cerr << __func__ << " ctor" << std::endl;
}
Junction::~Junction(){
  std::cerr << __func__ << " dtor" << std::endl;
}
void Junction::set_x(const double  x) { x_ = x;};
double Junction::get_x() const { return x_; };
/**
 * \brief set the length of the junction
 */
void Junction::set_len(const double  len) { len_ = len; set_len_dependency(); };
double Junction::get_len() const { return len_; };
/**
 *
 */
void Junction::set_len_dependency(){
  lower_boundary_ = x_ - len_;
  upper_boundary_ = x_ + len_;
}
/**
 * set the occupation status of the junction
 */
void Junction::set_occupation(const short pos){
  if(pos >=0 && pos<max_occupation_){
    occupation_list_status_[pos] = true;
  }else {
    std::cerr << "sorry position can not be more than " << max_occupation_ -1 << " ...exit "<< std::endl;
    exit(0);
  }
}
void Junction::reset_occupation(const short pos){
  if(pos >=0 && pos<max_occupation_){
    occupation_list_status_[pos] = false;
  }else {
    std::cerr << "sorry position can not be more than " << max_occupation_ -1 << " ...exit "<< std::endl;
    exit(0);
  }
}
int Junction::get_counter() const{
  return counter_;
}
void Junction::set_label(const int label){
  id_ = label;
}
int Junction::get_label() const{
  return id_;
}
void Junction::set_is_blocked(const bool is_blocked){
  is_blocked_ = is_blocked;
}
bool Junction::get_is_blocked() const {
  return is_blocked_;
}
/**
* Occupy the junction with a particle.
*
*/
void Junction::occupy(Particle* p,bool is_forward){///If forward then queue the list in the back or in the front
  if(is_forward)
    this->occupations_.emplace_back(p);
  else
    this->occupations_.emplace_front(p);
}
/**
 *
 *
 */
void Junction::unoccupy(bool is_forward ){
  if(occupations_.size() > 0){
    if(is_forward)
      this->occupations_.pop_front();
    else
      this->occupations_.pop_back();
  }else {
    std::cerr << "sorry can not unoccupy the junction as there is no particle" << std::endl;
  }
}
/**
 *
 *
 */
void Junction::print_occupation() const{
  for(const auto& it:occupations_){
    it->print();
  }
}
int Junction::get_junction_count() {
  return counter_;
}
void Junction::set_length_fixed_size(const double& l){
}


#ifdef MAKE_TEST
int main(int argc, char** argv){
  return 0;
}
#endif

#endif
