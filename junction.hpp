#ifndef JUNCTION_HPP
#define JUNCTION_HPP
/**
 *A Junction  class to implement quasi one dimensional models
 */
#include <list>
#include <iostream>
#include <string>
#include <vector>
class Junction{
  public:
    Junction();
    Junction(const double& x);
    Junction( double&& x);
    ~Junction();
    int get_counter() const;
    void set_label(const int label);

int get_label() const;
    inline bool does_belong_to(const double& x) { return (x>= lower_boundary_ && x <= upper_boundary_); }
    void set_x(const double  x);
    void set_len(const double  len);

    double get_x() const;
    unsigned get_id() const;
    double get_len() const;
    bool is_occupied() const;
    void set_occupation(const short pos);
    void reset_occupation(const short pos);

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
    std::vector<bool> occupation_list_status_; ///Store the list of ids occupying the junction
    static unsigned  counter_;
    static int max_occupation_;
};

/**
 *return whether or not the junction is occupied
 */
bool Junction::is_occupied() const {
bool is_occupied = {false};
for(const auto it:occupation_list_status_) {
is_occupied |= it;
}
return is_occupied;
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

void Junction::set_len(const double  len) { len_ = len; set_len_dependency(); };
double Junction::get_len() const { return len_; };

void Junction::set_len_dependency(){
  len_half_ = len_/2.0;
  lower_boundary_ = x_ - len_half_;
  upper_boundary_ = x_ + len_half_;

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
#endif
