#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include<iostream>
#include<fstream>
#include <thread>
#include <cmath>
#include <string>
#include <mutex>
using USH = unsigned short int;
class Particle{
  public:
    Particle();
    Particle(const double& x);
    Particle(double&& x);
    inline double mod_len(const double x) {
      if(x>=0)
        return std::fmod(x,system_len_);
      return (system_len_ - std::fmod(std::fabs(x),system_len_));

    }

    void save();
    ~Particle();
    decltype(auto) get_x() const;
    double get_system_len() const;
    double get_len() const;
    void set_x(const double& x);
    void set_x(const double&& x);
    void set_X(const double& x);
    void set_X(const double&& x);
    double get_X() const;
    void set_len(const double x) ;
    void set_system_len(const double x) ;
    unsigned int get_label() const;
/****** THIS IS  A VERY IMPORTANT FUNCTION ****************/
    inline bool does_belong_to(const double& x) {//Check whether or not a point below to the particle geometrical boundaries
      auto lower_boundary_ =  mod_len(x_ - len_half_);
      auto upper_boundary_ =  mod_len(x_ + len_half_);
      std::cout << x << " " << len_ << " " << lower_boundary_ << " " << upper_boundary_ << std::endl;
      if(lower_boundary_ > upper_boundary_) { ///This is a special case

        return (mod_len(x)>= lower_boundary_ || mod_len(x)<= upper_boundary_);
      }

      return (mod_len(x)>= lower_boundary_ && mod_len(x) <= upper_boundary_);
    }
    inline  void move(const double& dx) {x_ += dx; X_ += dx;}
    inline  void move(const double&& dx) {x_ += dx; X_ += dx;}
    void set_filename(std::string& filename);
    void time_evolve(const int& n);
    void print() const;

    void set_state_label(const USH&);
    USH get_state_label() const;

    void set_state_label_id(const USH&);
    USH get_state_label_id()const;



  private:
    void  ctor_helper();
    double x_ = {0};
    double x_last_ = {0};
    double X_ = {0};
    double X_last_ = {0};
    unsigned id_ = {0};
    static unsigned counter_;;
    static std::string filename_;
    static std::ofstream file_;
    static std::once_flag seed_initializer_;
    static double len_;
    static double len_half_;
    static double system_len_;
    /*** STATE LABEL ****/
    USH state_label_;
    USH state_label_id_;
};
double Particle::len_ = {0};
double Particle::len_half_ = {0};
double Particle::system_len_ = {0};
unsigned Particle::counter_ = {0};
std::once_flag Particle::seed_initializer_;
std::ofstream Particle::file_;
std::string Particle::filename_;
void Particle::ctor_helper(){
  id_ = counter_++;
  std::cerr << "creating particle with id:" << id_ << std::endl;
}
Particle::Particle():x_(0){
  ctor_helper();
  std::cerr << __func__ << " ctor" << std::endl;
}
Particle::Particle(const double& x):x_(x){
  ctor_helper();
  std::cerr << __func__ << " ctor" << std::endl;
}
Particle::~Particle(){
  std::cerr << __func__ << " dtor" << std::endl;
}
decltype(auto) Particle::get_x() const{
  return x_;
}
void Particle::set_x(const double& x){
  x_ = x;
}
void Particle::set_x(const double&& x){
  x_ = x;
}
void Particle::set_X(const double& x){
  X_ = x;
}
void Particle::set_X(const double&& x){
  X_ = x;
}
double Particle::get_X() const{///Uncorrected coordinate
  this->X_;
}
void Particle::time_evolve(const int& n){
  std::call_once(seed_initializer_,[](){
      std::srand(std::time(NULL));
      });
  for(auto i=0;i<n;++i){
    double dx = static_cast<double>(std::rand())/RAND_MAX - 0.5;
    move(dx);
    save();
  }
}
void Particle::set_filename(std::string& filename){
  filename_ = filename;
}
void Particle::save(){
  if(!file_.is_open()){
    file_.open(filename_.c_str());
    if(!file_.good()){
      std::cerr << "could not open file " << filename_ << std::endl;
    }
  }
  file_ << x_ << std::endl;
}
void Particle::set_len(const double x) {
  len_ = x;
  len_half_ = len_/2.0;
}
void Particle::set_system_len(const double x) {
  system_len_ = x;
}
double Particle::get_system_len() const{
  return system_len_;
}
double Particle::get_len() const{
  return len_;
}

void Particle::print() const{
  std::cout << "id:" << id_ << " x:" << x_ << " " << " len:" << len_<< std::endl;
}

unsigned int Particle::get_label() const{
  return this->id_;
}

void Particle::set_state_label(const USH& label)  {
  this->state_label_ = label;
}
USH Particle::get_state_label() const {
  return this->state_label_;
}

void Particle::set_state_label_id(const USH& id)  {
  this->state_label_id_ = id;
}

USH Particle::get_state_label_id() const {
  return this->state_label_id_;
}



#endif
