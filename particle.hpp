#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include<iostream>
#include<fstream>
#include <thread>
#include <mutex>
#include <string>
class Particle{
  public:
    Particle();
    Particle(const double& x);
    Particle(double&& x);
    void save();
    ~Particle();
    double get_x() const;
    inline  void update(const double& dx) {x_ += dx; }
    inline  void update(const double&& dx) {x_ += dx; }
    void set_filename(std::string& filename);
    void time_evolve(const int& n);
  private:
    void  ctor_helper();
    double x_;
    unsigned id_;
    static unsigned counter_;
    static std::string filename_;
    static std::ofstream file_;
    static std::once_flag seed_initializer_;
};
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
double Particle::get_x() const{
  return x_;
}
void Particle::time_evolve(const int& n){
  std::call_once(seed_initializer_,[](){
  std::srand(std::time(NULL));
  });

  for(auto i=0;i<n;++i){
    double dx = static_cast<double>(std::rand())/RAND_MAX - 0.5;
    update(dx);
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
#endif
