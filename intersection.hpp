#include "junction.hpp"
/**
 *Intersection consists of two junctions
 */
class Intersection{
  public:
    Intersection();
    ~Intersection();
    int get_counter() const;
    void initialize();
    int get_label() const;
      void print_intersection_labels() const ;
    void print_intersection_coordinates()  const;
    void allocate_labels(const int pos,const int label) ;
    auto get_intersection_ptr() const;
    static int get_junction_count() ;
  private:
    void ctor_helper();

    std::vector<Junction> intersection_;
    unsigned id_;
    static unsigned counter_;
    static unsigned num_junctions_per_intersection_;
};

//std::vector<Junction> Intersection::intersection_;
//In case of 3D it will be 3
//          -----
//         |     |
//         |     |
//          -----
//
unsigned Intersection::num_junctions_per_intersection_ = {2};
unsigned Intersection::counter_ = {0};
void Intersection::allocate_labels(const int pos,const int label) {
  if(pos >=0 && pos < num_junctions_per_intersection_){
    std::cerr << "intersection " << this->get_label() << " settting label " << label  << " to position " << pos << std::endl;
    intersection_[pos].set_label(label);
  }else {
    std::cerr << "position " << pos  << " is not acceptable" << std::endl;
    std::cerr << "allowed range:>=0 && " <<  " <" << num_junctions_per_intersection_ << std::endl;
  }
}
Intersection::Intersection(){
  std::cerr << __func__ << " ctor" << std::endl;
  ctor_helper();
  initialize();
}
Intersection::~Intersection(){
  std::cerr << __func__ << " dtor" << std::endl;
}
void Intersection::ctor_helper() {
  id_ = counter_++;
}
int Intersection::get_counter() const{
  return counter_;
}
void Intersection::initialize() {
  intersection_.resize(num_junctions_per_intersection_);
}
int Intersection::get_label() const {
  return id_;
}

void Intersection::print_intersection_coordinates()   const{
  for(auto& it:intersection_) {
    std::cout << it.get_x() << " ";
  }
  std::cout << std::endl;
}


void Intersection::print_intersection_labels()   const{
  for(auto& it:intersection_) {
    std::cout << it.get_label() << " ";
  }
  std::cout << std::endl;
}


auto Intersection::get_intersection_ptr() const {
  return &intersection_;
}
int Intersection::get_junction_count(){
  return num_junctions_per_intersection_;
}

#ifdef MAKE_TEST
int main(int argc, char** argv){
  return 0;
}
#endif
