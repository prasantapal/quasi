#include "junction.hpp"

class Intersection{
  public:
    Intersection();
    ~Intersection();
  private:
    void ctor_helper();
    std::vector<Junction> intersection_;
    unsigned id_;
    static unsigned counter_;
    static unsigned num_junctions_per_intersection_;
};
//In case of 3D it will be 3
//          -----
//         |     |
//         |     |
//          -----
//
unsigned Intersection::num_junctions_per_intersection_ = {2};
unsigned Intersection::counter_ = {0};


Intersection::Intersection(){
  std::cerr << __func__ << " ctor" << std::endl;
}

Intersection::~Intersection(){

  std::cerr << __func__ << " dtor" << std::endl;

}

void Intersection::ctor_helper() {
  id_ = counter_++;
}

