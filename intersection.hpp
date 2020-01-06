#include "junction.hpp"
#include <map>
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
    auto& get_intersection_ref() const;
    void set_junction_conjugate();
    static int get_junction_count() ;
    bool is_blocked(const int) const;
    void print_conjugate() const;

const std::multimap<int,int>& get_junction_conjugate() const;
  private:
    void ctor_helper();
    std::vector<Junction> intersection_;
    std::multimap<int,int> junction_conjugate_;
    unsigned id_;
    static unsigned counter_;
    static unsigned num_junctions_per_intersection_;
};

const std::multimap<int,int>& Intersection::get_junction_conjugate() const{
return this->junction_conjugate_;
}
/**
 * set conjugate junctions for each junction index
 * By the nature of it it is a multimap from a single juncion index to several other (in theory for higher dimensions)
 */
void Intersection::set_junction_conjugate() {
  for(auto& junction:intersection_){
    auto label = junction.get_label();
    std::vector<Junction> other_junctions;
    std::vector<Junction>::iterator iter = intersection_.begin();
    while ((iter = std::find_if(iter, intersection_.end(), [=](const auto& val){ return label != val.get_label() ;})) != intersection_.end()) {
      junction_conjugate_.emplace(label, iter->get_label());
      iter++;
    }
  }
}
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
auto& Intersection::get_intersection_ref() const {
  return intersection_;
}
bool Intersection::is_blocked(const int index) const{
//this->junction_.at(index);
}

void Intersection::print_conjugate() const{
  std::cout << "printing junction conjugate for intersection " << id_ << std::endl;
  for(auto& it:junction_conjugate_){

    std::cout << it.first << " " << it.second << std::endl;
  }
   std::cout << std::endl;

}
#ifdef MAKE_TEST
int main(int argc, char** argv){
  return 0;
}
#endif
