/**
 *A Junction  class to implement quasi one dimensional models
 */
class Junction{
  public:
    Junction();
    Junction(const double& x);
    Junction( double&& x);
    ~Junction();
    void set_x(const double  x);
    double get_x() const;
    unsigned get_id() const;
    double get_len() const;
    inline void update(const double& dx) { x_ += dx;}
  private:
    void ctor_helper();
    double x_;
    unsigned id_;
    static double len_;
    static unsigned  counter_;
};

unsigned Junction::counter_ = {0};
double Junction::len_ = {0.0};

void Junction::ctor_helper(){
  id_ = counter_++;
  std::cerr << "setting id:" << id_ << std::endl;
}

double Junction::get_len() const {
};
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
