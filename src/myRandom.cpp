
#include "../include/myRandom.h"

unsigned long int myRandom::seed = std::chrono::system_clock::now().time_since_epoch().count();;
mt19937 myRandom::m_g(myRandom::seed);
uniform_real_distribution<double> myRandom::double_dist;
uniform_int_distribution<int> myRandom::int_dist;

