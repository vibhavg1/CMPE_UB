
#ifndef CMPE_BB_MCKP_H
#define CMPE_BB_MCKP_H

#include <vector>

using namespace std;

struct Item {
    long double profit;
    long double cost;
    int pos_in_bin;
    int bin_id;
};

typedef vector<Item> Bin;

typedef vector<Bin> MCKP;

struct MCKP_HELPER {
    static long double
    LP_solve(vector<vector<long double> > &weights, vector<vector<long double> > &profits, long double max_cost,
             vector<vector<long double> > &solution);
};

#endif //CMPE_BB_MCKP_H
