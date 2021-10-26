#ifndef CMPE_BB_CG_H
#define CMPE_BB_CG_H

#include <MN.h>
#include <cfloat>

typedef enum {
    MINI_BUCKET,
    MINI_BUCKET_MM,
    JOIN_GRAPH
} PROP_TYPE;

struct CG {
    MN *mn1;
    MN *mn2;
    // Each cluster consists of a potential
    vector<LogPotential> clusters;
    // Given an ordering of variables the following maps clusters to order
    // Example if a variable j
    vector<vector<int> > order2clusters;
    vector<int> to_edge;
    vector<int> order;
    vector<int> mn1_pot2clusters;
    vector<int> mn2_pot2clusters;

    CG(MN *mn1_, MN *mn2_, vector<int> &order_, int i_bound);

    void initialize(long double lambda);

    static pair<long double, long double>
    ub_search(MN &mn1_, MN &mn2_, vector<int> &order_, int i_bound, long double q, PROP_TYPE propType = MINI_BUCKET);

    long double miniBucketBound();

    long double miniBucketBoundMM();

    long double joinGraphBound(int iter_max = 10);

    void propagate();

    long double getSolution() {
        for (int i = order.size() - 1; i > -1; i--) {
            Variable *var = mn1->variables[order[i]];
            vector<long double> f_max(var->domain_size, 0.0);
            for (int j = 0; j < order2clusters[i].size(); j++) {
                int c_id = order2clusters[i][j];
                for (int k = 0; k < var->domain_size; k++) {
                    var->value = k;
                    f_max[k] += clusters[c_id].getValue();
                }
            }
            int maxElementIndex = max_element(f_max.begin(), f_max.end()) - f_max.begin();
            var->value = maxElementIndex;
        }
        return mn1->getValue();
    }
};

#endif //CMPE_BB_CG_H
