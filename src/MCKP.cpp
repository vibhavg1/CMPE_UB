
#include <MCKP.h>
#include <algorithm>
#include <limits>
#include <iostream>

struct KP_Item {
    long double profit;
    long double cost;
    int i;
    int j;
};

bool Item_sorter_efficiency(KP_Item const &lhs, KP_Item &rhs) {
    long double eff1 = lhs.profit / lhs.cost;
    long double eff2 = rhs.profit / rhs.cost;
    return eff1 > eff2;
}

void sortBinEfficiency(vector<KP_Item> &bin) {
    sort(bin.begin(), bin.end(), &Item_sorter_efficiency);
}

bool Item_sorter(Item const &lhs, Item &rhs) {
    return lhs.cost < rhs.cost;
}

bool dominateIP(Item const &lhs, Item &rhs) {
    return (lhs.cost <= rhs.cost && lhs.profit >= rhs.profit);
}

bool dominateLP(Item const &it1, Item &it2, Item &it3) {
    if (it1.cost < it2.cost && it2.cost < it3.cost &&
        it1.profit < it2.profit && it2.profit < it3.profit) {
        double term1 = (it3.profit - it2.profit) / (it3.cost - it2.cost);
        double term2 = (it2.profit - it1.profit) / (it2.cost - it1.cost);
        return (term1 >= term2);
    } else return false;
}

void sortBin(Bin &bin) {
    sort(bin.begin(), bin.end(), &Item_sorter);
}

void RemoveDominatedItemsLP(Bin &bin_) {
    if (bin_.empty()) return;
    Bin bin = bin_;
    sortBin(bin);
    bin_ = vector<Item>();

    bin_.emplace_back(bin[0]);
    for (int k = 1; k < bin.size(); k++) {
        bool estDomine = false;
        int j = 0;
        while (!estDomine && j < k) {
            estDomine = dominateIP(bin[j], bin[k]);

            if (k < bin.size() - 1) {
                int l = k + 1;
                while (!estDomine && l < bin.size()) {
                    estDomine = dominateLP(bin[j], bin[k], bin[l]);
                    l++;
                }
            }
            j++;
        }
        if (!estDomine) bin_.push_back(bin[k]);
    }
}


void print_mckp(MCKP &mckp) {
    cout << "Num bins = " << mckp.size() << endl;
    for (int i = 0; i < mckp.size(); i++) {
        cout << "Bin = " << i + 1 << ": ";
        for (int j = 0; j < mckp[i].size(); j++) {
            cout << "(" << mckp[i][j].cost << "," << mckp[i][j].profit << ") ";
        }
        cout << endl;
    }
}

long double make_commensurate(vector<vector<long double> > &c) {
    long double min_value = std::numeric_limits<double>::max();
    for (int i = 0; i < c.size(); i++) {
        for (int j = 0; j < c[i].size(); j++) {
            if (min_value > c[i][j]) min_value = c[i][j];
        }
    }
    if (min_value < 0) {
        min_value -= 0.1;
        for (int i = 0; i < c.size(); i++) {
            for (int j = 0; j < c[i].size(); j++) {
                c[i][j] += -min_value;
            }
        }
        return (min_value) * c.size();
    }
    return 0.0;
}


long double MCKP_HELPER::LP_solve(vector<vector<long double> > &weights, vector<vector<long double> > &profits,
                                  long double max_cost,
                                  vector<vector<long double> > &solution) {
    // Begin: Construct the MCKP from profits and weights
    if (weights.size() != profits.size()) {
        cerr << "Mismatch in the number of Bins\n";
        exit(-1);
    }
    long double obj_add = make_commensurate(profits);
    long double constaint_add = make_commensurate(weights);
    max_cost += -constaint_add;
    MCKP mckp = vector<Bin>(weights.size());
    for (int i = 0; i < weights.size(); i++) {
        mckp[i] = vector<Item>(weights[i].size());
        if (weights[i].size() != profits[i].size()) {
            cerr << "Mismatch in the number of items in Bin " << i << "\n";
            exit(-1);
        }
        for (int j = 0; j < weights[i].size(); j++) {
            mckp[i][j].profit = profits[i][j];
            mckp[i][j].cost = weights[i][j];
            mckp[i][j].pos_in_bin = j;
            mckp[i][j].bin_id = i;
        }
    }
    //print_mckp(mckp);
    // End: Cosntruct MCKP

    // Construct Greedy solution to MCKP
    int num_bins = weights.size();
    // Step 1. Remove Dominated items in each bin


    vector<bool> is_multi_item_bin(num_bins, false);
    vector<vector<long double> > lp_solution(num_bins);
    for (int i = 0; i < num_bins; i++)
        lp_solution[i] = vector<long double>(mckp[i].size(), 0.0);
    for (int i = 0; i < num_bins; i++) {
        RemoveDominatedItemsLP(mckp[i]);
        if ((int) mckp[i].size() > 1) {
            is_multi_item_bin[i] = true;
        } else {
            if ((int) mckp[i].size() == 1) lp_solution[i][mckp[i][0].pos_in_bin] = 1.0;
        }
    }

    // Create an instance of knapsack

    long double knapsack_instance_max_cost = max_cost;
    long double lp_profit = 0.0;


    vector<KP_Item> knapsack_instance;

    for (int i = 0; i < num_bins; i++) {
        knapsack_instance_max_cost -= mckp[i][0].cost;
        lp_profit += mckp[i][0].profit;
        for (int j = 1; j < mckp[i].size(); j++) {
            KP_Item item;
            item.profit = mckp[i][j].profit - mckp[i][j - 1].profit;
            item.cost = mckp[i][j].cost - mckp[i][j - 1].cost;
            item.i = i;
            item.j = j;
            knapsack_instance.push_back(item);
        }
    }
    long double residual_capacity = knapsack_instance_max_cost;
    // Get LP bound on KP
    sortBinEfficiency(knapsack_instance);
    vector<long double> kp_x(knapsack_instance.size(), 0.0);
    long double kp_cost = 0.0, kp_profit = 0.0;
    for (int i = 0; i < knapsack_instance.size(); i++) {
        if (knapsack_instance[i].cost + kp_cost <= knapsack_instance_max_cost) {
            kp_profit += knapsack_instance[i].profit;
            kp_cost += knapsack_instance[i].cost;
            kp_x[i] = 1.0;

            lp_profit = lp_profit + knapsack_instance[i].profit;
            residual_capacity -= knapsack_instance[i].cost;
            int a = knapsack_instance[i].i;
            int b = knapsack_instance[i].j;
            lp_solution[mckp[a][b].bin_id][mckp[a][b].pos_in_bin] = 1.0;
            lp_solution[mckp[a][b - 1].bin_id][mckp[a][b - 1].pos_in_bin] = 0.0;

        } else {
            kp_x[i] = (knapsack_instance_max_cost - kp_cost) / knapsack_instance[i].cost;
            kp_profit += knapsack_instance[i].profit * kp_x[i];
            int a = knapsack_instance[i].i;
            int b = knapsack_instance[i].j;
            lp_solution[mckp[a][b].bin_id][mckp[a][b].pos_in_bin] = kp_x[i];
            lp_solution[mckp[a][b - 1].bin_id][mckp[a][b - 1].pos_in_bin] = 1.0 - kp_x[i];
            lp_profit = lp_profit + knapsack_instance[i].profit * kp_x[i];
            break;
        }
    }
    solution = lp_solution;
    //cout<<"Upper bound is "<<lp_profit<<endl;
    return lp_profit + obj_add;
}
