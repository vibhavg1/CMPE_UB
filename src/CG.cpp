
#include <CG.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>

CG::CG(MN *mn1_, MN *mn2_, vector<int> &order_, int i_bound) : mn1(mn1_), mn2(mn2_), order(order_) {

    vector<vector<LogPotential *> > buckets(order.size());
    vector<int> var_pos(order.size());
    for (int i = 0; i < order.size(); i++) {
        var_pos[order[i]] = i;
    }
    vector<vector<int> > origin(order.size());
    mn1_pot2clusters = vector<int>(mn1->lpotentials.size());
    // origin[i][j]=-1 if originated from mn1
    for (auto &lpotential : mn1->lpotentials) {
        int bucket_func = order.size();
        for (auto var : lpotential->variables) {
            if (bucket_func > var_pos[var->id]) {
                bucket_func = var_pos[var->id];
            }
        }
        buckets[bucket_func].push_back(lpotential);
        origin[bucket_func].push_back(-1);
    }
    mn2_pot2clusters = vector<int>(mn2->lpotentials.size());
    // origin[i][j]=-2 if originated from mn2
    for (auto &lpotential : mn2->lpotentials) {
        int bucket_func = order.size();
        for (auto var : lpotential->variables) {
            if (bucket_func > var_pos[var->id]) {
                bucket_func = var_pos[var->id];
            }
        }
        buckets[bucket_func].push_back(lpotential);
        origin[bucket_func].push_back(-2);
    }
    order2clusters = vector<vector<int> >(order.size());
    vector<LogPotential *> to_delete;

    // Process
    to_edge = vector<int>();
    for (int i = 0; i < order.size(); i++) {
        //cout<<"{rocessing "<<i<<endl;
        //vector<set<int> >mini_bucket_func_ids;
        //vector<vector<Variable*> > mini_bucket_vars;
        for (int j = 0; j < buckets[i].size(); j++) {
            // check if the function can be added to a minibucket
            //buckets[i][j]->print();
            int min_cluster_id = -1;
            size_t min_count = std::numeric_limits<size_t>::max();
            for (int k = 0; k < order2clusters[i].size(); k++) {
                int cluster_id = order2clusters[i][k];
                vector<Variable *> tmp;
                Variable::set_union_variables(clusters[cluster_id].variables, buckets[i][j]->variables, tmp);
                if ((int) tmp.size() <= i_bound) {
                    if ((tmp.size() - clusters[cluster_id].variables.size()) < min_count) {
                        min_cluster_id = cluster_id;
                        min_count = tmp.size() - clusters[cluster_id].variables.size();
                    }
                }
            }
            if (min_cluster_id == -1) {
                to_edge.push_back(-1);
                order2clusters[i].push_back(clusters.size());
                if (origin[i][j] > -1) {
                    to_edge[origin[i][j]] = clusters.size();
                    //edges.emplace_back(origin[i][j], clusters.size());
                } else if (origin[i][j] == -1) {
                    mn1_pot2clusters[buckets[i][j]->id] = clusters.size();
                } else if (origin[i][j] == -2) {
                    mn2_pot2clusters[buckets[i][j]->id] = clusters.size();
                }
                clusters.emplace_back();
                clusters[clusters.size() - 1].add(buckets[i][j]);
            } else {
                clusters[min_cluster_id].add(buckets[i][j]);
                if (origin[i][j] > -1) {
                    to_edge[origin[i][j]] = min_cluster_id;
                    //edges.emplace_back(origin[i][j], min_cluster_id);
                } else if (origin[i][j] == -1) {
                    mn1_pot2clusters[buckets[i][j]->id] = min_cluster_id;
                } else if (origin[i][j] == -2) {
                    mn2_pot2clusters[buckets[i][j]->id] = min_cluster_id;
                }
            }
        }
        for (int j = 0; j < order2clusters[i].size(); j++) {
            int cluster_id = order2clusters[i][j];
            vector<Variable *> tmp;
            tmp.push_back(mn1->variables[order[i]]);
            Variable::set_difference_variables(clusters[cluster_id].variables, tmp, tmp);
            auto *out = new LogPotential(tmp);
            to_delete.push_back(out);
            if (!out->variables.empty()) {
                int bucket_func = order.size();
                for (auto var : out->variables) {
                    if (bucket_func > var_pos[var->id]) {
                        bucket_func = var_pos[var->id];
                    }
                }
                buckets[bucket_func].push_back(out);
                origin[bucket_func].push_back(cluster_id);
            }
        }
    }
    //cout<<"Cleaning up\n";
    // Clean up
    for (auto &i : to_delete) {
        delete i;
    }
}

void CG::initialize(long double lambda) {
    for (auto &cluster : clusters) {
        size_t numvals = cluster.table.size();
        cluster.table = vector<long double>(numvals, 0.0);
    }
    vector<int> var_pos(order.size());
    for (int i = 0; i < order.size(); i++) {
        var_pos[order[i]] = i;
    }
    for (auto &lpotential : mn1->lpotentials) {
        int c_id = mn1_pot2clusters[lpotential->id];
        clusters[c_id].add(lpotential);
    }
    for (auto &lpotential : mn2->lpotentials) {
        int c_id = mn2_pot2clusters[lpotential->id];
        LogPotential pot(*lpotential);
        for (long double &tab_entry : pot.table)
            tab_entry *= -lambda;
        clusters[c_id].add(&pot);
    }
}

long double CG::miniBucketBound() {
    long double ub = 0.0;
    for (int i = 0; i < order.size(); i++) {
        for (int j = 0; j < order2clusters[i].size(); j++) {
            int c = order2clusters[i][j];
            LogPotential pot;
            clusters[c].maxOutVariable(pot, mn1->variables[order[i]]);
            //if (i>order.size()-3){
            //  pot.print();
            //}
            if (to_edge[c] > -1) {
                if (pot.variables.empty()) {
                    cerr << "Something wrong minibucket-bound 1\n";
                }
                clusters[to_edge[c]].add(&pot);
            } else {
                if (!pot.variables.empty()) {
                    cerr << "Something wrong minibucket-bound 2\n";
                }
//cout<<pot.table[0]<<" ";
                ub += pot.table[0];
//
            }
        }
    }
//cout<<endl;
    return ub;
}

long double CG::miniBucketBoundMM() {
    long double ub = 0.0;
    for (int i = 0; i < order.size(); i++) {
        // Find the common variables in all minibuckets
        vector<Variable *> common_variables;
        //cout<<"Mini buckets before: "<<mini_buckets.size()<<"\n";
        for (int j = 0; j < order2clusters[i].size(); j++) {
            int c_id = order2clusters[i][j];
            if (j == 0) common_variables = clusters[c_id].variables;
            else
                Variable::set_intersection_variables(common_variables, clusters[c_id].variables, common_variables);
        }

        // Find the max-out function for each minibucket
        LogPotential avg;
        vector<LogPotential> max_pots(order2clusters[i].size());
        //cout<<"Printing fmax\n";
        for (int j = 0; j < order2clusters[i].size(); j++) {
            int c_id = order2clusters[i][j];
            clusters[c_id].maxOutVariables(max_pots[j], common_variables);
            //mini_buckets[j].f_max.print();
            if (avg.variables.empty()) avg = max_pots[j];
            else avg.add(&max_pots[j]);
        }

        for (int j = 0; j < avg.table.size(); j++) {
            avg.table[j] /= (long double) order2clusters[i].size();
        }
        for (int j = 0; j < order2clusters[i].size(); j++) {
            int c = order2clusters[i][j];
            clusters[c].subtract(&max_pots[j]);
            clusters[c].add(&avg);
            LogPotential pot;
            clusters[c].maxOutVariable(pot, mn1->variables[order[i]]);
            if (to_edge[c] > -1) {
                if (pot.variables.empty()) {
                    cerr << "Something wrong minibucket-bound 1\n";
                }
                clusters[to_edge[c]].add(&pot);
            } else {
                if (!pot.variables.empty()) {
                    cerr << "Something wrong minibucket-bound 2\n";
                }
//cout<<pot.table[0]<<" ";
                ub += pot.table[0];
//
            }
        }
    }
    return ub;
}

void CG::propagate() {
    for (int iter = 0; iter < 10; iter++) {
        long double error = 0.0;
        for (int i = 0; i < order2clusters.size(); i++) {
            for (int j = 0; j < order2clusters[i].size(); j++) {
                LogPotential f1, f2;
                int c1 = order2clusters[i][j];
                int c2 = to_edge[c1];
                if (c2 > -1) {
                    vector<Variable *> common_variables;
                    Variable::set_intersection_variables(clusters[c1].variables, clusters[c2].variables,
                                                         common_variables);
                    clusters[c1].maxOutVariables(f1, common_variables);
                    clusters[c2].maxOutVariables(f2, common_variables);
                    f1.subtract(&f2);
                    for (int k = 0; k < f1.table.size(); k++) {
                        error += fabs(f1.table[k]);
                        f1.table[k] *= 0.5;
                    }
                    clusters[c1].subtract(&f1);
                    clusters[c2].add(&f1);
                }
            }
            for (int j = 1; j < order2clusters[i].size(); j++) {
                LogPotential f1, f2;
                int c1 = order2clusters[i][j];
                int c2 = order2clusters[i][j - 1];
                vector<Variable *> common_variables;
                Variable::set_intersection_variables(clusters[c1].variables, clusters[c2].variables, common_variables);
                clusters[c1].maxOutVariables(f1, common_variables);
                clusters[c2].maxOutVariables(f2, common_variables);
                f1.subtract(&f2);
                for (int k = 0; k < f1.table.size(); k++) {
                    error += fabs(f1.table[k]);
                    f1.table[k] *= 0.5;
                }
                clusters[c1].subtract(&f1);
                clusters[c2].add(&f1);
            }
        }
        //cout<<"Error = "<<error<<endl;
        if (error < 1e-3) break;
    }
}

long double CG::joinGraphBound(int iter_max) {
    for (int iter = 0; iter < iter_max; iter++) {
        long double error = 0.0;
        for (int i = 0; i < order2clusters.size(); i++) {
            for (int j = 0; j < order2clusters[i].size(); j++) {
                LogPotential f1, f2;
                int c1 = order2clusters[i][j];
                int c2 = to_edge[c1];
                if (c2 > -1) {
                    vector<Variable *> common_variables;
                    Variable::set_intersection_variables(clusters[c1].variables, clusters[c2].variables,
                                                         common_variables);
                    clusters[c1].maxOutVariables(f1, common_variables);
                    clusters[c2].maxOutVariables(f2, common_variables);
                    f1.subtract(&f2);
                    for (int k = 0; k < f1.table.size(); k++) {
                        error += fabs(f1.table[k]);
                        f1.table[k] *= 0.5;
                    }
                    clusters[c1].subtract(&f1);
                    clusters[c2].add(&f1);
                }
            }
            for (int j = 1; j < order2clusters[i].size(); j++) {
                LogPotential f1, f2;
                int c1 = order2clusters[i][j];
                int c2 = order2clusters[i][j - 1];
                vector<Variable *> common_variables;
                Variable::set_intersection_variables(clusters[c1].variables, clusters[c2].variables, common_variables);
                clusters[c1].maxOutVariables(f1, common_variables);
                clusters[c2].maxOutVariables(f2, common_variables);
                f1.subtract(&f2);
                for (int k = 0; k < f1.table.size(); k++) {
                    error += fabs(f1.table[k]);
                    f1.table[k] *= 0.5;
                }
                clusters[c1].subtract(&f1);
                clusters[c2].add(&f1);
            }
        }
        //cout<<"Error = "<<error<<endl;
        //if (error <1e-3) break;
    }
    return miniBucketBound();
}

pair<long double, long double>
CG::ub_search(MN &mn1_, MN &mn2_, vector<int> &order_, int i_bound, long double q, PROP_TYPE propType) {
    //long double MN::minibucket_ub_search(MN &mn1, MN &mn2, vector<int> &order, int i_bound, long double q,bool mm) {
    CG cg(&mn1_, &mn2_, order_, i_bound);
    long double lambda = 1.0;
    //long double ucmpe = MN::minibucket_ub(mn1, mn2, order,0.0, i_bound,q);
    long double best_l = 1.0;
    long double eta = 0.1;
    cg.initialize(0.0);
    long double ucmpe = cg.miniBucketBound();
    int j = 0;
    //cout<<"UMPE via CG = "<<ucmpe<<endl;
    ucmpe = DBL_MAX;
    if (propType == MINI_BUCKET) {
        cout << "Running the Mini Bucket Algorithm\n";
    } else if (propType == MINI_BUCKET_MM) {
        cout << "Running the Mini Bucket MM Algorithm\n";
    } else if (propType == JOIN_GRAPH) {
        cout << "Running the Join Graph Algorithm\n";
    }
    std::time_t start_time = std::time(nullptr);
    int jg_iter = 40;
    for (int iter = 0; iter < 100; iter++) {
        long double v1 = DBL_MAX;
        long double v2 = DBL_MAX;
        if (lambda + eta >= 0.0) {
            cg.initialize(lambda + eta);
            if (propType == MINI_BUCKET) {
                v1 = cg.miniBucketBound() + (lambda + eta) * q;
            } else if (propType == MINI_BUCKET_MM) {
                v1 = cg.miniBucketBoundMM() + (lambda + eta) * q;
            } else if (propType == JOIN_GRAPH) {
                v1 = cg.joinGraphBound(jg_iter) + (lambda + eta) * q;
            } else {
                cerr << "Invalid propagation type specified\n";
                exit(-1);
            }
        }
        if (lambda - eta >= 0.0) {
            cg.initialize(lambda - eta);
            if (propType == MINI_BUCKET) {
                v2 = cg.miniBucketBound() + (lambda - eta) * q;
            } else if (propType == MINI_BUCKET_MM) {
                v2 = cg.miniBucketBoundMM() + (lambda - eta) * q;
            } else if (propType == JOIN_GRAPH) {
                v2 = cg.joinGraphBound(jg_iter) + (lambda - eta) * q;
            } else {
                cerr << "Invalid propagation type specified\n";
                exit(-1);
            }

        }
        long double profit = 0.0;
        if (v1 < v2) {
            lambda = lambda + eta;
            profit = v1;
        } else {
            lambda = lambda - eta;
            profit = v2;
        }
        if (profit < ucmpe) {
            ucmpe = profit;
            best_l = lambda;
            j = 0;
        } else {
            j = j + 1;
            lambda = best_l;
            eta /= 2.0;
            if (eta < 1e-7) {
                break;
            }
        }
        std::time_t curr_time = std::time(nullptr);
        cout << std::setprecision(20) << iter + 1 << "," << curr_time - start_time << "," << ucmpe << endl;
        //cout<<"--UCMPE = "<<ucmpe<<" lambda = "<<lambda<<" Eta = "<<eta<<" v1, v2 "<<v1<<" "<<v2<<endl;
    }
    //cout<<best_l<<endl;
    //cout<<"Best lambda = "<<best_l<<endl;
    //return ucmpe;
    cg.getSolution();
    return pair<long double, long double>(ucmpe, best_l);
}