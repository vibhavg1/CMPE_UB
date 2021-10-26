
#include <MN.h>
#include <myRandom.h>
#include <list>
#include <ortools/linear_solver/linear_solver.h>
#include <MCKP.h>

using namespace operations_research;


struct Lambda {
    int vid;
    LogPotential *f1;
    int v1;
    LogPotential *f2;
    int v2;
    vector<long double> values;

    long double getValue(LogPotential *f, int index) {
        //cout<<f->id<<" "<<f1->id<<" "<<f2->id<<" index = "<<index<<" values = "<<values.size()<<endl;
        if (f->id == f1->id) return values[index];
        if (f->id == f2->id) return -values[index];
        else {
            cout << "Something wrong\n";
            return -1.0;
        }
    }

    void updateValue(long double epsilon, vector<vector<int> > &solution) {
        values[solution[f1->id][v1]] -= epsilon;
        values[solution[f2->id][v2]] += epsilon;
    }

    void print() {
        cout << vid << " " << f1->id << " " << f2->id << " " << v1 << " " << v2 << " " << values.size() << endl;
    }
};


long double solve_using_LP(vector<LogPotential *> &f, vector<LogPotential *> &g,
                           vector<vector<Lambda *> > &f2lambdas, long double logp,
                           vector<vector<int> > &solution) {


    vector<vector<long double> > weights(g.size());
    vector<vector<long double> > profits(g.size());
    vector<vector<long double> > x;

    for (int i = 0; i < weights.size(); i++) {
        int domain_size = Variable::getDomainSize(g[i]->variables);
        weights[i] = vector<long double>(domain_size);
        for (int j = 0; j < domain_size; j++) {
            Variable::setAddress(g[i]->variables, j);
            int entry = Variable::getAddress(g[i]->variables);
            weights[i][j] = g[i]->table[entry];
        }
    }


    // Create the objective function.

    for (int i = 0; i < f.size(); i++) {
        int domain_size = Variable::getDomainSize(f[i]->variables);
        profits[i] = vector<long double>(domain_size);
        for (int j = 0; j < domain_size; j++) {
            Variable::setAddress(f[i]->variables, j);
            int entry = Variable::getAddress(f[i]->variables);
            long double val = f[i]->table[entry];
            for (int k = 0; k < f2lambdas[i].size(); k++) {
                int curr_value = f2lambdas[i][k]->f1->variables[f2lambdas[i][k]->v1]->value;
                val += f2lambdas[i][k]->getValue(f[i], curr_value);

            }
            profits[i][j] = val;
        }
    }

    long double obj_value = MCKP_HELPER::LP_solve(weights, profits, logp, x);
    solution = vector<vector<int> >(f.size());
    for (int i = 0; i < f.size(); i++) {
        solution[i] = vector<int>(f[i]->variables.size());
        int domain_size = Variable::getDomainSize(f[i]->variables);
        for (int j = 0; j < domain_size; j++) {
            if (x[i][j] > 0.5) {
                Variable::setAddress(f[i]->variables, j);
                for (int k = 0; k < f[i]->variables.size(); k++)
                    solution[i][k] = f[i]->variables[k]->value;
                break;
            }
        }
    }

    return obj_value;

}

long double solve_optimally(vector<LogPotential *> &f, vector<LogPotential *> &g,
                            vector<vector<Lambda *> > &f2lambdas, long double logp,
                            vector<vector<int> > &solution) {
    //MPSolver solver("simple_mip_program", MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);

    MPSolver solver("simple_mip_program", MPSolver::CLP_LINEAR_PROGRAMMING);


    const double infinity = solver.infinity();
    // x[j] is an array of non-negative, integer variables.
    vector<vector<const MPVariable *> > x(g.size());
    for (int i = 0; i < g.size(); ++i) {
        x[i] = vector<const MPVariable *>(g[i]->table.size());
        for (int j = 0; j < x[i].size(); j++) {
            x[i][j] = solver.MakeBoolVar("");
        }
    }

    // Write the constraint that \sum_{i,j} g[i][j]*x[i][j] should be <=logp
    MPConstraint *constraint1 = solver.MakeRowConstraint(-infinity, logp, "");
    for (int i = 0; i < g.size(); i++) {
        int domain_size = Variable::getDomainSize(g[i]->variables);
        for (int j = 0; j < domain_size; j++) {
            Variable::setAddress(g[i]->variables, j);
            int entry = Variable::getAddress(g[i]->variables);
            constraint1->SetCoefficient(x[i][j], g[i]->table[entry]);
        }
    }
    // Hard constraint to make sure that exactly one value is chosen from each potential
    // Write the constraint that \sum_j x[i][j]=1 for each i
    for (int i = 0; i < g.size(); i++) {
        MPConstraint *constraint = solver.MakeRowConstraint(1.0, 1.0, "");
        for (int j = 0; j < g[i]->table.size(); j++) {
            constraint->SetCoefficient(x[i][j], 1.0);
        }
    }

    // Create the objective function.
    MPObjective *const objective = solver.MutableObjective();
    for (int i = 0; i < f.size(); i++) {
        int domain_size = Variable::getDomainSize(f[i]->variables);
        for (int j = 0; j < domain_size; j++) {
            Variable::setAddress(f[i]->variables, j);
            int entry = Variable::getAddress(f[i]->variables);
            long double val = f[i]->table[entry];
            for (int k = 0; k < f2lambdas[i].size(); k++) {
                int curr_value = f2lambdas[i][k]->f1->variables[f2lambdas[i][k]->v1]->value;
                val += f2lambdas[i][k]->getValue(f[i], curr_value);

            }
            objective->SetCoefficient(x[i][j], val);
        }
    }
    objective->SetMaximization();
    const MPSolver::ResultStatus result_status = solver.Solve();

    solution = vector<vector<int> >(f.size());
    for (int i = 0; i < f.size(); i++) {
        solution[i] = vector<int>(f[i]->variables.size());
        int domain_size = Variable::getDomainSize(f[i]->variables);
        for (int j = 0; j < domain_size; j++) {
            if (x[i][j]->solution_value() > 0.99) {
                Variable::setAddress(f[i]->variables, j);
                for (int k = 0; k < f[i]->variables.size(); k++)
                    solution[i][k] = f[i]->variables[k]->value;
                break;
            }
        }
    }

    return objective->Value();

}


long double MN::knapsack_ub(MN &mn1, MN &mn2, vector<int> &order, int i_bound, long double q) {
    vector<LogPotential *> new_pots1, new_pots2;
    // Step 1: Cluster the functions so that there are "i" or fewer variables in each cluster
    // Store these new functions in new_pots1 and new_pots2 for MN1 and MN2 respectively
    vector<vector<LogPotential *> > buckets(order.size());
    vector<vector<LogPotential *> > buckets1(order.size());
    vector<vector<LogPotential *> > buckets2(order.size());
    vector<int> var_pos(order.size());
    for (int i = 0; i < order.size(); i++) {
        var_pos[order[i]] = i;
    }
    for (auto &lpotential : mn1.lpotentials) {
        int bucket_func = order.size();
        for (auto var : lpotential->variables) {
            if (bucket_func > var_pos[var->id]) {
                bucket_func = var_pos[var->id];
            }
        }
        buckets1[bucket_func].emplace_back(lpotential);
        buckets[bucket_func].emplace_back(lpotential);
    }
    for (auto &lpotential : mn2.lpotentials) {
        int bucket_func = order.size();
        for (auto var : lpotential->variables) {
            if (bucket_func > var_pos[var->id]) {
                bucket_func = var_pos[var->id];
            }
        }
        buckets2[bucket_func].emplace_back(lpotential);
        buckets[bucket_func].emplace_back(lpotential);
    }

    for (int i = 0; i < order.size(); i++) {
        //cout<<"UB = "<<ub<<endl;
        vector<set<int> > mini_bucket_func_ids;
        vector<vector<Variable *> > mini_bucket_vars;
        for (int j = 0; j < buckets[i].size(); j++) {
            // check if the function can be added to a minibucket
            //buckets[i][j]->print();
            int bucket_id = -1;
            size_t min_count = std::numeric_limits<size_t>::max();
            for (int k = 0; k < mini_bucket_vars.size(); k++) {
                vector<Variable *> tmp;
                Variable::set_union_variables(mini_bucket_vars[k], buckets[i][j]->variables, tmp);

                if ((int) tmp.size() <= i_bound) {
                    if ((tmp.size() - mini_bucket_vars[k].size()) < min_count) {
                        bucket_id = k;
                        min_count = tmp.size() - mini_bucket_vars[k].size();
                    }
                }
            }
            if (bucket_id == -1) {
                mini_bucket_vars.emplace_back(buckets[i][j]->variables);
                mini_bucket_func_ids.emplace_back(set<int>());
                mini_bucket_func_ids[mini_bucket_func_ids.size() - 1].insert(j);
            } else {
                vector<Variable *> tmp;
                Variable::set_union_variables(mini_bucket_vars[bucket_id], buckets[i][j]->variables, tmp);
                mini_bucket_vars[bucket_id] = tmp;
                mini_bucket_func_ids[bucket_id].insert(j);
                //cout<<tmp.size()<<" "<<i_bound<<endl;
            }
        }
        //cout<<mini_bucket_vars.size()<<endl;
        vector<LogPotential *> curr_pots1(mini_bucket_vars.size()), curr_pots2(mini_bucket_vars.size());
        for (int k = 0; k < mini_bucket_vars.size(); k++) {
            curr_pots1[k] = new LogPotential(mini_bucket_vars[k]);
            curr_pots2[k] = new LogPotential(mini_bucket_vars[k]);
            new_pots1.push_back(curr_pots1[k]);
            new_pots2.push_back(curr_pots2[k]);
        }
        for (int j = 0; j < buckets1[i].size(); j++) {
            bool found = false;
            for (int k = 0; k < mini_bucket_vars.size(); k++) {
                if (Variable::is_included(mini_bucket_vars[k], buckets1[i][j]->variables)) {
                    curr_pots1[k]->add(buckets1[i][j]);
                    found = true;
                    break;
                }
            }
            if (!found) cerr << "Bucket elim mistake\n";
        }
        for (int j = 0; j < buckets2[i].size(); j++) {
            bool found = false;
            for (int k = 0; k < mini_bucket_vars.size(); k++) {
                if (Variable::is_included(mini_bucket_vars[k], buckets2[i][j]->variables)) {
                    curr_pots2[k]->add(buckets2[i][j]);
                    found = true;
                    break;
                }
            }
            if (!found) cerr << "Bucket elim mistake2\n";
        }
    }
    // Based on new_pots1 and new_pots2 construct the lambdas
    vector<vector<LogPotential *> > var2pots(mn1.variables.size());
    vector<vector<int> > var2pots_index(mn1.variables.size());
    vector<vector<Lambda *> > f2lambdas(new_pots1.size());
    vector<Lambda *> all_lambdas;

    for (int i = 0; i < new_pots1.size(); i++) {
        new_pots1[i]->id = i;
        for (int j = 0; j < new_pots1[i]->variables.size(); j++) {
            var2pots[new_pots1[i]->variables[j]->id].push_back(new_pots1[i]);
            var2pots_index[new_pots1[i]->variables[j]->id].push_back(j);
        }
    }

    int k = 0;
    for (int i = 0; i < var2pots.size(); i++) {
        for (int j = 0; j < var2pots[i].size(); j++) {
            for (int e = j + 1; e < var2pots[i].size(); e++) {
                all_lambdas.push_back(new Lambda());
                all_lambdas[k]->vid = mn1.variables[i]->id;
                all_lambdas[k]->f1 = var2pots[i][j];
                all_lambdas[k]->f2 = var2pots[i][e];
                all_lambdas[k]->v1 = var2pots_index[i][j];
                all_lambdas[k]->v2 = var2pots_index[i][e];
                all_lambdas[k]->values = vector<long double>(mn1.variables[i]->domain_size);
                for (int a = 0; a < all_lambdas[k]->values.size(); a++)
                    all_lambdas[k]->values[a] = myRandom::getDouble() * (rand() % 2 > 0 ? 1.0 : -1.0);;
                //all_lambdas[k]->print();
                f2lambdas[var2pots[i][j]->id].push_back(all_lambdas[k]);
                f2lambdas[var2pots[i][e]->id].push_back(all_lambdas[k]);
                k++;
            }
        }
    }

    long double epsilon = 0.1;
    vector<vector<int> > solution;
    long double ub = std::numeric_limits<long double>::max();
    long double best_ub = std::numeric_limits<long double>::max();
    int best_ub_counter = 0;
    cout << "Running the Knapsack Based Bounding Algorithm\n";
    std::time_t start_time = std::time(nullptr);
    for (int iter = 0; iter < 1000; iter++) {
        long double iter_ub = solve_optimally(new_pots1, new_pots2, f2lambdas, q, solution);
        // Uncomment the following to use fast linear time LP algorithms
        //long double iter_ub=solve_using_LP(new_pots1,new_pots2,f2lambdas,q,solution);
        if (ub < iter_ub) epsilon /= 2.0;
        for (int i = 0; i < all_lambdas.size(); i++)
            all_lambdas[i]->updateValue(epsilon, solution);
        //if(iter%10==0) cout<<"UB at iter = "<<iter<< " is "<<iter_ub<<endl;
        ub = iter_ub;
        if (best_ub > ub) {
            best_ub = ub;
            best_ub_counter = 0;
        }
        else {
            best_ub_counter++;
            if (best_ub_counter > 20) break;
        }
        std::time_t curr_time = std::time(nullptr);
        cout << std::setprecision(20) << iter + 1 << "," << curr_time - start_time << "," << best_ub << endl;
    }
    return best_ub;
}

