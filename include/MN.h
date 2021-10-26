

#ifndef CMPE_BB_MN_H
#define CMPE_BB_MN_H

#include <set>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <cfloat>
#include <cmath>


using namespace std;

/*
 * struct Variable
 * id:			the variable number (numbering starts from 0)
 * 				Default=-1 means that the variable is invalid
 * domain_size: the number of values in the domain of the variable.
 * 				For example, for binary variables domain_size=2
 * 				Default=-1 means that variable is in invalid state
 * 				and needs to be initialized
 * value:		the value currently assigned to the variable.
 * 				value should be in the set {0,...,domain_size}
 * 				Default=-1 means that variable is not assigned a value
 */
struct Variable {
    int id;
    int domain_size;
    int value;

    /*
     * Functions
     */
    Variable() :
            id(-1), domain_size(-1), value(-1) {
    }

    Variable(int id_, int domain_size_) :
            id(id_), domain_size(domain_size_), value(-1) {
    }

    ~Variable() {}

    /*
     * Get the domain size of a set of variables. This is the maximum size of the 1-D array
     * where you will store the function that maps each variable assignment to a real number.
     * If the domain size of each variable is d_i the maximum domain size equals
     * \prod_{i=1}^{k} d_i
     * In the function given below: variables.size()=k; d_i = variables[i]->domain_size
     */
    inline static int getDomainSize(const vector<Variable *> &variables) {
        int domain_size = 1;
        for (auto variable : variables)
            domain_size *= variable->domain_size;
        return domain_size;
    }

    /*
     * Get a 1-dimensional index into an array from a k-dimensional value assignment
     * If the domain size of each variable is d_i and the current value assigned
     * to a variable is v_i then address equals \sum_{i=1}^{k} v_i * \prod_{j=1}^{i-1} d_j
     * In the function given below: variables.size()=k; d_i = variables[i]->domain_size
     * and v_i = variables[i]->value
     */
    inline static int getAddress(const vector<Variable *> &variables) {
        int add_ress = 0;
        int multiplier = 1;
        for (auto variable : variables) {
            add_ress += (multiplier * variable->value);
            multiplier *= variable->domain_size;
        }
        return add_ress;
    }

    /*
     * Reverse of the function getAddress
     * Convert a 1-dimensional index into an array to a k-dimensional value assignment
     * If the domain size of each variable is d_i and the 1-D address is a then
     * v_i = Floor(a/\prod_{j=1}^{i-1} d_j) mod d_i
     * In the function given below: variables.size()=k; d_i = variables[i]->domain_size
     * and v_i = variables[i]->value
     */
    inline static void setAddress(const vector<Variable *> &variables, const int add_ress_) {
        int add_ress = add_ress_;
        for (auto variable : variables) {
            variable->value = add_ress % variable->domain_size;
            add_ress /= variable->domain_size;
        }
    }

    static bool compareVariablePtr(const Variable *a, const Variable *b) {
        return a->id < b->id;
    }

    inline static void sort_variables(vector<Variable *> &variables) {
        sort(variables.begin(), variables.end(), compareVariablePtr);
    }

    inline static void set_union_variables(const vector<Variable *> &var_set1,
                                           const vector<Variable *> &var_set2,
                                           vector<Variable *> &out_) {
        vector<Variable *> out;
        set_union(var_set1.begin(), var_set1.end(), var_set2.begin(),
                  var_set2.end(), back_inserter(out), compareVariablePtr);
        out_ = out;
    }

    inline static void set_intersection_variables(const vector<Variable *> &var_set1,
                                                  const vector<Variable *> &var_set2,
                                                  vector<Variable *> &out_) {
        vector<Variable *> out;
        set_intersection(var_set1.begin(), var_set1.end(), var_set2.begin(),
                         var_set2.end(), back_inserter(out), compareVariablePtr);
        out_ = out;
    }

    inline static bool is_included(const vector<Variable *> &var_set1, const vector<Variable *> &var_set2) {
        return includes(var_set1.begin(), var_set1.end(), var_set2.begin(),
                        var_set2.end(), compareVariablePtr);
    }

    inline static void set_difference_variables(const vector<Variable *> &var_set1,
                                                const vector<Variable *> &var_set2,
                                                vector<Variable *> &out_) {
        vector<Variable *> out;
        set_difference(var_set1.begin(), var_set1.end(), var_set2.begin(),
                       var_set2.end(), back_inserter(out), compareVariablePtr);
        out_ = out;
    }
};

/*
 * Note that all potentials are log-potentials.
 * A log-potential is a pair <X,table> where
 *       X is a the scope of the potential and
 *       table is a set of weights, one for each possible assignment to X
 *       We assume that variables are sorted. This helps fast max_out and min_out operations
 */
struct LogPotential {
    vector<Variable *> variables;
    vector<long double> table;
    int id;

    LogPotential() : id(-1) {
        table = vector<long double>(1, 0.0);
    }

    LogPotential(const LogPotential &l) = default;

    explicit LogPotential(vector<Variable *> &variables_) : variables(variables_), id(-1) {
        table = vector<long double>(Variable::getDomainSize(variables), 0.0);
    }

    ~LogPotential() {
        table.clear();
        variables.clear();
    }

    inline long double getValue() const {
        return table[Variable::getAddress(variables)];
    }

    void instantiate_evidence(LogPotential *l, Variable *var, int value) {
        vector<Variable *> tmp(1);
        tmp[0] = var;
        Variable::set_difference_variables(variables, tmp, l->variables);
        if (l->variables.size() == this->variables.size()) {
            l->table = this->table;
        } else {
            var->value = value;
            l->table = vector<long double>(Variable::getDomainSize(l->variables));
            for (int i = 0; i < l->table.size(); i++) {
                Variable::setAddress(l->variables, i);
                l->table[i] = this->getValue();
            }
        }
    }

    void sort() {
        vector<Variable *> variables_ = variables;
        vector<long double> table_(table.size());
        Variable::sort_variables(variables_);
        for (int i = 0; i < table.size(); i++) {
            Variable::setAddress(variables_, i);
            table_[i] = getValue();
        }
        variables = variables_;
        table = table_;
    }

    void print() {
        cout << "n= " << variables.size() << endl;
        for (int i = 0; i < variables.size(); i++) {
            cout << variables[i]->id << " ";
        }
        cout << endl;
        cout << "Table entries = " << table.size() << endl;
        for (int i = 0; i < table.size(); i++) {
            cout << table[i] << " ";
        }
        cout << endl;
    }

    // In place add Potential
    void add(LogPotential *l) {
        vector<Variable *> variables_;
        Variable::set_union_variables(variables, l->variables, variables_);
        vector<long double> table_(Variable::getDomainSize(variables_), 0.0);
        for (int i = 0; i < table_.size(); i++) {
            Variable::setAddress(variables_, i);
            table_[i] = this->getValue() + l->getValue();
        }
        variables = variables_;
        table = table_;
    }

    // In place add Potential
    void subtract(LogPotential *l) {
        vector<Variable *> variables_;
        Variable::set_union_variables(variables, l->variables, variables_);
        vector<long double> table_(Variable::getDomainSize(variables_), 0.0);
        for (int i = 0; i < table_.size(); i++) {
            Variable::setAddress(variables_, i);
            table_[i] = this->getValue() - l->getValue();
        }
        variables = variables_;
        table = table_;
    }

    void maxOutVariables(LogPotential &f_max, vector<Variable *> &common_variables) const {
        f_max = LogPotential();
        if (common_variables.empty()) {
            return;
        }
        if (Variable::is_included(this->variables, common_variables)) {
            f_max.variables = common_variables;
            f_max.table = vector<long double>(Variable::getDomainSize(f_max.variables), -DBL_MAX);
            for (int i = 0; i < this->table.size(); i++) {
                Variable::setAddress(this->variables, i);
                int j = Variable::getAddress(f_max.variables);
                f_max.table[j] = max(f_max.table[j], this->table[i]);
            }
        } else {
            cerr << "Maxing out error: variables\n";
        }
    }

    void maxOutVariable(LogPotential &l, Variable *var) const {
        vector<Variable *> tmp;
        tmp.push_back(var);
        if (Variable::is_included(this->variables, tmp)) {
            Variable::set_difference_variables(this->variables, tmp, l.variables);
            l.table = vector<long double>(Variable::getDomainSize(l.variables), -DBL_MAX);
            for (int i = 0; i < l.table.size(); i++) {
                Variable::setAddress(l.variables, i);
                for (int v = 0; v < var->domain_size; v++) {
                    var->value = v;
                    long double curr_value = this->getValue();
                    if (l.table[i] < curr_value) l.table[i] = curr_value;
                }
            }
        } else {
            cerr << "Maxing out error\n";
        }
    }

    //Useful static functions
    static void add(const vector<LogPotential *> &lpotentials, LogPotential &out) {
        out = LogPotential();
        for (const auto &lpotential : lpotentials) {
            Variable::set_union_variables(lpotential->variables, out.variables, out.variables);
        }
        out.table = vector<long double>(Variable::getDomainSize(out.variables), 0.0);
        for (int i = 0; i < out.table.size(); i++) {
            Variable::setAddress(out.variables, i);
            for (const auto &lpotential : lpotentials) {
                out.table[i] += lpotential->getValue();
            }
        }
    }

    static void subtract(LogPotential *lpot1, LogPotential *lpot2, LogPotential &out) {
        vector<Variable *> tmp;
        if (!(Variable::is_included(lpot1->variables, lpot2->variables))) {
            return;
        }
        out = *lpot1;
        for (int i = 0; i < out.table.size(); i++) {
            Variable::setAddress(out.variables, i);
            out.table[i] -= lpot2->getValue();
        }
    }

    static void MaxOutVariable(const vector<LogPotential *> &lpotentials, Variable *var, LogPotential &out) {
        out = LogPotential();
        for (const auto &lpotential : lpotentials) {
            vector<Variable *> tmp;
            Variable::set_union_variables(lpotential->variables, out.variables, tmp);
            out.variables = tmp;
        }
        vector<Variable *> curr_variable(1), tmp;
        curr_variable[0] = var;
        Variable::set_difference_variables(out.variables, curr_variable, tmp);
        if (out.variables.size() == tmp.size()) {
            cerr << "Something wrong\n";
            exit(-1);
        }
        out.variables = tmp;
        out.table = vector<long double>(Variable::getDomainSize(out.variables), 0.0);
        //cout<<"Status: "<<endl;
        //cout<<var->id<<" "<<out.table.size()<<endl;
        for (int i = 0; i < out.table.size(); i++) {
            Variable::setAddress(out.variables, i);
            for (int v = 0; v < var->domain_size; v++) {
                var->value = v;
                long double curr_tab_value = 0.0;
                for (const auto &lpotential : lpotentials) {
                    curr_tab_value += lpotential->getValue();
                    if (isnan(curr_tab_value)) {
                        cerr << "Nan generated \n";
                        cerr << lpotential->getValue() << endl;
                        cerr << "------\n";
                        lpotential->print();
                        exit(-1);
                    }
                }
                if (v == 0) out.table[i] = curr_tab_value;
                else if (curr_tab_value > out.table[i]) out.table[i] = curr_tab_value;
            }
        }
        //if (out.table.size()>=1) {
        //  cout << "Table entry " << out.table[0] << endl;
        //}
    }

    static void MinOutVariable(const vector<LogPotential *> &lpotentials, Variable *var, LogPotential &out) {
        out = LogPotential();
        for (const auto &lpotential : lpotentials) {
            vector<Variable *> tmp;
            Variable::set_union_variables(lpotential->variables, out.variables, tmp);
            out.variables = tmp;
        }
        vector<Variable *> curr_variable(1), tmp;
        curr_variable[0] = var;
        Variable::set_difference_variables(out.variables, curr_variable, tmp);
        out.variables = tmp;
        out.table = vector<long double>(Variable::getDomainSize(out.variables), 0.0);

        for (int i = 0; i < out.table.size(); i++) {
            Variable::setAddress(out.variables, i);
            for (int v = 0; v < var->domain_size; v++) {
                var->value = v;
                long double curr_tab_value = 0.0;
                for (const auto &lpotential : lpotentials) {
                    curr_tab_value += lpotential->getValue();
                }
                if (v == 0) out.table[i] = curr_tab_value;
                else if (curr_tab_value < out.table[i]) out.table[i] = curr_tab_value;
            }
        }
    }
};

class Utils {
public:
    static void getMinDegreeOrder(vector<Variable *> &variables, vector<LogPotential *> &functions, vector<int> &order);

    static void getMinFillOrder(vector<Variable *> &variables, vector<LogPotential *> &functions, vector<int> &order);
};

struct MN {
    vector<Variable *> variables;
    vector<LogPotential *> lpotentials;

    MN() = default;

    ~MN() {
        for (int i = 0; i < lpotentials.size(); i++)
            delete (lpotentials[i]);
    }

    void copy(MN &mn) {
        mn.variables = variables;
        mn.lpotentials = vector<LogPotential *>(lpotentials.size());
        for (int i = 0; i < lpotentials.size(); i++) {
            mn.lpotentials[i] = new LogPotential();
            mn.lpotentials[i]->variables = lpotentials[i]->variables;
            mn.lpotentials[i]->table = lpotentials[i]->table;
        }
    }

    void readMN(string filename_);

    void readMN2(string filename_, MN &mn1);

    inline long double getValue() {
        long double logp = 0.0;
        for (auto &lpotential : lpotentials)
            logp += lpotential->getValue();
        return logp;
    }

    long double instantiate_evidence(MN &mn, Variable *var, int value) {
        long double ret_value = 0.0;
        mn.variables = this->variables;
        var->value = value;
        mn.lpotentials = vector<LogPotential *>();
        for (int i = 0; i < lpotentials.size(); i++) {
            LogPotential *l = new LogPotential();
            lpotentials[i]->instantiate_evidence(l, var, value);
            if ((int) l->table.size() == 1) {
                ret_value += l->table[0];
                delete (l);
            } else {
                mn.lpotentials.push_back(l);
            }
        }
        return ret_value;
    }

    void ILP_Bound(MN &mn_c, long double logp);

    static long double knapsack_ub(MN &mn1, MN &mn2, vector<int> &order, int i_bound, long double q);


};


#endif //CMPE_BB_MN_H
