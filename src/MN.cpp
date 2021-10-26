
#include <MN.h>
#include <myRandom.h>
#include <fstream>
#include <list>
#include <cfloat>


void Utils::getMinDegreeOrder(vector<Variable *> &variables, vector<LogPotential *> &functions, vector<int> &order) {
    vector<set<int> > degree(variables.size());
    vector<set<int> > graph(variables.size());

    for (auto &function : functions) {
        for (int j = 0; j < function->variables.size(); j++) {
            for (int k = j + 1; k < function->variables.size(); k++) {
                int a = function->variables[j]->id;
                int b = function->variables[k]->id;
                graph[a].insert(b);
                graph[b].insert(a);
            }
        }
    }
    int min_degree = variables.size();
    int max_degree = -1;
    for (int i = 0; i < variables.size(); i++) {
        degree[graph[i].size()].insert(i);
        if (graph[i].size() < min_degree) {
            min_degree = graph[i].size();
        }
        if (graph.size() > max_degree) {
            max_degree = graph.size();
        }
    }
    order = vector<int>();
    while (order.size() != variables.size()) {
        while (degree[min_degree].empty()) {
            min_degree++;
        }
        int curr = *(degree[min_degree].begin());
        degree[min_degree].erase(curr);
        for (auto i = graph[curr].begin(); i != graph[curr].end(); i++) {
            degree[graph[*i].size()].erase(*i);
            graph[*i].erase(curr);
            degree[graph[*i].size()].insert(*i);
            if (graph[*i].size() < min_degree) { min_degree = graph[*i].size(); }
        }
        order.push_back(curr);
    }
}

void Utils::getMinFillOrder(vector<Variable *> &variables, vector<LogPotential *> &functions, vector<int> &order) {
    double estimate = 0.0;
    int max_cluster_size = 0;
    order = vector<int>(variables.size());
    vector<set<int> > clusters(variables.size());
    vector<vector<bool> > adj_matrix(variables.size());

    // Create the interaction graph of the functions in this graphical model - i.e.
    // create a graph structure such that an edge is drawn between variables in the
    // model that appear in the same function
    for (int i = 0; i < variables.size(); i++) {
        adj_matrix[i] = vector<bool>(variables.size());
    }
    vector<set<int> > graph(variables.size());
    vector<bool> processed(variables.size());
    for (auto &function : functions) {
        for (int j = 0; j < function->variables.size(); j++) {
            for (int k = j + 1; k < function->variables.size(); k++) {
                int a = function->variables[j]->id;
                int b = function->variables[k]->id;
                graph[a].insert(b);
                graph[b].insert(a);
                adj_matrix[a][b] = true;
                adj_matrix[b][a] = true;
            }
        }
    }
    list<int> zero_list;

    // For i = 1 to number of variables in the model
    // 1) Identify the variables that if deleted would add the fewest number of edges to the
    //    interaction graph
    // 2) Choose a variable, pi(i), from among this set
    // 3) Add an edge between every pair of non-adjacent neighbors of pi(i)
    // 4) Delete pi(i) from the interaction graph
    for (int i = 0; i < variables.size(); i++) {
        // Find variables with the minimum number of edges added
        double min = DBL_MAX;
        int min_id = -1;
        bool first = true;

        // Flag indicating whether the variable to be removed is from the
        // zero list - i.e. adds no edges to interaction graph when deleted
        bool fromZeroList = false;
        // Vector to keep track of the ID of each minimum fill variable
        vector<int> minFillIDs;

        // If there are no variables that, when deleted, add no edges...
        if (zero_list.empty()) {

            // For each unprocessed (non-deleted) variable
            for (int j = 0; j < variables.size(); j++) {
                if (processed[j])
                    continue;
                double curr_min = 0.0;
                for (auto a = graph[j].begin();
                     a != graph[j].end(); a++) {
                    auto b = a;
                    b++;
                    for (; b != graph[j].end(); b++) {
                        if (!adj_matrix[*a][*b]) {
                            curr_min += (variables[*a]->domain_size
                                         * variables[*b]->domain_size);
                            if (curr_min > min)
                                break;
                        }
                    }
                    if (curr_min > min)
                        break;
                }

                // Store the first non-deleted variable as a potential minimum
                if (first) {
                    minFillIDs.push_back(j);
                    min = curr_min;
                    first = false;
                } else {
                    // If this is a new minimum...
                    if (min > curr_min) {
                        min = curr_min;
                        minFillIDs.clear();
                        minFillIDs.push_back(j);
                    }
                        // Otherwise, if the number of edges removed is also a minimum, but
                        // the minimum is zero
                    else if (curr_min < DBL_MIN) {
                        zero_list.push_back(j);
                    }
                        // Else if this is another potential min_fill
                    else if (min == curr_min) {
                        minFillIDs.push_back(j);
                    }
                }
            }
        }
            // Else...delete variables from graph that don't add any edges
        else {
            min_id = zero_list.front();
            zero_list.pop_front();
            fromZeroList = true;
        }

        // If not from zero_list, choose one of the variables at random
        // from the set of min fill variables
        if (!fromZeroList) {
            int indexInVector;
            indexInVector = myRandom::getInt(minFillIDs.size());
            min_id = minFillIDs[indexInVector];
        }

        //cout<<"order["<<i<<"]= "<<min_id<<" "<<flush;
        assert(min_id != -1);
        order[i] = min_id;
        // Now form the cluster
        clusters[i] = graph[min_id];
        clusters[i].insert(min_id);

        // Trinagulate min id and remove it from the graph
        for (auto a = graph[min_id].begin();
             a != graph[min_id].end(); a++) {
            auto b = a;
            b++;
            for (; b != graph[min_id].end(); b++) {
                if (!adj_matrix[*a][*b]) {
                    adj_matrix[*a][*b] = true;
                    adj_matrix[*b][*a] = true;
                    graph[*a].insert(*b);
                    graph[*b].insert(*a);
                }
            }
        }
        for (auto a = graph[min_id].begin();
             a != graph[min_id].end(); a++) {
            graph[*a].erase(min_id);
            adj_matrix[*a][min_id] = false;
            adj_matrix[min_id][*a] = false;
        }
        graph[min_id].clear();
        processed[min_id] = true;
    }
    /*
    // compute the estimate
    for (auto & cluster : clusters) {
        if ((int) cluster.size() > max_cluster_size)
            max_cluster_size = (int) cluster.size();
        double curr_estimate = 1.0;
        for (std::__1::__tree_const_iterator<int, std::__1::__tree_node<int, void *> *, long>::value_type j : cluster) {
            curr_estimate *= (double) variables[j]->d;
        }
        estimate += curr_estimate;
    }
    cout<<"Max cluster size = "<<max_cluster_size<<endl;
    */
}


// Read the Markov network
void MN::readMN(string filename) {
    ifstream infile(filename.c_str());
    int num_variables;
    string tmp_string;
    infile >> tmp_string;
    if (tmp_string.compare("MARKOV") != 0) {
        cerr << "Not a Markov network\n";
        //exit(-1);
        //return;
    }
    infile >> num_variables;
    // Read domains
    variables = vector<Variable *>(num_variables);
    for (int i = 0; i < num_variables; i++) {
        int domain_size;
        infile >> domain_size;
        variables[i] = new Variable(i, domain_size);
    }
    int num_functions;
    infile >> num_functions;
    vector<vector<Variable *> > scope(num_functions);
    for (int i = 0; i < num_functions; i++) {
        // Read parents of variables
        int num_vars_in_func;
        infile >> num_vars_in_func;
        scope[i] = vector<Variable *>(num_vars_in_func);
        for (int j = 0; j < num_vars_in_func; j++) {
            int temp;
            infile >> temp;
            scope[i][j] = variables[temp];
        }
    }
    lpotentials = vector<LogPotential *>(num_functions);
    srand(100000000L);
    for (int i = 0; i < num_functions; i++) {
        int num_entries;
        infile >> num_entries;
        lpotentials[i] = new LogPotential();
        lpotentials[i]->variables = scope[i];
        lpotentials[i]->id = i;
        int num_values = Variable::getDomainSize(scope[i]);
        lpotentials[i]->table = vector<long double>(num_values);
        for (int j = 0; j < num_values; j++) {
            Variable::setAddress(scope[i], j);
            long double value;
            infile >> value;
            int entry = Variable::getAddress(lpotentials[i]->variables);
            if (value > 0.0)
                lpotentials[i]->table[entry] = log(value);
            else {
                lpotentials[i]->table[entry] =
                        ((double) rand() / ((double) RAND_MAX + 1.0)) * (rand() % 2 > 0 ? 1.0 : -1.0);
                //cerr << "Cannot handle zeros: Log-potentials\n";
                //exit(-1);
            }
        }
        lpotentials[i]->sort();
    }
    infile.close();
}

void MN::readMN2(string filename, MN &mn1) {
    ifstream infile(filename);
    int num_variables;
    string tmp_string;
    infile >> tmp_string;
    if (tmp_string.compare("MARKOV") != 0) {
        cerr << "Not a Markov network\n";
        //exit(-1);
        //return;
    }
    infile >> num_variables;
    if (num_variables != mn1.variables.size()) {
        cerr << "Markov networks do not match in number of variables\n";
        exit(-1);
        return;
    }
    // Read domains
    variables = mn1.variables;
    for (int i = 0; i < num_variables; i++) {
        int domain_size;
        infile >> domain_size;
        if (variables[i]->domain_size != domain_size) {
            cerr << "Variables in Markov networks do not match; different domains\n";
            exit(-1);
            return;
        }
    }
    int num_functions;
    infile >> num_functions;
    vector<vector<Variable *> > scope(num_functions);
    for (int i = 0; i < num_functions; i++) {
        // Read parents of variables
        int num_vars_in_func;
        infile >> num_vars_in_func;
        scope[i] = vector<Variable *>(num_vars_in_func);
        for (int j = 0; j < num_vars_in_func; j++) {
            int temp;
            infile >> temp;
            scope[i][j] = variables[temp];
        }
    }
    lpotentials = vector<LogPotential *>(num_functions);
    srand(100000000L);
    for (int i = 0; i < num_functions; i++) {
        int num_entries;
        infile >> num_entries;
        lpotentials[i] = new LogPotential();
        lpotentials[i]->variables = scope[i];
        lpotentials[i]->id = i;
        int num_values = Variable::getDomainSize(scope[i]);
        lpotentials[i]->table = vector<long double>(num_values);
        for (int j = 0; j < num_values; j++) {
            Variable::setAddress(scope[i], j);
            long double value;
            infile >> value;
            int entry = Variable::getAddress(lpotentials[i]->variables);
            if (value > 0.0)
                lpotentials[i]->table[entry] =
                        log(value) + ((double) rand() / ((double) RAND_MAX + 1.0)) * (rand() % 2 > 0 ? 1.0 : -1.0);
            else {
                lpotentials[i]->table[entry] =
                        ((double) rand() / ((double) RAND_MAX + 1.0)) * (rand() % 2 > 0 ? 1.0 : -1.0);
                //cerr << "Cannot handle zeros: Log-potentials\n";
                //exit(-1);
            }
        }
        lpotentials[i]->sort();
    }
    infile.close();
}

