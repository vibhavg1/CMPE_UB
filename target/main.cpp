

#include <vector>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <sstream>


#include "MN.h"
#include "CG.h"


using namespace std;

void print_help(const string &program_name) {
    cerr << "Usage: " << program_name << " -m1 <uaifilename1> -m2 <uaifilename2> -q <q-value>\n";
    cerr << "\t Other Options: [-i -l -s ]\n";
    cerr << "-------------------------------------------------------------------------\n";
    cerr << "\t\t Details on Required Option\n";
    cerr << "\t\t\t uaifilename1 and uaifilename2: are evidence instantiated Markov networks in UAI format\n";
    cerr << "\t\t\t q-value: (Real number): constraint on weight of the assignment in CMPE\n";
    cerr << "-------------------------------------------------------------------------\n";
    cerr << "\t\t Details on Other Options and Default values\n";
    cerr << "\t\t\t -l     [double]: lambda value (default 0.1)\n";
    cerr << "\t\t\t -i     [int]: i-bound to use (default 10)\n";
    cerr << "\t\t\t -s     [int]: Seed for Repeatability; default 1000000L\n";

}

void print_help_mnist(const string &program_name) {
    cerr << "Usage: " << program_name << " -m <lr-filename> -t <test_data> -o <modified-test_data> -u <int>\n";
    cerr << "\t Other Options: [-i -l -s ]\n";
    cerr << "-------------------------------------------------------------------------\n";
    cerr << "\t\t Details on Required Option\n";
    cerr << "\t\t\t lr-filename: Logistic Regression model\n";
    cerr << "\t\t\t modified-test-data: Results of experiments will be stored here\n";
    cerr << "\t\t\t -u[int]: use 1 for univariate and 2 for grid\n";
    cerr << "-------------------------------------------------------------------------\n";
    cerr << "\t\t Details on Other Options and Default values\n";
    cerr << "\t\t\t -l     [double]: tol value (default 0.1)\n";
    cerr << "\t\t\t -i     [int]: i-bound to use (default 10)\n";
    cerr << "\t\t\t -s     [int]: Seed for Repeatability; default 1000000L\n";
    //cerr << "\t\t\t -si    [int]: print status every integer seconds; default 1\n";
    //cerr << "\t\t\t -w  [string]: Write file in MPS format and store it in string\n";
}
/*
 * The following four functions are for the MNIST dataset
 */

// Read the SVM weights from a file. Use Scikit learn to learn the linear SVM
void readWeights(vector<long double> &weights, string filename) {
    ifstream in(filename.c_str());
    if (!in.good()) {
        cerr
                << "Either the training, validation or test dataset file is missing\n";
        return;
    }
    weights.clear();
    string line;
    while (getline(in, line)) {
        stringstream inss(line);
        long double m;
        while (inss >> m) {
            weights.push_back(m);
            if (inss.peek() == ' ')
                inss.ignore();
        }
    }
    in.close();
}

// Read the MNIST data
void readData(vector<vector<int> > &data, string filename) {
    ifstream in(filename);
    if (!in.good())
        return;
    string line;
    data.clear();
    while (getline(in, line)) {
        stringstream inss(line);
        vector<int> row;
        int m;
        while (inss >> m) {
            row.push_back(m);
            if (inss.peek() == ' ')
                inss.ignore();
        }
        data.push_back(row);
    }
    int numvars = data[0].size();
    cout << "numvars = " << numvars << endl;
    in.close();
}

// Univariate Markov network: Distance function for MNIST
void constructUniVariateMN(vector<int> &example, MN &mn, vector<Variable *> &variables) {
    if (example.size() != variables.size())
        cerr << "Something wrong: in function main; size of variable and examples do not match\n";
    int nvars = variables.size();
    mn.variables = variables;
    mn.lpotentials = vector<LogPotential *>(nvars);
    for (int i = 0; i < nvars; i++) {
        vector<Variable *> v;
        v.push_back(mn.variables[i]);
        mn.lpotentials[i] = new LogPotential(v);
        mn.lpotentials[i]->id = i;
        mn.variables[i]->value = example[i];
        mn.lpotentials[i]->table[0] = 0.0;
        mn.lpotentials[i]->table[1] = 0.0;
        mn.lpotentials[i]->table[Variable::getAddress(v)] += 1.0;
    }
}

// Grid Markov network given the test example. Increase tolerance to increase the affinity in the Grid
void constructGridMN(vector<int> &example, MN &mn, vector<Variable *> &variables, long double tol = 1.0) {
    if (example.size() != variables.size())
        cerr << "Something wrong: in function main; size of variable and examples do not match\n";
    int nvars = variables.size();
    mn.variables = variables;
    // Construct Grid from variables
    int rows = sqrt(nvars);
    int cols = rows;

    cout << "Number of rows = " << rows << endl;
    vector<vector<Variable *> > grid(rows);
    for (int i = 0; i < rows; i++)
        grid[i] = vector<Variable *>(cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            grid[i][j] = variables[i * rows + j];
        }
    }
    int k = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols - 1; j++) {
            int v1 = i * rows + j;
            int v2 = i * rows + j + 1;
            LogPotential *lpot = new LogPotential();
            lpot->variables.push_back(variables[v1]);
            lpot->variables.push_back(variables[v2]);
            lpot->table = vector<long double>(4);
            lpot->table[0] = tol;
            lpot->table[1] = 0.0;
            lpot->table[2] = 0.0;
            lpot->table[3] = tol;
            variables[v1]->value = example[v1];
            variables[v2]->value = example[v2];
            lpot->table[Variable::getAddress(lpot->variables)] += tol * tol;
            lpot->id = k++;
            mn.lpotentials.push_back(lpot);
        }
    }
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows - 1; j++) {
            int v1 = i + j * rows;
            int v2 = i + (j + 1) * rows;
            LogPotential *lpot = new LogPotential();
            lpot->variables.push_back(variables[v1]);
            lpot->variables.push_back(variables[v2]);
            lpot->table = vector<long double>(4);
            lpot->table[0] = tol;
            lpot->table[1] = 0.0;
            lpot->table[2] = 0.0;
            lpot->table[3] = tol;
            variables[v1]->value = example[v1];
            variables[v2]->value = example[v2];
            lpot->table[Variable::getAddress(lpot->variables)] += tol * tol;
            lpot->id = k++;
            mn.lpotentials.push_back(lpot);
        }
    }
}

int main(int argc, char *argv[]) {
    srand(1000000L); // for repeatability
    string uai_filename1;
    string uai_filename2;
    int i_bound = 10;
    long double q;
    long double lambda;
    bool uaioption1 = false, uaioption2 = false, qoption = false;
    if (argc == 1) {
        print_help(argv[0]);
        exit(-1);
    }
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-m1") == 0) {
            uai_filename1 = argv[i + 1];
            uaioption1 = true;
        } else if (strcmp(argv[i], "-m2") == 0) {
            uai_filename2 = argv[i + 1];
            uaioption2 = true;
        } else if (strcmp(argv[i], "-l") == 0) {
            lambda = atof(argv[i + 1]);
        } else if (strcmp(argv[i], "-i") == 0) {
            i_bound = atoi(argv[i + 1]);
        } else if (strcmp(argv[i], "-s") == 0) {
            int seed = atoi(argv[i + 1]);
            srand(seed);
        } else if (strcmp(argv[i], "-si") == 0) {
            //GlobalSearchOptions::print_interval = atoi(argv[i + 1]);
        } else if (strcmp(argv[i], "-q") == 0) {
            q = atof(argv[i + 1]);
            qoption = true;
        } else if (strcmp(argv[i], "-h") == 0) {
            print_help(argv[0]);
            exit(-1);
        }
    }
    if (!uaioption1) {
        cerr << "UAI file1 not specified\n";
        print_help(argv[0]);
        exit(-1);
    }
    if (!uaioption2) {
        cerr << "UAI file2 not specified\n";
        print_help(argv[0]);
        exit(-1);
    }
    if (!qoption) {
        cerr << "Q not specified\n";
        print_help(argv[0]);
        exit(-1);
    }
    cout << "Starting Experiment on UAI File " << uai_filename1 << endl;
    cerr << "i-bound=" << i_bound << endl;
    cout << "i-bound=" << i_bound << endl;
    cout << "q=" << q << endl;
    MN mn1, mn2;
    mn1.readMN(uai_filename1);
    mn2.readMN2(uai_filename2, mn1);

    if (mn1.variables.size() != mn2.variables.size()) {
        cerr << "Variable size mismatch\n";
        cerr << "Code requires the two Markov networks be defined over the same set of variables\n";
        exit(-1);
    }

    vector<int> order;
    Utils::getMinFillOrder(mn1.variables, mn1.lpotentials, order);
    long double ub;


    ub = CG::ub_search(mn1, mn2, order, i_bound, q, MINI_BUCKET).first;
    cout << "UB MBE is " << ub << endl;
    ub = CG::ub_search(mn1, mn2, order, i_bound, q, MINI_BUCKET_MM).first;
    cout << "UB MM is " << ub << endl;
    ub = CG::ub_search(mn1, mn2, order, i_bound, q, JOIN_GRAPH).first;
    cout << "JG Bound  = " << ub << endl;
    ub = MN::knapsack_ub(mn1, mn2, order, i_bound, q);
    cout << "UB Knapsack method is " << ub << endl;
    cout << "End Experiment on UAI File " << uai_filename1 << endl;
    cerr << "Exiting\n";

    exit(0);
}

// Main file for running the MNIST experiments
int main_mnist(int argc, char *argv[]) {
    srand(1000000L);
    string lr_filename;
    string test_filename;
    string out_filename;
    int max_time = 1200;
    int i_bound = 10;
    int sampling_number = 1000;
    bool univariate = true;
    long double q;
    long double tol = 1.0;
    bool uaioption1 = false, uaioption2 = false, outoption = false, qoption = false;
    if (argc == 1) {
        print_help(argv[0]);
        exit(-1);
    }
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-m") == 0) {
            lr_filename = argv[i + 1];
            uaioption1 = true;
        } else if (strcmp(argv[i], "-t") == 0) {
            test_filename = argv[i + 1];
            uaioption2 = true;
        } else if (strcmp(argv[i], "-o") == 0) {
            out_filename = argv[i + 1];
            outoption = true;
        } else if (strcmp(argv[i], "-l") == 0) {
            tol = atof(argv[i + 1]);
        } else if (strcmp(argv[i], "-u") == 0) {
            int c = atoi(argv[i + 1]);
            if (c == 1) univariate = true;
            else if (c == 2) univariate = false;
            else exit(-1);
        } else if (strcmp(argv[i], "-i") == 0) {
            i_bound = atoi(argv[i + 1]);
        } else if (strcmp(argv[i], "-s") == 0) {
            int seed = atoi(argv[i + 1]);
            srand(seed);
        } else if (strcmp(argv[i], "-h") == 0) {
            print_help(argv[0]);
            exit(-1);
        }
    }
    if (!uaioption1) {
        cerr << "UAI file1 not specified\n";
        print_help(argv[0]);
        exit(-1);
    }
    if (!uaioption2) {
        cerr << "UAI file2 not specified\n";
        print_help(argv[0]);
        exit(-1);
    }
    if (!outoption) {
        cerr << "Output file not specified\n";
        print_help(argv[0]);
        exit(-1);
    }

    // Construct the first Markov network using LR
    MN mn_lr;
    vector<long double> weights;
    readWeights(weights, lr_filename);
    cout << "Number of variables =" << weights.size() - 1 << endl;
    mn_lr.variables = vector<Variable *>(weights.size() - 1);
    for (int i = 0; i < mn_lr.variables.size(); i++) {
        mn_lr.variables[i] = new Variable(i, 2);
    }
    mn_lr.lpotentials = vector<LogPotential *>(mn_lr.variables.size());
    for (int i = 0; i < mn_lr.variables.size(); i++) {
        vector<Variable *> v;
        v.push_back(mn_lr.variables[i]);
        mn_lr.lpotentials[i] = new LogPotential(v);
        mn_lr.lpotentials[i]->id = i;
        mn_lr.lpotentials[i]->table[0] = 0.0;
        mn_lr.lpotentials[i]->table[1] = weights[i];
    }
    q = -weights[weights.size() - 1];
    vector<vector<int> > data;
    readData(data, test_filename);
    ofstream out(out_filename, ofstream::out);
    for (int e = 0; e < data.size(); e++) {
        // Construct the MN1
        MN mn1;
        if (univariate) {
            constructUniVariateMN(data[e], mn1, mn_lr.variables);
        } else {
            constructGridMN(data[e], mn1, mn_lr.variables, tol);
        }
        vector<int> order;
        Utils::getMinFillOrder(mn1.variables, mn1.lpotentials, order);


        long double ub = CG::ub_search(mn1, mn_lr, order, i_bound, q, MINI_BUCKET).first;
        cout << "UB MBE = " << ub << endl;

        ub = CG::ub_search(mn1, mn_lr, order, i_bound, q, MINI_BUCKET_MM).first;
        cout << "UB MBE MM= " << ub << endl;

        ub = CG::ub_search(mn1, mn_lr, order, i_bound, q, JOIN_GRAPH).first;
        cout << "UB JG = " << ub << endl;
        ub = MN::knapsack_ub(mn1, mn_lr, order, i_bound, q);
        cout << "Knapsack UB = " << ub << endl;
    }
    return 1;
}