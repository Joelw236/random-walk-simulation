#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <string>
using namespace std;

mt19937 gen(random_device{}());

vector<int> random_walk(int n, double p){

    int S = 0;
    vector<int> vals;

    bernoulli_distribution step_up(1.0 - p);

    vals.push_back(S);

    for(int i = 0; i < n; i++){

        if(step_up(gen)){
            S += 1;
        }
        else{
            S -= 1;
        }

        vals.push_back(S);
    }

    return vals;
}

vector<int> Sn(const vector<int>& val){
    vector<int> S;
    for(int i = 1; i < val.size(); i++){
        if(val[i-1] < val[i]){
            S.push_back(1);
        }
        else{
            S.push_back(-1);
        }
    }
    return S;
}

double Ev(const vector<int>& val){
    double E = 0;
    int n = val.size();
    if (n == 0) return 0.0;
    for(int i = 0; i < n; i++){
        E += static_cast<double>(val[i]) / n;
    }
    return E;
}

double Var(const vector<int>& val){
    double mean = Ev(val);
    double var = 0;
    for(int i = 0; i < val.size(); i++ ){
        var += (val[i] - mean)*(val[i] - mean);
    }
    var = var/(val.size()-1);
    return var;
}

double drift_check(const vector<vector<int>>& data, int t) {
    double sum = 0;
    int count = 0;

    for (const auto& walk : data) {
        if (t + 1 < walk.size()) {
            sum += (walk[t + 1] - walk[t]);
            count++;
        }
    }

    return sum / count;
}


vector<double> drift_over_time(const vector<vector<int>>& data) {
    if (data.empty()) return {};

    size_t T = data[0].size();

    for (const auto& walk : data) {
        T = min(T, walk.size());
    }

    vector<double> drift(T - 1, 0.0);
    vector<int> count(T - 1, 0);

    for (const auto& walk : data) {
        for (size_t t = 0; t < T - 1; t++) {
            drift[t] += walk[t + 1] - walk[t];
            count[t]++;
        }
    }

    for (size_t t = 0; t < T - 1; t++) {
        drift[t] /= count[t];
    }

    return drift;
}


vector<double> quad_variation(const vector<vector<int>>& data) {
    if (data.empty()) return {};

    size_t T = data[0].size();

    for (const auto& walk : data) {
        T = min(T, walk.size());
    }

    vector<double> qv(T - 1, 0.0);
    vector<int> count(T - 1, 0);

    for (const auto& walk : data) {
        for (size_t t = 0; t < T - 1; t++) {

            int x = walk[t];
            int x_next = walk[t + 1];

            int increment = x_next - x;

            qv[t] += increment * increment;

            count[t]++;
        }
    }

    for (size_t t = 0; t < T - 1; t++) {
        qv[t] /= count[t];
    }

    return qv;
}

// empirical variance over time across all walks.
vector<double> variance_over_time(const vector<vector<int>>& data) {
    if (data.empty()) return {};

    size_t T = data[0].size();
    for (const auto& walk : data) {
        T = min(T, walk.size());
    }

    vector<double> mean(T, 0.0);
    vector<double> mean_sq(T, 0.0);
    vector<double> var(T, 0.0);

    for (const auto& walk : data) {
        for (size_t t = 0; t < T; t++) {
            mean[t] += walk[t];
            mean_sq[t] += static_cast<double>(walk[t]) * static_cast<double>(walk[t]);
        }
    }

    for (size_t t = 0; t < T; t++) {
        mean[t] /= data.size();
        mean_sq[t] /= data.size();
        var[t] = mean_sq[t] - mean[t] * mean[t];
    }

    return var;
}

//  write a matrix-like CSV (one walk per row).
void write_matrix_csv(const vector<vector<int>>& data, const string& filename) {
    ofstream file(filename);
    for (const vector<int>& v : data) {
        for (size_t i = 0; i < v.size(); i++) {
            file << v[i];
            if (i != v.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
}

//  write terminal positions for histogram / normal fit.
void write_terminal_csv(const vector<vector<int>>& data, const string& filename) {
    ofstream file(filename);
    file << "walk_index,terminal_position,scaled_terminal\n";
    if (data.empty()) return;

    double n = static_cast<double>(data[0].size() - 1);
    double scale = sqrt(n);

    for (size_t i = 0; i < data.size(); i++) {
        int terminal = data[i].back();
        file << i << "," << terminal << "," << (terminal / scale) << "\n";
    }
}

// write variance convergence.
void write_variance_convergence_csv(const vector<vector<int>>& data, const string& filename) {
    ofstream file(filename);
    file << "t,variance,variance_over_t,theory_variance\n";
    if (data.empty()) return;

    vector<double> var = variance_over_time(data);
    for (size_t t = 0; t < var.size(); t++) {
        double theory = static_cast<double>(t);
        double ratio = (t == 0 ? 0.0 : var[t] / static_cast<double>(t));
        file << t << "," << var[t] << "," << ratio << "," << theory << "\n";
    }
}

// parameter study for p-values.
void write_parameter_study_csv(int steps, int walks_per_p, const vector<double>& ps, const string& filename) {
    ofstream file(filename);
    file << "p,empirical_step_mean,theoretical_step_mean,empirical_terminal_mean,theoretical_terminal_mean\n";

    for (double p : ps) {
        vector<vector<int>> local_data;
        local_data.reserve(walks_per_p);

        for (int i = 0; i < walks_per_p; i++) {
            local_data.push_back(random_walk(steps, p));
        }

        double step_sum = 0.0;
        long long step_count = 0;
        double terminal_sum = 0.0;

        for (const auto& walk : local_data) {
            terminal_sum += walk.back();
            for (size_t t = 0; t + 1 < walk.size(); t++) {
                step_sum += static_cast<double>(walk[t + 1] - walk[t]);
                step_count++;
            }
        }

        double empirical_step_mean = step_sum / static_cast<double>(step_count);
        double theoretical_step_mean = 1.0 - 2.0 * p;
        double empirical_terminal_mean = terminal_sum / static_cast<double>(walks_per_p);
        double theoretical_terminal_mean = steps * theoretical_step_mean;

        file << p << ","
             << empirical_step_mean << ","
             << theoretical_step_mean << ","
             << empirical_terminal_mean << ","
             << theoretical_terminal_mean << "\n";
    }
}

int main() {
    // Baseline experiment: symmetric random walk.
    int num_walks = 5000;
    int steps = 1000;
    double p = 0.5;

    vector<vector<int>> data;
    data.reserve(num_walks);

    for(int i = 0; i < num_walks; i++){
        data.push_back(random_walk(steps, p));
    }

    // File 1: all walks, one row per walk.
    write_matrix_csv(data, "data.csv");

    // File 2: drift and quadratic variation over time.
    vector<double> drift = drift_over_time(data);
    vector<double> qv = quad_variation(data);

    ofstream stats("random_walk_stats.csv");
    stats << "t,drift,quadratic_variation\n";
    for(size_t t = 0; t < drift.size(); t++){
        stats << t << "," << drift[t] << "," << qv[t] << "\n";
    }
    stats.close();

    // File 3: terminal distribution for histogram / Gaussian overlay.
    write_terminal_csv(data, "terminal_distribution.csv");

    // File 4: empirical variance over time.
    write_variance_convergence_csv(data, "variance_convergence.csv");

    // File 5: drift study for different p.
    vector<double> ps = {0.40, 0.45, 0.50, 0.55, 0.60};
    write_parameter_study_csv(steps, 2000, ps, "parameter_study.csv");

    cout << "CSV files erstellt\n";
    return 0;
}


// g++ main.cpp -O2 -o main.exe
// .\main.exe
