#include "MSMS_Cluster_Sim.h"
#include <time.h>
#include <stdlib.h>
#include <cmath>

void connect_network (Network* net) {
    vector<double> dist;
    double deg_array[] = {0, 3, 45, 160, 393, 715, 883, 989, 897, 752,
                          697, 755, 856, 1085, 1224, 1452, 1578, 1622, 1711, 1584,
                          1514, 1355, 1209, 990, 805, 597, 477, 353, 232, 177,
                          126, 90, 69, 54, 36, 29, 21, 16, 10, 5,
                          8, 13, 7, 9, 3, 1, 3, 3, 2, 5,
                          0, 1, 0, 0, 0, 1};
    dist.assign(deg_array,deg_array+56);

    dist = normalize_dist(dist, sum(dist));
    net->rand_connect_user(dist);
}


int main(int argc, char* argv[]) {
    Network* net = new Network("EpiNet", Network::Undirected);
    const int N = 1e5;
    net->populate(N);
    //net->erdos_renyi(9);
    connect_network(net);

    const double CHI = 0.7;
    const int j_max = 1;
    const int num_years = 50;
    const double r_zero = 5;
    const int heterotypic_immunity = 7;
    const vector<int> initial_exposed = {5, 5};
    const double mu = 0.01333;
    const double cluster_jump = 0.2;

    MSMS_Cluster_Sim* sim = new MSMS_Cluster_Sim(net, CHI, r_zero, heterotypic_immunity, initial_exposed, num_years, mu, cluster_jump);
    for ( int j = 0; j < j_max; j++) {
        cerr << "Repetition: " << j << endl;
        sim->run_simulation();
        vector< vector<int> > final_sizes = sim->epidemic_sizes();
        for (auto v: final_sizes) {
            cout << v[0] << " + " << v[1] << " = " << ((float) v[0] + v[1])/N << endl;
        }
        //sim->reset();
    }
    cerr << "derp\n";

    delete sim;
    delete net;
}
