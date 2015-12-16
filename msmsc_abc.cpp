#include "MSMS_Cluster_Sim.h"
#include "AbcSmc.h"
#include <time.h>
#include <stdlib.h>
#include <cmath>
#include <gsl/gsl_rng.h>


const vector<int> POP_SIZES    = {(int) 1e5};
const int BURNIN               = 50;
const int FITTED_LENGTH        = 2016-1984;
const int RUN_LENGTH           = BURNIN + FITTED_LENGTH;
const int HETEROTYPIC_IMMUNITY = 7;
const double MU                = 0.02;

const gsl_rng* RNG = gsl_rng_alloc(gsl_rng_taus2);

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

vector<long double> simulator(vector<long double> args, const unsigned long int rng_seed, const ABC::MPI_par*) {
    vector<long double> metrics;
    const double r_zero               = (double) args[0];
    const double chi                  = (double) args[1];
    const vector<int> initial_exposed = {(int) args[2], (int) args[2]};
    const double cluster_jump         = (double) args[3];

    for ( int N: POP_SIZES) {
        Network* net = new Network("EpiNet", Network::Undirected);
        net->get_rng()->seed(rng_seed);
        net->populate(N);
        connect_network(net);

        MSMS_Cluster_Sim* sim = new MSMS_Cluster_Sim(net, chi, r_zero, HETEROTYPIC_IMMUNITY, initial_exposed, RUN_LENGTH, MU, cluster_jump);
        sim->run_simulation();
        vector< vector<int> > final_sizes = sim->epidemic_sizes();
        vector<vector<int> > (final_sizes.begin()+BURNIN, final_sizes.end()).swap(final_sizes); // discard burn-in data
        for (auto v: final_sizes) {
//            cout << v[0] << " + " << v[1] << " = " << ((float) v[0] + v[1])/N << endl;
            metrics.push_back(v[0]);
        }
        for (auto v: final_sizes) metrics.push_back(v[1]);

        delete sim;
        delete net;
    }
    return metrics;
}


void usage() {
    cerr << "\n\tUsage: ./abc_sql abc_config_sql.json --process\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate -n <number of simulations per database write>\n\n";
}


int main(int argc, char* argv[]) {

    if (not (argc == 3 or argc == 5 or argc == 6) ) {
        usage();
        exit(100);
    }

    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = -1;

    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) {
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            buffer_size = atoi(argv[++i]);
        } else {
            usage();
            exit(101);
        }
    }

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (process_db) {
        gsl_rng_set(RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id
        abc->process_database(RNG);
    }

    if (simulate_db) {
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }

    return 0;
}

