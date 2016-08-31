#include "MSMS_Cluster_Sim.h"
#include "AbcSmc.h"
#include <time.h>
#include <stdlib.h>
#include <cmath>
#include <gsl/gsl_rng.h>


const int BURNIN               = 50;
const int FITTED_LENGTH        = 2014-1985;
const int RUN_LENGTH           = BURNIN + FITTED_LENGTH;
const int HETEROTYPIC_IMMUNITY = 7;
const double MU                = 1.0/79.0;

const gsl_rng* RNG = gsl_rng_alloc(gsl_rng_taus2);

const vector<int> POP_SIZES = {267051, 218948, 141004, 153813, 212494,
                               137046, 191325, 122308, 2147857, 229055,
                               127498, 398423, 191164, 117157, 193259,
                               277728, 139210, 807071, 453187};

double metric(vector<double> &incidence, int idx) {
    double val = 0;
    ABC::Map<ABC::Col> col(&incidence[0], incidence.size()); // copy data from vector into Col
    switch (idx){
        case 0: val = ABC::mean(col);                           break;
        case 1: val = quantile(incidence, 0.0);                 break;
        case 2: val = quantile(incidence, 0.25);                break;
        case 3: val = quantile(incidence, 0.5);                 break;
        case 4: val = quantile(incidence, 0.75);                break;
        case 5: val = quantile(incidence, 1.0);                 break;
        case 6: val = sqrt(ABC::variance(col, ABC::mean(col))); break;
        case 7: val = ABC::skewness(col);                       break;
        case 8: val = ABC::median_crossings(col);               break;
        default: cerr << "Unknown metric requested\n"; exit(-132);
    }
    return val;
}


void connect_network (Network* net) {
    vector<double> dist = {0, 3, 45, 160, 393, 715, 883, 989, 897, 752,
                          697, 755, 856, 1085, 1224, 1452, 1578, 1622, 1711, 1584,
                          1514, 1355, 1209, 990, 805, 597, 477, 353, 232, 177,
                          126, 90, 69, 54, 36, 29, 21, 16, 10, 5,
                          8, 13, 7, 9, 3, 1, 3, 3, 2, 5,
                          0, 1, 0, 0, 0, 1};

    dist = normalize_dist(dist, sum(dist));
    net->rand_connect_user(dist);
}

vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par*) {
    gsl_rng_set(RNG, rng_seed); // We're using two different RNGs for different things, unfortunately
    vector<double> metrics;
    const double r_zero               = (double) args[0];
    const double chi                  = (double) args[1];
    const vector<double> initial_exposed_per_100k = {pow(10,args[2]), pow(10,args[2])};
    const double cluster_jump         = (double) args[3];
    const double h                    = (double) args[4]; // reported fraction
    const double noise                = (double) args[5]; // per cap annual probability of non-flu ili
    const int per_cap                 = 100000;

    vector< vector<double> > case_mat(POP_SIZES.size());

    cerr << "serial: " << serial << endl;
    cerr << "pars:"; for (double par: args) cerr << " " << par; cerr << endl;
    for ( unsigned int i = 0; i < POP_SIZES.size(); ++i) {
        cerr << i << "\t:";
        const int N = POP_SIZES[i];
        const vector<int> initial_exposed = {(int) (0.5 + N*initial_exposed_per_100k[0]/per_cap), (int) (0.5 + N*initial_exposed_per_100k[1]/per_cap)};
        Network* net = new Network("EpiNet", Network::Undirected);
        net->get_rng()->seed(rng_seed);
        net->populate(N);
        connect_network(net);

        MSMS_Cluster_Sim* sim = new MSMS_Cluster_Sim(net, chi, r_zero, HETEROTYPIC_IMMUNITY, initial_exposed, RUN_LENGTH, MU, cluster_jump);
        sim->run_simulation();
        vector< vector<int> > final_sizes = sim->epidemic_sizes();
        vector<vector<int> > (final_sizes.begin()+BURNIN, final_sizes.end()).swap(final_sizes); // discard burn-in data
        for (auto v: final_sizes) { // {H1, H3} pair for each year
            const int fs = (int) (0.5 + per_cap*((float) v[0] + v[1])/N); // H1+H3 final size per 100k
            cerr << " " << fs;
            assert(fs >= 0);
            const double reported_flu = fs == 0 ? 0 : gsl_ran_binomial(RNG,h,fs);
            const double reported_ili = reported_flu + gsl_ran_binomial(RNG,noise,per_cap); 
            case_mat[i].push_back(reported_ili);
//            cout << v[0] << " + " << v[1] << " = " << ((float) v[0] + v[1])/N << endl;
            //metrics.push_back(v[0]);
        }
        cerr << endl;

        delete sim;
        delete net;
    }

    for (unsigned int met_idx = 0; met_idx < 9; ++met_idx) {
        for (vector<double> incidence: case_mat) {
            metrics.push_back(metric(incidence, met_idx));
        }
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

