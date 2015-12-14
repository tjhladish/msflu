#ifndef MSMS_CLUSTER_SIM_H
#define MSMS_CLUSTER_SIM_H
#include "Simulator.h"
#include <unordered_set>

//#include <assert.h>
//#include <queue>

const int NEVER = INT_MAX;

enum StrainType { H1N1,
                  H3N2,
                  NumStrains }; // must be last


struct Infection { // strain type is the container index, i.e. position in InfectionHistory
    Infection() { clear(); }
    //Infection(int c, int y, int t): cluster(c), year(y), timestep(t) {}

    Infection* update(int c, int y, int t) { cluster = c; year = y; timestep = t; return this; }
    void clear() { cluster = -1; year = NEVER; timestep = NEVER; }
    int cluster;
    int year;
    int timestep;
};


struct InfectionHistory { // needs to persist across years
    vector<Infection*> infection_history;
    Infection* last;

    InfectionHistory() { 
        for (int i = 0; i<(int) NumStrains; ++i) {
            infection_history.push_back(new Infection());
        }
        last = nullptr;
    }

    ~InfectionHistory() { 
        for (auto& ih: infection_history) { 
            delete ih;
        }
        last = nullptr;
    }

    Infection* infect( const StrainType s, const int c, const int y, const int t ) { last = infection_history[(int) s]; return last->update(c, y, t); }
    Infection* strain( const StrainType s ) { return infection_history[(int) s]; }
    void clear_history() { for (auto& ih: infection_history) ih->clear(); }
};


class MSMS_Cluster_Sim: public Simulator {
    int HeterotypicImmunityDuration;         // heterotypic immunity lasts 7 infectious periods, (i.e. ~3 weeks)
    double CI;                               // cross immunity between homotypic clusters
    double Mu;                               // annual birth/death probability
    double ClusterJump;                      // probability that a cluster jump will occur
    vector<int> Clusters;                    // current cluster number for each strain
    vector<int> InitialNumExposed;           // how many nodes for each strain to expose to spark epidemic
    vector<InfectionHistory*> NodeHist;      // infection history for all nodes
    vector<Node*> Infected;                  // list of infected nodes, so we don't have to check all nodes when transmitting
    double T;                                // naive transmission probability
    vector<Node*> Recovered;                 // we don't need this for the simulation per se, but we do need it
                                             // in order to quickly reset the network to completely susceptible
    vector< vector<int> > EpidemicSizes;     // annual time series of epidemic sizes, by strain
    int Year;
    int NumYears;
    
    public:
        MSMS_Cluster_Sim():Simulator(), CI(0), Clusters({0,0}), T(1.0), HeterotypicImmunityDuration(7), NumYears(1), ClusterJump(1.0) {};
        MSMS_Cluster_Sim(Network* net, double cluster_crossimmunity, double r_zero, int hid, vector<int> initial_exposed, int num_years, double mu, double cj)
            :CI(cluster_crossimmunity), 
             Clusters({0,0}), 
             HeterotypicImmunityDuration(hid), 
             InitialNumExposed(initial_exposed),
             NumYears(num_years),
             Mu(mu),
             ClusterJump(cj) {
                set_network(net);
                calc_naive_transmissibility(r_zero);
            }

        ~MSMS_Cluster_Sim() {
            clear_node_histories();
        }

        void initialize_node_histories() {
            assert(net->size() > 0); 
            NodeHist = vector<InfectionHistory*>(net->size());
            for (auto &e: NodeHist) e = new InfectionHistory();
        }
        void clear_node_histories() {
            cerr << "node hist size: " << NodeHist.size() << endl;
            for(auto e: NodeHist) {
                delete e;
            }
            NodeHist.clear();
        }
        
        void set_network( Network* net ) { 
            Simulator::set_network( net );
            initialize_node_histories();
        };

        void set_cluster(StrainType s, int c) { Clusters[s] = c; }
        void increment_cluster(StrainType s) { ++Clusters[s]; }

        double calc_naive_transmissibility(double R_zero) { 
            const double temp_T = R_zero * calc_critical_transmissibility(); 
            this->T = temp_T < 0 ? 0 :
                      temp_T > 1 ? 1 :
                      temp_T;
            cerr << "T: " << T << endl; return T;
        }
//double calc_naive_transmissibility(double R_zero) { this->T = 0.4; cerr << "T: " << T << endl; return T; }
        double get_naive_transmissibility() { return T; }

        double get_transmissibility(Node* end, const int year, const int timestep, const StrainType strain) {
            return get_naive_transmissibility() * get_susceptibility(end, year, timestep, strain); // Transmission probability given immunity
        }


        Infection* infect(Node* node, StrainType strain, const int cluster, const int year, const int timestep) {
//            if (node->get_state() != -1) cerr << "double infection\n";
            node->set_state(strain);
            NodeHist[node->get_id()]->infect(strain, cluster, year, timestep);
        }


        inline double get_susceptibility(Node* node, const int year, const int timestep, const StrainType challenging_strain) { 
            const int node_id = node->get_id();
            assert(node_id < NodeHist.size());
            Infection* last_infection = NodeHist[node_id]->last;

            double susceptibility = 1.0; // default is no protection / fully susceptible

            if (last_infection) {
                // node has already been infected during this epidemic year
                if (year == last_infection->year and timestep - last_infection->timestep <= HeterotypicImmunityDuration) {
                    // perfect short-term heterotypic immunity
                    susceptibility = 0.0;
                } else {
                    // no recent infection
                    Infection* last_homotypic_infection = NodeHist[node_id]->strain(challenging_strain);
                    if (last_homotypic_infection->year != NEVER) {
                        // node has been infected with this strain; check clusters
                        const int delta = Clusters[(int) challenging_strain] - last_homotypic_infection->cluster;
                        susceptibility = 1.0 - pow(CI, delta);
                    }
                }
            }
            return susceptibility;
        }


        vector<Node*> rand_expose (vector<int> initial_exposures_by_strain, const int year, const int timestep) {
            vector<StrainType> exposure_types; 
            int total_exposures = 0;
            for ( unsigned int s = 0; s < initial_exposures_by_strain.size(); ++s ) {
                total_exposures += initial_exposures_by_strain[s];
                exposure_types.resize(total_exposures, (StrainType) s);
            }
            shuffle(exposure_types, mtrand); // randomize exposures so there is no player 1 advantage

            vector<Node*> patients_zero;
            vector<Node*> exposed_nodes = rand_choose_nodes(total_exposures);
            for ( unsigned int i = 0; i < exposed_nodes.size(); ++i) {
                const StrainType strain = exposure_types[i];
                Node* node = exposed_nodes[i];
                //const int node_id = exposed_nodes[i]->get_id();
                const double susceptibility = get_susceptibility(node, year, timestep, strain); 
                if (mtrand->rand() < susceptibility) {
                    infect(node, strain, Clusters[strain], year, timestep);
                    Infected.push_back(node);
                    patients_zero.push_back(node);
                }
            };
            return patients_zero;
        }


        vector<Node*> rand_infect (vector<int>) {
            cerr << "rand_expose() should be used instead of rand_infect()\n";
            assert(false);
        }


        void age_network(const double mu) { // birth rate == death rate
            // should be called just before an epidemic is simulated, after resetting (if epidemics have already been simualted)
            assert(Infected.size() == 0);
            assert(Recovered.size() == 0);
            assert(time == 0);

            const int num_affected_nodes = rand_binomial(net->size(), mu, mtrand);
            assert(num_affected_nodes != -1); // indicates bad parameters were passed to rand_binomial
            vector<Node*> phoenix_nodes = rand_choose_nodes(num_affected_nodes);
            for (auto& n: phoenix_nodes) {
                NodeHist[n->get_id()]->clear_history();
            }
        }


        void step_simulation () {
            time++;
            assert(Infected.size() > 0);
            vector<Node*> new_infected;
shuffle(Infected, mtrand);
            for (Node* inode: Infected) {
                StrainType strain = (StrainType) inode->get_state();
                vector<Node*> neighbors = inode->get_neighbors();
                for (Node* test: neighbors) {
                    if ( mtrand->rand() < get_transmissibility(test, Year, time, strain) ) {
                        infect(test, strain, Clusters[strain], Year, time);
                        new_infected.push_back( test );
                    }
                }
                Recovered.push_back(inode);
            }
            Infected = new_infected;
        }


        vector<int> final_size_by_strain() {
            vector< unordered_set<Node*> > recovered_sets(NumStrains);
cerr << "tallying recovered: " << Recovered.size() << endl;
            for (Node* node: Recovered) { 
                InfectionHistory* ih = NodeHist[node->get_id()];
                for (int s = 0; s < NumStrains; ++s) {
                    if (ih->strain((StrainType) s)->year == Year) recovered_sets[s].insert(node);
                }
            } 
            vector<int> tally(NumStrains, 0);
            for (int s = 0; s < NumStrains; ++s) tally[s] = recovered_sets[s].size();
            return tally;
        }

        int epidemic_size() { return 0; } // dummy implementation, required by base class
        vector< vector<int> > epidemic_sizes() { return EpidemicSizes; } // dummy implementation, required by base class

        void calculate_average_transmissibility(vector<double> &average_tk, double &average_t, StrainType strain) {
            average_t = 0;
            vector<int> deg_dist = net->get_deg_dist();
            vector<double> tk( deg_dist.size() );    // total transmissibility by degree -- used to calc mean
		    vector<Node*> nodes = net->get_nodes();

            for (int i = 0; i < (signed) nodes.size(); i++) {
                int deg = nodes[i]->deg();
                tk[deg] += get_transmissibility(nodes[i], Year, time, strain);
            }
            
            average_tk.resize( tk.size() ); // average transmissibility of degree k nodes
            int deg_sum = sum(net->get_deg_series()); // sum of all degrees in the network, i.e. the total number of stubs
            for (int deg = 0; deg < (signed) tk.size(); deg++) {
                average_tk[deg] = tk[deg]/ ((double) deg_dist[deg]);
                average_t += tk[deg] * deg / ((double) deg_sum);
            }
        }


        void simulate_epidemic() { // simulate a single multi-strain epidemic
            while (Infected.size() > 0) {
                //cerr << ".";
                step_simulation();
            }
        }


        void run_simulation() {
            EpidemicSizes = vector< vector<int> > (NumYears, vector<int>(NumStrains, 0));
            for (int y = 0; y<NumYears; ++y) {
                Year = y;
                reset_time();
                set_these_nodes_to_state(net->get_nodes(), -1); // state indicates strain of current infection

                if (y > 0) {
                    age_network(Mu);
                    for (int s = 0; s < (int) NumStrains; ++s) {
                        if (mtrand->rand() < ClusterJump) increment_cluster((StrainType) s); 
                    }     
                }

                vector<Node*> p_zeroes = rand_expose(InitialNumExposed, y, time);
                cerr << "initally infected: " << p_zeroes.size() << endl;
                simulate_epidemic();
                EpidemicSizes[y] = final_size_by_strain();
                Recovered.clear();
            }
        }

        void reset() {
            reset_time();
            
            set_these_nodes_to_state(net->get_nodes(), -1);
            clear_node_histories();
            Infected.clear();
            Recovered.clear();
        }
        
};

#endif
