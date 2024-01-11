#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Synapse.hpp"
#include "../engines/analytic/Analytic_engine.hpp"
#include "../engines/stochastic/Gillespie_engine.hpp"

using namespace std;

#define N_FORKS 1 // Note that it is (2^N_FORKS - 1)*3 compartments (if 2 synapses on each dend seg)!
void fork_dendrite(Dendritic_segment* ds, size_t depth=0) {
  if (depth < N_FORKS) {
    auto ds1 = new Dendritic_segment(*ds, ds->get_name() + "-1");
    new Synapse(*ds1, "s_" + ds1->get_name() + "_1");
    new Synapse(*ds1, "s_" + ds1->get_name() + "_2", .6, 6);
    fork_dendrite(ds1, depth+1);

    auto ds2 = new Dendritic_segment(*ds, ds->get_name() + "-2");
    new Synapse(*ds2, "s_" + ds2->get_name() + "_1");
    new Synapse(*ds2, "s_" + ds2->get_name() + "_2", .6, 6);
    fork_dendrite(ds2, depth+1);
  }
}

int main() {

  Soma soma("soma" /*,Parameters of the soma*/);
  
  ///// Linear neuron
  // Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1");
  // new Synapse(*p_ds, "s_1_1");
  // new Synapse(*p_ds, "s_1_2");
  
  // for(unsigned int i=0; i<100; ++i) {
  //   p_ds = new Dendritic_segment(*p_ds, "d_" + to_string(i+2));
  //   new Synapse(*p_ds, "s_" + to_string(i+2) + "_1");
  //   new Synapse(*p_ds, "s_" + to_string(i+2) + "_2");
  // }


  ///// Branching neuron
  Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1");
  new Synapse(*p_ds, "s_1_1", .6, 6);
  new Synapse(*p_ds, "s_1_2", .6, 6);

  fork_dendrite(p_ds);


  Neuron neuron(soma, "Test_neuron");
  
  cout << neuron << endl;


  // cout << "----------------- ANALYTIC ENGINE -----------------\n";
  Analytic_engine ae(neuron);

  arma::mat covariances(15,15);
  arma::vec expectations(15), variances(15);

  ofstream ofs_expectations("expectations"),
    ofs_covariances("covariances"),
    ofs_variances("variances"),
    ofs_correlations("correlations");

  double dt = .1;
  double t_max = 5e3;
  for(double t=0; t<t_max; t+=dt) {
    ofs_expectations << t << ',' << (expectations = ae.get_expectations()).t() << endl;

    ofs_covariances << "t=" << t << endl
                    << "covariances:\n" << (covariances = ae.get_covariances()) << endl;

    
    ofs_correlations << "t=" << t << std::endl;
    for(size_t i=0; i<15; ++i) {
      for(size_t j=0; j<15; ++j)
        ofs_correlations << covariances(i,j) - expectations(i)*expectations(j) << ',';
      ofs_correlations << std::endl;
    }        

    ofs_variances << t << ',';
    for(size_t i=0; i<15; ++i)
      ofs_variances << sqrt(covariances(i,i) - expectations(i)*expectations(i)) << ',';
    ofs_variances << endl;
    
    if(t==0)
      ae.nonstationary_covariances_direct_ODE_solver_step(dt, true);
    else
      ae.nonstationary_covariances_direct_ODE_solver_step(dt);

    ae.nonstationary_expectations(t, true, false);
      
    // if(t==0)
    //   ae.sem_nonstationary_expectations(t, true, false);
    // else
    //   ae.sem_nonstationary_expectations(t);
  }

  ae.stationary_expectations().stationary_covariances();

  

  // double t_max = 5e3;
  // double o1_dt = .01, o2_dt = o1_dt*10;
  
  // for(double t=0; t<t_max; t+=o1_dt) {
  //   // cout << "t=" << t << endl
  //   //      << "expextations: " << (expectations = ae.get_expectations()).t() << endl;
  //   ofs_expectations << t << ',' << (expectations = ae.get_expectations()) << endl;
  //   if(t==0)
  //     ae.sem_nonstationary_expectations_direct_ODE_solver_step(o1_dt, true);
  //   else
  //     ae.sem_nonstationary_expectations_direct_ODE_solver_step(o1_dt);


  //   if(t)
    
  //   ofs_covariances << "t=" << t << endl
  //                   << "covariances:\n" << (covariances = ae.get_covariances()) << endl;

  //   ofs_variances << t << ',';
  //   for(size_t i=0; i<15; ++i)
  //     ofs_variances << covariances(i,i) - expectations(i)*expectations(i) << ',';
  //   ofs_variances << endl;

  //   ae.sem_nonstationary_covariances_direct_ODE_solver_step(o2_dt);

  // }


  // double 
  // for(double t=0; t<t_max; t+=o2_dt) {
  //   // cout << "t=" << t << endl
  //   //      << "covariances:\n" << (covariances = ae.get_covariances()) << endl;
  //   // cout << "Variances:" << endl;

  //   ofs_covariances << "t=" << t << endl
  //        << "covariances:\n" << (covariances = ae.get_covariances()) << endl;

  //   ofs_variances << t << ',';
  //   for(size_t i=0; i<15; ++i)
  //     ofs_variances << covariances(i,i) - expectations(i)*expectations(i) << ',';
  //   ofs_variances << endl;

  //   ae.sem_nonstationary_covariances_direct_ODE_solver_step(o2_dt);
  // }

  ofs_expectations.close();
  ofs_covariances.close();
  ofs_variances.close();  

  
  
  // for(double t=0; t<10000; t+=dt) {
  //   cout << "t=" << t << endl
  //        << "expextations: " << (expectations = ae.get_expectations()).t() << endl
  //        << "covariances:\n" << (covariances = ae.get_covariances()) << endl;
  //   cout << "Variances:" << endl;
  //   for(size_t i=0; i<15; ++i)
  //     cout << covariances(i,i) - expectations(i)*expectations(i) << ',';
  //   cout << endl;

  //   if(t==0)
  //     ae.sem_nonstationary_expectations_direct_ODE_solver_step(dt, true);
  //   //      ae.sem_nonstationary_covariances_direct_ODE_solver_step(dt, true);
  //   else
  //     //      ae.sem_nonstationary_covariances_direct_ODE_solver_step(dt);
  //     ae.sem_nonstationary_expectations_direct_ODE_solver_step(dt);
  //   ae.sem_nonstationary_covariances_direct_ODE_solver_step(dt);
  // }
  
  return 0;
}
