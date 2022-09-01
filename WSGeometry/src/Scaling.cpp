#include "RcppArmadillo.h"
#include "divergence.h"
#include "timer.h"



//[[Rcpp::export]]
Rcpp::List StabilizedScaling_Rcpp(const arma::mat &cost_matrix,
                                  const arma::vec &supply,
                                  const arma::vec &demand,
                                  int supply_div_type,
                                  int demand_div_type,
                                  const arma::vec supply_div_parameters,
                                  const arma::vec demand_div_parameters,
                                  int iter_max,
                                  const arma::vec &epsvec,
                                  double tol,
                                  int check_interval = 1,
                                  double indicator_slack = 1e-6,
                                  double clamp = 1e100) {
    
    divergence_function<double>::indicator_slack = indicator_slack;
    
    timer time({"pre", "supply_workers", "demand_workers", "loop_rest", "post"});
    time.start_section(0);
    
    divergence_function<double> *div_supply = 0;
    divergence_function<double> *div_demand = 0;
    
    try {
        div_supply = make_divergence<double>((divergence_type)supply_div_type, supply_div_parameters);
        div_demand = make_divergence<double>((divergence_type)demand_div_type, demand_div_parameters);
    } catch (const std::invalid_argument &ex) {
        Rcpp::Rcout << "ERROR: " << ex.what() << std::endl;
        
        if (div_supply)
            delete div_supply;
        if (div_demand)
            delete div_demand;
        
        return Rcpp::List();
    }
    
    
    int epsind = 0;
    double epsi = epsvec[epsind];
    int last_epsi_change = 0;
    arma::ivec epsi_iterations(epsvec.size());
    
    // linear part of iterates
    arma::vec supply_linear(supply.size(), arma::fill::ones);
    arma::vec demand_linear(demand.size(), arma::fill::ones);
    
    arma::vec supply_linear_old(supply.size(), arma::fill::ones);
    arma::vec demand_linear_old(demand.size(), arma::fill::ones);
    
    // logarithmic part of iterates used to keep linear part low
    arma::vec supply_log(supply.size(), arma::fill::zeros);
    arma::vec demand_log(demand.size(), arma::fill::zeros);
    
    auto get_kernel = [&]()->arma::mat{
        double epsi_inv = 1. / epsi;
        arma::mat res(supply.size(), demand.size());
        for (int i = 0; i < (int)supply.size(); i++)
            for (int j = 0; j < (int)demand.size(); j++)
                res(i, j) = std::min(clamp, std::exp(epsi_inv * (supply_log[i] + demand_log[j] - cost_matrix(i, j))));
        // Rcpp::Rcerr << "kernel values " << res.min() << " to " << res.max() << std::endl;
        return res;
    };
    
    arma::mat kernel = get_kernel();
    arma::mat kernel_t = kernel.t();
    
    auto get_primal = [&]()->arma::mat{
        return get_kernel() % (supply_linear * demand_linear.t());
    };
    
    auto primal_cost = [&]()->double{
        double res = 0, epsi_inv = 1.0 / epsi, sd, c, t0, t1;
        
        for (int i = 0; i < (int)supply.size(); i++)
            for (int j = 0; j < (int)demand.size(); j++) {
                sd = supply_log[i] + demand_log[j];
                c = cost_matrix(i, j);
                t0 = std::exp((sd - c) * epsi_inv) * (sd - epsi) + epsi * std::exp(-c * epsi_inv);
                t1 = supply[i] * demand[j];
                res += (t0 && t1) ? t0 : 0;
            }
        
        auto pi = get_primal();
        arma::vec supply_marginal = arma::sum(pi, 1).as_col();
        arma::vec demand_marginal = arma::sum(pi, 0).as_col();
        res += div_supply->eval(supply_marginal, supply);
        res += div_demand->eval(demand_marginal, demand);
        
        return res;
    };
    auto dual_cost = [&]()->double{
        double res = 0, epsi_inv = 1.0 / epsi, sd, c, t0, t1;
        
        for (int i = 0; i < (int)supply.size(); i++)
            for (int j = 0; j < (int)demand.size(); j++) {
                sd = (supply_log[i] + demand_log[j]);
                c = cost_matrix(i, j);
                t0 = std::exp((sd - c) * epsi_inv) - std::exp(-c * epsi_inv);
                t1 = supply[i] * demand[j];
                res += (t0 && t1) ? epsi * t0 : 0;
            }
            
        res += div_supply->eval_conjugate(-supply_log, supply);
        res += div_demand->eval_conjugate(-demand_log, demand);
        
        return -res;
    };
    
    double scaling_limit = 1e30;
    bool rescale = 0, run = 1;
    double last_step_size = std::numeric_limits<double>::infinity();
    
    time.start_section(3);
    
    int k = 1;
    for (; k <= iter_max && run; k++) {
        
        try {
            Rcpp::checkUserInterrupt();
        } catch (...) {
            break; 
        }
        
        supply_linear_old = supply_linear;
        demand_linear_old = demand_linear;
        
        time.start_section(1);
        supply_linear = kernel * demand_linear;
        supply_linear = div_supply->eval_proxdiv(supply_linear, supply_log, supply, epsi);
        
        for (auto x = supply_linear.begin(); x < supply_linear.end(); x++)
            if (*x > clamp)
                *x = clamp;
        
        time.start_section(2);
        
        demand_linear = kernel_t * supply_linear;
        demand_linear = div_demand->eval_proxdiv(demand_linear, demand_log, demand, epsi);
        
        for (auto x = demand_linear.begin(); x < demand_linear.end(); x++)
            if (*x > clamp)
                *x = clamp;
        
        time.start_section(3);
        
        if (supply_linear.has_inf() || demand_linear.has_inf()) {
            Rcpp::Rcout << "WARNING: numeric overflow in single step!" << std::endl;
            break;
        }
        if (supply_linear.has_nan() || demand_linear.has_nan()) {
            Rcpp::Rcout << "WARNING: NaN in dual!" << std::endl;
            break;
        }
        
        if (k % check_interval == 0 || k == iter_max) {
            
            last_step_size = std::max(arma::max(arma::abs(supply_linear_old - supply_linear)),
                                      arma::max(arma::abs(demand_linear_old - demand_linear)));
            
            if (last_step_size < tol || 
                k / (double)iter_max > (epsind + 1) / (double)epsvec.size() ||
                k == iter_max) {
                
                epsi_iterations[epsind] = k - last_epsi_change;
                last_epsi_change = k;
                
                if (epsind == (int)epsvec.size() - 1) {
                    run = 0;
                } else {
                    epsi = epsvec[++epsind];
                }
                
                rescale = 1;
            }
        }
        
        if (arma::max(arma::abs(supply_linear)) > scaling_limit || 
            arma::max(arma::abs(demand_linear)) > scaling_limit ||
            rescale) {
            
            // Rcpp::Rcerr << "rescaling: " << k << " " << arma::max(arma::abs(supply_linear)) << " " << arma::max(arma::abs(demand_linear)) << std::endl;
            
            for (auto *slin = supply_linear.begin(), *slog = supply_log.begin(); slin != supply_linear.end(); slin++, slog++) {
                if (*slin) {
                    *slog += epsi * std::log(*slin);
                    *slin = 1;
                }
            }
            for (auto *dlin = demand_linear.begin(), *dlog = demand_log.begin(); dlin != demand_linear.end(); dlin++, dlog++) {
                if (*dlin) {
                    *dlog += epsi * std::log(*dlin);
                    *dlin = 1;
                }
            }
            
            kernel = get_kernel();
            kernel_t = kernel.t();
            
            rescale = 0;
        }
    }
    time.start_section(4);
    
    if (k > iter_max)
        k--;
    
    auto pi = get_primal();
    double pcost = primal_cost();
    double dcost = dual_cost();
    
    delete div_demand;
    delete div_supply;
    
    // Rcpp::Rcerr << time.format_times() << std::endl;
    
    return Rcpp::List::create(Rcpp::Named("transportPlan") = pi,
                              Rcpp::Named("primalCost") = pcost,
                              Rcpp::Named("dualCost") = dcost,
                              Rcpp::Named("supplyDual") = supply_log.t(),
                              Rcpp::Named("demandDual") = demand_log.t(),
                              Rcpp::Named("finalStepSize") = last_step_size, 
                              Rcpp::Named("iterations") = k,
                              Rcpp::Named("epsiIterations") = epsi_iterations.t());
}