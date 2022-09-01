#include "RcppArmadillo.h"
#include "divergence.h"
#include <thread>
#include "timer.h"
#include "threading.h"

class sinkhorn_work : public work
{
protected:
    int istart, iend;
    bool transpose;
    arma::mat c;
    arma::vec log_measure;
    const double &epsi;
    arma::vec &dual_in;
    arma::subview_col<double> dual_out;
    divergence_function<double> *div;
    
public:
    sinkhorn_work(int istart, int iend, 
                    bool transpose,
                    const arma::mat &cost_matrix,
                    const arma::vec &log_measure,
                    const double &epsi,
                    arma::vec &dual_in,
                    arma::vec &dual_out,
                    divergence_function<double> *div) : 
    istart(istart), iend(iend), transpose(transpose),
    c(transpose ? arma::mat(cost_matrix.cols(istart, iend - 1).t()) : cost_matrix.rows(istart, iend - 1)),
    log_measure(log_measure), epsi(epsi),
    dual_in(dual_in), dual_out(dual_out.subvec(istart, iend - 1)),
    div(div) { }
    
    ~sinkhorn_work() { }
    
    void execute() override {
        double epsi_inv = 1. / epsi;
        
        // make copy of other dual to prevent false sharing
        arma::vec d(dual_in);
        
        for (int i = 0; i < iend - istart; i++) {
            
            // extract row of c to use (just for simplicity)
            auto cc = c.row(i);
            
            // calculate MSE
            double p = -1e100, q = 0, a;
            for (int j = 0; j < (int)d.size(); j++) {
                a = (d[j] - cc[j]) * epsi_inv + log_measure[j];
                
                if (a > p) {
                    q = q * std::exp(p - a) + 1;
                    p = a;
                } else {
                    q += std::exp(a - p);
                }
            }
            
            // calculate aprox and set output
            dual_out[i] = -div->eval_aprox(epsi * (p + std::log(q)), epsi);
            // dual_out[i] = -epsi * (p + std::log(q)); // exact marginals for debugging
        }
    }
    
};


// divergences:
// KL: type = 1, parameters = (lambda,)
// TV: type = 2, parameters = (lambda,) or (lambda_neg, lambda_pos)
// RG: type = 3, parameters = (lower, upper)
// Power: type = 4, parameters = (lambda, p), special case Berg: p = 0, KL: p = 1
// LSR: type = 5, parameters = (lower, upper, lambda_neg, lambda_pos)

//[[Rcpp::export]]
Rcpp::List Sinkhorn_Rcpp(const arma::mat &cost_matrix,
                         const arma::vec &supply,
                         const arma::vec &demand,
                         int supply_div_type,
                         int demand_div_type,
                         const arma::vec supply_div_parameters,
                         const arma::vec demand_div_parameters,
                         int iter_max,
                         const arma::vec &epsvec,
                         double tol,
                         int thread_cnt = -1,
                         int max_lines_per_work = 5,
                         double indicator_slack = 1e-6) {
    
    divergence_function<double>::indicator_slack = indicator_slack;
    
    timer time({"pre", "supply_workers", "demand_workers", "loop_rest", "post"});
    time.start_section(0);
    
    
    int min_lines_to_multithread = 30;
    if (thread_cnt < 1) {
        if ((int)supply.size() >= min_lines_to_multithread || 
            (int)demand.size() >= min_lines_to_multithread)
            thread_cnt = (int)std::thread::hardware_concurrency();
        else
            thread_cnt = 1;
    }
    
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
    
    arma::vec supply_log = arma::log(supply);
    arma::vec demand_log = arma::log(demand);
    
    arma::vec supply_dual(supply.size(), arma::fill::zeros);
    arma::vec demand_dual(demand.size(), arma::fill::zeros);
    arma::vec supply_dual_old(supply.size(), arma::fill::zeros);
    arma::vec demand_dual_old(supply.size(), arma::fill::zeros);
    
    auto primal = [&]()->arma::mat{
        arma::mat pi = -cost_matrix;
        pi.each_row() += demand_dual.t();
        pi.each_col() += supply_dual;
        pi = arma::exp(pi / epsi);
        pi.each_row() %= demand.t();
        pi.each_col() %= supply;
        return pi;
    };
    
    auto primal_cost = [&]()->double{
        double res = 0, epsi_inv = 1.0 / epsi, sd, c, t0, t1;
        
        for (int i = 0; i < (int)supply.size(); i++)
            for (int j = 0; j < (int)demand.size(); j++) {
                sd = supply_dual[i] + demand_dual[j];
                c = cost_matrix(i, j);
                t0 = std::exp((sd - c) * epsi_inv) * (sd - epsi) + epsi * std::exp(-c * epsi_inv);
                t1 = supply[i] * demand[j];
                res += (t0 && t1) ? t0 * t1 : 0;
            }
        
        auto pi = primal();
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
                sd = (supply_dual[i] + demand_dual[j]);
                c = cost_matrix(i, j);
                t0 = std::exp((sd - c) * epsi_inv) - std::exp(-c * epsi_inv);
                t1 = supply[i] * demand[j];
                res += (t0 && t1) ? epsi * t0 * t1 : 0;
            }
            
        res += div_supply->eval_conjugate(-supply_dual, supply);
        res += div_demand->eval_conjugate(-demand_dual, demand);
        
        return -res;
    };
    
    double last_step_size = std::numeric_limits<double>::infinity();
    
    thread_pool pool(thread_cnt);
    
    // create work
    int supply_work_cnt = std::max(1, std::min((int)supply.size(), (int)(supply.size() + max_lines_per_work - 1) / max_lines_per_work));
    int demand_work_cnt = std::max(1, std::min((int)demand.size(), (int)(demand.size() + max_lines_per_work - 1) / max_lines_per_work));
    std::vector<work*> supply_work, demand_work;
    for (int i = 0; i < supply_work_cnt; i++) {
        int istart = supply.size() * i / supply_work_cnt;
        int iend = supply.size() * (i + 1) / supply_work_cnt;
        supply_work.push_back(
            new sinkhorn_work(istart, iend, false, cost_matrix, demand_log,
                              epsi, demand_dual, supply_dual, div_supply)
        );
    }
    for (int i = 0; i < demand_work_cnt; i++) {
        int istart = demand.size() * i / demand_work_cnt;
        int iend = demand.size() * (i + 1) / demand_work_cnt;
        demand_work.push_back(
            new sinkhorn_work(istart, iend, true, cost_matrix, supply_log,
                              epsi, supply_dual, demand_dual, div_demand)
        );
    }
    
    time.start_section(3);
    int k = 1;
    for (; k <= iter_max; k++) {
        
        try {
            Rcpp::checkUserInterrupt();
        } catch (...) {
            break; 
        }
        
        supply_dual_old = supply_dual;
        demand_dual_old = demand_dual;
        
        time.start_section(1);
        pool.do_work(supply_work);
        time.start_section(2);
        pool.do_work(demand_work);
        time.start_section(3);
        
        last_step_size = std::max(arma::max(arma::abs(supply_dual_old - supply_dual)),
                                  arma::max(arma::abs(demand_dual_old - demand_dual)));
        
        if (last_step_size < tol || 
            k / (double)iter_max > (epsind + 1) / (double)epsvec.size() ||
            k == iter_max) {
            
            epsi_iterations[epsind] = k - last_epsi_change;
            last_epsi_change = k;
            
            if (epsind == (int)epsvec.size() - 1) {
                break;
            } else {
                epsi = epsvec[++epsind];
            }
        }
    }
    
    time.start_section(4);
    
    double pcost = primal_cost();
    double dcost = dual_cost();
    arma::mat pi = primal();
    
    delete div_demand;
    delete div_supply;
    for (auto x : supply_work)
        delete x;
    for (auto x : demand_work)
        delete x;
    
    time.end_section();
    // Rcpp::Rcout << time.format_times() << std::endl;
    
    return Rcpp::List::create(Rcpp::Named("transportPlan") = pi,
                              Rcpp::Named("primalCost") = pcost,
                              Rcpp::Named("dualCost") = dcost,
                              Rcpp::Named("supplyDual") = supply_dual.t(),
                              Rcpp::Named("demandDual") = demand_dual.t(),
                              Rcpp::Named("finalStepSize") = last_step_size, 
                              Rcpp::Named("iterations") = k,
                              Rcpp::Named("epsiIterations") = epsi_iterations.t());
}