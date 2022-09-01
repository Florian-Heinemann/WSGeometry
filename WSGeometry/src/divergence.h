// 
// This file contains various definitions and implementations to work with divergence functions 
// as well as some helping utility functions required.
// 

#ifndef DIVERGENCE_H
#define DIVERGENCE_H

#include "RcppArmadillo.h"

enum divergence_type {
  KL = 1,
  TV = 2,
  RG = 3,
  Power = 4,
  LSR = 5
};

// lambert W-function (defined by W(x) exp(W(x)) = x)
template<typename T>
arma::Mat<T> lambertW(const arma::Mat<T> &x);
template<typename T>
T lambertW(T x);

// W(exp(x))
template<typename T>
arma::Mat<T> lambertWexp(const arma::Mat<T> &x);
template<typename T>
T lambertWexp(T x);


template<typename T>
arma::Mat<T> zeros_like(const arma::Mat<T> &x) {
    arma::Mat<T> res(x);
    res.fill(0);
    return res;
}

// a class for entropy-divergence functions
template<typename T>
class divergence_function 
{
public:
    virtual ~divergence_function() {}
    
    // allowed violation of indicator functions
    static T indicator_slack;
    
    // evaluate divergece D(a | b)
    virtual T eval(const arma::Mat<T> &a, const arma::Mat<T> &b) = 0;
    
    // evaluate convex conjugate D*(a | b)
    virtual T eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) = 0;
    
    // evaluate aprox^epsi_phi*(x) = argmin_y epsi * exp((x - y) / epsi) + phi*(x) where phi is the entropy function
    virtual T eval_aprox(const T &x, const T &epsi) = 0;
    virtual arma::Mat<T> eval_aprox(const arma::Mat<T> &x, const T &epsi) = 0;
    
    // evaluate proxdiv^epsi_D
    virtual arma::Mat<T> eval_proxdiv(const arma::Mat<T> &s, // main argument
                                      const arma::Mat<T> &u, // log-stabilization factor
                                      const arma::Mat<T> &b, // reference measure
                                      const T &epsi) = 0;
};

template<typename T>
T divergence_function<T>::indicator_slack = 1e-6;

// create a divergence_function from given type and parameters; 
// returns a pointer with ownership transferred to the caller: the caller should call delete once it is done
template<typename T> 
divergence_function<T>* make_divergence(divergence_type type, arma::vec parameters);

// Kullback-Leibler divergence
template<typename T>
class KL_divergence : public divergence_function<T>
{
protected:
    T lambda;
public:
    KL_divergence(T lambda);
    ~KL_divergence() override {}
    
    T eval(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_aprox(const T &x, const T &epsi) override;
    arma::Mat<T> eval_aprox(const arma::Mat<T> &x, const T &epsi) override;
    arma::Mat<T> eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) override;
};

// (Asymmetric) total variation divergence
template<typename T>
class TV_divergence : public divergence_function<T>
{
protected:
    // slope for > 1 and < 1 densities
    T lambda_pos, lambda_neg;
public:
    TV_divergence(T lambda);
    TV_divergence(T lambda_pos, T lambda_neg);
    ~TV_divergence() override {}
    
    T eval(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_aprox(const T &x, const T &epsi) override;
    arma::Mat<T> eval_aprox(const arma::Mat<T> &x, const T &epsi) override;
    arma::Mat<T> eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) override;
};

// Range divergence
template<typename T>
class RG_divergence : public divergence_function<T>
{
protected:
    // lower and upper bound of range
    T lower, upper;
public:
    RG_divergence(T lower, T upper);
    ~RG_divergence() override {}
    
    T eval(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_aprox(const T &x, const T &epsi) override;
    arma::Mat<T> eval_aprox(const arma::Mat<T> &x, const T &epsi) override;
    arma::Mat<T> eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) override;
};

// Linear Soft Range divergence (combination of RG and TV)
template<typename T>
class LSR_divergence : public divergence_function<T>
{
protected:
    // lower and upper bound of range
    T lower, upper;
    // positive and negative slope
    T lambda_pos, lambda_neg;
public:
    LSR_divergence(T lower, T upper, T lambda_neg, T lambda_pos);
    ~LSR_divergence() override {}
    
    T eval(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_aprox(const T &x, const T &epsi) override;
    arma::Mat<T> eval_aprox(const arma::Mat<T> &x, const T &epsi) override;
    arma::Mat<T> eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) override;
};

// Berg divergence
template<typename T>
class berg_divergence : public divergence_function<T>
{
protected:
    // strength
    T lambda;
public:
    berg_divergence(T lambda);
    ~berg_divergence() override {}
    
    T eval(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_aprox(const T &x, const T &epsi) override;
    arma::Mat<T> eval_aprox(const arma::Mat<T> &x, const T &epsi) override;
    arma::Mat<T> eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) override;
};

// Power divergence (for p != 0, 1)
template<typename T>
class power_divergence : public divergence_function<T>
{
protected:
    // strength and power (!= 0, 1)
    T lambda, p;
public:
    power_divergence(T lambda, T p);
    ~power_divergence() override {}

    T eval(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) override;
    T eval_aprox(const T &x, const T &epsi) override;
    arma::Mat<T> eval_aprox(const arma::Mat<T> &x, const T &epsi) override;
    arma::Mat<T> eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) override;
};



// lambert W-function
template<typename T>
arma::Mat<T> lambertW(const arma::Mat<T> &x) {
    
    // find initial approximation (values optimized empirically; small values approximated by sqrt, medium by rational, large by log - loglog)
    const static T p[] = {0, 15,
                          1.79611733, 0.99924628, 0.36787942,
                          485.22388221, 9520.53522247, 204.35198655, 13085.94588194, 5368.86916901,
                          0.31044538};
    
    auto init_f = [](const T &x)->T{ 
        if (x < p[0])
            return p[2] * std::sqrt(std::max(0.0, x + p[4])) - p[3];
        if (x < p[1])
            return (p[5] + x * p[6] + x * x * p[7]) / (p[8] + x * p[9]);
        else {
            T l = std::log(x);
            return l - std::log(l) + p[10];
        }
    };
    
    auto step = [](T &y, const T &x)->void{
        T e = std::exp(y);
        T t = (y * e - x);
        y -= t / (e * (1 + y) - t * (2 + y) / (2 + 2 * y));
    };
    
    // initialize result matrix
    arma::Mat<T> y(x);
    y.transform(init_f);
    
    // apply Halley's method (3 steps are sufficient to reach double precision)
    for (auto xi = x.begin(), yi = y.begin(); xi != x.end(); xi++, yi++) {
        for (int i = 0; i < 3; i++)
            step(*yi, *xi);
    }
    
    return y;
}

template<typename T>
T lambertW(T x) {
    // find initial approximation (values optimized empirically; small values approximated by sqrt, medium by rational, large by log - loglog)
    const static T p[] = {0, 15,
                          1.79611733, 0.99924628, 0.36787942,
                          485.22388221, 9520.53522247, 204.35198655, 13085.94588194, 5368.86916901,
                          0.31044538};
    
    auto init_f = [](const T &x)->T{ 
        if (x < p[0])
            return p[2] * std::sqrt(std::max(0.0, x + p[4])) - p[3];
        if (x < p[1])
            return (p[5] + x * p[6] + x * x * p[7]) / (p[8] + x * p[9]);
        else {
            T l = std::log(x);
            return l - std::log(l) + p[10];
        }
    };
    
    auto step = [](T &y, const T &x)->void{
        T e = std::exp(y);
        T t = (y * e - x);
        y -= t / (e * (1 + y) - t * (2 + y) / (2 + 2 * y));
    };
    
    // initialize result matrix
    T y = init_f(x);
    
    // apply Halley's method (3 steps are sufficient to reach double precision)
    for (int i = 0; i < 3; i++)
        step(y, x);
    
    return y;
}

template<typename T>
arma::Mat<T> lambertWexp(const arma::Mat<T> &x) {
    
    // find initial approximation (values optimized empirically; small values approximated by sqrt, medium by rational, large by log - loglog)
    const static T p[] = {-4, 8,
                          0.66887365, 0.35780384, 0.04420618,
                          0.1316505};
    
    auto init_f = [](const T &x)->T{ 
        if (x < p[0])
            return 0;
        if (x < p[1])
            return p[2] + (p[3] + p[4] * x) * x;
        else {
            return x - std::log(x) + p[5];
        }
    };
    
    auto step = [](T &y, const T &x)->void{
        T e = y - std::exp(x - y);
        y -= e / (1 + y - e * (2 + y) / (2 + 2 * y));
    };
    
    // initialize result matrix
    arma::Mat<T> y(x);
    y.transform(init_f);
    
    // apply Halley's method (3 steps are sufficient to reach double precision)
    for (auto xi = x.begin(), yi = y.begin(); xi != x.end(); xi++, yi++) {
        for (int i = 0; i < 3; i++)
            step(*yi, *xi);
    }
    
    return y;
}

template<typename T>
T lambertWexp(T x) {
    // find initial approximation (values optimized empirically; small values approximated by sqrt, medium by rational, large by log - loglog)
    const static T p[] = {-4, 8,
                          0.66887365, 0.35780384, 0.04420618,
                          0.1316505};
    
    auto init_f = [](const T &x)->T{ 
        if (x < p[0])
            return 0;
        if (x < p[1])
            return p[2] + (p[3] + p[4] * x) * x;
        else {
            return x - std::log(x) + p[5];
        }
    };
    
    auto step = [](T &y, const T &x)->void{
        T e = y - std::exp(x - y);
        y -= e / (1 + y - e * (2 + y) / (2 + 2 * y));
    };
    
    // initialize result matrix
    T y = init_f(x);
    
    // apply Halley's method (3 steps are sufficient to reach double precision)
    for (int i = 0; i < 3; i++)
        step(y, x);
    
    return y;
}

template<typename T> 
divergence_function<T>* make_divergence(divergence_type type, arma::vec parameters) {
    divergence_function<T> *res = 0;
    
    auto param_len_error = [&](std::string type, std::string correct_num)->void{
        std::stringstream ss;
        ss << type << " divergence needs " << correct_num << " parameter(s) but got " << parameters.size();
        throw new std::invalid_argument(ss.str());
    };
    
    switch (type) {
    case KL: {
        if (parameters.size() == 1) {
            res = new KL_divergence<double>(parameters[0]);
        } else {
            param_len_error("KL", "1");
        }
    } break;
    case TV: {
        if (parameters.size() == 1) {
            res = new TV_divergence<double>(parameters[0]);
        } else if (parameters.size() == 2) {
            res = new TV_divergence<double>(parameters[0], parameters[1]);
        } else {
            param_len_error("TV", "1 or 2");
        }
    } break;
    case RG: {
        if (parameters.size() == 2) {
            res = new RG_divergence<double>(parameters[0], parameters[1]);
        } else {
            param_len_error("RG", "2");
        }
    } break;
    case Power: {
        if (parameters.size() == 2) {
            if (parameters[1] == 0) { // Berg
                res = new berg_divergence<double>(parameters[0]);
            } else if (parameters[1] == 1) { // KL
                res = new KL_divergence<double>(parameters[0]);
            } else { // generic Power
                res = new power_divergence<double>(parameters[0], parameters[1]);
            }
        } else {
            param_len_error("Power", "2");
        }
    } break;
    case LSR: {
        if (parameters.size() == 4) {
            res = new LSR_divergence<double>(parameters[0], parameters[1], parameters[2], parameters[3]);
        } else {
            param_len_error("LSR", "4");
        }
    } break;
    default: {
        throw new std::invalid_argument("the supplied divergence type is not valid");
    }
    }
    return res;
}



// #################### KL ########################

template<typename T>
KL_divergence<T>::KL_divergence(T lambda) : lambda(lambda) {
    if (lambda <= 0)
        throw new std::invalid_argument("lambda has to be positive");
};
template<typename T>
T KL_divergence<T>::eval(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T res = 0;
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end() && bi != b.end(); ai++, bi++) {
        if (*ai < 0 || (*ai > 0 && *bi == 0))
            return std::numeric_limits<T>::infinity();
        res += (*ai ? *ai * std::log(*ai / *bi) : 0) - *ai + *bi;
    }
    return lambda * res;
}
template<typename T>
T KL_divergence<T>::eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T res = 0;
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ai++, bi++)
        res += *bi ? (std::exp(*ai / lambda) - 1) * *bi : 0;
    return lambda * res;
}
template<typename T>
T KL_divergence<T>::eval_aprox(const T &x, const T &epsi) {
    return lambda / (epsi + lambda) * x;
}
template<typename T>
arma::Mat<T> KL_divergence<T>::eval_aprox(const arma::Mat<T> &x, const T &epsi) {
    return lambda / (epsi + lambda) * x;
}
template<typename T>
arma::Mat<T> KL_divergence<T>::eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) {
    arma::Mat<T> res = zeros_like(s);
    T ile = 1. / (lambda + epsi);
    for (auto si = s.begin(), ui = u.begin(), bi = b.begin(), resi = res.begin(); si != s.end(); si++, ui++, bi++, resi++)
        if (*bi) {
            if (!*si)
                *resi = std::numeric_limits<T>::infinity();
            else
                *resi = std::exp((lambda * std::log(*bi / *si) - *ui) * ile);
        }
    return res;
}


// #################### TV ########################

template<typename T>
TV_divergence<T>::TV_divergence(T lambda) : lambda_pos(lambda), lambda_neg(lambda) { }
template<typename T>
TV_divergence<T>::TV_divergence(T lambda_pos, T lambda_neg) : lambda_pos(lambda_pos), lambda_neg(lambda_neg) { }
template<typename T>
T TV_divergence<T>::eval(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T pos = 0, neg = 0, d;
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end() && bi != b.end(); ai++, bi++) {
        d = *ai - *bi;
        if (d > 0)
            neg += d;
        else
            pos -= d;
    }
    return pos * lambda_pos + neg * lambda_neg;
}
template<typename T>
T TV_divergence<T>::eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T res = 0;
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ai++, bi++) {
        if (*ai < -lambda_neg - divergence_function<T>::indicator_slack || *ai > lambda_pos + divergence_function<T>::indicator_slack)
            return std::numeric_limits<T>::infinity();
        res += *ai * *bi;
    }
    return res;
}
template<typename T>
T TV_divergence<T>::eval_aprox(const T &x, const T &epsi) {
    return std::min(std::max(x, -lambda_neg), lambda_pos);
}
template<typename T>
arma::Mat<T> TV_divergence<T>::eval_aprox(const arma::Mat<T> &x, const T &epsi) {
    return arma::clamp(x, -lambda_neg, lambda_pos);
}
template<typename T>
arma::Mat<T> TV_divergence<T>::eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) {
    arma::Mat<T> res = zeros_like(s);
    T epsi_inv = 1. / epsi;
    for (auto si = s.begin(), ui = u.begin(), bi = b.begin(), resi = res.begin(); si != s.end(); si++, ui++, bi++, resi++)
        if (*bi) {
            if (!*si) 
                *resi = std::exp(epsi_inv * (lambda_pos - *ui));
            else
                *resi = std::max(std::exp(epsi_inv * (-lambda_pos - *ui)), std::min(std::exp(epsi_inv * (lambda_neg - *ui)), *bi / *si));
        } else {
            *resi = std::exp(epsi_inv * (-lambda_pos - *ui));
        }
    return res;
}


// #################### RG ########################

template<typename T>
RG_divergence<T>::RG_divergence(T lower, T upper) : lower(lower), upper(upper) {
    if (lower < 0)
        throw new std::invalid_argument("lower boundary may not be negative");
    if (lower > upper)
        throw new std::invalid_argument("lower boundary may not be larger than upper");
}
template<typename T>
T RG_divergence<T>::eval(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ai++, bi++)
        if (*ai < lower * *bi - divergence_function<T>::indicator_slack || *ai > upper * *bi + divergence_function<T>::indicator_slack)
            return std::numeric_limits<T>::infinity();
    return 0;
}
template<typename T>
T RG_divergence<T>::eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T res = 0;
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ai++, bi++)
        res += *bi ? std::max(*ai * lower, *ai * upper) * *bi : 0;
    return res;
}
template<typename T>
T RG_divergence<T>::eval_aprox(const T &x, const T &epsi) {
    return std::min(std::max(0.0, x - epsi * std::log(upper)), x - epsi * std::log(lower));
}
template<typename T>
arma::Mat<T> RG_divergence<T>::eval_aprox(const arma::Mat<T> &x, const T &epsi) {
    return arma::min(arma::max(zeros_like(x), x - epsi * std::log(upper)), x - epsi * std::log(lower));
}
template<typename T>
arma::Mat<T> RG_divergence<T>::eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) {
    arma::Mat<T> res = zeros_like(s);
    T epsi_inv = 1. / epsi, lb, ub;
    for (auto si = s.begin(), ui = u.begin(), bi = b.begin(), resi = res.begin(); si != s.end(); si++, ui++, bi++, resi++)
        if (*bi) {
            if (!*si)
                *resi = lower ? std::numeric_limits<T>::infinity() : std::exp(-*ui * epsi_inv);
            else
                *resi = std::max(lower * *bi / *si, std::min(upper * *bi / *si, std::exp(-*ui * epsi_inv)));
        }
    return res;
}


// #################### LSR ########################

template<typename T>
LSR_divergence<T>::LSR_divergence(T lower, T upper, T lambda_neg, T lambda_pos) : 
    lower(lower), upper(upper), lambda_pos(lambda_pos), lambda_neg(lambda_neg) {
    if (lambda_neg < 0 || lambda_pos < 0)
        throw new std::invalid_argument("lambda has to be non-negative");
    if (lower < 0)
        throw new std::invalid_argument("lower boundary may not be negative");
    if (lower > upper)
        throw new std::invalid_argument("lower boundary may not be larger than upper");
}
template<typename T>
T LSR_divergence<T>::eval(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T pos = 0, neg = 0;
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end() && bi != b.end(); ai++, bi++) {
        neg += std::max(0.0, lower * *bi - *ai);
        pos += std::max(0.0, *ai - upper * *bi);
    }
    return pos ? pos * lambda_pos : 0 + neg ? neg * lambda_neg : 0;
}
template<typename T>
T LSR_divergence<T>::eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T res = 0;
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ai++, bi++) {
        if (*ai < -lambda_neg - divergence_function<T>::indicator_slack || *ai > lambda_pos + divergence_function<T>::indicator_slack)
            return std::numeric_limits<T>::infinity();
        res += std::max(*ai * lower, *ai * upper) * *bi;
    }
    return res;
}
template<typename T>
T LSR_divergence<T>::eval_aprox(const T &x, const T &epsi) {
    return std::max(-lambda_neg, std::min(lambda_pos, std::min(std::max((T)0, x - epsi * std::log(upper)), x - epsi * std::log(lower))));
}
template<typename T>
arma::Mat<T> LSR_divergence<T>::eval_aprox(const arma::Mat<T> &x, const T &epsi) {
    return arma::clamp(arma::min(arma::max(zeros_like(x), x - epsi * std::log(upper)), x - epsi * std::log(lower)), -lambda_neg, lambda_pos);
}
template<typename T>
arma::Mat<T> LSR_divergence<T>::eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) {
    arma::Mat<T> res = zeros_like(s);
    T epsi_inv = 1. / epsi, t, lb, ub;
    for (auto si = s.begin(), ui = u.begin(), bi = b.begin(), resi = res.begin(); si != s.end(); si++, ui++, bi++, resi++) {
        if (*bi) {
            if (!*si) {
                *resi = lower ? std::exp(epsi_inv * (lambda_pos - *ui)) : std::exp(- epsi_inv * *ui);
            } else {
                t = std::max(lower * *bi / *si, std::min(upper * *bi / *si, std::exp(-*ui * epsi_inv)));
                *resi = std::max(std::exp(epsi_inv * (-lambda_pos - *ui)), std::min(std::exp(epsi_inv * (lambda_neg - *ui)), t));
            }
        } else {
            *resi = std::exp(-epsi_inv * (lambda_pos + *ui));
        }
    }   
    return res;
}


// #################### Berg ########################

template<typename T>
berg_divergence<T>::berg_divergence(T lambda) : lambda(lambda) {
    if (lambda <= 0)
        throw new std::invalid_argument("lambda has to be positive");
}
template<typename T>
T berg_divergence<T>::eval(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T res = 0;
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ai++, bi++) {
        if (*ai < 0 || (*ai == 0 && *bi > 0))
            return std::numeric_limits<T>::infinity();
        else if (*ai >= 0 && *bi == 0)
            res += *ai;
        else 
            res += *ai - *bi - std::log(*ai / *bi) * *bi;
    }
    return lambda * res;
}
template<typename T>
T berg_divergence<T>::eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T res = 0;
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ai++, bi++) {
        if (*bi) {
            if (*ai >= lambda)
                return std::numeric_limits<T>::infinity();
            res -= std::log(1 - *ai / lambda) * *bi;
        } else {
            if (*ai > lambda + divergence_function<T>::indicator_slack)
                return std::numeric_limits<T>::infinity();
        }
    }
    return lambda * res;
}
template<typename T>
T berg_divergence<T>::eval_aprox(const T &x, const T &epsi) {
    return lambda - epsi * lambertWexp(std::log(lambda / epsi) + (lambda - x) / epsi);
}
template<typename T>
arma::Mat<T> berg_divergence<T>::eval_aprox(const arma::Mat<T> &x, const T &epsi) {
    return lambda - epsi * lambertWexp(arma::Mat<T>(std::log(lambda / epsi) + (lambda - x) / epsi));
}
template<typename T>
arma::Mat<T> berg_divergence<T>::eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) {
    arma::Mat<T> res = zeros_like(s);
    T epsi_inv = 1. / epsi, t;
    for (auto si = s.begin(), ui = u.begin(), bi = b.begin(), resi = res.begin(); si != s.end(); si++, ui++, bi++, resi++) {
        if (*bi) {
            if (!*si) {
                *resi = std::exp(-epsi_inv * (lambda + *ui));
            } else {
                t = lambda * *bi / (epsi * *si);
                *resi = t / lambertWexp(std::log(t) + (lambda + *ui) * epsi_inv);
            }
        } else {
            *resi = std::exp(-epsi_inv * (lambda + *ui));
        }
    }   
    return res;
}


// #################### Power ########################

template<typename T>
power_divergence<T>::power_divergence(T lambda, T p) : lambda(lambda), p(p) {
    if (lambda <= 0)
        throw new std::invalid_argument("lambda has to be positive");
    if (p == 0)
        throw new std::invalid_argument("p may not be 0, use Berg instead");
    if (p == 1)
        throw new std::invalid_argument("p may not be 1, use KL instead");
}
template<typename T>
T power_divergence<T>::eval(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T res = 0, t0 = 1. / p, t1 = 1. / (p - 1.);
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ai++, bi++) {
        if (*ai != 0 && *bi == 0) {
            if (p > 1)
                return std::numeric_limits<T>::infinity();
            else
                res += *ai * t1;
        } else {
            if (*ai < 0)
                return std::numeric_limits<T>::infinity();
            else if (*ai > 0) 
                res += std::pow(*ai / *bi, p) * *bi * t0 * t1 - *ai * t1 + *bi * t0;
            else
                res += *bi * t0;
        }
    }
    return lambda * res;
}
template<typename T>
T power_divergence<T>::eval_conjugate(const arma::Mat<T> &a, const arma::Mat<T> &b) {
    T res = 0, t0 = 1. / p, t1 = 1. / (p - 1.);
    for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ai++, bi++) {
        if (*ai != 0 && *bi == 0) {
            if (p < 1 && *ai > lambda / (1 - p))
                return std::numeric_limits<T>::infinity();
        } else {
            if ((1 - p) * *ai <= lambda)
                res += lambda * t0 * (std::pow(1 - (1 - p) / lambda * *ai, p * t1) - 1) * *bi;
            else if (p > 1)
                res -= lambda * t0 * *bi;
            else
                return std::numeric_limits<T>::infinity();
        }
    }
    return lambda * res;
}
template<typename T>
T power_divergence<T>::eval_aprox(const T &x, const T &epsi) {
    return lambda / (1 - p) - epsi / (1 - p) * lambertWexp(std::log(lambda / epsi) + (lambda - (1 - p) * x) / epsi);
}
template<typename T>
arma::Mat<T> power_divergence<T>::eval_aprox(const arma::Mat<T> &x, const T &epsi) {
    return lambda / (1 - p) - epsi / (1 - p) * lambertWexp(arma::Mat<T>(std::log(lambda / epsi) + (lambda - (1 - p) * x) / epsi));
}
template<typename T>
arma::Mat<T> power_divergence<T>::eval_proxdiv(const arma::Mat<T> &s, const arma::Mat<T> &u, const arma::Mat<T> &b, const T &epsi) {
    arma::Mat<T> res = zeros_like(s);
    T epsi_inv = 1. / epsi;
    for (auto si = s.begin(), ui = u.begin(), bi = b.begin(), resi = res.begin(); si != s.end(); si++, ui++, bi++, resi++) {
        if (*bi) {
            if (!*si) {
                *resi = std::exp(-epsi_inv * (lambda / (1 - p) + *ui));
            } else {
                *resi = *bi / *si * std::pow(epsi / lambda * lambertWexp(std::log(lambda * epsi_inv) 
                                                                             + (p - 1) * std::log(*si / *bi) 
                                                                             + (lambda + *ui * (1 - p)) * epsi_inv), 1 / (p - 1));
            }   
        } else {
            if (p > 1)
                *resi = 0;
            else
                *resi = std::exp(-epsi_inv * (lambda / (1 - p) + *ui));
        }
    }   
    return res;
}


#endif