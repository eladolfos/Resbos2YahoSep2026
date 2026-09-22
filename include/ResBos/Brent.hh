#ifndef CROSSING_HH
#define CROSSING_HH

#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <vector>

// Class to hold the needed information for preforming the Brent algorithm
class Brent {
    public:
        Brent(const std::function<double(double)>& func_, const double& tol_=5E-8) 
            : func(func_), tol(tol_) {};
        void SetArgs(const std::vector<double> &args_) { args = args_; }
        double CalcRoot(double, double);
        double Minimize(const double&, const double&, const double &b=1E99) const;

    private:
        static inline void shift(double &a, double &b, double &c, const double &d) {
            a = b; b = c; c = d;
        }
        void GetBracket(double &a, double &c, double &b,
                        double &fa, double &fb, double &fc) const;        
        static inline void swap(double &fa, double &fb, double &a, double &b) {
            if(std::abs(fa) < std::abs(fb)) {
                double tmp = a;
                a = b;
                b = tmp;
                tmp = fa;
                fa = fb;
                fb = tmp;
            }
            return;
        };
        std::function<double(double)> func;
        double tol;
        bool mflag;
        std::vector<double> args;
        static const size_t itmax = 100;
        static constexpr double tiny = 1.0E-20;
        static constexpr double glimit = 100.0;
        static constexpr double gold = 1.618034;
        static constexpr double cgold = 0.3819660;
        static constexpr double eps = std::numeric_limits<double>::epsilon()*1.0E-3;
};

#endif
