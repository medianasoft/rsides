#include <exception>
#include <memory>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <string>
#include <unistd.h>
#include <numeric>
#include <cmath>
#include <map>
#include <algorithm>
#include <utility>
#include <iomanip>

#include <Rcpp.h>
#include <RcppEigen.h>

#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

using namespace rapidxml;
using namespace std;
using namespace Rcpp;

using Eigen::Map;                       
using Eigen::MatrixXd;                  
using Eigen::VectorXd; 

namespace ublas = boost::numeric::ublas; 

#define asv(x) (as<vector<double>>(x))

typedef std::pair<double,double> ddpair;
typedef std::pair<double,int> dipair;

Rcpp::Environment base("package:base");
Function UniquePreserveOrder = base["unique"];

// Global parameters
double precision = 0.001;
int max_iter = 10;
vector<int> biomarker_level;

struct ScoreResult {
    Rcpp::NumericMatrix score;
    Rcpp::NumericMatrix information;
};

// Structure to store all characteristic of a patient subgroup
struct ModelCovariates {
    vector<double> cov1;
    vector<double> cov2;
    vector<double> cov3;
    vector<double> cov4;
    vector<double> cov_class;
};

// Structure to store all characteristic of a patient subgroup
struct SingleSubgroup {
    // Splitting criterion
    double criterion;
    // Splitting criterion on p-value scale
    double criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    double adjusted_criterion_pvalue;
    double test_statistic;
    double pvalue;
    double prom_estimate;
    double prom_sderr;
    double prom_sd;
    double adjusted_pvalue;
    // Vector of biomarker values to define the current patient subgroup
    std::vector<double> value;
    // 1 if <=, 2 if >, 3 if =
    int sign;
    // Size of the current patient subgroup (control, treatment)
    int size;    
    int size_control;
    int size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    std::vector<int> subgroup_rows;
    // Index of biomarker used in the current patient subgroup
    int biomarker_index;
    // Level of the current patient subgroup
    int level;
    // Parent index for the current patient subgroup
    int parent_index;
    // Indexes of child subgroups for the current patient subgroup
    std::vector<int> child_index;     
    // Error code for testing (0 if no errors)
    int code;
    // Is the current patient subgroup terminal (0/1)
    int terminal_subgroup;

    int signat;

    std::vector<SingleSubgroup> subgroups;
};

struct LogRankdata {
    double t;
    double cens;
    double id;
    bool operator < (const LogRankdata& str) const
    {
        return (t < str.t);
    }
    bool operator > (const LogRankdata& str) const
    {
        return (t >= str.t);
    }
};


bool DIPairDown (const dipair& l, const dipair& r) { 
    return l.first > r.first; 
}

void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// Function to randomize permutations
// https://gallery.rcpp.org/articles/stl-random-shuffle/
template<typename VectorType>
void shuffle_vector(vector<VectorType> &vec) {
    int n = vec.size();
    int j;

    // Fisher-Yates Shuffle Algorithm
    for (int i = 0; i < n - 1; i++) {
        j = i + R::runif(0, n - i);
        std::swap(vec[i], vec[j]);
    }
}

double rcpp_pnorm(const double &x) {
    NumericVector vec_input(1), vec_output(1);
    vec_input[0] = x;
    vec_output = Rcpp::pnorm(vec_input);
    return vec_output[0];
}

double rcpp_pt(const double &x, const double &df) {
    NumericVector vec_input(1), vec_output(1);
    vec_input[0] = x;
    vec_output = Rcpp::pt(vec_input, df);
    return vec_output[0];
}

// Supportive function used in logrank test  
void TupleSort(std::vector<double> &in1, std::vector<double> &in2, std::vector<double> &in3,vector<LogRankdata> &vec) {
    vec.resize(in1.size());
    for (int i = 0; i < in1.size();++i) {
        vec[i].t = in1[i];
        vec[i].cens = in2[i];
        vec[i].id = in3[i];
    }
    sort(vec.begin(),vec.end());
}

double vecsum(const vector<double> &vec) {
    int i, m = vec.size();
    double sum = 0.0;
    for (i = 0; i < m; ++i) sum += vec[i];
    return sum;
}

// Hazard ratio estimate based on an exponential distribution assumption
double HazardRatio(const std::vector<double> &outcome, const std::vector<double> &censor, std::vector<double> &treatment, const int &direction) {

    double rate_control, rate_treatment, hazard_ratio, outcome_control = 0.0, outcome_treatment = 0.0, n_control = 0.0, n_treatment = 0.0, censor_control = 0.0, censor_treatment = 0.0;    
    int i, m = outcome.size();

    for (i = 0; i < m; i++) {

        if (treatment[i] == 0.0) {
            n_control++;
            outcome_control += outcome[i];
            censor_control += censor[i];
        } else {
            n_treatment++;
            outcome_treatment += outcome[i];
            censor_treatment += censor[i];            
        }

    }

    // Hazard rates
    rate_control = (n_control + 0.0 - censor_control) / outcome_control;
    rate_treatment = (n_treatment + 0.0 - censor_treatment) / outcome_treatment;

    // Hazard ratio
    if (direction == 1) hazard_ratio = rate_treatment / rate_control; else hazard_ratio = rate_control / rate_treatment;

    return hazard_ratio;  

}


// Treatment effect test for time-to-event outcome variable (logrank test)
double LRTest(std::vector<double> &xy, std::vector<double> &cxy, std::vector<double> &ixy, const int &direction)
{

    std::vector<LogRankdata> vec;

    TupleSort(xy,cxy,ixy,vec);
    std::vector<double> t;
    t.reserve(vec.size());
    std::vector<int> m1;
    m1.reserve(vec.size());
    std::vector<int> m2;
    m2.reserve(vec.size());
    std::vector<int> cens1,cens2;
    cens1.reserve(vec.size());cens2.reserve(vec.size());
    int curit=0;
    for (int i =0; i <vec.size();++i) {
        t.push_back(vec[i].t);
        m1.push_back(0);
        m2.push_back(0);
        cens1.push_back(0);
        cens2.push_back(0);
        if(vec[i].id == 1)
        {
            m1[curit] += 1-vec[i].cens;
            cens1[curit]+=vec[i].cens;
        }
        else
        {
            m2[curit] += 1-vec[i].cens;
            cens2[curit]+=vec[i].cens;
        }
        if(i+1<vec.size())
            while (i+1<vec.size() && vec[i].t==vec[i+1].t)  // VIC: Additional check "i+1<vec.size() &&" )
            {
                ++i;
                if(vec[i].id == 1)
                {
                    m1[curit] += 1-vec[i].cens;
                    cens1[curit]+=vec[i].cens;
                }
                else
                {
                    m2[curit] += 1-vec[i].cens;
                    cens2[curit] += (int)vec[i].cens;
                }
                if(i+1==vec.size())
                    break;
            }
        ++curit;
    }

    int n1 = 0, n2;
    for(int i =0; i <vec.size();++i)
        if (vec[i].id ==1)
            ++n1;
    n2 = vec.size() - n1;
    std::vector<double> oe1,voe;
    oe1.resize(curit);voe.resize(curit);
    double s1 = 0, s2 = 0, oe = 0, v = 0;

    for(int i = 0; i<curit;++i)
    {
        oe = m1[i]-((double)n1)/(n1+n2)*(m1[i]+m2[i]);
        oe1[i] = oe;
        v = (m1[i]+m2[i])*((double)n1/(n1+n2))*(1-(double)n1/(n1+n2))*(n1+n2-m1[i]-m2[i])/(n1+n2-1);
        v = (v!=v)?0:v;
        voe[i] = v;
        s1 += oe;
        s2 += v;
        n1 -= m1[i]+cens1[i];
        n2 -= m2[i]+cens2[i];
    }

    s1 = s1 * (-direction + 0.0);

    return (s2>std::numeric_limits<double>::epsilon())?s1/sqrt(s2):s1/std::numeric_limits<double>::epsilon();

}

// Compute the size of a subgroup tree
void TreeSize(vector<SingleSubgroup> sub, int &length, int &width){
    int l=0,w=0, totl=0,totw=0;
    for(int i = 0; i<sub.size();++i){
        if (!sub[i].subgroups.empty()){
            TreeSize(sub[i].subgroups,l,w);
            totl += l;
            totw = std::max(w,totw);
        }
    }
    totl += sub.size();
    ++totw;
    length = totl;
    width =totw;
}

double Quantile(const vector<double> &vec, const double &prob) {

    int i, m = vec.size();

    Rcpp::NumericVector rvec(m);

    for(i = 0; i<m; i++) rvec[i] = vec[i];

  // Obtain environment containing function
  Rcpp::Environment base("package:stats"); 

  // Make function callable from C++
  Rcpp::Function rquantile = base["quantile"];    

  // Call the function and receive its list output
  Rcpp::NumericVector res = rquantile(Rcpp::_["x"] = rvec,
                                      Rcpp::_["probs"]  = prob); 

  // Return test object in list structure
  return res[0];

}

vector<double> QuantileTransform(const vector<double> &vec, const int &nperc) {

    int i, j, m = vec.size();
    vector<double> complete_vec, res, transform(m);

    double prob, temp;

    for (i = 0; i < m; i++) {
        if (!::isnan(vec[i])) complete_vec.push_back(vec[i]);
    }

    for (i = 0; i < nperc - 1; i++) {
        prob = (i + 1.0) / (nperc + 0.0);
        temp = round(1000.0 * Quantile(complete_vec, prob)) / 1000.0;
        res.push_back(temp);
    }

    for (i = 0; i < m; i++) {
        if (!::isnan(vec[i])) {    
            if (vec[i] <= res[0]) transform[i] = res[0];
            if (vec[i] > res[nperc - 2]) transform[i] = res[nperc - 2];
            for (j = 0; j < nperc - 2; j++) {
                if (vec[i] > res[j] && vec[i] <= res[j + 1]) transform[i] = res[j + 1];
            }
        } else {
            transform[i] = vec[i];
        }

    }

    return transform;

}


// Compute an approximate penalized criterion on log scale: log(crit)+log(G)
double LogerrApprx(const double &zcrit, const double &g, const int &sides) {
  
  double x = zcrit / sqrt(2.0);
  double res = x * x + log(377.0 / 324.0 * x + sqrt(1 + 314.0 / 847.0 * x * x)) - log(g) + (2 - sides) * log(2.0);
  res = -res;
  return res;

}

// Compute an approximate computing of -log(pvalue), pvalue is 2-sided based on 2*(1-F(x))
double Log1minusP(const double &pval, const double &zcrit, const double &g, const int &sides) {
  
  // Computes LN(2*(1-F(z_crit)), F(x) is normal cdf
  double res;

  if (pval > 1E-16) {
    res = log(1.0 - pow(1.0 - pval, g));
  }
  else {
    res= LogerrApprx(zcrit, g, 2);
  }
  
  return res;

}

// Compute the number of unique values or levels for a biomarker (missing values are excluded)
int CountUniqueValues(const vector<double> &vec) {

    int length = vec.size();
    int count = 1;
    int found, i, j;
    double value1, value2;

    for (i = 1; i < length; i++) {

        found = 0;

        value1 = vec[i];

        // Compare with other values if the current value is non-missing
        if (::isnan(value1) < 1) {

            for (j = 0; j < i; j++) {

                value2 = vec[j];
                if (::isnan(value2) < 1)   
                  if (value1 == value2) found = found + 1;
                
                }

            if (found == 0) count = count + 1;    
        }
    }

    return count;
 
}

// Compute parameters used in local multiplicity adjustment for Criteria 1 and 2
double AdjParameterCriteria(const int &biomarker_type, const int &n_unique, const int &criterion_type) {

    double adjustment = 1.0;
    int n_levels;
    
    // Continuous biomarker
    if (biomarker_type == 1) {

        if (n_unique > 400) {
            n_levels = 400;
        }
        else {
           n_levels = n_unique;
        }

        if (n_levels ==  2 && criterion_type == 1) adjustment =  0; 
        if (n_levels ==  2 && criterion_type == 2) adjustment =  0; 
        if (n_levels ==  3 && criterion_type == 1) adjustment =  0.3062132; 
        if (n_levels ==  3 && criterion_type == 2) adjustment =  0.1386845; 
        if (n_levels ==  4 && criterion_type == 1) adjustment =  0.3456545; 
        if (n_levels ==  4 && criterion_type == 2) adjustment =  0.1585821; 
        if (n_levels ==  5 && criterion_type == 1) adjustment =  0.3690571; 
        if (n_levels ==  5 && criterion_type == 2) adjustment =  0.1648235; 
        if (n_levels ==  6 && criterion_type == 1) adjustment =  0.3844297; 
        if (n_levels ==  6 && criterion_type == 2) adjustment =  0.1673766; 
        if (n_levels ==  7 && criterion_type == 1) adjustment =  0.3952477; 
        if (n_levels ==  7 && criterion_type == 2) adjustment =  0.1685605; 
        if (n_levels ==  8 && criterion_type == 1) adjustment =  0.403248; 
        if (n_levels ==  8 && criterion_type == 2) adjustment =  0.1691396; 
        if (n_levels ==  9 && criterion_type == 1) adjustment =  0.4093902; 
        if (n_levels ==  9 && criterion_type == 2) adjustment =  0.1694212; 
        if (n_levels ==  10 && criterion_type == 1) adjustment =  0.4142456; 
        if (n_levels ==  10 && criterion_type == 2) adjustment =  0.1695453; 
        if (n_levels ==  11 && criterion_type == 1) adjustment =  0.4181749; 
        if (n_levels ==  11 && criterion_type == 2) adjustment =  0.1695823; 
        if (n_levels ==  12 && criterion_type == 1) adjustment =  0.4214165; 
        if (n_levels ==  12 && criterion_type == 2) adjustment =  0.1695697; 
        if (n_levels ==  13 && criterion_type == 1) adjustment =  0.4241338; 
        if (n_levels ==  13 && criterion_type == 2) adjustment =  0.1695287; 
        if (n_levels ==  14 && criterion_type == 1) adjustment =  0.4264429; 
        if (n_levels ==  14 && criterion_type == 2) adjustment =  0.1694715; 
        if (n_levels ==  15 && criterion_type == 1) adjustment =  0.428428; 
        if (n_levels ==  15 && criterion_type == 2) adjustment =  0.1694055; 
        if (n_levels ==  16 && criterion_type == 1) adjustment =  0.430152; 
        if (n_levels ==  16 && criterion_type == 2) adjustment =  0.1693352; 
        if (n_levels ==  17 && criterion_type == 1) adjustment =  0.4316625; 
        if (n_levels ==  17 && criterion_type == 2) adjustment =  0.1692635; 
        if (n_levels ==  18 && criterion_type == 1) adjustment =  0.4329963; 
        if (n_levels ==  18 && criterion_type == 2) adjustment =  0.169192; 
        if (n_levels ==  19 && criterion_type == 1) adjustment =  0.4341823; 
        if (n_levels ==  19 && criterion_type == 2) adjustment =  0.1691219; 
        if (n_levels ==  20 && criterion_type == 1) adjustment =  0.4352434; 
        if (n_levels ==  20 && criterion_type == 2) adjustment =  0.1690538; 
        if (n_levels ==  21 && criterion_type == 1) adjustment =  0.4361981; 
        if (n_levels ==  21 && criterion_type == 2) adjustment =  0.1689881; 
        if (n_levels ==  22 && criterion_type == 1) adjustment =  0.4370614; 
        if (n_levels ==  22 && criterion_type == 2) adjustment =  0.1689248; 
        if (n_levels ==  23 && criterion_type == 1) adjustment =  0.4378457; 
        if (n_levels ==  23 && criterion_type == 2) adjustment =  0.1688642; 
        if (n_levels ==  24 && criterion_type == 1) adjustment =  0.4385613; 
        if (n_levels ==  24 && criterion_type == 2) adjustment =  0.1688062; 
        if (n_levels ==  25 && criterion_type == 1) adjustment =  0.4392165; 
        if (n_levels ==  25 && criterion_type == 2) adjustment =  0.1687507; 
        if (n_levels ==  26 && criterion_type == 1) adjustment =  0.4398188; 
        if (n_levels ==  26 && criterion_type == 2) adjustment =  0.1686977; 
        if (n_levels ==  27 && criterion_type == 1) adjustment =  0.4403741; 
        if (n_levels ==  27 && criterion_type == 2) adjustment =  0.1686471; 
        if (n_levels ==  28 && criterion_type == 1) adjustment =  0.4408876; 
        if (n_levels ==  28 && criterion_type == 2) adjustment =  0.1685987; 
        if (n_levels ==  29 && criterion_type == 1) adjustment =  0.4413639; 
        if (n_levels ==  29 && criterion_type == 2) adjustment =  0.1685525; 
        if (n_levels ==  30 && criterion_type == 1) adjustment =  0.4418069; 
        if (n_levels ==  30 && criterion_type == 2) adjustment =  0.1685083; 
        if (n_levels ==  31 && criterion_type == 1) adjustment =  0.4422197; 
        if (n_levels ==  31 && criterion_type == 2) adjustment =  0.1684661; 
        if (n_levels ==  32 && criterion_type == 1) adjustment =  0.4426055; 
        if (n_levels ==  32 && criterion_type == 2) adjustment =  0.1684257; 
        if (n_levels ==  33 && criterion_type == 1) adjustment =  0.4429666; 
        if (n_levels ==  33 && criterion_type == 2) adjustment =  0.168387; 
        if (n_levels ==  34 && criterion_type == 1) adjustment =  0.4433055; 
        if (n_levels ==  34 && criterion_type == 2) adjustment =  0.1683499; 
        if (n_levels ==  35 && criterion_type == 1) adjustment =  0.443624; 
        if (n_levels ==  35 && criterion_type == 2) adjustment =  0.1683143; 
        if (n_levels ==  36 && criterion_type == 1) adjustment =  0.4439239; 
        if (n_levels ==  36 && criterion_type == 2) adjustment =  0.1682802; 
        if (n_levels ==  37 && criterion_type == 1) adjustment =  0.4442067; 
        if (n_levels ==  37 && criterion_type == 2) adjustment =  0.1682475; 
        if (n_levels ==  38 && criterion_type == 1) adjustment =  0.444474; 
        if (n_levels ==  38 && criterion_type == 2) adjustment =  0.1682161; 
        if (n_levels ==  39 && criterion_type == 1) adjustment =  0.4447268; 
        if (n_levels ==  39 && criterion_type == 2) adjustment =  0.1681858; 
        if (n_levels ==  40 && criterion_type == 1) adjustment =  0.4449664; 
        if (n_levels ==  40 && criterion_type == 2) adjustment =  0.1681568; 
        if (n_levels ==  41 && criterion_type == 1) adjustment =  0.4451937; 
        if (n_levels ==  41 && criterion_type == 2) adjustment =  0.1681288; 
        if (n_levels ==  42 && criterion_type == 1) adjustment =  0.4454096; 
        if (n_levels ==  42 && criterion_type == 2) adjustment =  0.1681019; 
        if (n_levels ==  43 && criterion_type == 1) adjustment =  0.445615; 
        if (n_levels ==  43 && criterion_type == 2) adjustment =  0.1680759; 
        if (n_levels ==  44 && criterion_type == 1) adjustment =  0.4458107; 
        if (n_levels ==  44 && criterion_type == 2) adjustment =  0.1680509; 
        if (n_levels ==  45 && criterion_type == 1) adjustment =  0.4459971; 
        if (n_levels ==  45 && criterion_type == 2) adjustment =  0.1680267; 
        if (n_levels ==  46 && criterion_type == 1) adjustment =  0.4461751; 
        if (n_levels ==  46 && criterion_type == 2) adjustment =  0.1680034; 
        if (n_levels ==  47 && criterion_type == 1) adjustment =  0.4463452; 
        if (n_levels ==  47 && criterion_type == 2) adjustment =  0.1679809; 
        if (n_levels ==  48 && criterion_type == 1) adjustment =  0.4465078; 
        if (n_levels ==  48 && criterion_type == 2) adjustment =  0.1679591; 
        if (n_levels ==  49 && criterion_type == 1) adjustment =  0.4466635; 
        if (n_levels ==  49 && criterion_type == 2) adjustment =  0.1679381; 
        if (n_levels ==  50 && criterion_type == 1) adjustment =  0.4468126; 
        if (n_levels ==  50 && criterion_type == 2) adjustment =  0.1679177; 
        if (n_levels ==  51 && criterion_type == 1) adjustment =  0.4469556; 
        if (n_levels ==  51 && criterion_type == 2) adjustment =  0.167898; 
        if (n_levels ==  52 && criterion_type == 1) adjustment =  0.4470928; 
        if (n_levels ==  52 && criterion_type == 2) adjustment =  0.1678789; 
        if (n_levels ==  53 && criterion_type == 1) adjustment =  0.4472246; 
        if (n_levels ==  53 && criterion_type == 2) adjustment =  0.1678604; 
        if (n_levels ==  54 && criterion_type == 1) adjustment =  0.4473513; 
        if (n_levels ==  54 && criterion_type == 2) adjustment =  0.1678424; 
        if (n_levels ==  55 && criterion_type == 1) adjustment =  0.4474732; 
        if (n_levels ==  55 && criterion_type == 2) adjustment =  0.1678251; 
        if (n_levels ==  56 && criterion_type == 1) adjustment =  0.4475905; 
        if (n_levels ==  56 && criterion_type == 2) adjustment =  0.1678082; 
        if (n_levels ==  57 && criterion_type == 1) adjustment =  0.4477035; 
        if (n_levels ==  57 && criterion_type == 2) adjustment =  0.1677918; 
        if (n_levels ==  58 && criterion_type == 1) adjustment =  0.4478125; 
        if (n_levels ==  58 && criterion_type == 2) adjustment =  0.1677759; 
        if (n_levels ==  59 && criterion_type == 1) adjustment =  0.4479175; 
        if (n_levels ==  59 && criterion_type == 2) adjustment =  0.1677604; 
        if (n_levels ==  60 && criterion_type == 1) adjustment =  0.4480189; 
        if (n_levels ==  60 && criterion_type == 2) adjustment =  0.1677454; 
        if (n_levels ==  61 && criterion_type == 1) adjustment =  0.4481169; 
        if (n_levels ==  61 && criterion_type == 2) adjustment =  0.1677308; 
        if (n_levels ==  62 && criterion_type == 1) adjustment =  0.4482115; 
        if (n_levels ==  62 && criterion_type == 2) adjustment =  0.1677166; 
        if (n_levels ==  63 && criterion_type == 1) adjustment =  0.4483029; 
        if (n_levels ==  63 && criterion_type == 2) adjustment =  0.1677027; 
        if (n_levels ==  64 && criterion_type == 1) adjustment =  0.4483914; 
        if (n_levels ==  64 && criterion_type == 2) adjustment =  0.1676893; 
        if (n_levels ==  65 && criterion_type == 1) adjustment =  0.4484771; 
        if (n_levels ==  65 && criterion_type == 2) adjustment =  0.1676761; 
        if (n_levels ==  66 && criterion_type == 1) adjustment =  0.44856; 
        if (n_levels ==  66 && criterion_type == 2) adjustment =  0.1676634; 
        if (n_levels ==  67 && criterion_type == 1) adjustment =  0.4486403; 
        if (n_levels ==  67 && criterion_type == 2) adjustment =  0.1676509; 
        if (n_levels ==  68 && criterion_type == 1) adjustment =  0.4487182; 
        if (n_levels ==  68 && criterion_type == 2) adjustment =  0.1676388; 
        if (n_levels ==  69 && criterion_type == 1) adjustment =  0.4487937; 
        if (n_levels ==  69 && criterion_type == 2) adjustment =  0.167627; 
        if (n_levels ==  70 && criterion_type == 1) adjustment =  0.448867; 
        if (n_levels ==  70 && criterion_type == 2) adjustment =  0.1676154; 
        if (n_levels ==  71 && criterion_type == 1) adjustment =  0.448938; 
        if (n_levels ==  71 && criterion_type == 2) adjustment =  0.1676042; 
        if (n_levels ==  72 && criterion_type == 1) adjustment =  0.4490071; 
        if (n_levels ==  72 && criterion_type == 2) adjustment =  0.1675932; 
        if (n_levels ==  73 && criterion_type == 1) adjustment =  0.4490741; 
        if (n_levels ==  73 && criterion_type == 2) adjustment =  0.1675824; 
        if (n_levels ==  74 && criterion_type == 1) adjustment =  0.4491393; 
        if (n_levels ==  74 && criterion_type == 2) adjustment =  0.167572; 
        if (n_levels ==  75 && criterion_type == 1) adjustment =  0.4492026; 
        if (n_levels ==  75 && criterion_type == 2) adjustment =  0.1675617; 
        if (n_levels ==  76 && criterion_type == 1) adjustment =  0.4492643; 
        if (n_levels ==  76 && criterion_type == 2) adjustment =  0.1675517; 
        if (n_levels ==  77 && criterion_type == 1) adjustment =  0.4493242; 
        if (n_levels ==  77 && criterion_type == 2) adjustment =  0.1675419; 
        if (n_levels ==  78 && criterion_type == 1) adjustment =  0.4493825; 
        if (n_levels ==  78 && criterion_type == 2) adjustment =  0.1675324; 
        if (n_levels ==  79 && criterion_type == 1) adjustment =  0.4494393; 
        if (n_levels ==  79 && criterion_type == 2) adjustment =  0.167523; 
        if (n_levels ==  80 && criterion_type == 1) adjustment =  0.4494946; 
        if (n_levels ==  80 && criterion_type == 2) adjustment =  0.1675139; 
        if (n_levels ==  81 && criterion_type == 1) adjustment =  0.4495485; 
        if (n_levels ==  81 && criterion_type == 2) adjustment =  0.1675049; 
        if (n_levels ==  82 && criterion_type == 1) adjustment =  0.449601; 
        if (n_levels ==  82 && criterion_type == 2) adjustment =  0.1674962; 
        if (n_levels ==  83 && criterion_type == 1) adjustment =  0.4496522; 
        if (n_levels ==  83 && criterion_type == 2) adjustment =  0.1674876; 
        if (n_levels ==  84 && criterion_type == 1) adjustment =  0.4497021; 
        if (n_levels ==  84 && criterion_type == 2) adjustment =  0.1674792; 
        if (n_levels ==  85 && criterion_type == 1) adjustment =  0.4497508; 
        if (n_levels ==  85 && criterion_type == 2) adjustment =  0.167471; 
        if (n_levels ==  86 && criterion_type == 1) adjustment =  0.4497983; 
        if (n_levels ==  86 && criterion_type == 2) adjustment =  0.167463; 
        if (n_levels ==  87 && criterion_type == 1) adjustment =  0.4498447; 
        if (n_levels ==  87 && criterion_type == 2) adjustment =  0.1674551; 
        if (n_levels ==  88 && criterion_type == 1) adjustment =  0.44989; 
        if (n_levels ==  88 && criterion_type == 2) adjustment =  0.1674473; 
        if (n_levels ==  89 && criterion_type == 1) adjustment =  0.4499342; 
        if (n_levels ==  89 && criterion_type == 2) adjustment =  0.1674398; 
        if (n_levels ==  90 && criterion_type == 1) adjustment =  0.4499774; 
        if (n_levels ==  90 && criterion_type == 2) adjustment =  0.1674323; 
        if (n_levels ==  91 && criterion_type == 1) adjustment =  0.4500196; 
        if (n_levels ==  91 && criterion_type == 2) adjustment =  0.1674251; 
        if (n_levels ==  92 && criterion_type == 1) adjustment =  0.4500608; 
        if (n_levels ==  92 && criterion_type == 2) adjustment =  0.1674179; 
        if (n_levels ==  93 && criterion_type == 1) adjustment =  0.4501012; 
        if (n_levels ==  93 && criterion_type == 2) adjustment =  0.1674109; 
        if (n_levels ==  94 && criterion_type == 1) adjustment =  0.4501406; 
        if (n_levels ==  94 && criterion_type == 2) adjustment =  0.167404; 
        if (n_levels ==  95 && criterion_type == 1) adjustment =  0.4501792; 
        if (n_levels ==  95 && criterion_type == 2) adjustment =  0.1673973; 
        if (n_levels ==  96 && criterion_type == 1) adjustment =  0.4502169; 
        if (n_levels ==  96 && criterion_type == 2) adjustment =  0.1673907; 
        if (n_levels ==  97 && criterion_type == 1) adjustment =  0.4502538; 
        if (n_levels ==  97 && criterion_type == 2) adjustment =  0.1673842; 
        if (n_levels ==  98 && criterion_type == 1) adjustment =  0.4502899; 
        if (n_levels ==  98 && criterion_type == 2) adjustment =  0.1673778; 
        if (n_levels ==  99 && criterion_type == 1) adjustment =  0.4503253; 
        if (n_levels ==  99 && criterion_type == 2) adjustment =  0.1673715; 
        if (n_levels ==  100 && criterion_type == 1) adjustment =  0.4503599; 
        if (n_levels ==  100 && criterion_type == 2) adjustment =  0.1673654; 
        if (n_levels ==  101 && criterion_type == 1) adjustment =  0.4503938; 
        if (n_levels ==  101 && criterion_type == 2) adjustment =  0.1673593; 
        if (n_levels ==  102 && criterion_type == 1) adjustment =  0.4504271; 
        if (n_levels ==  102 && criterion_type == 2) adjustment =  0.1673534; 
        if (n_levels ==  103 && criterion_type == 1) adjustment =  0.4504596; 
        if (n_levels ==  103 && criterion_type == 2) adjustment =  0.1673476; 
        if (n_levels ==  104 && criterion_type == 1) adjustment =  0.4504915; 
        if (n_levels ==  104 && criterion_type == 2) adjustment =  0.1673418; 
        if (n_levels ==  105 && criterion_type == 1) adjustment =  0.4505228; 
        if (n_levels ==  105 && criterion_type == 2) adjustment =  0.1673362; 
        if (n_levels ==  106 && criterion_type == 1) adjustment =  0.4505534; 
        if (n_levels ==  106 && criterion_type == 2) adjustment =  0.1673307; 
        if (n_levels ==  107 && criterion_type == 1) adjustment =  0.4505835; 
        if (n_levels ==  107 && criterion_type == 2) adjustment =  0.1673252; 
        if (n_levels ==  108 && criterion_type == 1) adjustment =  0.450613; 
        if (n_levels ==  108 && criterion_type == 2) adjustment =  0.1673199; 
        if (n_levels ==  109 && criterion_type == 1) adjustment =  0.4506419; 
        if (n_levels ==  109 && criterion_type == 2) adjustment =  0.1673146; 
        if (n_levels ==  110 && criterion_type == 1) adjustment =  0.4506703; 
        if (n_levels ==  110 && criterion_type == 2) adjustment =  0.1673094; 
        if (n_levels ==  111 && criterion_type == 1) adjustment =  0.4506981; 
        if (n_levels ==  111 && criterion_type == 2) adjustment =  0.1673043; 
        if (n_levels ==  112 && criterion_type == 1) adjustment =  0.4507254; 
        if (n_levels ==  112 && criterion_type == 2) adjustment =  0.1672993; 
        if (n_levels ==  113 && criterion_type == 1) adjustment =  0.4507522; 
        if (n_levels ==  113 && criterion_type == 2) adjustment =  0.1672944; 
        if (n_levels ==  114 && criterion_type == 1) adjustment =  0.4507786; 
        if (n_levels ==  114 && criterion_type == 2) adjustment =  0.1672895; 
        if (n_levels ==  115 && criterion_type == 1) adjustment =  0.4508044; 
        if (n_levels ==  115 && criterion_type == 2) adjustment =  0.1672847; 
        if (n_levels ==  116 && criterion_type == 1) adjustment =  0.4508298; 
        if (n_levels ==  116 && criterion_type == 2) adjustment =  0.16728; 
        if (n_levels ==  117 && criterion_type == 1) adjustment =  0.4508548; 
        if (n_levels ==  117 && criterion_type == 2) adjustment =  0.1672754; 
        if (n_levels ==  118 && criterion_type == 1) adjustment =  0.4508793; 
        if (n_levels ==  118 && criterion_type == 2) adjustment =  0.1672708; 
        if (n_levels ==  119 && criterion_type == 1) adjustment =  0.4509033; 
        if (n_levels ==  119 && criterion_type == 2) adjustment =  0.1672663; 
        if (n_levels ==  120 && criterion_type == 1) adjustment =  0.450927; 
        if (n_levels ==  120 && criterion_type == 2) adjustment =  0.1672619; 
        if (n_levels ==  121 && criterion_type == 1) adjustment =  0.4509502; 
        if (n_levels ==  121 && criterion_type == 2) adjustment =  0.1672575; 
        if (n_levels ==  122 && criterion_type == 1) adjustment =  0.4509731; 
        if (n_levels ==  122 && criterion_type == 2) adjustment =  0.1672532; 
        if (n_levels ==  123 && criterion_type == 1) adjustment =  0.4509955; 
        if (n_levels ==  123 && criterion_type == 2) adjustment =  0.167249; 
        if (n_levels ==  124 && criterion_type == 1) adjustment =  0.4510176; 
        if (n_levels ==  124 && criterion_type == 2) adjustment =  0.1672448; 
        if (n_levels ==  125 && criterion_type == 1) adjustment =  0.4510393; 
        if (n_levels ==  125 && criterion_type == 2) adjustment =  0.1672407; 
        if (n_levels ==  126 && criterion_type == 1) adjustment =  0.4510607; 
        if (n_levels ==  126 && criterion_type == 2) adjustment =  0.1672367; 
        if (n_levels ==  127 && criterion_type == 1) adjustment =  0.4510817; 
        if (n_levels ==  127 && criterion_type == 2) adjustment =  0.1672327; 
        if (n_levels ==  128 && criterion_type == 1) adjustment =  0.4511024; 
        if (n_levels ==  128 && criterion_type == 2) adjustment =  0.1672288; 
        if (n_levels ==  129 && criterion_type == 1) adjustment =  0.4511227; 
        if (n_levels ==  129 && criterion_type == 2) adjustment =  0.1672249; 
        if (n_levels ==  130 && criterion_type == 1) adjustment =  0.4511427; 
        if (n_levels ==  130 && criterion_type == 2) adjustment =  0.1672211; 
        if (n_levels ==  131 && criterion_type == 1) adjustment =  0.4511624; 
        if (n_levels ==  131 && criterion_type == 2) adjustment =  0.1672173; 
        if (n_levels ==  132 && criterion_type == 1) adjustment =  0.4511818; 
        if (n_levels ==  132 && criterion_type == 2) adjustment =  0.1672136; 
        if (n_levels ==  133 && criterion_type == 1) adjustment =  0.4512009; 
        if (n_levels ==  133 && criterion_type == 2) adjustment =  0.1672099; 
        if (n_levels ==  134 && criterion_type == 1) adjustment =  0.4512197; 
        if (n_levels ==  134 && criterion_type == 2) adjustment =  0.1672063; 
        if (n_levels ==  135 && criterion_type == 1) adjustment =  0.4512382; 
        if (n_levels ==  135 && criterion_type == 2) adjustment =  0.1672027; 
        if (n_levels ==  136 && criterion_type == 1) adjustment =  0.4512564; 
        if (n_levels ==  136 && criterion_type == 2) adjustment =  0.1671992; 
        if (n_levels ==  137 && criterion_type == 1) adjustment =  0.4512743; 
        if (n_levels ==  137 && criterion_type == 2) adjustment =  0.1671957; 
        if (n_levels ==  138 && criterion_type == 1) adjustment =  0.451292; 
        if (n_levels ==  138 && criterion_type == 2) adjustment =  0.1671923; 
        if (n_levels ==  139 && criterion_type == 1) adjustment =  0.4513094; 
        if (n_levels ==  139 && criterion_type == 2) adjustment =  0.1671889; 
        if (n_levels ==  140 && criterion_type == 1) adjustment =  0.4513266; 
        if (n_levels ==  140 && criterion_type == 2) adjustment =  0.1671855; 
        if (n_levels ==  141 && criterion_type == 1) adjustment =  0.4513435; 
        if (n_levels ==  141 && criterion_type == 2) adjustment =  0.1671823; 
        if (n_levels ==  142 && criterion_type == 1) adjustment =  0.4513601; 
        if (n_levels ==  142 && criterion_type == 2) adjustment =  0.167179; 
        if (n_levels ==  143 && criterion_type == 1) adjustment =  0.4513765; 
        if (n_levels ==  143 && criterion_type == 2) adjustment =  0.1671758; 
        if (n_levels ==  144 && criterion_type == 1) adjustment =  0.4513927; 
        if (n_levels ==  144 && criterion_type == 2) adjustment =  0.1671726; 
        if (n_levels ==  145 && criterion_type == 1) adjustment =  0.4514086; 
        if (n_levels ==  145 && criterion_type == 2) adjustment =  0.1671695; 
        if (n_levels ==  146 && criterion_type == 1) adjustment =  0.4514243; 
        if (n_levels ==  146 && criterion_type == 2) adjustment =  0.1671664; 
        if (n_levels ==  147 && criterion_type == 1) adjustment =  0.4514398; 
        if (n_levels ==  147 && criterion_type == 2) adjustment =  0.1671633; 
        if (n_levels ==  148 && criterion_type == 1) adjustment =  0.4514551; 
        if (n_levels ==  148 && criterion_type == 2) adjustment =  0.1671603; 
        if (n_levels ==  149 && criterion_type == 1) adjustment =  0.4514702; 
        if (n_levels ==  149 && criterion_type == 2) adjustment =  0.1671573; 
        if (n_levels ==  150 && criterion_type == 1) adjustment =  0.451485; 
        if (n_levels ==  150 && criterion_type == 2) adjustment =  0.1671544; 
        if (n_levels ==  151 && criterion_type == 1) adjustment =  0.4514997; 
        if (n_levels ==  151 && criterion_type == 2) adjustment =  0.1671515; 
        if (n_levels ==  152 && criterion_type == 1) adjustment =  0.4515141; 
        if (n_levels ==  152 && criterion_type == 2) adjustment =  0.1671486; 
        if (n_levels ==  153 && criterion_type == 1) adjustment =  0.4515284; 
        if (n_levels ==  153 && criterion_type == 2) adjustment =  0.1671458; 
        if (n_levels ==  154 && criterion_type == 1) adjustment =  0.4515424; 
        if (n_levels ==  154 && criterion_type == 2) adjustment =  0.167143; 
        if (n_levels ==  155 && criterion_type == 1) adjustment =  0.4515563; 
        if (n_levels ==  155 && criterion_type == 2) adjustment =  0.1671402; 
        if (n_levels ==  156 && criterion_type == 1) adjustment =  0.45157; 
        if (n_levels ==  156 && criterion_type == 2) adjustment =  0.1671375; 
        if (n_levels ==  157 && criterion_type == 1) adjustment =  0.4515835; 
        if (n_levels ==  157 && criterion_type == 2) adjustment =  0.1671348; 
        if (n_levels ==  158 && criterion_type == 1) adjustment =  0.4515968; 
        if (n_levels ==  158 && criterion_type == 2) adjustment =  0.1671321; 
        if (n_levels ==  159 && criterion_type == 1) adjustment =  0.45161; 
        if (n_levels ==  159 && criterion_type == 2) adjustment =  0.1671295; 
        if (n_levels ==  160 && criterion_type == 1) adjustment =  0.451623; 
        if (n_levels ==  160 && criterion_type == 2) adjustment =  0.1671268; 
        if (n_levels ==  161 && criterion_type == 1) adjustment =  0.4516358; 
        if (n_levels ==  161 && criterion_type == 2) adjustment =  0.1671243; 
        if (n_levels ==  162 && criterion_type == 1) adjustment =  0.4516484; 
        if (n_levels ==  162 && criterion_type == 2) adjustment =  0.1671217; 
        if (n_levels ==  163 && criterion_type == 1) adjustment =  0.4516609; 
        if (n_levels ==  163 && criterion_type == 2) adjustment =  0.1671192; 
        if (n_levels ==  164 && criterion_type == 1) adjustment =  0.4516733; 
        if (n_levels ==  164 && criterion_type == 2) adjustment =  0.1671167; 
        if (n_levels ==  165 && criterion_type == 1) adjustment =  0.4516854; 
        if (n_levels ==  165 && criterion_type == 2) adjustment =  0.1671142; 
        if (n_levels ==  166 && criterion_type == 1) adjustment =  0.4516975; 
        if (n_levels ==  166 && criterion_type == 2) adjustment =  0.1671118; 
        if (n_levels ==  167 && criterion_type == 1) adjustment =  0.4517094; 
        if (n_levels ==  167 && criterion_type == 2) adjustment =  0.1671094; 
        if (n_levels ==  168 && criterion_type == 1) adjustment =  0.4517211; 
        if (n_levels ==  168 && criterion_type == 2) adjustment =  0.167107; 
        if (n_levels ==  169 && criterion_type == 1) adjustment =  0.4517327; 
        if (n_levels ==  169 && criterion_type == 2) adjustment =  0.1671046; 
        if (n_levels ==  170 && criterion_type == 1) adjustment =  0.4517441; 
        if (n_levels ==  170 && criterion_type == 2) adjustment =  0.1671023; 
        if (n_levels ==  171 && criterion_type == 1) adjustment =  0.4517554; 
        if (n_levels ==  171 && criterion_type == 2) adjustment =  0.1671; 
        if (n_levels ==  172 && criterion_type == 1) adjustment =  0.4517666; 
        if (n_levels ==  172 && criterion_type == 2) adjustment =  0.1670977; 
        if (n_levels ==  173 && criterion_type == 1) adjustment =  0.4517776; 
        if (n_levels ==  173 && criterion_type == 2) adjustment =  0.1670954; 
        if (n_levels ==  174 && criterion_type == 1) adjustment =  0.4517885; 
        if (n_levels ==  174 && criterion_type == 2) adjustment =  0.1670932; 
        if (n_levels ==  175 && criterion_type == 1) adjustment =  0.4517993; 
        if (n_levels ==  175 && criterion_type == 2) adjustment =  0.167091; 
        if (n_levels ==  176 && criterion_type == 1) adjustment =  0.45181; 
        if (n_levels ==  176 && criterion_type == 2) adjustment =  0.1670888; 
        if (n_levels ==  177 && criterion_type == 1) adjustment =  0.4518205; 
        if (n_levels ==  177 && criterion_type == 2) adjustment =  0.1670866; 
        if (n_levels ==  178 && criterion_type == 1) adjustment =  0.4518309; 
        if (n_levels ==  178 && criterion_type == 2) adjustment =  0.1670845; 
        if (n_levels ==  179 && criterion_type == 1) adjustment =  0.4518412; 
        if (n_levels ==  179 && criterion_type == 2) adjustment =  0.1670824; 
        if (n_levels ==  180 && criterion_type == 1) adjustment =  0.4518514; 
        if (n_levels ==  180 && criterion_type == 2) adjustment =  0.1670803; 
        if (n_levels ==  181 && criterion_type == 1) adjustment =  0.4518614; 
        if (n_levels ==  181 && criterion_type == 2) adjustment =  0.1670782; 
        if (n_levels ==  182 && criterion_type == 1) adjustment =  0.4518713; 
        if (n_levels ==  182 && criterion_type == 2) adjustment =  0.1670761; 
        if (n_levels ==  183 && criterion_type == 1) adjustment =  0.4518812; 
        if (n_levels ==  183 && criterion_type == 2) adjustment =  0.1670741; 
        if (n_levels ==  184 && criterion_type == 1) adjustment =  0.4518909; 
        if (n_levels ==  184 && criterion_type == 2) adjustment =  0.1670721; 
        if (n_levels ==  185 && criterion_type == 1) adjustment =  0.4519005; 
        if (n_levels ==  185 && criterion_type == 2) adjustment =  0.1670701; 
        if (n_levels ==  186 && criterion_type == 1) adjustment =  0.45191; 
        if (n_levels ==  186 && criterion_type == 2) adjustment =  0.1670681; 
        if (n_levels ==  187 && criterion_type == 1) adjustment =  0.4519194; 
        if (n_levels ==  187 && criterion_type == 2) adjustment =  0.1670662; 
        if (n_levels ==  188 && criterion_type == 1) adjustment =  0.4519287; 
        if (n_levels ==  188 && criterion_type == 2) adjustment =  0.1670642; 
        if (n_levels ==  189 && criterion_type == 1) adjustment =  0.4519378; 
        if (n_levels ==  189 && criterion_type == 2) adjustment =  0.1670623; 
        if (n_levels ==  190 && criterion_type == 1) adjustment =  0.4519469; 
        if (n_levels ==  190 && criterion_type == 2) adjustment =  0.1670604; 
        if (n_levels ==  191 && criterion_type == 1) adjustment =  0.4519559; 
        if (n_levels ==  191 && criterion_type == 2) adjustment =  0.1670585; 
        if (n_levels ==  192 && criterion_type == 1) adjustment =  0.4519648; 
        if (n_levels ==  192 && criterion_type == 2) adjustment =  0.1670567; 
        if (n_levels ==  193 && criterion_type == 1) adjustment =  0.4519736; 
        if (n_levels ==  193 && criterion_type == 2) adjustment =  0.1670548; 
        if (n_levels ==  194 && criterion_type == 1) adjustment =  0.4519823; 
        if (n_levels ==  194 && criterion_type == 2) adjustment =  0.167053; 
        if (n_levels ==  195 && criterion_type == 1) adjustment =  0.4519909; 
        if (n_levels ==  195 && criterion_type == 2) adjustment =  0.1670512; 
        if (n_levels ==  196 && criterion_type == 1) adjustment =  0.4519994; 
        if (n_levels ==  196 && criterion_type == 2) adjustment =  0.1670494; 
        if (n_levels ==  197 && criterion_type == 1) adjustment =  0.4520079; 
        if (n_levels ==  197 && criterion_type == 2) adjustment =  0.1670476; 
        if (n_levels ==  198 && criterion_type == 1) adjustment =  0.4520162; 
        if (n_levels ==  198 && criterion_type == 2) adjustment =  0.1670459; 
        if (n_levels ==  199 && criterion_type == 1) adjustment =  0.4520245; 
        if (n_levels ==  199 && criterion_type == 2) adjustment =  0.1670441; 
        if (n_levels ==  200 && criterion_type == 1) adjustment =  0.4520326; 
        if (n_levels ==  200 && criterion_type == 2) adjustment =  0.1670424; 
        if (n_levels ==  201 && criterion_type == 1) adjustment =  0.4520407; 
        if (n_levels ==  201 && criterion_type == 2) adjustment =  0.1670407; 
        if (n_levels ==  202 && criterion_type == 1) adjustment =  0.4520487; 
        if (n_levels ==  202 && criterion_type == 2) adjustment =  0.167039; 
        if (n_levels ==  203 && criterion_type == 1) adjustment =  0.4520566; 
        if (n_levels ==  203 && criterion_type == 2) adjustment =  0.1670373; 
        if (n_levels ==  204 && criterion_type == 1) adjustment =  0.4520645; 
        if (n_levels ==  204 && criterion_type == 2) adjustment =  0.1670356; 
        if (n_levels ==  205 && criterion_type == 1) adjustment =  0.4520722; 
        if (n_levels ==  205 && criterion_type == 2) adjustment =  0.167034; 
        if (n_levels ==  206 && criterion_type == 1) adjustment =  0.4520799; 
        if (n_levels ==  206 && criterion_type == 2) adjustment =  0.1670324; 
        if (n_levels ==  207 && criterion_type == 1) adjustment =  0.4520875; 
        if (n_levels ==  207 && criterion_type == 2) adjustment =  0.1670307; 
        if (n_levels ==  208 && criterion_type == 1) adjustment =  0.4520951; 
        if (n_levels ==  208 && criterion_type == 2) adjustment =  0.1670291; 
        if (n_levels ==  209 && criterion_type == 1) adjustment =  0.4521025; 
        if (n_levels ==  209 && criterion_type == 2) adjustment =  0.1670276; 
        if (n_levels ==  210 && criterion_type == 1) adjustment =  0.4521099; 
        if (n_levels ==  210 && criterion_type == 2) adjustment =  0.167026; 
        if (n_levels ==  211 && criterion_type == 1) adjustment =  0.4521172; 
        if (n_levels ==  211 && criterion_type == 2) adjustment =  0.1670244; 
        if (n_levels ==  212 && criterion_type == 1) adjustment =  0.4521244; 
        if (n_levels ==  212 && criterion_type == 2) adjustment =  0.1670229; 
        if (n_levels ==  213 && criterion_type == 1) adjustment =  0.4521316; 
        if (n_levels ==  213 && criterion_type == 2) adjustment =  0.1670213; 
        if (n_levels ==  214 && criterion_type == 1) adjustment =  0.4521387; 
        if (n_levels ==  214 && criterion_type == 2) adjustment =  0.1670198; 
        if (n_levels ==  215 && criterion_type == 1) adjustment =  0.4521457; 
        if (n_levels ==  215 && criterion_type == 2) adjustment =  0.1670183; 
        if (n_levels ==  216 && criterion_type == 1) adjustment =  0.4521527; 
        if (n_levels ==  216 && criterion_type == 2) adjustment =  0.1670168; 
        if (n_levels ==  217 && criterion_type == 1) adjustment =  0.4521596; 
        if (n_levels ==  217 && criterion_type == 2) adjustment =  0.1670153; 
        if (n_levels ==  218 && criterion_type == 1) adjustment =  0.4521664; 
        if (n_levels ==  218 && criterion_type == 2) adjustment =  0.1670139; 
        if (n_levels ==  219 && criterion_type == 1) adjustment =  0.4521732; 
        if (n_levels ==  219 && criterion_type == 2) adjustment =  0.1670124; 
        if (n_levels ==  220 && criterion_type == 1) adjustment =  0.4521799; 
        if (n_levels ==  220 && criterion_type == 2) adjustment =  0.167011; 
        if (n_levels ==  221 && criterion_type == 1) adjustment =  0.4521866; 
        if (n_levels ==  221 && criterion_type == 2) adjustment =  0.1670095; 
        if (n_levels ==  222 && criterion_type == 1) adjustment =  0.4521931; 
        if (n_levels ==  222 && criterion_type == 2) adjustment =  0.1670081; 
        if (n_levels ==  223 && criterion_type == 1) adjustment =  0.4521997; 
        if (n_levels ==  223 && criterion_type == 2) adjustment =  0.1670067; 
        if (n_levels ==  224 && criterion_type == 1) adjustment =  0.4522061; 
        if (n_levels ==  224 && criterion_type == 2) adjustment =  0.1670053; 
        if (n_levels ==  225 && criterion_type == 1) adjustment =  0.4522125; 
        if (n_levels ==  225 && criterion_type == 2) adjustment =  0.1670039; 
        if (n_levels ==  226 && criterion_type == 1) adjustment =  0.4522189; 
        if (n_levels ==  226 && criterion_type == 2) adjustment =  0.1670026; 
        if (n_levels ==  227 && criterion_type == 1) adjustment =  0.4522251; 
        if (n_levels ==  227 && criterion_type == 2) adjustment =  0.1670012; 
        if (n_levels ==  228 && criterion_type == 1) adjustment =  0.4522314; 
        if (n_levels ==  228 && criterion_type == 2) adjustment =  0.1669998; 
        if (n_levels ==  229 && criterion_type == 1) adjustment =  0.4522375; 
        if (n_levels ==  229 && criterion_type == 2) adjustment =  0.1669985; 
        if (n_levels ==  230 && criterion_type == 1) adjustment =  0.4522436; 
        if (n_levels ==  230 && criterion_type == 2) adjustment =  0.1669972; 
        if (n_levels ==  231 && criterion_type == 1) adjustment =  0.4522497; 
        if (n_levels ==  231 && criterion_type == 2) adjustment =  0.1669959; 
        if (n_levels ==  232 && criterion_type == 1) adjustment =  0.4522557; 
        if (n_levels ==  232 && criterion_type == 2) adjustment =  0.1669945; 
        if (n_levels ==  233 && criterion_type == 1) adjustment =  0.4522617; 
        if (n_levels ==  233 && criterion_type == 2) adjustment =  0.1669933; 
        if (n_levels ==  234 && criterion_type == 1) adjustment =  0.4522676; 
        if (n_levels ==  234 && criterion_type == 2) adjustment =  0.166992; 
        if (n_levels ==  235 && criterion_type == 1) adjustment =  0.4522734; 
        if (n_levels ==  235 && criterion_type == 2) adjustment =  0.1669907; 
        if (n_levels ==  236 && criterion_type == 1) adjustment =  0.4522792; 
        if (n_levels ==  236 && criterion_type == 2) adjustment =  0.1669894; 
        if (n_levels ==  237 && criterion_type == 1) adjustment =  0.452285; 
        if (n_levels ==  237 && criterion_type == 2) adjustment =  0.1669882; 
        if (n_levels ==  238 && criterion_type == 1) adjustment =  0.4522907; 
        if (n_levels ==  238 && criterion_type == 2) adjustment =  0.1669869; 
        if (n_levels ==  239 && criterion_type == 1) adjustment =  0.4522963; 
        if (n_levels ==  239 && criterion_type == 2) adjustment =  0.1669857; 
        if (n_levels ==  240 && criterion_type == 1) adjustment =  0.4523019; 
        if (n_levels ==  240 && criterion_type == 2) adjustment =  0.1669845; 
        if (n_levels ==  241 && criterion_type == 1) adjustment =  0.4523075; 
        if (n_levels ==  241 && criterion_type == 2) adjustment =  0.1669832; 
        if (n_levels ==  242 && criterion_type == 1) adjustment =  0.452313; 
        if (n_levels ==  242 && criterion_type == 2) adjustment =  0.166982; 
        if (n_levels ==  243 && criterion_type == 1) adjustment =  0.4523184; 
        if (n_levels ==  243 && criterion_type == 2) adjustment =  0.1669808; 
        if (n_levels ==  244 && criterion_type == 1) adjustment =  0.4523238; 
        if (n_levels ==  244 && criterion_type == 2) adjustment =  0.1669797; 
        if (n_levels ==  245 && criterion_type == 1) adjustment =  0.4523292; 
        if (n_levels ==  245 && criterion_type == 2) adjustment =  0.1669785; 
        if (n_levels ==  246 && criterion_type == 1) adjustment =  0.4523345; 
        if (n_levels ==  246 && criterion_type == 2) adjustment =  0.1669773; 
        if (n_levels ==  247 && criterion_type == 1) adjustment =  0.4523398; 
        if (n_levels ==  247 && criterion_type == 2) adjustment =  0.1669761; 
        if (n_levels ==  248 && criterion_type == 1) adjustment =  0.452345; 
        if (n_levels ==  248 && criterion_type == 2) adjustment =  0.166975; 
        if (n_levels ==  249 && criterion_type == 1) adjustment =  0.4523502; 
        if (n_levels ==  249 && criterion_type == 2) adjustment =  0.1669738; 
        if (n_levels ==  250 && criterion_type == 1) adjustment =  0.4523554; 
        if (n_levels ==  250 && criterion_type == 2) adjustment =  0.1669727; 
        if (n_levels ==  251 && criterion_type == 1) adjustment =  0.4523605; 
        if (n_levels ==  251 && criterion_type == 2) adjustment =  0.1669716; 
        if (n_levels ==  252 && criterion_type == 1) adjustment =  0.4523655; 
        if (n_levels ==  252 && criterion_type == 2) adjustment =  0.1669705; 
        if (n_levels ==  253 && criterion_type == 1) adjustment =  0.4523705; 
        if (n_levels ==  253 && criterion_type == 2) adjustment =  0.1669694; 
        if (n_levels ==  254 && criterion_type == 1) adjustment =  0.4523755; 
        if (n_levels ==  254 && criterion_type == 2) adjustment =  0.1669683; 
        if (n_levels ==  255 && criterion_type == 1) adjustment =  0.4523805; 
        if (n_levels ==  255 && criterion_type == 2) adjustment =  0.1669672; 
        if (n_levels ==  256 && criterion_type == 1) adjustment =  0.4523854; 
        if (n_levels ==  256 && criterion_type == 2) adjustment =  0.1669661; 
        if (n_levels ==  257 && criterion_type == 1) adjustment =  0.4523902; 
        if (n_levels ==  257 && criterion_type == 2) adjustment =  0.166965; 
        if (n_levels ==  258 && criterion_type == 1) adjustment =  0.452395; 
        if (n_levels ==  258 && criterion_type == 2) adjustment =  0.1669639; 
        if (n_levels ==  259 && criterion_type == 1) adjustment =  0.4523998; 
        if (n_levels ==  259 && criterion_type == 2) adjustment =  0.1669629; 
        if (n_levels ==  260 && criterion_type == 1) adjustment =  0.4524046; 
        if (n_levels ==  260 && criterion_type == 2) adjustment =  0.1669618; 
        if (n_levels ==  261 && criterion_type == 1) adjustment =  0.4524093; 
        if (n_levels ==  261 && criterion_type == 2) adjustment =  0.1669608; 
        if (n_levels ==  262 && criterion_type == 1) adjustment =  0.4524139; 
        if (n_levels ==  262 && criterion_type == 2) adjustment =  0.1669597; 
        if (n_levels ==  263 && criterion_type == 1) adjustment =  0.4524186; 
        if (n_levels ==  263 && criterion_type == 2) adjustment =  0.1669587; 
        if (n_levels ==  264 && criterion_type == 1) adjustment =  0.4524232; 
        if (n_levels ==  264 && criterion_type == 2) adjustment =  0.1669577; 
        if (n_levels ==  265 && criterion_type == 1) adjustment =  0.4524277; 
        if (n_levels ==  265 && criterion_type == 2) adjustment =  0.1669567; 
        if (n_levels ==  266 && criterion_type == 1) adjustment =  0.4524323; 
        if (n_levels ==  266 && criterion_type == 2) adjustment =  0.1669556; 
        if (n_levels ==  267 && criterion_type == 1) adjustment =  0.4524368; 
        if (n_levels ==  267 && criterion_type == 2) adjustment =  0.1669546; 
        if (n_levels ==  268 && criterion_type == 1) adjustment =  0.4524412; 
        if (n_levels ==  268 && criterion_type == 2) adjustment =  0.1669536; 
        if (n_levels ==  269 && criterion_type == 1) adjustment =  0.4524456; 
        if (n_levels ==  269 && criterion_type == 2) adjustment =  0.1669527; 
        if (n_levels ==  270 && criterion_type == 1) adjustment =  0.45245; 
        if (n_levels ==  270 && criterion_type == 2) adjustment =  0.1669517; 
        if (n_levels ==  271 && criterion_type == 1) adjustment =  0.4524544; 
        if (n_levels ==  271 && criterion_type == 2) adjustment =  0.1669507; 
        if (n_levels ==  272 && criterion_type == 1) adjustment =  0.4524587; 
        if (n_levels ==  272 && criterion_type == 2) adjustment =  0.1669497; 
        if (n_levels ==  273 && criterion_type == 1) adjustment =  0.452463; 
        if (n_levels ==  273 && criterion_type == 2) adjustment =  0.1669488; 
        if (n_levels ==  274 && criterion_type == 1) adjustment =  0.4524672; 
        if (n_levels ==  274 && criterion_type == 2) adjustment =  0.1669478; 
        if (n_levels ==  275 && criterion_type == 1) adjustment =  0.4524715; 
        if (n_levels ==  275 && criterion_type == 2) adjustment =  0.1669469; 
        if (n_levels ==  276 && criterion_type == 1) adjustment =  0.4524757; 
        if (n_levels ==  276 && criterion_type == 2) adjustment =  0.1669459; 
        if (n_levels ==  277 && criterion_type == 1) adjustment =  0.4524798; 
        if (n_levels ==  277 && criterion_type == 2) adjustment =  0.166945; 
        if (n_levels ==  278 && criterion_type == 1) adjustment =  0.452484; 
        if (n_levels ==  278 && criterion_type == 2) adjustment =  0.1669441; 
        if (n_levels ==  279 && criterion_type == 1) adjustment =  0.4524881; 
        if (n_levels ==  279 && criterion_type == 2) adjustment =  0.1669431; 
        if (n_levels ==  280 && criterion_type == 1) adjustment =  0.4524921; 
        if (n_levels ==  280 && criterion_type == 2) adjustment =  0.1669422; 
        if (n_levels ==  281 && criterion_type == 1) adjustment =  0.4524962; 
        if (n_levels ==  281 && criterion_type == 2) adjustment =  0.1669413; 
        if (n_levels ==  282 && criterion_type == 1) adjustment =  0.4525002; 
        if (n_levels ==  282 && criterion_type == 2) adjustment =  0.1669404; 
        if (n_levels ==  283 && criterion_type == 1) adjustment =  0.4525042; 
        if (n_levels ==  283 && criterion_type == 2) adjustment =  0.1669395; 
        if (n_levels ==  284 && criterion_type == 1) adjustment =  0.4525081; 
        if (n_levels ==  284 && criterion_type == 2) adjustment =  0.1669386; 
        if (n_levels ==  285 && criterion_type == 1) adjustment =  0.452512; 
        if (n_levels ==  285 && criterion_type == 2) adjustment =  0.1669377; 
        if (n_levels ==  286 && criterion_type == 1) adjustment =  0.4525159; 
        if (n_levels ==  286 && criterion_type == 2) adjustment =  0.1669368; 
        if (n_levels ==  287 && criterion_type == 1) adjustment =  0.4525198; 
        if (n_levels ==  287 && criterion_type == 2) adjustment =  0.166936; 
        if (n_levels ==  288 && criterion_type == 1) adjustment =  0.4525236; 
        if (n_levels ==  288 && criterion_type == 2) adjustment =  0.1669351; 
        if (n_levels ==  289 && criterion_type == 1) adjustment =  0.4525275; 
        if (n_levels ==  289 && criterion_type == 2) adjustment =  0.1669342; 
        if (n_levels ==  290 && criterion_type == 1) adjustment =  0.4525312; 
        if (n_levels ==  290 && criterion_type == 2) adjustment =  0.1669334; 
        if (n_levels ==  291 && criterion_type == 1) adjustment =  0.452535; 
        if (n_levels ==  291 && criterion_type == 2) adjustment =  0.1669325; 
        if (n_levels ==  292 && criterion_type == 1) adjustment =  0.4525387; 
        if (n_levels ==  292 && criterion_type == 2) adjustment =  0.1669317; 
        if (n_levels ==  293 && criterion_type == 1) adjustment =  0.4525424; 
        if (n_levels ==  293 && criterion_type == 2) adjustment =  0.1669308; 
        if (n_levels ==  294 && criterion_type == 1) adjustment =  0.4525461; 
        if (n_levels ==  294 && criterion_type == 2) adjustment =  0.16693; 
        if (n_levels ==  295 && criterion_type == 1) adjustment =  0.4525498; 
        if (n_levels ==  295 && criterion_type == 2) adjustment =  0.1669292; 
        if (n_levels ==  296 && criterion_type == 1) adjustment =  0.4525534; 
        if (n_levels ==  296 && criterion_type == 2) adjustment =  0.1669283; 
        if (n_levels ==  297 && criterion_type == 1) adjustment =  0.452557; 
        if (n_levels ==  297 && criterion_type == 2) adjustment =  0.1669275; 
        if (n_levels ==  298 && criterion_type == 1) adjustment =  0.4525606; 
        if (n_levels ==  298 && criterion_type == 2) adjustment =  0.1669267; 
        if (n_levels ==  299 && criterion_type == 1) adjustment =  0.4525642; 
        if (n_levels ==  299 && criterion_type == 2) adjustment =  0.1669259; 
        if (n_levels ==  300 && criterion_type == 1) adjustment =  0.4525677; 
        if (n_levels ==  300 && criterion_type == 2) adjustment =  0.1669251; 
        if (n_levels ==  301 && criterion_type == 1) adjustment =  0.4525712; 
        if (n_levels ==  301 && criterion_type == 2) adjustment =  0.1669243; 
        if (n_levels ==  302 && criterion_type == 1) adjustment =  0.4525747; 
        if (n_levels ==  302 && criterion_type == 2) adjustment =  0.1669235; 
        if (n_levels ==  303 && criterion_type == 1) adjustment =  0.4525781; 
        if (n_levels ==  303 && criterion_type == 2) adjustment =  0.1669227; 
        if (n_levels ==  304 && criterion_type == 1) adjustment =  0.4525816; 
        if (n_levels ==  304 && criterion_type == 2) adjustment =  0.1669219; 
        if (n_levels ==  305 && criterion_type == 1) adjustment =  0.452585; 
        if (n_levels ==  305 && criterion_type == 2) adjustment =  0.1669211; 
        if (n_levels ==  306 && criterion_type == 1) adjustment =  0.4525884; 
        if (n_levels ==  306 && criterion_type == 2) adjustment =  0.1669204; 
        if (n_levels ==  307 && criterion_type == 1) adjustment =  0.4525917; 
        if (n_levels ==  307 && criterion_type == 2) adjustment =  0.1669196; 
        if (n_levels ==  308 && criterion_type == 1) adjustment =  0.4525951; 
        if (n_levels ==  308 && criterion_type == 2) adjustment =  0.1669188; 
        if (n_levels ==  309 && criterion_type == 1) adjustment =  0.4525984; 
        if (n_levels ==  309 && criterion_type == 2) adjustment =  0.1669181; 
        if (n_levels ==  310 && criterion_type == 1) adjustment =  0.4526017; 
        if (n_levels ==  310 && criterion_type == 2) adjustment =  0.1669173; 
        if (n_levels ==  311 && criterion_type == 1) adjustment =  0.452605; 
        if (n_levels ==  311 && criterion_type == 2) adjustment =  0.1669166; 
        if (n_levels ==  312 && criterion_type == 1) adjustment =  0.4526082; 
        if (n_levels ==  312 && criterion_type == 2) adjustment =  0.1669158; 
        if (n_levels ==  313 && criterion_type == 1) adjustment =  0.4526115; 
        if (n_levels ==  313 && criterion_type == 2) adjustment =  0.1669151; 
        if (n_levels ==  314 && criterion_type == 1) adjustment =  0.4526147; 
        if (n_levels ==  314 && criterion_type == 2) adjustment =  0.1669143; 
        if (n_levels ==  315 && criterion_type == 1) adjustment =  0.4526179; 
        if (n_levels ==  315 && criterion_type == 2) adjustment =  0.1669136; 
        if (n_levels ==  316 && criterion_type == 1) adjustment =  0.452621; 
        if (n_levels ==  316 && criterion_type == 2) adjustment =  0.1669129; 
        if (n_levels ==  317 && criterion_type == 1) adjustment =  0.4526242; 
        if (n_levels ==  317 && criterion_type == 2) adjustment =  0.1669122; 
        if (n_levels ==  318 && criterion_type == 1) adjustment =  0.4526273; 
        if (n_levels ==  318 && criterion_type == 2) adjustment =  0.1669114; 
        if (n_levels ==  319 && criterion_type == 1) adjustment =  0.4526304; 
        if (n_levels ==  319 && criterion_type == 2) adjustment =  0.1669107; 
        if (n_levels ==  320 && criterion_type == 1) adjustment =  0.4526335; 
        if (n_levels ==  320 && criterion_type == 2) adjustment =  0.16691; 
        if (n_levels ==  321 && criterion_type == 1) adjustment =  0.4526366; 
        if (n_levels ==  321 && criterion_type == 2) adjustment =  0.1669093; 
        if (n_levels ==  322 && criterion_type == 1) adjustment =  0.4526396; 
        if (n_levels ==  322 && criterion_type == 2) adjustment =  0.1669086; 
        if (n_levels ==  323 && criterion_type == 1) adjustment =  0.4526427; 
        if (n_levels ==  323 && criterion_type == 2) adjustment =  0.1669079; 
        if (n_levels ==  324 && criterion_type == 1) adjustment =  0.4526457; 
        if (n_levels ==  324 && criterion_type == 2) adjustment =  0.1669072; 
        if (n_levels ==  325 && criterion_type == 1) adjustment =  0.4526487; 
        if (n_levels ==  325 && criterion_type == 2) adjustment =  0.1669065; 
        if (n_levels ==  326 && criterion_type == 1) adjustment =  0.4526517; 
        if (n_levels ==  326 && criterion_type == 2) adjustment =  0.1669058; 
        if (n_levels ==  327 && criterion_type == 1) adjustment =  0.4526546; 
        if (n_levels ==  327 && criterion_type == 2) adjustment =  0.1669051; 
        if (n_levels ==  328 && criterion_type == 1) adjustment =  0.4526575; 
        if (n_levels ==  328 && criterion_type == 2) adjustment =  0.1669045; 
        if (n_levels ==  329 && criterion_type == 1) adjustment =  0.4526605; 
        if (n_levels ==  329 && criterion_type == 2) adjustment =  0.1669038; 
        if (n_levels ==  330 && criterion_type == 1) adjustment =  0.4526634; 
        if (n_levels ==  330 && criterion_type == 2) adjustment =  0.1669031; 
        if (n_levels ==  331 && criterion_type == 1) adjustment =  0.4526662; 
        if (n_levels ==  331 && criterion_type == 2) adjustment =  0.1669024; 
        if (n_levels ==  332 && criterion_type == 1) adjustment =  0.4526691; 
        if (n_levels ==  332 && criterion_type == 2) adjustment =  0.1669018; 
        if (n_levels ==  333 && criterion_type == 1) adjustment =  0.452672; 
        if (n_levels ==  333 && criterion_type == 2) adjustment =  0.1669011; 
        if (n_levels ==  334 && criterion_type == 1) adjustment =  0.4526748; 
        if (n_levels ==  334 && criterion_type == 2) adjustment =  0.1669005; 
        if (n_levels ==  335 && criterion_type == 1) adjustment =  0.4526776; 
        if (n_levels ==  335 && criterion_type == 2) adjustment =  0.1668998; 
        if (n_levels ==  336 && criterion_type == 1) adjustment =  0.4526804; 
        if (n_levels ==  336 && criterion_type == 2) adjustment =  0.1668992; 
        if (n_levels ==  337 && criterion_type == 1) adjustment =  0.4526832; 
        if (n_levels ==  337 && criterion_type == 2) adjustment =  0.1668985; 
        if (n_levels ==  338 && criterion_type == 1) adjustment =  0.4526859; 
        if (n_levels ==  338 && criterion_type == 2) adjustment =  0.1668979; 
        if (n_levels ==  339 && criterion_type == 1) adjustment =  0.4526887; 
        if (n_levels ==  339 && criterion_type == 2) adjustment =  0.1668972; 
        if (n_levels ==  340 && criterion_type == 1) adjustment =  0.4526914; 
        if (n_levels ==  340 && criterion_type == 2) adjustment =  0.1668966; 
        if (n_levels ==  341 && criterion_type == 1) adjustment =  0.4526941; 
        if (n_levels ==  341 && criterion_type == 2) adjustment =  0.166896; 
        if (n_levels ==  342 && criterion_type == 1) adjustment =  0.4526968; 
        if (n_levels ==  342 && criterion_type == 2) adjustment =  0.1668953; 
        if (n_levels ==  343 && criterion_type == 1) adjustment =  0.4526995; 
        if (n_levels ==  343 && criterion_type == 2) adjustment =  0.1668947; 
        if (n_levels ==  344 && criterion_type == 1) adjustment =  0.4527022; 
        if (n_levels ==  344 && criterion_type == 2) adjustment =  0.1668941; 
        if (n_levels ==  345 && criterion_type == 1) adjustment =  0.4527048; 
        if (n_levels ==  345 && criterion_type == 2) adjustment =  0.1668935; 
        if (n_levels ==  346 && criterion_type == 1) adjustment =  0.4527074; 
        if (n_levels ==  346 && criterion_type == 2) adjustment =  0.1668929; 
        if (n_levels ==  347 && criterion_type == 1) adjustment =  0.45271; 
        if (n_levels ==  347 && criterion_type == 2) adjustment =  0.1668923; 
        if (n_levels ==  348 && criterion_type == 1) adjustment =  0.4527126; 
        if (n_levels ==  348 && criterion_type == 2) adjustment =  0.1668917; 
        if (n_levels ==  349 && criterion_type == 1) adjustment =  0.4527152; 
        if (n_levels ==  349 && criterion_type == 2) adjustment =  0.166891; 
        if (n_levels ==  350 && criterion_type == 1) adjustment =  0.4527178; 
        if (n_levels ==  350 && criterion_type == 2) adjustment =  0.1668904; 
        if (n_levels ==  351 && criterion_type == 1) adjustment =  0.4527203; 
        if (n_levels ==  351 && criterion_type == 2) adjustment =  0.1668899; 
        if (n_levels ==  352 && criterion_type == 1) adjustment =  0.4527229; 
        if (n_levels ==  352 && criterion_type == 2) adjustment =  0.1668893; 
        if (n_levels ==  353 && criterion_type == 1) adjustment =  0.4527254; 
        if (n_levels ==  353 && criterion_type == 2) adjustment =  0.1668887; 
        if (n_levels ==  354 && criterion_type == 1) adjustment =  0.4527279; 
        if (n_levels ==  354 && criterion_type == 2) adjustment =  0.1668881; 
        if (n_levels ==  355 && criterion_type == 1) adjustment =  0.4527304; 
        if (n_levels ==  355 && criterion_type == 2) adjustment =  0.1668875; 
        if (n_levels ==  356 && criterion_type == 1) adjustment =  0.4527329; 
        if (n_levels ==  356 && criterion_type == 2) adjustment =  0.1668869; 
        if (n_levels ==  357 && criterion_type == 1) adjustment =  0.4527354; 
        if (n_levels ==  357 && criterion_type == 2) adjustment =  0.1668863; 
        if (n_levels ==  358 && criterion_type == 1) adjustment =  0.4527378; 
        if (n_levels ==  358 && criterion_type == 2) adjustment =  0.1668858; 
        if (n_levels ==  359 && criterion_type == 1) adjustment =  0.4527403; 
        if (n_levels ==  359 && criterion_type == 2) adjustment =  0.1668852; 
        if (n_levels ==  360 && criterion_type == 1) adjustment =  0.4527427; 
        if (n_levels ==  360 && criterion_type == 2) adjustment =  0.1668846; 
        if (n_levels ==  361 && criterion_type == 1) adjustment =  0.4527451; 
        if (n_levels ==  361 && criterion_type == 2) adjustment =  0.1668841; 
        if (n_levels ==  362 && criterion_type == 1) adjustment =  0.4527475; 
        if (n_levels ==  362 && criterion_type == 2) adjustment =  0.1668835; 
        if (n_levels ==  363 && criterion_type == 1) adjustment =  0.4527499; 
        if (n_levels ==  363 && criterion_type == 2) adjustment =  0.1668829; 
        if (n_levels ==  364 && criterion_type == 1) adjustment =  0.4527522; 
        if (n_levels ==  364 && criterion_type == 2) adjustment =  0.1668824; 
        if (n_levels ==  365 && criterion_type == 1) adjustment =  0.4527546; 
        if (n_levels ==  365 && criterion_type == 2) adjustment =  0.1668818; 
        if (n_levels ==  366 && criterion_type == 1) adjustment =  0.4527569; 
        if (n_levels ==  366 && criterion_type == 2) adjustment =  0.1668813; 
        if (n_levels ==  367 && criterion_type == 1) adjustment =  0.4527593; 
        if (n_levels ==  367 && criterion_type == 2) adjustment =  0.1668807; 
        if (n_levels ==  368 && criterion_type == 1) adjustment =  0.4527616; 
        if (n_levels ==  368 && criterion_type == 2) adjustment =  0.1668802; 
        if (n_levels ==  369 && criterion_type == 1) adjustment =  0.4527639; 
        if (n_levels ==  369 && criterion_type == 2) adjustment =  0.1668796; 
        if (n_levels ==  370 && criterion_type == 1) adjustment =  0.4527662; 
        if (n_levels ==  370 && criterion_type == 2) adjustment =  0.1668791; 
        if (n_levels ==  371 && criterion_type == 1) adjustment =  0.4527685; 
        if (n_levels ==  371 && criterion_type == 2) adjustment =  0.1668786; 
        if (n_levels ==  372 && criterion_type == 1) adjustment =  0.4527707; 
        if (n_levels ==  372 && criterion_type == 2) adjustment =  0.166878; 
        if (n_levels ==  373 && criterion_type == 1) adjustment =  0.452773; 
        if (n_levels ==  373 && criterion_type == 2) adjustment =  0.1668775; 
        if (n_levels ==  374 && criterion_type == 1) adjustment =  0.4527752; 
        if (n_levels ==  374 && criterion_type == 2) adjustment =  0.166877; 
        if (n_levels ==  375 && criterion_type == 1) adjustment =  0.4527774; 
        if (n_levels ==  375 && criterion_type == 2) adjustment =  0.1668764; 
        if (n_levels ==  376 && criterion_type == 1) adjustment =  0.4527797; 
        if (n_levels ==  376 && criterion_type == 2) adjustment =  0.1668759; 
        if (n_levels ==  377 && criterion_type == 1) adjustment =  0.4527819; 
        if (n_levels ==  377 && criterion_type == 2) adjustment =  0.1668754; 
        if (n_levels ==  378 && criterion_type == 1) adjustment =  0.4527841; 
        if (n_levels ==  378 && criterion_type == 2) adjustment =  0.1668749; 
        if (n_levels ==  379 && criterion_type == 1) adjustment =  0.4527862; 
        if (n_levels ==  379 && criterion_type == 2) adjustment =  0.1668744; 
        if (n_levels ==  380 && criterion_type == 1) adjustment =  0.4527884; 
        if (n_levels ==  380 && criterion_type == 2) adjustment =  0.1668738; 
        if (n_levels ==  381 && criterion_type == 1) adjustment =  0.4527906; 
        if (n_levels ==  381 && criterion_type == 2) adjustment =  0.1668733; 
        if (n_levels ==  382 && criterion_type == 1) adjustment =  0.4527927; 
        if (n_levels ==  382 && criterion_type == 2) adjustment =  0.1668728; 
        if (n_levels ==  383 && criterion_type == 1) adjustment =  0.4527948; 
        if (n_levels ==  383 && criterion_type == 2) adjustment =  0.1668723; 
        if (n_levels ==  384 && criterion_type == 1) adjustment =  0.452797; 
        if (n_levels ==  384 && criterion_type == 2) adjustment =  0.1668718; 
        if (n_levels ==  385 && criterion_type == 1) adjustment =  0.4527991; 
        if (n_levels ==  385 && criterion_type == 2) adjustment =  0.1668713; 
        if (n_levels ==  386 && criterion_type == 1) adjustment =  0.4528012; 
        if (n_levels ==  386 && criterion_type == 2) adjustment =  0.1668708; 
        if (n_levels ==  387 && criterion_type == 1) adjustment =  0.4528033; 
        if (n_levels ==  387 && criterion_type == 2) adjustment =  0.1668703; 
        if (n_levels ==  388 && criterion_type == 1) adjustment =  0.4528053; 
        if (n_levels ==  388 && criterion_type == 2) adjustment =  0.1668698; 
        if (n_levels ==  389 && criterion_type == 1) adjustment =  0.4528074; 
        if (n_levels ==  389 && criterion_type == 2) adjustment =  0.1668693; 
        if (n_levels ==  390 && criterion_type == 1) adjustment =  0.4528095; 
        if (n_levels ==  390 && criterion_type == 2) adjustment =  0.1668688; 
        if (n_levels ==  391 && criterion_type == 1) adjustment =  0.4528115; 
        if (n_levels ==  391 && criterion_type == 2) adjustment =  0.1668684; 
        if (n_levels ==  392 && criterion_type == 1) adjustment =  0.4528135; 
        if (n_levels ==  392 && criterion_type == 2) adjustment =  0.1668679; 
        if (n_levels ==  393 && criterion_type == 1) adjustment =  0.4528156; 
        if (n_levels ==  393 && criterion_type == 2) adjustment =  0.1668674; 
        if (n_levels ==  394 && criterion_type == 1) adjustment =  0.4528176; 
        if (n_levels ==  394 && criterion_type == 2) adjustment =  0.1668669; 
        if (n_levels ==  395 && criterion_type == 1) adjustment =  0.4528196; 
        if (n_levels ==  395 && criterion_type == 2) adjustment =  0.1668664; 
        if (n_levels ==  396 && criterion_type == 1) adjustment =  0.4528216; 
        if (n_levels ==  396 && criterion_type == 2) adjustment =  0.166866; 
        if (n_levels ==  397 && criterion_type == 1) adjustment =  0.4528236; 
        if (n_levels ==  397 && criterion_type == 2) adjustment =  0.1668655; 
        if (n_levels ==  398 && criterion_type == 1) adjustment =  0.4528255; 
        if (n_levels ==  398 && criterion_type == 2) adjustment =  0.166865; 
        if (n_levels ==  399 && criterion_type == 1) adjustment =  0.4528275; 
        if (n_levels ==  399 && criterion_type == 2) adjustment =  0.1668646; 
        if (n_levels ==  400 && criterion_type == 1) adjustment =  0.4528294; 
        if (n_levels ==  400 && criterion_type == 2) adjustment =  0.1668641; 

    }
      
    // Nominal biomarker
    if (biomarker_type == 2) {

        if (n_unique > 20) {
            n_levels = 20;
        }
        else {
           n_levels = n_unique;
        }

        if (n_levels ==  2 && criterion_type == 1) adjustment =  0; 
        if (n_levels ==  2 && criterion_type == 2) adjustment =  0; 
        if (n_levels ==  3 && criterion_type == 1) adjustment =  0.1864386; 
        if (n_levels ==  3 && criterion_type == 2) adjustment =  0.1914214; 
        if (n_levels ==  4 && criterion_type == 1) adjustment =  0.1751265; 
        if (n_levels ==  4 && criterion_type == 2) adjustment =  0.2203869; 
        if (n_levels ==  5 && criterion_type == 1) adjustment =  0.1598419; 
        if (n_levels ==  5 && criterion_type == 2) adjustment =  0.2308939; 
        if (n_levels ==  6 && criterion_type == 1) adjustment =  0.1438832; 
        if (n_levels ==  6 && criterion_type == 2) adjustment =  0.2359842; 
        if (n_levels ==  7 && criterion_type == 1) adjustment =  0.1288371; 
        if (n_levels ==  7 && criterion_type == 2) adjustment =  0.2388894; 
        if (n_levels ==  8 && criterion_type == 1) adjustment =  0.1153799; 
        if (n_levels ==  8 && criterion_type == 2) adjustment =  0.2407454; 
        if (n_levels ==  9 && criterion_type == 1) adjustment =  0.1036859; 
        if (n_levels ==  9 && criterion_type == 2) adjustment =  0.2420337; 
        if (n_levels ==  10 && criterion_type == 1) adjustment =  0.09367185; 
        if (n_levels ==  10 && criterion_type == 2) adjustment =  0.2429851; 
        if (n_levels ==  11 && criterion_type == 1) adjustment =  0.0851444; 
        if (n_levels ==  11 && criterion_type == 2) adjustment =  0.2437209; 
        if (n_levels ==  12 && criterion_type == 1) adjustment =  0.07788102; 
        if (n_levels ==  12 && criterion_type == 2) adjustment =  0.24431; 
        if (n_levels ==  13 && criterion_type == 1) adjustment =  0.07167023; 
        if (n_levels ==  13 && criterion_type == 2) adjustment =  0.2447942; 
        if (n_levels ==  14 && criterion_type == 1) adjustment =  0.06632762; 
        if (n_levels ==  14 && criterion_type == 2) adjustment =  0.2452002; 
        if (n_levels ==  15 && criterion_type == 1) adjustment =  0.06169949; 
        if (n_levels ==  15 && criterion_type == 2) adjustment =  0.2455463; 
        if (n_levels ==  16 && criterion_type == 1) adjustment =  0.05766069; 
        if (n_levels ==  16 && criterion_type == 2) adjustment =  0.2458451; 
        if (n_levels ==  17 && criterion_type == 1) adjustment =  0.05411047; 
        if (n_levels ==  17 && criterion_type == 2) adjustment =  0.2461058; 
        if (n_levels ==  18 && criterion_type == 1) adjustment =  0.050968; 
        if (n_levels ==  18 && criterion_type == 2) adjustment =  0.2463355; 
        if (n_levels ==  19 && criterion_type == 1) adjustment =  0.04816834; 
        if (n_levels ==  19 && criterion_type == 2) adjustment =  0.2465394; 
        if (n_levels ==  20 && criterion_type == 1) adjustment =  0.0456591; 
        if (n_levels ==  20 && criterion_type == 2) adjustment =  0.2467217; 


    }

    return adjustment;

}

// Treatment effect test for binary outcome variable (comparison of proportions test)
double PropTestStatistic(const double &treatment_event_count, const double &control_event_count,
                         const double &treatment_count, const double &control_count,
                         const int &direction, int &error_flag) {

    double pooled_prop, variance;
    double statistic = 0;

    // Set the flag to error
    error_flag = 1;

    if (treatment_count >1 && control_count >1) {
        pooled_prop = (treatment_event_count + control_event_count)/(treatment_count + control_count);
        variance = pooled_prop * (1 - pooled_prop) * (treatment_count + control_count) / (treatment_count * control_count);
        if (variance > 1.0e-10) {
            statistic = direction * (treatment_event_count / treatment_count - control_event_count / control_count) / sqrt(variance);
            // No error
            error_flag = 0;
        }
    }
    return statistic;
}


// Treatment effect test in overall group of patients when the outcome variable is binary
SingleSubgroup BinOutOverallAnalysis(const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const int &direction) {

    SingleSubgroup res;

    int i;
    int n = treatment.size();

    double treatment_count = 0; 
    double control_count = 0;

    double treatment_event_count = 0;
    double control_event_count = 0;

    double test_statistic;
    int error_flag;

    for(i = 0; i < n; i++) {
        if (treatment[i] == 0) {
            control_event_count += outcome[i];
            control_count ++;
        }
        if (treatment[i] == 1) {
            treatment_event_count += outcome[i];
            treatment_count ++;
        }
    }

    test_statistic = PropTestStatistic(treatment_event_count, control_event_count, 
                                  treatment_count, control_count, direction, error_flag);

    double pvalue = 1.0 - rcpp_pnorm(test_statistic);

    // Save the results

    // Splitting criterion
    res.criterion = 0;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = 0;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = 0;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    res.adjusted_pvalue = -1.0;    
    res.prom_estimate = (treatment_event_count/treatment_count-control_event_count/control_count);
    res.prom_sderr = -1.0;   
    res.prom_sd = -1.0;   
    res.size_control = control_count;
    res.size_treatment = treatment_count;    

    // Vector of biomarker values to define the current patient subgroup
    vector<double> value;
    value.clear();
    value.push_back(0);
    res.value = value; 
    // 1 if <=, 2 if >, 3 if =
    res.sign = 0;
    // Size of the current patient subgroup
    res.size = control_count + treatment_count;
    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows;
    subgroup_rows.clear();
    subgroup_rows.push_back(0);
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup 
    res.biomarker_index = 0;
    // Level of the current patient subgroup  
    res.level = 0;
    // Parent index for the current patient subgroup 
    res.parent_index = -1;
    // Indexes of child subgroups for the current patient subgroup 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(0);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors)  
    res.code = 0;
    // Is the current patient subgroup terminal (0/1) 
    res.terminal_subgroup = 0;

    return res;    

}


// Compute parameters used in local multiplicity adjustment
double AdjParameter(const int &biomarker_type, const int &n_unique, const int &criterion_type) {

    double adjustment, adjustment1, adjustment2;
    double ave_rho1, ave_rho2;

    if (criterion_type == 1 || criterion_type == 2) {
        ave_rho1 = AdjParameterCriteria(biomarker_type, n_unique, criterion_type);
        if (biomarker_type == 1) adjustment = pow(n_unique - 1, 1 - ave_rho1);
        if (biomarker_type == 2) adjustment = pow(pow(2, n_unique-1) - 1, 1 - ave_rho1);
    }

    if (criterion_type == 3) {
        ave_rho1 = AdjParameterCriteria(biomarker_type, n_unique, 1);
        ave_rho2 = AdjParameterCriteria(biomarker_type, n_unique, 2);

        if (biomarker_type == 1) {
            adjustment1 = pow(n_unique - 1, 1 - ave_rho1);
            adjustment2 = pow(n_unique - 1, 1 - ave_rho2);
        }

        if (biomarker_type == 2) {
            adjustment1 = pow(pow(2, n_unique-1) - 1, 1 - ave_rho1);
            adjustment2 = pow(pow(2, n_unique-1) - 1, 1 - ave_rho2);
        }

        adjustment = (adjustment1 + 2.0 * adjustment2) / 3.0;
      
    }

    return adjustment;

}


void ExtractPvalues(const vector<SingleSubgroup> &single_level, string par_info, int &iterator, int par_index, vector<int> &signat, vector<double> &pvalue){

    int n_subgroups = single_level.size();

    for(int i = 0; i < n_subgroups; i++) {

        bool skip = false;
        for(int j =0; j < iterator;++j){
            if (single_level[i].signat == signat[j]){
                skip = true;
                break;
            }
        }
        signat[iterator] = single_level[i].signat;
        if(skip){
            continue;
        }

        string cur_info(par_info);
 
        pvalue.push_back(single_level[i].pvalue);
 
        int index = par_index*100 + i+1;

        ++iterator;
        if(!single_level[i].subgroups.empty()){
            ExtractPvalues(single_level[i].subgroups,cur_info,iterator,index,signat, pvalue);
        }
    }
}


void IterateSubgroupSummaryCSV(const vector<SingleSubgroup> &single_level, ofstream &file, string par_info, int &iterator, int par_index, vector<int> &signat, const vector<int> &biomarker_index, const vector<double> &adj_pvalue){

    int n_subgroups = single_level.size();
    ostringstream vals;

    for(int i = 0; i < n_subgroups; i++) {

        bool skip = false;
        for(int j =0; j < iterator;++j){
            if (single_level[i].signat == signat[j]){
                skip = true;
                break;
            }
        }
        signat[iterator] = single_level[i].signat;
        if(skip){
            continue;
        }

        // Biomarker index
        vals.str("");
        vals << "   <component description=\"\" biomarker=\"" << biomarker_index[single_level[i].biomarker_index - 1]<<"\" ";
        // 1 if <=, 2 if >, 3 if =
        if (single_level[i].sign == 1)
            vals << "sign=\"&lt;=\" ";
        if (single_level[i].sign == 2)
            vals << "sign=\"&gt;\" ";
        if (single_level[i].sign == 3)
            vals << "sign=\"=\" ";
        // Numeric biomarker
        if (single_level[i].sign != 3) {
        	vals << "value=\"" << single_level[i].value[0] << "\"/> \n";
        } 
        else {
        	// Nominal biomarker
        	vals << "value=\"";
        	for (int j = 0; j <single_level[i].value.size();++j){
            	vals <<  single_level[i].value[j] << " ";
        	}
        	vals <<"\"/> \n";
        }
        
        string cur_info(par_info + vals.str());
        file<< " <subgroup> \n" << "  <definition> \n" << cur_info<< "  </definition> \n";

        file<< "  <parameters size_control=\""<<single_level[i].size_control<<"\" ";
        file<< " size_treatment=\""<<single_level[i].size_treatment <<"\" ";
        file<< "splitting_criterion=\""<< single_level[i].criterion <<"\" ";
        file<< "splitting_criterion_log_p_value=\""<< - single_level[i].adjusted_criterion_pvalue <<"\" ";
        file<< "p_value=\""<< single_level[i].pvalue <<"\" ";
        file<< "prom_estimate=\""<< single_level[i].prom_estimate <<"\" ";
        file<< "prom_sderr=\""<< single_level[i].prom_sderr <<"\" ";
        file<< "prom_sd=\""<< single_level[i].prom_sd <<"\" ";
        file<< "adjusted_p_value=\""<< adj_pvalue[iterator] <<"\"/> \n";
        file<< " </subgroup> \n";

        int index = par_index*100 + i+1;

        ++iterator;
        if(!single_level[i].subgroups.empty()){
            IterateSubgroupSummaryCSV(single_level[i].subgroups,file,cur_info,iterator,index,signat, biomarker_index, adj_pvalue);
        }
    }
}

struct ANCOVAResult{
    
    double estimate;
    double sderr;
    double test_stat;
    double pvalue;

};

// Compute an approximate penalized criterion on log scale: log(crit)+log(G)
double LogerrApprxSupp(const double &zcrit, const double &g, const int &sides) {
  
  double x = zcrit / sqrt(2.0);
  double res = x * x + log(377.0 / 324.0 * x + sqrt(1 + 314.0 / 847.0 * x * x)) - log(g) + (2 - sides) * log(2.0);
  res = -res;

  return res;

}

// Compute an approximate computing of -log(pvalue), pvalue is 2-sided based on 2*(1-F(x))
double Log1minusPSupp(const double &pval, const double &zcrit, const double &g, const int &sides) {
  
  // Computes LN(2*(1-F(z_crit)), F(x) is normal cdf
  double res;
  double threshold = 1e-16;

  if (pval > threshold) {
    res = log(1.0 - pow(1.0 - pval, g));
  }
  else {
    res= LogerrApprxSupp(zcrit, g, 2);
  }
 
  return res;

}

bool InvertMatrix (const ublas::matrix<double>& input, ublas::matrix<double>& inverse) { 
    using namespace boost::numeric::ublas; 
    typedef permutation_matrix<std::size_t> pmatrix; 
    // create a working copy of the input 
    matrix<double> A(input); 
    pmatrix pm(A.size1()); 
    // perform LU-factorization 
    int res = lu_factorize(A,pm); 
    if( res != 0 ) return false; 
    // create identity matrix of "inverse" 
    inverse.assign(ublas::identity_matrix<double>(A.size1())); 
    // back substitute to get the inverse 
    lu_substitute(A, pm, inverse); 

    return true; 
}

// Create a list of unique centers
vector<double> ListUniqueCenters(const vector<double> &vec) {

    int length = vec.size();
    int found, i, j;
    double value1, value2;
    vector<double> unique;

    for (i = 1; i < length; i++) {

        found = 0;

        value1 = vec[i];

        // Compare with other values if the current value is non-missing
        if (::isnan(value1) < 1) {

            for (j = 0; j < i; j++) {

                value2 = vec[j];
                if (::isnan(value2) < 1)   
                  if (value1 == value2) found = found + 1;
                
                }

            if (found == 0) unique.push_back(value1);    
        }
    }

    return unique;
 
}

// Create a list of unique values
vector<double> ListUniqueValues(const vector<double> &vec) {

    int length = vec.size();
    int found, i, j;
    double value1, value2;
    vector<double> unique;

    for (i = 0; i < length; i++) {

        found = 0;

        value1 = vec[i]; 

        // Compare with other values if the current value is non-missing
        if (::isnan(value1) < 1) {

            for (j = 0; j < i; j++) {

                value2 = vec[j];
                if (::isnan(value2) < 1)   
                  if (value1 == value2) found = found + 1;
                
                }

            if (found == 0) unique.push_back(value1);    
        }
    }

    return unique;
 
}


// Analyze treatment effect in the overall group of patients using an ANCOVA model with cov1, cov2 and center
ANCOVAResult ContANCOVA(const std::vector<double> &ancova_treatment, const std::vector<double> &ancova_outcome, const ModelCovariates &model_covariates, const int &analysis_method, const int &direction)  {

    ANCOVAResult ancova_result;
    int i, j, nrow = ancova_treatment.size(), dim, ncenters;
    std::vector<double> ancova_cov1, ancova_cov2, ancova_cov3, ancova_cov4, ancova_center;

    ancova_cov1 = model_covariates.cov1; 
    ancova_cov2 = model_covariates.cov2; 
    ancova_cov3 = model_covariates.cov3; 
    ancova_cov4 = model_covariates.cov4; 
    ancova_center = model_covariates.cov_class; 

    ublas::matrix<double> x(1, 1);

    // Y matrix
    ublas::matrix<double> y(nrow, 1);
    for (j = 0; j < nrow; j++) y(j, 0) = ancova_outcome[j];

    // X matrix
    if (analysis_method == 2) {    

        dim = 3;
        x.resize(nrow, dim);

        for(i = 0; i < dim; i++) {
          for (j = 0; j < nrow; j++) {
            if (i == 0) x(j, i) = 1.0;
            if (i == 1) x(j, i) = ancova_treatment[j];
            if (i == 2) x(j, i) = ancova_cov1[j];
          }
        }

    }

    if (analysis_method == 3) {    

        dim = 4;
        x.resize(nrow, dim);

        for(i = 0; i < dim; i++) {
          for (j = 0; j < nrow; j++) {
            if (i == 0) x(j, i) = 1.0;
            if (i == 1) x(j, i) = ancova_treatment[j];
            if (i == 2) x(j, i) = ancova_cov1[j];
            if (i == 3) x(j, i) = ancova_cov2[j];
          }
        }

    }

    // List of centers within this subset
    if (analysis_method == 4) {

        std::vector<double> unique_centers = ListUniqueCenters(ancova_center);
        sort(unique_centers.begin(), unique_centers.end());
        ncenters = unique_centers.size();

        dim = ncenters + 4;
        x.resize(nrow, dim);

        // Design matrix (without the first center)
        ublas::matrix<double> center_matrix(nrow, ncenters);
        for(j = 0; j < nrow; j++) {
          for (i = 0; i < ncenters; i++) {
            center_matrix(j, i) = 0.0;
            if (ancova_center[j] == unique_centers[i]) center_matrix(j, i) = 1.0;
          }
        }

        for(i = 0; i < dim; i++) {
          for (j = 0; j < nrow; j++) {
            if (i == 0) x(j, i) = 1.0;
            if (i == 1) x(j, i) = ancova_treatment[j];
            if (i == 2) x(j, i) = ancova_cov1[j];
            if (i == 3) x(j, i) = ancova_cov2[j];
            if (i > 3) x(j, i) = center_matrix(j, i - 4);
          }
        }

    }

    // List of centers within this subset
    if (analysis_method == 5) {

        std::vector<double> unique_centers = ListUniqueCenters(ancova_center);
        sort(unique_centers.begin(), unique_centers.end());
        ncenters = unique_centers.size();

        dim = ncenters + 5;
        x.resize(nrow, dim);

        // Design matrix (without the first center)
        ublas::matrix<double> center_matrix(nrow, ncenters);
        for(j = 0; j < nrow; j++) {
          for (i = 0; i < ncenters; i++) {
            center_matrix(j, i) = 0.0;
            if (ancova_center[j] == unique_centers[i]) center_matrix(j, i) = 1.0;
          }
        }

        for(i = 0; i < dim; i++) {
          for (j = 0; j < nrow; j++) {
            if (i == 0) x(j, i) = 1.0;
            if (i == 1) x(j, i) = ancova_treatment[j];
            if (i == 2) x(j, i) = ancova_cov1[j];
            if (i == 3) x(j, i) = ancova_cov2[j];
            if (i == 4) x(j, i) = ancova_cov3[j];
            if (i > 4) x(j, i) = center_matrix(j, i - 5);
          }
        }

    }

    // List of centers within this subset
    if (analysis_method == 6) {

        std::vector<double> unique_centers = ListUniqueCenters(ancova_center);
        sort(unique_centers.begin(), unique_centers.end());
        ncenters = unique_centers.size();

        dim = ncenters + 6;
        x.resize(nrow, dim);

        // Design matrix (without the first center)
        ublas::matrix<double> center_matrix(nrow, ncenters);
        for(j = 0; j < nrow; j++) {
          for (i = 0; i < ncenters; i++) {
            center_matrix(j, i) = 0.0;
            if (ancova_center[j] == unique_centers[i]) center_matrix(j, i) = 1.0;
          }
        }

        for(i = 0; i < dim; i++) {
          for (j = 0; j < nrow; j++) {
            if (i == 0) x(j, i) = 1.0;
            if (i == 1) x(j, i) = ancova_treatment[j];
            if (i == 2) x(j, i) = ancova_cov1[j];
            if (i == 3) x(j, i) = ancova_cov2[j];
            if (i == 4) x(j, i) = ancova_cov3[j];
            if (i == 5) x(j, i) = ancova_cov4[j];
            if (i > 5) x(j, i) = center_matrix(j, i - 6);
          }
        }

    }

    // Inverse of X tr(X)  
    ublas::matrix<double> trxx = prod(trans(x), x);
    ublas::matrix<double> trxx_inv(dim, dim); 

    // Estimate
    ublas::matrix<double> txy = prod(trans(x), y);
    ublas::matrix<double> est = prod(trxx_inv, txy);

    // Residuals
    ublas::matrix<double> xest = prod(x, est);

    // SSE
    double sse = 0.0;
    for (j = 0; j < nrow; j++) sse += (y(j, 0) - xest(j, 0)) * (y(j, 0) - xest(j, 0));
    int p = dim;
    double sigma = sqrt(sse / (nrow - p));
    vector<double> sterr(dim);
    for (i = 0; i < dim; i++) sterr[i] = sigma * sqrt(trxx_inv(i, i));

    double t_stat = -5.0;

    if (abs(sterr[1]) >= 0.00000001) t_stat = direction * est(1, 0) / sterr[1];
    
    if (std::isnan(t_stat) == 1) t_stat = -5.0;

    // Treatment effect estimate
    double estimate = est(1, 0);

    int df;

    // Test statistic 
    if (nrow - p > 0) {
        df = nrow - p;
    } else {
        df = 10;
    }

    // One-sided p-value
    double pvalue = 1.0 - rcpp_pt(t_stat, df);

    ancova_result.estimate = estimate;
    ancova_result.sderr = sterr[1];
    ancova_result.test_stat = t_stat;
    ancova_result.pvalue = pvalue;

    return ancova_result;

}

// Analyze treatment effect in the overall group of patients using an ANCOVA model
SingleSubgroup OverallAnalysisCont(const std::vector<double> &treatment, const std::vector<double> &outcome, const ModelCovariates &model_covariates, const int &analysis_method, const int &direction)  {
  
    SingleSubgroup res;

    int i;
    int n = treatment.size();

    std::vector<double> ancova_treatment = treatment;
    std::vector<double> ancova_outcome = outcome;

    // Extract the covariates
    ANCOVAResult ancova_result = ContANCOVA(ancova_treatment, ancova_outcome, model_covariates, analysis_method, direction);

    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control++; else size_treatment++;
    }

    // Save the results

    // Splitting criterion
    res.criterion = 0;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = 0;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = 0;
    res.test_statistic = ancova_result.test_stat;
    res.pvalue = ancova_result.pvalue;
    res.prom_estimate = ancova_result.estimate;
    res.prom_sderr = ancova_result.sderr;   
    res.prom_sd = -1.0;   

    res.adjusted_pvalue = -1.0;
    // Vector of biomarker values to define the current patient subgroup
    vector<double> value;
    value.clear();
    value.push_back(0);
    res.value = value; 
    // 1 if <=, 2 if >, 3 if =
    res.sign = 0;
    // Size of the current patient subgroup
    res.size_control = size_control;
    res.size_treatment = size_treatment;
    res.size = size_control + size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows;
    subgroup_rows.clear();
    subgroup_rows.push_back(0);
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup 
    res.biomarker_index = 0;
    // Level of the current patient subgroup  
    res.level = 0;
    // Parent index for the current patient subgroup 
    res.parent_index = -1;
    // Indexes of child subgroups for the current patient subgroup 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(0);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors)  
    res.code = 0;
    // Is the current patient subgroup terminal (0/1) 
    res.terminal_subgroup = 0;

    return res;    

}

NumericVector VectorPower(const NumericVector vec, double power) {
    NumericVector res = NumericVector(vec.size());
    for (int i = 0; i < vec.size(); i++) {
        res[i] = pow(vec[i], power);
    }
    return res;
}

double EuclideanDistance(NumericVector x, NumericVector y) {
    return sqrt(Rcpp::sum(VectorPower((x - y), 2)));
}

NumericVector LRAlphaSteps(const double start, const double end, const double step, const double value_to_skip=0.0) {
    double val = start;
    const int arr_len = (int)((end - start) / step);
    NumericVector arr(arr_len);
    int index = 0;
    while (val <= end) {
        if (val == value_to_skip) {
            val += step;
            index += 1;
            continue;
        }
        arr[index] = val;
        val += step;
        index += 1;
    }
    return arr;
}

NumericVector LRLogitProb(NumericMatrix X, NumericVector beta) {
    //   # compute vector p of probabilities for logistic regression with logit link
    NumericVector p;
    Eigen::Map<Eigen::MatrixXd> X_eigen = as<Eigen::Map<Eigen::MatrixXd> >(X);
    Eigen::Map<Eigen::MatrixXd> beta_eigen = as<Eigen::Map<Eigen::MatrixXd> >(beta);
    p = 1 / (1 + exp(NumericVector(wrap(-1 * (X_eigen * beta_eigen)))));
    return p;
}

double LRLogLikelihood(NumericVector y, NumericVector m, NumericVector p) {
    //   # binomial log likelihood function
    //   # input:   vectors: y = counts; m = sample sizes; p = probabilities
    //   # output: log-likelihood l, a scalar

    NumericVector log_p = log(p);
    NumericVector log_p_conj = log(1 - p);

    Eigen::Map<Eigen::MatrixXd> y_eigen = as<Eigen::Map<Eigen::MatrixXd> >(y);
    Eigen::Map<Eigen::MatrixXd> m_eigen = as<Eigen::Map<Eigen::MatrixXd> >(m);
    Eigen::Map<Eigen::MatrixXd> log_p_eigen = as<Eigen::Map<Eigen::MatrixXd> >(log_p);
    Eigen::Map<Eigen::MatrixXd> log_p_conj_eigen = as<Eigen::Map<Eigen::MatrixXd> >(log_p_conj);

    Eigen::MatrixXd y_trans = y_eigen.transpose();
    Eigen::MatrixXd my_trans = (m_eigen - y_eigen).transpose();

    double ll = as<double>(wrap((y_trans * log_p_eigen) + (my_trans * log_p_conj_eigen)));
    return ll;
}

ANCOVAResult BinANCOVA(const std::vector<double> &ancova_treatment, const std::vector<double> &ancova_outcome, const NumericMatrix &cont_covariates, const NumericMatrix &class_covariates, const int &direction, const int &cont_covariates_flag, const int &class_covariates_flag) {

        int i, j, k, l,
                  n_cont_covariates = 0, 
                  n_class_covariates = 0,
                  nobs = ancova_treatment.size(),
                  n_total_covariates,
                  nlevels = 0;

        vector<double> current_class_covariate, unique_levels;                  
        // Outcome
        NumericVector outcome(nobs), temp_vec; 

        if (cont_covariates_flag >= 1) n_cont_covariates = cont_covariates.ncol(); 

        if (class_covariates_flag >= 1) {   

            n_class_covariates = class_covariates.ncol();

            // Compute the number of columns for categorical covariates
            for (i = 0; i < n_class_covariates; i++) {

                temp_vec = class_covariates(_, i);
                current_class_covariate = as<vector<double>>(temp_vec);
                nlevels += (CountUniqueValues(current_class_covariate) - 1);

            }

        }

        // Total number of covariates in the model, including the intercept and treatment
        n_total_covariates = n_cont_covariates + nlevels + 2;

        // Matrix of covariates
        NumericMatrix covariates(nobs, n_total_covariates); 
        NumericVector intercept(nobs);

        for (i = 0; i < nobs; i++) {

            outcome[i] = ancova_outcome[i];
            covariates(i, 0) = 1.0;
            covariates(i, 1) = ancova_treatment[i];
            intercept[i] = 1.0;            

            // Add continuous covariates
            if (cont_covariates_flag >= 1) {
                for (j = 0; j < n_cont_covariates; j++) covariates(i, j + 2) = cont_covariates(i, j);
            }

        }

        if (class_covariates_flag >= 1) { 

            // Add categorical covariates
            l = n_cont_covariates + 1;
            for (j = 0; j < n_class_covariates; j++) {

                // List of unique levels for the current categorical covariate
                temp_vec = class_covariates(_ , j);
                current_class_covariate = as<vector<double>>(temp_vec);

                unique_levels = ListUniqueValues(current_class_covariate);
                sort(unique_levels.begin(), unique_levels.end());
                nlevels = unique_levels.size();

                // Design matrix for the current categorical covariate (without the first level)
                for (i = 0; i < nobs; i++) {
                    for (k = 1; k < nlevels; k++) {
                        covariates(i, l + k) = 0.0;
                        if (current_class_covariate[i] == unique_levels[k]) covariates(i, l + k) = 1.0;
                    }

                }

                l += nlevels - 1;

            }
        }

        // Initial values of model parameters and standard errors
        NumericVector beta1(n_total_covariates), se(n_total_covariates);

        double eps_beta = precision;
        double eps_ll = precision;

        NumericVector beta2 = NumericVector(beta1.size()) + (-1 * pow(10.0, 15));
        double beta_diff = EuclideanDistance(beta1, beta2);

        NumericVector pp1 = LRLogitProb(covariates, beta1);
        NumericVector pp2 = LRLogitProb(covariates, beta2);
        double llike1 = LRLogLikelihood(outcome, intercept, pp1);
        double llike2 = LRLogLikelihood(outcome, intercept, pp2);
        double diff_like = abs(llike1 - llike2);

        if(std::isnan(diff_like)) {
            diff_like = pow(10.0, 9);
        }

        NumericVector alpha_steps = LRAlphaSteps(-1.0, 2.0, 0.1);

        int iter = 0;
        while ((iter < max_iter) && (beta_diff > eps_beta) && (diff_like > eps_ll)) {
            iter += 1;
            NumericVector mu2;

            beta2 = clone(beta1);
            NumericVector pp2 = LRLogitProb(covariates, beta2);
            mu2 = intercept * pp2; // mean

            NumericMatrix v2;
            NumericVector multiplied1 = intercept * pp2 * (1 - pp2);
            Eigen::Map<Eigen::VectorXd> multiplied1_eigen = as<Eigen::Map<Eigen::VectorXd> >(multiplied1);
            Eigen::MatrixXd v_eigen = multiplied1_eigen.asDiagonal();
            v2 = NumericMatrix(wrap(v_eigen));

            Eigen::Map<Eigen::MatrixXd> covariates_eigen = as<Eigen::Map<Eigen::MatrixXd> >(covariates);
            NumericVector ymu2 = outcome - mu2;
            Eigen::Map<Eigen::MatrixXd> ymu2_eigen = as<Eigen::Map<Eigen::MatrixXd> >(ymu2);
            Eigen::MatrixXd score2 = covariates_eigen.transpose() * ymu2_eigen;

            Eigen::MatrixXd to_solve = covariates_eigen.transpose() * v_eigen * covariates_eigen;
            Eigen::VectorXd increm_eigen = to_solve.colPivHouseholderQr().solve(score2);
            NumericVector increm = NumericVector(wrap(increm_eigen));

            const int len_alpha_steps = alpha_steps.size();
            NumericVector llike_alpha_step(len_alpha_steps);
            for(int i = 0; i < len_alpha_steps; i++) {
                llike_alpha_step[i] = - pow(10.0, 15);
            }

            for(int i = 0; i < len_alpha_steps; i++) {
                llike_alpha_step[i] = LRLogLikelihood(
                    outcome, intercept, LRLogitProb(covariates, beta2 + alpha_steps[i] * increm)
                );
            }

            int which_max_increment = which_max(llike_alpha_step);
            beta1 = beta2 + alpha_steps[which_max_increment] * increm;
            
            beta_diff = EuclideanDistance(beta1, beta2);

            llike2 = llike1;
            llike1 = LRLogLikelihood(outcome, intercept, LRLogitProb(covariates, beta1));
            diff_like = abs(llike1 - llike2);
        }

        // Standard errors
        NumericVector pp = LRLogitProb(covariates, beta1);
        NumericVector multiplied1 = intercept * pp * (1 - pp);
        Eigen::Map<Eigen::VectorXd> multiplied1_eigen = as<Eigen::Map<Eigen::VectorXd> >(multiplied1);
        Eigen::Map<Eigen::MatrixXd> covariates_eigen = as<Eigen::Map<Eigen::MatrixXd> >(covariates);
        Eigen::MatrixXd v1_eigen = multiplied1_eigen.asDiagonal();
        se = sqrt(NumericVector(wrap(
            (covariates_eigen.transpose() * v1_eigen * covariates_eigen).inverse().diagonal()
        )));

    ANCOVAResult ancova_result;

    // Treatment effect estimate and test
    double test_stat = (direction + 0.0) * beta1[1] / se[1];
    // double test_stat = 1.0;
    double pvalue = 1.0 - rcpp_pnorm(test_stat);

    ancova_result.estimate = exp((direction + 0.0) * beta1[1]);
    ancova_result.sderr = se[1];
    ancova_result.test_stat = test_stat;
    ancova_result.pvalue = pvalue;

    return ancova_result;


}

// qq


Rcpp::NumericMatrix VectorToMatrix(Rcpp::NumericVector row, const string type) {
    Rcpp::NumericMatrix row_matrix;
    if (type == "row") {
        row_matrix = Rcpp::NumericMatrix(1, row.size());
        for(int i = 0; i < row.size(); i++) {
            row_matrix(0, i) = row(i);
        }
    } else if (type == "column") {
        row_matrix = Rcpp::NumericMatrix(row.size(), 1);
        for(int i = 0; i < row.size(); i++) {
            row_matrix(i, 0) = row(i);
        }
    } else {
        throw;
    }
    return row_matrix;
}

NumericVector ConvertToNumericVector(const vector<double> &vec) {

    int i, m = vec.size();
    NumericVector num_vec(m);

    for (i = 0; i < m; i++) num_vec[i] = vec[i];

    return num_vec;

}

Rcpp::NumericVector Sequence(const int first, const int last) {
    NumericVector y(abs(last - first) + 1);
    if (first < last) {
        std::iota(y.begin(), y.end(), first);
    } else {
        std::iota(y.begin(), y.end(), last);
        std::reverse(y.begin(), y.end());
    }
    return y;
}

Rcpp::IntegerVector ROrder(NumericVector x) {
    if (is_true(any(duplicated(x)))) {
        Rcout << "There are duplicates in 'x'; order not guaranteed to match that of R's base::order" << endl;
        throw;
    }
    NumericVector sorted = clone(x).sort();
    return match(sorted, x);
}


Rcpp::List PrepareSurvdata(const vector<double> &time, const vector<double> &censor) {

    NumericVector time_v, cens_v;

    time_v = ConvertToNumericVector(time);
    cens_v = ConvertToNumericVector(censor);
    
    NumericVector index_obs = Sequence(0, time_v.length() - 1);

    LogicalVector sel = (cens_v == 0);
    NumericVector event_to_index = index_obs[sel];
    
    NumericVector selected_times = time_v[sel];
    NumericVector event_times = UniquePreserveOrder(selected_times);

    IntegerVector ord = ROrder(event_times) - 1;
    NumericVector event_times_ordered = event_times[ord]; //rearrange_by_order(event_times, ord);
    NumericVector ord_in_full = event_to_index[ord]; //rearrange_by_order(event_to_index, ord);

    int n = time_v.length();
    int nev = event_times.length();
    NumericMatrix risk_sets = NumericMatrix(nev, n);
    IntegerVector nrisk = IntegerVector(nev);

    IntegerVector ord_times = ROrder(time_v);

    int i1 = nev - 1;
    for (int j = 0; j < n; j++) {
        int j1 = n - j - 1;  
        while (time_v[ord_times[j1] - 1] < event_times_ordered[i1] && i1 > 0) {
            i1 = i1 - 1;
        }
        if (time_v[ord_times[j1] - 1] >= event_times_ordered[i1]) {
            nrisk[i1] = nrisk[i1] + 1;
            risk_sets(i1, nrisk[i1] - 1) = ord_times[j1];
        }
    }

    return Rcpp::List::create(
        Named("times_ind") = ord_in_full,
        Named("risk_sets") = risk_sets,
        Named("nrisk") = nrisk
    );
};

ScoreResult ComputeScoreInformation(
    Rcpp::NumericVector beta, Rcpp::NumericMatrix covariates, Rcpp::NumericVector event_time_indices, 
    Rcpp::NumericVector num_event_times, Rcpp::NumericMatrix risk_sets
){

    Rcpp::NumericMatrix sum = Rcpp::NumericMatrix(1, covariates.ncol());
    Rcpp::NumericMatrix v = Rcpp::NumericMatrix(covariates.ncol(), covariates.ncol());
    // store covariate-specific weighted averages for each risk set
    Rcpp::NumericMatrix xmeans = Rcpp::NumericMatrix(num_event_times.size(), covariates.ncol());
    // store provisional means
    Rcpp::NumericMatrix xm = Rcpp::NumericMatrix(1, covariates.ncol());

    for (int i = 0; i < covariates.ncol(); i++) { 
        xm[i] = Rcpp::mean(covariates(_, i));
    }

    Rcpp::NumericVector num = Rcpp::NumericVector(covariates.ncol());
    Rcpp::NumericMatrix num2 = Rcpp::NumericMatrix(covariates.ncol(),covariates.ncol());
    Rcpp::NumericMatrix num_cum = Rcpp::NumericMatrix(covariates.ncol(),covariates.ncol());
    double den = 0;
    for (int i = 1; i <= num_event_times.size(); i++) {
        int i1 = (num_event_times.size() - i); // will cycle from the last risk set to first to use accumulated data

        for (int j = 0; j < num_event_times[i1]; j++) {
            // Rcpp::NumericVector row_in_risks = covariates(risk_sets[i1, j], _);
            // Rcpp::NumericVector mult = row_in_risks * beta;
            double w = exp(Rcpp::sum(covariates(risk_sets(i1, j) - 1, _) * beta));
            num += covariates(risk_sets(i1, j) - 1, _) * w;

            Rcpp::NumericVector retrieved1 = covariates(risk_sets(i1, j) - 1, _);
            Rcpp::NumericMatrix transformed_to_matrix1 = VectorToMatrix(retrieved1, "row");
            Eigen::Map<Eigen::MatrixXd> transformed_to_matrix1_eigen = as<Eigen::Map<Eigen::MatrixXd> >(transformed_to_matrix1);
            Eigen::Map<Eigen::MatrixXd> xm_eigen1 = as<Eigen::Map<Eigen::MatrixXd> >(xm);
            Rcpp::NumericMatrix normalized = Rcpp::NumericMatrix(Rcpp::wrap(transformed_to_matrix1_eigen - xm_eigen1));

            Eigen::Map<Eigen::MatrixXd> normalized_eigen = as<Eigen::Map<Eigen::MatrixXd> >(normalized);
            Eigen::Map<Eigen::MatrixXd> normalized_transp_eigen = as<Eigen::Map<Eigen::MatrixXd> >(Rcpp::transpose(normalized));
            Eigen::MatrixXd multiplied = normalized_transp_eigen * normalized_eigen;
            num2 += Rcpp::NumericMatrix(Rcpp::wrap(multiplied * w));
            den += w;
        }

        Eigen::Map<Eigen::VectorXd> num_eigen = as<Eigen::Map<Eigen::VectorXd> >(num);
        Rcpp::NumericVector division_res = Rcpp::NumericVector(Rcpp::wrap(num_eigen / den));
        xmeans(i1,_) = division_res;

        Rcpp::NumericVector retrieved2 = covariates(event_time_indices(i1), _);
        Rcpp::NumericVector retrieved3 = xmeans(i1, _);
        sum += VectorToMatrix(retrieved2, "row") - VectorToMatrix(retrieved3, "row");
        
        num_cum = num2;

        Rcpp::NumericVector retrieved4 = xmeans(i1, _);
        Rcpp::NumericMatrix transformed_to_matrix2 = VectorToMatrix(retrieved4, "row");
        Eigen::Map<Eigen::MatrixXd> transformed_to_matrix2_eigen = as<Eigen::Map<Eigen::MatrixXd> >(transformed_to_matrix2);
        Eigen::Map<Eigen::MatrixXd> xm_eigen2 = as<Eigen::Map<Eigen::MatrixXd> >(xm);
        Rcpp::NumericMatrix normalized = Rcpp::NumericMatrix(Rcpp::wrap(transformed_to_matrix2_eigen - xm_eigen2));
        

        Eigen::Map<Eigen::MatrixXd> num_cum_eigen = as<Eigen::Map<Eigen::MatrixXd> >(num_cum);
        Eigen::Map<Eigen::MatrixXd> normalized_eigen = as<Eigen::Map<Eigen::MatrixXd> >(normalized);
        Eigen::Map<Eigen::MatrixXd> normalized_transp_eigen = as<Eigen::Map<Eigen::MatrixXd> >(Rcpp::transpose(normalized));
        Eigen::MatrixXd multiplied = normalized_transp_eigen * normalized_eigen;

        Rcpp::NumericMatrix num_corrected = NumericMatrix(Rcpp::wrap(num_cum_eigen - (multiplied * den)));

        v += num_corrected / den;
    }

    ScoreResult res = ScoreResult();
    res.score = sum;
    res.information = v;
    return(res);
}

ANCOVAResult SurvANCOVA(const std::vector<double> &ancova_treatment, const std::vector<double> &ancova_outcome, const std::vector<double> &ancova_censor, const NumericMatrix &cont_covariates, const NumericMatrix &class_covariates, const int &direction, const int &cont_covariates_flag, const int &class_covariates_flag) {

        int i, j, k, l,
                  n_cont_covariates = 0, 
                  n_class_covariates = 0,
                  nobs = ancova_treatment.size(),
                  n_total_covariates,
                  nlevels = 0, iter = 0;

        double diff = 1.0;
    
        vector<double> current_class_covariate, unique_levels;                  
        // Outcome
        NumericVector outcome(nobs), temp_vec; 

        if (cont_covariates_flag >= 1) n_cont_covariates = cont_covariates.ncol(); 

        if (class_covariates_flag >= 1) {   

            n_class_covariates = class_covariates.ncol();

            // Compute the number of columns for categorical covariates
            for (i = 0; i < n_class_covariates; i++) {

                temp_vec = class_covariates(_, i);
                current_class_covariate = as<vector<double>>(temp_vec);
                nlevels += (CountUniqueValues(current_class_covariate) - 1);

            }

        }

        // Total number of covariates in the model, including the intercept and treatment
        n_total_covariates = n_cont_covariates + nlevels + 1;

        // Matrix of covariates
        NumericMatrix covariates(nobs, n_total_covariates); 
        NumericVector intercept(nobs);

        for (i = 0; i < nobs; i++) {

            outcome[i] = ancova_outcome[i];
            covariates(i, 0) = ancova_treatment[i];
            intercept[i] = 1.0;            

            // Add continuous covariates
            if (cont_covariates_flag >= 1) {
                for (j = 0; j < n_cont_covariates; j++) covariates(i, j + 1) = cont_covariates(i, j);
            }

        }

        if (class_covariates_flag >= 1) { 

            // Add categorical covariates
            l = n_cont_covariates;
            for (j = 0; j < n_class_covariates; j++) {

                // List of unique levels for the current categorical covariate
                temp_vec = class_covariates(_ , j);
                current_class_covariate = as<vector<double>>(temp_vec);

                unique_levels = ListUniqueValues(current_class_covariate);
                sort(unique_levels.begin(), unique_levels.end());
                nlevels = unique_levels.size();

                // Design matrix for the current categorical covariate (without the first level)
                for (i = 0; i < nobs; i++) {
                    for (k = 1; k < nlevels; k++) {

                        covariates(i, l + k) = 0.0;
                        if (current_class_covariate[i] == unique_levels[k]) covariates(i, l + k) = 1.0;
                    }

                }

                l += nlevels - 1;

            }
        }

    Rcpp::List survdata = PrepareSurvdata(ancova_outcome, ancova_censor);

    Rcpp::NumericMatrix risk_sets = survdata["risk_sets"];
    Rcpp::NumericVector num_event_times = survdata["nrisk"];
    Rcpp::NumericVector event_time_indices = survdata["times_ind"];
    
    // Initial values of model parameters and standard errors
    NumericVector beta(n_total_covariates), se(n_total_covariates);

    while (diff > precision && iter < max_iter) {

        iter = iter + 1;

        Rcpp::NumericVector beta_prev = clone(beta);
        ScoreResult res_iter = ComputeScoreInformation(
            beta, covariates, event_time_indices, num_event_times, risk_sets
        );

        Eigen::Map<Eigen::MatrixXd> score_eigen = as<Eigen::Map<Eigen::MatrixXd> >(res_iter.score);
        Eigen::Map<Eigen::MatrixXd> information_eigen = as<Eigen::Map<Eigen::MatrixXd> >(res_iter.information);
        MatrixXd solved = information_eigen.inverse();
        Eigen::MatrixXd multiplied = score_eigen * solved;

        beta += NumericVector(Rcpp::wrap(multiplied));
        
        diff = Rcpp::sum(Rcpp::abs(beta - beta_prev));
    }

    ScoreResult res_iter = ComputeScoreInformation(
            beta, covariates, event_time_indices, num_event_times, risk_sets
    );

    Eigen::Map<Eigen::MatrixXd> information_eigen = as<Eigen::Map<Eigen::MatrixXd> >(res_iter.information);
    se = Rcpp::NumericVector(Rcpp::wrap(information_eigen.inverse().diagonal()));

    ANCOVAResult ancova_result;

    // Treatment effect estimate and test
    double test_stat = (-direction + 0.0) * beta[0] / sqrt(se[0]);
    // double test_stat = 1.0;
    double pvalue = 1.0 - rcpp_pnorm(test_stat);

    ancova_result.estimate = exp((direction + 0.0) * beta[0]);
    ancova_result.sderr = sqrt(se[0]);
    ancova_result.test_stat = test_stat;
    ancova_result.pvalue = pvalue;

    return ancova_result;

}




// Analyze treatment effect in the overall group of patients using a logistic regression model
SingleSubgroup OverallAnalysisBin(const std::vector<double> &treatment, const std::vector<double> &outcome, const NumericMatrix &cont_covariates, const NumericMatrix &class_covariates, const int &direction, const int &cont_covariates_flag, const int &class_covariates_flag)  {
  
    SingleSubgroup res;

    int i;
    int n = treatment.size();

    // Extract the covariates 
    ANCOVAResult ancova_result = BinANCOVA(treatment, outcome, cont_covariates, class_covariates, direction, cont_covariates_flag, class_covariates_flag);

    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control++; else size_treatment++;
    }

    // Save the results

    // Splitting criterion
    res.criterion = 0;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = 0;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = 0;
    res.test_statistic = ancova_result.test_stat;
    res.pvalue = ancova_result.pvalue;
    res.prom_estimate = ancova_result.estimate;
    res.prom_sderr = ancova_result.sderr;   
    res.prom_sd = -1.0;   

    res.adjusted_pvalue = -1.0;
    // Vector of biomarker values to define the current patient subgroup
    vector<double> value;
    value.clear();
    value.push_back(0);
    res.value = value; 
    // 1 if <=, 2 if >, 3 if =
    res.sign = 0;
    // Size of the current patient subgroup
    res.size_control = size_control;
    res.size_treatment = size_treatment;
    res.size = size_control + size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows;
    subgroup_rows.clear();
    subgroup_rows.push_back(0);
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup 
    res.biomarker_index = 0;
    // Level of the current patient subgroup  
    res.level = 0;
    // Parent index for the current patient subgroup 
    res.parent_index = -1;
    // Indexes of child subgroups for the current patient subgroup 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(0);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors)  
    res.code = 0;
    // Is the current patient subgroup terminal (0/1) 
    res.terminal_subgroup = 0;

    return res;    

}

// Analyze treatment effect in the overall group of patients using the Cox PH regression model
SingleSubgroup OverallAnalysisSurv(const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &censor, const NumericMatrix &cont_covariates, const NumericMatrix &class_covariates, const int &direction, const int &cont_covariates_flag, const int &class_covariates_flag)  {
  
    SingleSubgroup res;

    int i;
    int n = treatment.size();

    // Extract the covariates 
    ANCOVAResult ancova_result = SurvANCOVA(treatment, outcome, censor, cont_covariates, class_covariates, direction, cont_covariates_flag, class_covariates_flag);

    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control++; else size_treatment++;
    }

    // Save the results

    // Splitting criterion
    res.criterion = 0;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = 0;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = 0;
    res.test_statistic = ancova_result.test_stat;
    res.pvalue = ancova_result.pvalue;
    res.prom_estimate = ancova_result.estimate;
    res.prom_sderr = ancova_result.sderr;   
    res.prom_sd = -1.0;   

    res.adjusted_pvalue = -1.0;
    // Vector of biomarker values to define the current patient subgroup
    vector<double> value;
    value.clear();
    value.push_back(0);
    res.value = value; 
    // 1 if <=, 2 if >, 3 if =
    res.sign = 0;
    // Size of the current patient subgroup
    res.size_control = size_control;
    res.size_treatment = size_treatment;
    res.size = size_control + size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows;
    subgroup_rows.clear();
    subgroup_rows.push_back(0);
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup 
    res.biomarker_index = 0;
    // Level of the current patient subgroup  
    res.level = 0;
    // Parent index for the current patient subgroup 
    res.parent_index = -1;
    // Indexes of child subgroups for the current patient subgroup 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(0);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors)  
    res.code = 0;
    // Is the current patient subgroup terminal (0/1) 
    res.terminal_subgroup = 0;

    return res;    

}

// Supportive function for processing biomarker vectors with missing observations
void na_index_rows(vector<double> X, vector<int> &bio_rows)
{
    int n=X.size();
    for (int i = 0; i < n; i++) {
        if(X[i]!=X[i]) {
            bio_rows[i]=0;
        }
    }
}



// Supportive function for sorting vectors (from smallest to largest)
bool DDPairSortUp (const ddpair& l, const ddpair& r) { 
    return l.first < r.first; 
}

// Supportive function for sorting vectors (from largest to smallest)
bool DIPairSortDown (const dipair& l, const dipair& r) { 
    return l.first > r.first; 
}

// Supportive function for sorting vectors (from smallest to largest)
bool DIPairSortUp (const dipair& l, const dipair& r) { 
    return l.first < r.first; 
}


// Splitting criterion function
double CriterionFunction(const double &test_left, const double &test_right, const int &criterion_type) {

    double criterion;

    // Differential criterion function
    if (criterion_type == 1) criterion = std::abs(test_left - test_right) / sqrt(2);
    // Maximum criterion function
    if (criterion_type == 2) criterion = std::max(test_left, test_right);
    // Directional criterion function
    if (criterion_type == 3) {
        if (std::max(test_left, test_right) > 0) {
            if (std::min(test_left, test_right) > 0) {
                criterion = std::abs(test_left - test_right);
            }
            else {
                criterion = std::max(test_left, test_right);
            }
        } 
        else {
            criterion = -5.0;
        }
    }

    return criterion;
}

// Compute criterion p-value
double CriterionPvalue(const double &criterion, const int &criterion_type) {

    double pvalue, temp;

    // Differential criterion function
    if (criterion_type == 1) pvalue = 2.0 * rcpp_pnorm(-criterion); 

    // Maximum criterion function
    if (criterion_type == 2) pvalue = 2.0 * rcpp_pnorm(-criterion); 

    // Directional criterion function
    if (criterion_type == 3) {
        temp = rcpp_pnorm(criterion / sqrt(2.0)) - 1.0;
        pvalue = 1.0 - (1.0 - 4.0 * temp * temp + 2.0 * (2.0 * rcpp_pnorm(criterion) - 1.0)) / 3.0;
    }
    
    if (pvalue > 1.0) pvalue = 1.0; 

    return pvalue;

}


// Treatment effect test for continuous outcome variable (t test)
double TTestStatistic(const double &treatment_mean, const double &treatment_var,
                         const double &control_mean, const double &control_var,
                         const double &treatment_count, const double &control_count,
                         const int &direction, int &error_flag) {

    double pooled_var;
    double statistic = 0;

    // Set the flag to error
    error_flag = 1;

    if (treatment_count >1 && control_count >1) {
        pooled_var = (treatment_var * (treatment_count - 1) + control_var * (control_count - 1))/(treatment_count + control_count - 2);
        if (pooled_var > 1.0e-10) {
            statistic = direction * (treatment_mean - control_mean) / sqrt(pooled_var * (1/treatment_count + 1/control_count));
            // No error
            error_flag = 0;
        }
        
    }

    return statistic;
}



// Treatment effect test in overall group of patients when the outcome variable is continuous
SingleSubgroup ContOutOverallAnalysis(const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const int &direction) {

    SingleSubgroup res;

    int i;
    int n = treatment.size();

    double treatment_count = 0; 
    double control_count = 0;

    double treatment_sum = 0;
    double control_sum = 0;

    double treatment_sum_squares = 0;
    double control_sum_squares = 0;

    int error_flag;

    double treatment_mean, treatment_var, control_mean, control_var, test_statistic;

    for(i = 0; i < n; i++) {
        if (treatment[i] == 0) {
            control_sum += outcome[i];
            control_sum_squares += outcome[i] * outcome[i];
            control_count ++;
        }
        if (treatment[i] == 1) {
            treatment_sum += outcome[i];
            treatment_sum_squares += outcome[i] * outcome[i];
            treatment_count ++;
        }
    }

    treatment_mean = treatment_sum / treatment_count;
    treatment_var = (treatment_sum_squares / treatment_count - treatment_mean * treatment_mean) * treatment_count / (treatment_count - 1); 
    control_mean = control_sum / control_count;
    control_var = (control_sum_squares / control_count - control_mean * control_mean) * control_count / (control_count - 1); 

    test_statistic = TTestStatistic(treatment_mean, treatment_var,
                                        control_mean, control_var,
                                        treatment_count, control_count,
                                        direction, error_flag);

    double pvalue = 1.0 - rcpp_pnorm(test_statistic);

    double pooled_var = (treatment_var * (treatment_count - 1) + control_var * (control_count - 1))/(treatment_count + control_count - 2);

    // Save the results

    // Splitting criterion
    res.criterion = 0;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = 0;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = 0;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    res.adjusted_pvalue = -1.0;    
    res.prom_estimate = (treatment_mean - control_mean); 
    res.prom_sderr = -1.0;   
    res.prom_sd = sqrt(pooled_var);   
    res.size_control = control_count;
    res.size_treatment = treatment_count;

    // Vector of biomarker values to define the current patient subgroup
    vector<double> value;
    value.clear();
    value.push_back(0);
    res.value = value; 
    // 1 if <=, 2 if >, 3 if =
    res.sign = 0;
    // Size of the current patient subgroup
    res.size = control_count + treatment_count;
    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows;
    subgroup_rows.clear();
    subgroup_rows.push_back(0);
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup 
    res.biomarker_index = 0;
    // Level of the current patient subgroup  
    res.level = 0;
    // Parent index for the current patient subgroup 
    res.parent_index = -1;
    // Indexes of child subgroups for the current patient subgroup 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(0);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors)  
    res.code = 0;
    // Is the current patient subgroup terminal (0/1) 
    res.terminal_subgroup = 0;

    return res;    

}


// Treatment effect test in overall group of patients when the outcome variable is time-to-event
SingleSubgroup SurvOutOverallAnalysis(const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const int &direction) {

    SingleSubgroup res;

    vector<double> outcome1(outcome);
    vector<double> outcome_censor1(outcome_censor);
    vector<double> treatment1(treatment);

    double prom_estimate = HazardRatio(outcome1, outcome_censor1, treatment1, direction);

    double test_statistic = LRTest(outcome1, outcome_censor1, treatment1, direction); 

    double treatment_count = 0; 
    double control_count = 0;
    int i, n = treatment.size();

    for(i = 0; i < n; i++) {
        if (treatment[i] == 0) {
            control_count ++;
        }
        if (treatment[i] == 1) {
            treatment_count ++;
        }
    }

    double pvalue = 1.0 - rcpp_pnorm(test_statistic);

    // Save the results

    // Splitting criterion
    res.criterion = 0;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = 0;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = 0;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    res.adjusted_pvalue = -1.0;    
    res.prom_estimate = prom_estimate;
    res.prom_sderr = -1.0;   
    res.prom_sd = -1.0;   
    res.size_control = control_count; 
    res.size_treatment = treatment_count;    
    // Vector of biomarker values to define the current patient subgroup
    vector<double> value;
    value.clear();
    value.push_back(0);
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = 0;
    // Size of the current patient subgroup
    res.size = treatment.size();
    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows;
    subgroup_rows.clear();
    subgroup_rows.push_back(0);
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup
    res.biomarker_index = 0;
    // Level of the current patient subgroup
    res.level = 0;
    // Parent index for the current patient subgroup
    res.parent_index = -1;
    // Indexes of child subgroups for the current patient subgroup
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(0);
    res.child_index = child_index;
    // Error code for testing (0 if no errors)
    res.code = 0;
    // Is the current patient subgroup terminal (0/1)
    res.terminal_subgroup = 0;

    return res;

}


int Included(const std::vector<int> &vec, const int &value) {

    int i, outcome = 0, m = vec.size();

    for(i = 0; i < m; i++) {
        if (vec[i] == value) outcome = 1;
    }

    return outcome;

}


// Analyze treatment effect in the overall group of patients
SingleSubgroup OverallAnalysis(const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const int &outcome_type, const int &direction) {
    
    SingleSubgroup res;

    // Outcome variable type
    // Continuous outcome
    if (outcome_type == 1) {
        res = ContOutOverallAnalysis(treatment, outcome, outcome_censor, direction);
    }
    // Binary outcome
    if (outcome_type == 2) {
        res = BinOutOverallAnalysis(treatment, outcome, outcome_censor, direction);
    }
    // Time-to-event outcome
    if (outcome_type == 3) {
        res = SurvOutOverallAnalysis(treatment, outcome, outcome_censor, direction);
    }

    return res;    

    }


// Find all child subgroups for a numerical biomarker when the outcome variable is continuous 
SingleSubgroup ContOutContBio(const std::vector<double> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter) {

    SingleSubgroup res;

    int i;

    // Total number of observations
    int n = biomarker.size();

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<double> biomarker_selected;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            n_parent_rows++;    
        }
    }

    // Structures to sort the columns by biomarker values
    vector<ddpair> treatment_sorted;
    vector<ddpair> outcome_sorted;
    vector<ddpair> outcome_censor_sorted;

    for(i = 0; i < n_parent_rows; i++) {
        treatment_sorted.push_back(ddpair(biomarker_selected[i], treatment_selected[i]));
        outcome_sorted.push_back(ddpair(biomarker_selected[i], outcome_selected[i]));
        outcome_censor_sorted.push_back(ddpair(biomarker_selected[i], outcome_censor_selected[i]));
    }    

    sort(treatment_sorted.begin(), treatment_sorted.end(), DDPairSortUp);
    sort(outcome_sorted.begin(), outcome_sorted.end(), DDPairSortUp);
    sort(outcome_censor_sorted.begin(), outcome_censor_sorted.end(), DDPairSortUp);

    double criterion_max=numeric_limits<double>::quiet_NaN(), value_max=numeric_limits<double>::quiet_NaN();
    int code = 1;//no partition found
    int sign=-1;
    int max_found = 0;
    int error_flag_left, error_flag_right;

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)

    double treatment_count_left = 0; 
    double treatment_count_right = 0;
    double control_count_left = 0;
    double control_count_right = 0;

    double treatment_sum_left = 0;
    double treatment_sum_right = 0;
    double control_sum_left = 0;
    double control_sum_right = 0;

    double treatment_sum_squares_left = 0; 
    double treatment_sum_squares_right = 0;
    double control_sum_squares_left = 0;
    double control_sum_squares_right = 0;

    double treatment_mean_left,treatment_var_left, control_mean_left, control_var_left, treatment_mean_right, treatment_var_right, control_mean_right, control_var_right, test_left, test_right, test_statistic, criterion, prom_estimate, prom_sd;

    for(i = 0; i < nmin-1; i++) {
        // Observations (rows) are sorted by biomarker values
        if (treatment_sorted[i].second == 0) {
            control_sum_left += outcome_sorted[i].second;
            control_sum_squares_left += outcome_sorted[i].second * outcome_sorted[i].second;
            control_count_left ++;
        }
        if (treatment_sorted[i].second == 1) {
            treatment_sum_left += outcome_sorted[i].second;
            treatment_sum_squares_left += outcome_sorted[i].second * outcome_sorted[i].second;
            treatment_count_left ++;
        }
    }

    for(i = nmin-1; i < n_parent_rows; i++) {
        if (treatment_sorted[i].second == 0) {
            control_sum_right += outcome_sorted[i].second;
            control_sum_squares_right += outcome_sorted[i].second * outcome_sorted[i].second;
            control_count_right ++;
        }
        if (treatment_sorted[i].second == 1) {
            treatment_sum_right += outcome_sorted[i].second;
            treatment_sum_squares_right += outcome_sorted[i].second * outcome_sorted[i].second;
            treatment_count_right ++;
        }
    }

    for(i = nmin-1; i < n_parent_rows - nmin; i++) {
        // If the current biomarker value is equal to the next value, recompute statistics
            if (treatment_sorted[i].second == 0) {
                control_sum_left += outcome_sorted[i].second;
                control_sum_squares_left += outcome_sorted[i].second * outcome_sorted[i].second;
                control_sum_right -= outcome_sorted[i].second;
                control_sum_squares_right -= outcome_sorted[i].second * outcome_sorted[i].second;
                control_count_left ++;
                control_count_right --;
            }
            if (treatment_sorted[i].second == 1) {
                treatment_sum_left += outcome_sorted[i].second;
                treatment_sum_squares_left += outcome_sorted[i].second * outcome_sorted[i].second;
                treatment_sum_right -= outcome_sorted[i].second;
                treatment_sum_squares_right -= outcome_sorted[i].second * outcome_sorted[i].second;
                treatment_count_left++;
                treatment_count_right--;
            }

        // If the current biomarker value is not equal to the next value, recompute statistics, save the value of the splitting criterion and biomarker value
        if (treatment_sorted[i].first != treatment_sorted[i + 1].first) {
            code = 0;
            treatment_mean_left = treatment_sum_left / treatment_count_left;
            treatment_var_left = (treatment_sum_squares_left / treatment_count_left - treatment_mean_left * treatment_mean_left) * treatment_count_left / (treatment_count_left - 1); 
            treatment_mean_right = treatment_sum_right / treatment_count_right;
            treatment_var_right = (treatment_sum_squares_right / treatment_count_right - treatment_mean_right * treatment_mean_right) * treatment_count_right / (treatment_count_right - 1); 
            control_mean_left = control_sum_left / control_count_left;
            control_var_left = (control_sum_squares_left / control_count_left - control_mean_left * control_mean_left) * control_count_left / (control_count_left - 1);
            control_mean_right = control_sum_right / control_count_right;
            control_var_right = (control_sum_squares_right / control_count_right - control_mean_right * control_mean_right) * control_count_right / (control_count_right - 1);

            error_flag_left = 0;
            error_flag_right = 0;

            test_left = TTestStatistic(treatment_mean_left,
                                                treatment_var_left,
                                                control_mean_left,
                                                control_var_left,
                                                treatment_count_left,
                                                control_count_left,
                                                direction, error_flag_left);
            test_right = TTestStatistic(treatment_mean_right,
                                                treatment_var_right,
                                                control_mean_right,
                                                control_var_right,
                                                treatment_count_right,
                                                control_count_right,
                                                direction, error_flag_right);
            criterion = CriterionFunction(test_left, test_right, criterion_type);

            // Find the maximum value of the splitting criterion if no error are detected
            if (error_flag_left == 0 && error_flag_right == 0) {  
                if (max_found == 0) {
                    criterion_max = criterion;
                    // Sorted by biomarker values
                    value_max = treatment_sorted[i].first;
                    max_found = 1;
                    if (test_left > test_right) {
                        sign = 1; // <=
                        test_statistic = test_left;
                        prom_estimate = (treatment_mean_left - control_mean_left);
                        prom_sd = sqrt((treatment_var_left * (treatment_count_left - 1) + control_var_left * (control_count_left - 1))/(treatment_count_left + control_count_left- 2));
                    } else {
                        sign = 2; // >             
                        test_statistic = test_right;
                        prom_estimate = (treatment_mean_right - control_mean_right);
                        prom_sd = sqrt((treatment_var_right * (treatment_count_right - 1) + control_var_right * (control_count_right - 1))/(treatment_count_right + control_count_right- 2));
                    }

                } else {
                    if (criterion_max < criterion) { 
                        criterion_max = criterion;      
                        // Sorted by biomarker values
                        value_max = treatment_sorted[i].first;
                        if (test_left > test_right) {
                            sign = 1; // <=
                            test_statistic = test_left;
                            prom_estimate =  (treatment_mean_left - control_mean_left);
                            prom_sd = sqrt((treatment_var_left * (treatment_count_left - 1) + control_var_left * (control_count_left - 1))/(treatment_count_left + control_count_left- 2));
                        } else {
                            sign = 2; // >             
                            test_statistic = test_right;
                            prom_estimate =  (treatment_mean_right - control_mean_right);
                            prom_sd = sqrt((treatment_var_right * (treatment_count_right - 1) + control_var_right * (control_count_right - 1))/(treatment_count_right + control_count_right- 2));
                        }
                    }
                }
            }
        }
    }

    // Save the split which corresponds to the maximal spitting criterion
    vector<double> value;
    value.clear();
    value.push_back(value_max);

    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows(parent_rows);

    // Left subgroup (x <= c)
    if (sign == 1) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] > value_max) subgroup_rows[i] = 0;
        }
    }

    // Right subgroup (x > c)
    if (sign == 2) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] <= value_max) subgroup_rows[i] = 0;
        }    
    }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    double pvalue = 1.0;
    if (test_statistic == test_statistic) pvalue = 1.0 - rcpp_pnorm(test_statistic);

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    
    res.prom_estimate = prom_estimate; 
    res.prom_sderr = -1.0;   
    res.prom_sd = prom_sd;   
    res.size_control = size_control;
    res.size_treatment = size_treatment;

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = sign;
    // Size of the current patient subgroup
    res.size = size_control + size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change) 
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change) 
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change) 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors) (this temporary value will change) 
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change) 
    res.terminal_subgroup = 0;

    return res;

}

// Find all child subgroups for a numerical biomarker when the outcome variable is continuous 
SingleSubgroup ContOutContBioANCOVA(const std::vector<double> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter, const ModelCovariates &model_covariates, const int &analysis_method) {

    SingleSubgroup res;

    int i, j;

    // Total number of observations
    int n = biomarker.size();

    // Extract the covariates
    std::vector<double> ancova_cov1, ancova_cov2, ancova_cov3, ancova_cov4, ancova_center;

    ancova_cov1 = model_covariates.cov1; 
    ancova_cov2 = model_covariates.cov2; 
    ancova_cov3 = model_covariates.cov3; 
    ancova_cov4 = model_covariates.cov4; 
    ancova_center = model_covariates.cov_class; 

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<double> biomarker_selected;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    vector<double> ancova_cov1_selected;
    vector<double> ancova_cov2_selected;
    vector<double> ancova_cov3_selected;
    vector<double> ancova_cov4_selected;    
    vector<double> ancova_center_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            ancova_cov1_selected.push_back(ancova_cov1[i]);
            ancova_cov2_selected.push_back(ancova_cov2[i]);
            ancova_cov3_selected.push_back(ancova_cov3[i]);
            ancova_cov4_selected.push_back(ancova_cov4[i]);           
            ancova_center_selected.push_back(ancova_center[i]);
            n_parent_rows++;    
        }
    }

    // List of sorted unique values of biomarker
    vector<double> unique = ListUniqueValues(biomarker_selected);
    sort(unique.begin(),unique.end());

    double criterion_max=numeric_limits<double>::quiet_NaN(), value_max=numeric_limits<double>::quiet_NaN();
    int code = 1;//no partition found
    int sign=-1;
    int max_found = 0;
    int error_flag_left, error_flag_right;

    error_flag_left = 0; error_flag_right = 0;

    double pvalue, test_statistic, criterion, prom_estimate, prom_sderr;
    int count_low, count_high;

    vector<double> ancova_treatment_left, ancova_outcome_left, ancova_cov1_left, ancova_cov2_left, ancova_cov3_left, ancova_cov4_left, ancova_center_left, ancova_treatment_right, ancova_outcome_right, ancova_cov1_right, ancova_cov2_right, ancova_cov3_right, ancova_cov4_right, ancova_center_right;

    ANCOVAResult ancova_result_left, ancova_result_right;

    ModelCovariates model_covariates_local;

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)

    for(i = 0; i < unique.size(); i++) {

        // Number of observations in biomarker-low and biomarker-high groups
        count_low = 0; count_high = 0;
        for(j = 0; j < biomarker_selected.size(); j++) {

            // Non-missing
            if (::isnan(biomarker_selected[j]) < 1.0) { 

                if (biomarker_selected[j] <= unique[i]) count_low++;
                if (biomarker_selected[j] > unique[i]) count_high++;

            }
        
        }

        if (count_low >= nmin && count_high >= nmin) {

            code = 0;

        ancova_treatment_left.clear(); ancova_outcome_left.clear(); ancova_cov1_left.clear(); ancova_cov2_left.clear(); ancova_cov3_left.clear(); ancova_cov4_left.clear(); ancova_center_left.clear(); ancova_treatment_right.clear(); ancova_outcome_right.clear(); ancova_cov1_right.clear(); ancova_cov2_right.clear(); ancova_cov3_right.clear(); ancova_cov4_right.clear(); ancova_center_right.clear(); 

            for(j = 0; j < biomarker_selected.size(); j++) {

                // Non-missing
                if (::isnan(biomarker_selected[j]) < 1.0) { 

                    if (biomarker_selected[j] <= unique[i]) {

                        ancova_treatment_left.push_back(treatment_selected[j]);
                        ancova_outcome_left.push_back(outcome_selected[j]);
                        ancova_cov1_left.push_back(ancova_cov1_selected[j]);
                        ancova_cov2_left.push_back(ancova_cov2_selected[j]);
                        ancova_cov3_left.push_back(ancova_cov3_selected[j]);
                        ancova_cov4_left.push_back(ancova_cov4_selected[j]);
                        ancova_center_left.push_back(ancova_center_selected[j]);

                    }
                    if (biomarker_selected[j] > unique[i]) {

                        ancova_treatment_right.push_back(treatment_selected[j]);
                        ancova_outcome_right.push_back(outcome_selected[j]);
                        ancova_cov1_right.push_back(ancova_cov1_selected[j]);
                        ancova_cov2_right.push_back(ancova_cov2_selected[j]);
                        ancova_cov3_right.push_back(ancova_cov3_selected[j]);
                        ancova_cov4_right.push_back(ancova_cov4_selected[j]);
                        ancova_center_right.push_back(ancova_center_selected[j]);

                    }

                }
            
            }

            model_covariates_local.cov1 = ancova_cov1_left; 
            model_covariates_local.cov2 = ancova_cov2_left; 
            model_covariates_local.cov3 = ancova_cov3_left; 
            model_covariates_local.cov4 = ancova_cov4_left; 
            model_covariates_local.cov_class = ancova_center_left; 

            ancova_result_left = ContANCOVA(ancova_treatment_left, ancova_outcome_left, model_covariates_local, analysis_method,direction);

            model_covariates_local.cov1 = ancova_cov1_right; 
            model_covariates_local.cov2 = ancova_cov2_right; 
            model_covariates_local.cov3 = ancova_cov3_right; 
            model_covariates_local.cov4 = ancova_cov4_right; 
            model_covariates_local.cov_class = ancova_center_right; 

            ancova_result_right = ContANCOVA(ancova_treatment_right, ancova_outcome_right, model_covariates_local, analysis_method, direction);

            // ancova_result_left = ANCOVAOld(ancova_treatment_left, ancova_outcome_left, ancova_cov1_left, ancova_cov2_left, ancova_center_left, analysis_method);
            // ancova_result_right = ANCOVAOld(ancova_treatment_right, ancova_outcome_right, ancova_cov1_right, ancova_cov2_right, ancova_center_right, analysis_method);

            criterion = CriterionFunction(ancova_result_left.test_stat, ancova_result_right.test_stat, criterion_type);

            // Find the maximum value of the splitting criterion if no error are detected
            if (error_flag_left == 0 && error_flag_right == 0) { 
                if (max_found == 0) {
                    criterion_max = criterion;
                    // Sorted by biomarker values
                    value_max = unique[i];
                    max_found = 1;
                    if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                        sign = 1; // <=
                        test_statistic = ancova_result_left.test_stat;
                        pvalue = ancova_result_left.pvalue;
                        prom_estimate = ancova_result_left.estimate;
                        prom_sderr = ancova_result_left.sderr;
                    } else {
                        sign = 2; // >             
                        test_statistic = ancova_result_right.test_stat;
                        pvalue = ancova_result_right.pvalue;
                        prom_estimate = ancova_result_right.estimate;
                        prom_sderr = ancova_result_right.sderr;
                    }

                } else {
                    if (criterion_max < criterion) { 
                        criterion_max = criterion;      
                        // Sorted by biomarker values
                        value_max = unique[i];
                       if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                            sign = 1; // <=
                            test_statistic = ancova_result_left.test_stat;
                            pvalue = ancova_result_left.pvalue;
                            prom_estimate = ancova_result_left.estimate;
                            prom_sderr = ancova_result_left.sderr;
                        } else {
                            sign = 2; // >             
                            test_statistic = ancova_result_right.test_stat;
                            pvalue = ancova_result_right.pvalue;
                            prom_estimate = ancova_result_right.estimate;
                            prom_sderr = ancova_result_right.sderr;
                        }
                    }
                }
            }

        }

    }

    // Save the split which corresponds to the maximal spitting criterion
    vector<double> value;
    value.clear();
    value.push_back(value_max);

    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows(parent_rows);

    // Left subgroup (x <= c)
    if (sign == 1) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] > value_max) subgroup_rows[i] = 0;
        }
    }

    // Right subgroup (x > c)
    if (sign == 2) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] <= value_max) subgroup_rows[i] = 0;
        }    
    }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    res.prom_estimate = prom_estimate;
    res.prom_sderr = prom_sderr;
    res.prom_sd = -1.0;   
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = sign;
    // Size of the current patient subgroup
    res.size = size_control + size_treatment;
    res.size_control = size_control;
    res.size_treatment = size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change) 
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change) 
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change) 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors) (this temporary value will change) 
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change) 
    res.terminal_subgroup = 0;

    return res;

}

// Find all child subgroups for a nominal biomarker when the outcome variable is binary  
SingleSubgroup BinOutNomBioANCOVA(const std::vector<double> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter, const NumericMatrix &cont_covariates, const NumericMatrix &class_covariates, const int &cont_covariates_flag, const int &class_covariates_flag) {

    SingleSubgroup res;

    int i, j, k;

    double temp;

    // Total number of observations
    int n = biomarker.size();

    // Total number of covariates
    int n_cont_covariates = cont_covariates.ncol(), n_class_covariates = class_covariates.ncol();

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<int> biomarker_selected;
    vector<double> biomarker_selected_double;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            biomarker_selected_double.push_back(biomarker[i] + 0.0);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            n_parent_rows++;    
        }
    }

    // List of all possible partitions

    // Compute the number of unique values 
    int nlevels = CountUniqueValues(biomarker_selected_double);
    int ncombinations = pow(2, nlevels) - 2;

    ublas::matrix<int> h(ncombinations, nlevels);

    for (i = 0; i < nlevels; i++) {
        for (j = 0; j < ncombinations; j++) {
            temp=floor((j+1.0)/(pow(2.0, nlevels-i-1)));
            if (temp/2.0==floor(temp/2.0)) h(j,i)=1; else h(j,i)=0;
        }
    }

    double criterion_max=numeric_limits<double>::quiet_NaN();
    int code = 1;//no partition found
    int max_found = 0;
    int error_flag_left, error_flag_right;
    int sign = -1;

    error_flag_left = 0; error_flag_right = 0;

    double pvalue, test_statistic, criterion, prom_estimate, prom_sderr;
    int count_low, count_high, k_low, k_high;

    vector<double> ancova_treatment_left, ancova_outcome_left, ancova_treatment_right, ancova_outcome_right, value;

    ANCOVAResult ancova_result_left, ancova_result_right;

    NumericMatrix cont_covariates_selected(n_parent_rows, n_cont_covariates), class_covariates_selected(n_parent_rows, n_class_covariates);

    k = 0;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            cont_covariates_selected(k, _) = cont_covariates(i, _);
            class_covariates_selected(k, _) = class_covariates(i, _);
            k++;
        }
    }

    // Vectors of left and right indices
    vector<int> left_index, right_index, value_int;

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)

    for(i = 0; i < ncombinations; i++) {

        left_index.clear(); right_index.clear();

        ancova_treatment_left.clear(); ancova_outcome_left.clear(); ancova_treatment_right.clear(); ancova_outcome_right.clear();  

        // Left and right indices
        for(j = 0; j < nlevels; j++) {
            if (h(i, j) == 1) left_index.push_back(j + 1); else right_index.push_back(j + 1);
        }

        count_low = 0; count_high = 0;

        for(j = 0; j < biomarker_selected.size(); j++) {

            // Non-missing
            if (isnan(biomarker_selected[j]) < 1.0) { 

                if (Included(left_index, biomarker_selected[j]) == 1) {

                    ancova_treatment_left.push_back(treatment_selected[j]);
                    ancova_outcome_left.push_back(outcome_selected[j]);
                    count_low++;

                }
                if (Included(right_index, biomarker_selected[j]) == 1) {

                    ancova_treatment_right.push_back(treatment_selected[j]);
                    ancova_outcome_right.push_back(outcome_selected[j]);
                    count_high++;

                }

            }

        }

        if (count_low >= nmin && count_high >= nmin) {

            code = 0;

            NumericMatrix cont_covariates_left(count_low, n_cont_covariates), class_covariates_left(count_low, n_class_covariates), cont_covariates_right(count_high, n_cont_covariates), class_covariates_right(count_high, n_class_covariates);      

            ancova_treatment_left.clear(); ancova_outcome_left.clear(); ancova_treatment_right.clear(); ancova_outcome_right.clear(); 

            k_low = 0; k_high = 0;

            for(j = 0; j < biomarker_selected.size(); j++) {

                // Non-missing
                if (isnan(biomarker_selected[j]) < 1.0) { 

                if (Included(left_index, biomarker_selected[j]) == 1) {

                        ancova_treatment_left.push_back(treatment_selected[j]);
                        ancova_outcome_left.push_back(outcome_selected[j]);
                        cont_covariates_left(k_low, _) = cont_covariates_selected(j, _);
                        class_covariates_left(k_low, _) = class_covariates_selected(j, _);
                        k_low++;

                    }

                if (Included(right_index, biomarker_selected[j]) == 1) {


                        ancova_treatment_right.push_back(treatment_selected[j]);
                        ancova_outcome_right.push_back(outcome_selected[j]);
                        cont_covariates_right(k_high, _) = cont_covariates_selected(j, _);
                        class_covariates_right(k_high, _) = class_covariates_selected(j, _);
                        k_high++;

                    }

                }
            
            }

            ancova_result_left = BinANCOVA(ancova_treatment_left, ancova_outcome_left, cont_covariates_left, class_covariates_left, direction, cont_covariates_flag, class_covariates_flag);

            ancova_result_right = BinANCOVA(ancova_treatment_right, ancova_outcome_right, cont_covariates_right, class_covariates_right, direction, cont_covariates_flag, class_covariates_flag);

            criterion = CriterionFunction(ancova_result_left.test_stat, ancova_result_right.test_stat, criterion_type);

            // Find the maximum value of the splitting criterion if no error are detected
            if (error_flag_left == 0 && error_flag_right == 0) { 
                if (max_found == 0) {
                    criterion_max = criterion;
                    max_found = 1;
                    if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                        sign = 3; // ==
                        value_int = left_index;
                        test_statistic = ancova_result_left.test_stat;
                        pvalue = ancova_result_left.pvalue;
                        prom_estimate = ancova_result_left.estimate;
                        prom_sderr = ancova_result_left.sderr;
                    } else {
                        sign = 3; // ==           
                        value_int = right_index;  
                        test_statistic = ancova_result_right.test_stat;
                        pvalue = ancova_result_right.pvalue;
                        prom_estimate = ancova_result_right.estimate;
                        prom_sderr = ancova_result_right.sderr;
                    }

                } else {
                    if (criterion_max < criterion) { 
                        criterion_max = criterion;      
                       if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                            sign = 3; // ==
                            value_int = left_index;
                            test_statistic = ancova_result_left.test_stat;
                            pvalue = ancova_result_left.pvalue;
                            prom_estimate = ancova_result_left.estimate;
                            prom_sderr = ancova_result_left.sderr;
                        } else {
                            sign = 3; // ==             
                            value_int = right_index;
                            test_statistic = ancova_result_right.test_stat;
                            pvalue = ancova_result_right.pvalue;
                            prom_estimate = ancova_result_right.estimate;
                            prom_sderr = ancova_result_right.sderr;
                        }
                    }
                }
            }

        }

    }

    // Save the optimal partition
    value.clear();
    for (i = 0; i < value_int.size(); i++) value.push_back(value_int[i] + 0.0);

    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows(parent_rows);

    // Selected subgroup
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && Included(value_int, biomarker[i]) == 0) subgroup_rows[i] = 0;
        }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    res.prom_estimate = prom_estimate;
    res.prom_sderr = prom_sderr;   
    res.prom_sd = -1.0;   
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = 3;
    // Size of the current patient subgroup
    res.size_control = size_control;
    res.size_treatment = size_treatment;
    res.size = size_control + size_treatment;   
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change) 
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change) 
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change) 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors) (this temporary value will change) 
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change) 
    res.terminal_subgroup = 0;

    return res;

}

// Find all child subgroups for a nominal biomarker when the outcome variable is survival  
SingleSubgroup SurvOutNomBioANCOVA(const std::vector<double> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter, const NumericMatrix &cont_covariates, const NumericMatrix &class_covariates, const int &cont_covariates_flag, const int &class_covariates_flag) {

    SingleSubgroup res;

    int i, j, k;

    double temp;

    // Total number of observations
    int n = biomarker.size();

    // Total number of covariates
    int n_cont_covariates = cont_covariates.ncol(), n_class_covariates = class_covariates.ncol();

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<int> biomarker_selected;
    vector<double> biomarker_selected_double;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            biomarker_selected_double.push_back(biomarker[i] + 0.0);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            n_parent_rows++;    
        }
    }

    // List of all possible partitions

    // Compute the number of unique values 
    int nlevels = CountUniqueValues(biomarker_selected_double);
    int ncombinations = pow(2, nlevels) - 2;

    ublas::matrix<int> h(ncombinations, nlevels);

    for (i = 0; i < nlevels; i++) {
        for (j = 0; j < ncombinations; j++) {
            temp=floor((j+1.0)/(pow(2.0, nlevels-i-1)));
            if (temp/2.0==floor(temp/2.0)) h(j,i)=1; else h(j,i)=0;
        }
    }

    double criterion_max=numeric_limits<double>::quiet_NaN();
    int sign;
    int code = 1;//no partition found
    int max_found = 0;
    int error_flag_left, error_flag_right;

    error_flag_left = 0; error_flag_right = 0;

    double pvalue, test_statistic, criterion, prom_estimate, prom_sderr;
    int count_low, count_high, k_low, k_high;

    vector<double> ancova_treatment_left, ancova_outcome_left, ancova_treatment_right, ancova_outcome_right, ancova_censor_left, ancova_censor_right, value;

    ANCOVAResult ancova_result_left, ancova_result_right;

    NumericMatrix cont_covariates_selected(n_parent_rows, n_cont_covariates), class_covariates_selected(n_parent_rows, n_class_covariates);

    k = 0;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            cont_covariates_selected(k, _) = cont_covariates(i, _);
            class_covariates_selected(k, _) = class_covariates(i, _);
            k++;
        }
    }

    // Vectors of left and right indices
    vector<int> left_index, right_index, value_int;

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)

    for(i = 0; i < ncombinations; i++) {

        left_index.clear(); right_index.clear();

        ancova_treatment_left.clear(); ancova_outcome_left.clear(); ancova_treatment_right.clear(); ancova_outcome_right.clear(); ancova_censor_left.clear(), ancova_censor_right.clear(); 

        // Left and right indices
        for(j = 0; j < nlevels; j++) {
            if (h(i, j) == 1) left_index.push_back(j + 1); else right_index.push_back(j + 1);
        }

        count_low = 0; count_high = 0;

        for(j = 0; j < biomarker_selected.size(); j++) {

            // Non-missing
            if (isnan(biomarker_selected[j]) < 1.0) { 

                if (Included(left_index, biomarker_selected[j]) == 1) {

                    ancova_treatment_left.push_back(treatment_selected[j]);
                    ancova_outcome_left.push_back(outcome_selected[j]);
                    ancova_censor_left.push_back(outcome_censor_selected[j]);
                    count_low++;

                }
                if (Included(right_index, biomarker_selected[j]) == 1) {

                    ancova_treatment_right.push_back(treatment_selected[j]);
                    ancova_outcome_right.push_back(outcome_selected[j]);
                    ancova_censor_right.push_back(outcome_censor_selected[j]);
                    count_high++;

                }

            }

        }

        if (count_low >= nmin && count_high >= nmin) {

            code = 0;

            NumericMatrix cont_covariates_left(count_low, n_cont_covariates), class_covariates_left(count_low, n_class_covariates), cont_covariates_right(count_high, n_cont_covariates), class_covariates_right(count_high, n_class_covariates);      

            ancova_treatment_left.clear(); ancova_outcome_left.clear(); ancova_treatment_right.clear(); ancova_outcome_right.clear(); ancova_censor_left.clear(), ancova_censor_right.clear(); 

            k_low = 0; k_high = 0;

            for(j = 0; j < biomarker_selected.size(); j++) {

                // Non-missing
                if (isnan(biomarker_selected[j]) < 1.0) { 

                if (Included(left_index, biomarker_selected[j]) == 1) {

                        ancova_treatment_left.push_back(treatment_selected[j]);
                        ancova_outcome_left.push_back(outcome_selected[j]);
                        ancova_censor_left.push_back(outcome_censor_selected[j]);
                        cont_covariates_left(k_low, _) = cont_covariates_selected(j, _);
                        class_covariates_left(k_low, _) = class_covariates_selected(j, _);
                        k_low++;

                    }

                if (Included(right_index, biomarker_selected[j]) == 1) {

                        ancova_treatment_right.push_back(treatment_selected[j]);
                        ancova_outcome_right.push_back(outcome_selected[j]);
                        ancova_censor_right.push_back(outcome_censor_selected[j]);                        
                        cont_covariates_right(k_high, _) = cont_covariates_selected(j, _);
                        class_covariates_right(k_high, _) = class_covariates_selected(j, _);
                        k_high++;

                    }

                }
            
            }

            ancova_result_left = SurvANCOVA(ancova_treatment_left, ancova_outcome_left, ancova_censor_left, cont_covariates_left, class_covariates_left, direction, cont_covariates_flag, class_covariates_flag);

            ancova_result_right = SurvANCOVA(ancova_treatment_right, ancova_outcome_right, ancova_censor_right, cont_covariates_right, class_covariates_right, direction, cont_covariates_flag, class_covariates_flag);

            criterion = CriterionFunction(ancova_result_left.test_stat, ancova_result_right.test_stat, criterion_type);

            // Find the maximum value of the splitting criterion if no error are detected
            if (error_flag_left == 0 && error_flag_right == 0) { 
                if (max_found == 0) {
                    criterion_max = criterion;
                    max_found = 1;
                    if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                        sign = 3; // ==
                        value_int = left_index;
                        test_statistic = ancova_result_left.test_stat;
                        pvalue = ancova_result_left.pvalue;
                        prom_estimate = ancova_result_left.estimate;
                        prom_sderr = ancova_result_left.sderr;
                    } else {
                        sign = 3; // ==           
                        value_int = right_index;  
                        test_statistic = ancova_result_right.test_stat;
                        pvalue = ancova_result_right.pvalue;
                        prom_estimate = ancova_result_right.estimate;
                        prom_sderr = ancova_result_right.sderr;
                    }

                } else {
                    if (criterion_max < criterion) { 
                        criterion_max = criterion;      
                       if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                            sign = 3; // ==
                            value_int = left_index;
                            test_statistic = ancova_result_left.test_stat;
                            pvalue = ancova_result_left.pvalue;
                            prom_estimate = ancova_result_left.estimate;
                            prom_sderr = ancova_result_left.sderr;
                        } else {
                            sign = 3; // ==             
                            value_int = right_index;
                            test_statistic = ancova_result_right.test_stat;
                            pvalue = ancova_result_right.pvalue;
                            prom_estimate = ancova_result_right.estimate;
                            prom_sderr = ancova_result_right.sderr;
                        }
                    }
                }
            }

        }

    }

    // Save the optimal partition
    value.clear();
    for (i = 0; i < value_int.size(); i++) value.push_back(value_int[i] + 0.0);

    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows(parent_rows);

    // Selected subgroup
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && Included(value_int, biomarker[i]) == 0) subgroup_rows[i] = 0;
        }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    res.prom_estimate = prom_estimate;
    res.prom_sderr = prom_sderr;   
    res.prom_sd = -1.0;   
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = 3;
    // Size of the current patient subgroup
    res.size_control = size_control;
    res.size_treatment = size_treatment;
    res.size = size_control + size_treatment;   
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change) 
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change) 
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change) 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors) (this temporary value will change) 
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change) 
    res.terminal_subgroup = 0;

    return res;

}


// Find all child subgroups for a numerical biomarker when the outcome variable is binary 
SingleSubgroup BinOutContBioANCOVA(const std::vector<double> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter, const NumericMatrix &cont_covariates, const NumericMatrix &class_covariates, const int &cont_covariates_flag, const int &class_covariates_flag) {

    SingleSubgroup res;

    int i, j, k;

    // Total number of observations
    int n = biomarker.size();

    // Total number of covariates
    int n_cont_covariates = cont_covariates.ncol(), n_class_covariates = class_covariates.ncol();

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<double> biomarker_selected;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            n_parent_rows++;    
        }
    }

    // List of sorted unique values of biomarker
    vector<double> unique = ListUniqueValues(biomarker_selected);
    sort(unique.begin(),unique.end());

    double criterion_max=numeric_limits<double>::quiet_NaN(), value_max=numeric_limits<double>::quiet_NaN();
    int code = 1;//no partition found
    int sign=-1;
    int max_found = 0;
    int error_flag_left, error_flag_right;

    error_flag_left = 0; error_flag_right = 0;

    double pvalue, test_statistic, criterion, prom_estimate, prom_sderr;
    int count_low, count_high, k_low, k_high;

    vector<double> ancova_treatment_left, ancova_outcome_left, ancova_treatment_right, ancova_outcome_right;

    ANCOVAResult ancova_result_left, ancova_result_right;

    NumericMatrix cont_covariates_selected(n_parent_rows, n_cont_covariates), class_covariates_selected(n_parent_rows, n_class_covariates);

    k = 0;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            cont_covariates_selected(k, _) = cont_covariates(i, _);
            class_covariates_selected(k, _) = class_covariates(i, _);
            k++;
        }
    }

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)

    for(i = 0; i < unique.size(); i++) {

        // Number of observations in biomarker-low and biomarker-high groups
        count_low = 0; count_high = 0;
        for(j = 0; j < biomarker_selected.size(); j++) {

            // Non-missing
            if (::isnan(biomarker_selected[j]) < 1.0) { 

                if (biomarker_selected[j] <= unique[i]) count_low++;
                if (biomarker_selected[j] > unique[i]) count_high++;

            }
        
        }

        if (count_low >= nmin && count_high >= nmin) {

            code = 0;

            NumericMatrix cont_covariates_left(count_low, n_cont_covariates), class_covariates_left(count_low, n_class_covariates), cont_covariates_right(count_high, n_cont_covariates), class_covariates_right(count_high, n_class_covariates);      

            ancova_treatment_left.clear(); ancova_outcome_left.clear(); ancova_treatment_right.clear(); ancova_outcome_right.clear(); 

            k_low = 0; k_high = 0;

            for(j = 0; j < biomarker_selected.size(); j++) {

                // Non-missing
                if (::isnan(biomarker_selected[j]) < 1.0) { 

                    if (biomarker_selected[j] <= unique[i]) {

                        ancova_treatment_left.push_back(treatment_selected[j]);
                        ancova_outcome_left.push_back(outcome_selected[j]);
                        cont_covariates_left(k_low, _) = cont_covariates_selected(j, _);
                        class_covariates_left(k_low, _) = class_covariates_selected(j, _);
                        k_low++;

                    }
                    if (biomarker_selected[j] > unique[i]) {

                        ancova_treatment_right.push_back(treatment_selected[j]);
                        ancova_outcome_right.push_back(outcome_selected[j]);
                        cont_covariates_right(k_high, _) = cont_covariates_selected(j, _);
                        class_covariates_right(k_high, _) = class_covariates_selected(j, _);
                        k_high++;

                    }

                }
            
            }

            ancova_result_left = BinANCOVA(ancova_treatment_left, ancova_outcome_left, cont_covariates_left, class_covariates_left, direction, cont_covariates_flag, class_covariates_flag);

            ancova_result_right = BinANCOVA(ancova_treatment_right, ancova_outcome_right, cont_covariates_right, class_covariates_right, direction, cont_covariates_flag, class_covariates_flag);

            criterion = CriterionFunction(ancova_result_left.test_stat, ancova_result_right.test_stat, criterion_type);

            // Find the maximum value of the splitting criterion if no error are detected
            if (error_flag_left == 0 && error_flag_right == 0) { 
                if (max_found == 0) {
                    criterion_max = criterion;
                    // Sorted by biomarker values
                    value_max = unique[i];
                    max_found = 1;
                    if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                        sign = 1; // <=
                        test_statistic = ancova_result_left.test_stat;
                        pvalue = ancova_result_left.pvalue;
                        prom_estimate = ancova_result_left.estimate;
                        prom_sderr = ancova_result_left.sderr;
                    } else {
                        sign = 2; // >             
                        test_statistic = ancova_result_right.test_stat;
                        pvalue = ancova_result_right.pvalue;
                        prom_estimate = ancova_result_right.estimate;
                        prom_sderr = ancova_result_right.sderr;
                    }

                } else {
                    if (criterion_max < criterion) { 
                        criterion_max = criterion;      
                        // Sorted by biomarker values
                        value_max = unique[i];
                       if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                            sign = 1; // <=
                            test_statistic = ancova_result_left.test_stat;
                            pvalue = ancova_result_left.pvalue;
                            prom_estimate = ancova_result_left.estimate;
                            prom_sderr = ancova_result_left.sderr;
                        } else {
                            sign = 2; // >             
                            test_statistic = ancova_result_right.test_stat;
                            pvalue = ancova_result_right.pvalue;
                            prom_estimate = ancova_result_right.estimate;
                            prom_sderr = ancova_result_right.sderr;
                        }
                    }
                }
            }

        }

    }

    // Save the split which corresponds to the maximal spitting criterion
    vector<double> value;
    value.clear();
    value.push_back(value_max);

    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows(parent_rows);

    // Left subgroup (x <= c)
    if (sign == 1) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] > value_max) subgroup_rows[i] = 0;
        }
    }

    // Right subgroup (x > c)
    if (sign == 2) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] <= value_max) subgroup_rows[i] = 0;
        }    
    }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    res.prom_estimate = prom_estimate;
    res.prom_sderr = prom_sderr;
    res.prom_sd = -1.0;   
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = sign;
    // Size of the current patient subgroup
    res.size = size_control + size_treatment;
    res.size_control = size_control;
    res.size_treatment = size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change) 
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change) 
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change) 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors) (this temporary value will change) 
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change) 
    res.terminal_subgroup = 0;

    return res;

}

// Find all child subgroups for a numerical biomarker when the outcome variable is survival 
SingleSubgroup SurvOutContBioANCOVA(const std::vector<double> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter, const NumericMatrix &cont_covariates, const NumericMatrix &class_covariates, const int &cont_covariates_flag, const int &class_covariates_flag) {

    SingleSubgroup res;

    int i, j, k;

    // Total number of observations
    int n = biomarker.size();

    // Total number of covariates
    int n_cont_covariates = cont_covariates.ncol(), n_class_covariates = class_covariates.ncol();

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<double> biomarker_selected;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            n_parent_rows++;    
        }
    }

    // List of sorted unique values of biomarker
    vector<double> unique = ListUniqueValues(biomarker_selected);
    sort(unique.begin(),unique.end());

    double criterion_max=numeric_limits<double>::quiet_NaN(), value_max=numeric_limits<double>::quiet_NaN();
    int code = 1;//no partition found
    int sign=-1;
    int max_found = 0;
    int error_flag_left, error_flag_right;

    error_flag_left = 0; error_flag_right = 0;

    double pvalue, test_statistic, criterion, prom_estimate, prom_sderr;
    int count_low, count_high, k_low, k_high;

    vector<double> ancova_treatment_left, ancova_outcome_left, ancova_treatment_right, ancova_outcome_right, ancova_censor_left, ancova_censor_right;

    ANCOVAResult ancova_result_left, ancova_result_right;

    NumericMatrix cont_covariates_selected(n_parent_rows, n_cont_covariates), class_covariates_selected(n_parent_rows, n_class_covariates);

    k = 0;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            cont_covariates_selected(k, _) = cont_covariates(i, _);
            class_covariates_selected(k, _) = class_covariates(i, _);
            k++;
        }
    }

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)

    for(i = 0; i < unique.size(); i++) {

        // Number of observations in biomarker-low and biomarker-high groups
        count_low = 0; count_high = 0;
        for(j = 0; j < biomarker_selected.size(); j++) {

            // Non-missing
            if (::isnan(biomarker_selected[j]) < 1.0) { 

                if (biomarker_selected[j] <= unique[i]) count_low++;
                if (biomarker_selected[j] > unique[i]) count_high++;

            }
        
        }

        if (count_low >= nmin && count_high >= nmin) {

            code = 0;

            NumericMatrix cont_covariates_left(count_low, n_cont_covariates), class_covariates_left(count_low, n_class_covariates), cont_covariates_right(count_high, n_cont_covariates), class_covariates_right(count_high, n_class_covariates);      

            ancova_treatment_left.clear(); ancova_outcome_left.clear(); ancova_treatment_right.clear(); ancova_outcome_right.clear(); ancova_censor_left.clear(); ancova_censor_right.clear();

            k_low = 0; k_high = 0;

            for(j = 0; j < biomarker_selected.size(); j++) {

                // Non-missing
                if (::isnan(biomarker_selected[j]) < 1.0) { 

                    if (biomarker_selected[j] <= unique[i]) {

                        ancova_treatment_left.push_back(treatment_selected[j]);
                        ancova_outcome_left.push_back(outcome_selected[j]);
                        ancova_censor_left.push_back(outcome_censor_selected[j]);
                        cont_covariates_left(k_low, _) = cont_covariates_selected(j, _);
                        class_covariates_left(k_low, _) = class_covariates_selected(j, _);
                        k_low++;

                    }
                    if (biomarker_selected[j] > unique[i]) {

                        ancova_treatment_right.push_back(treatment_selected[j]);
                        ancova_outcome_right.push_back(outcome_selected[j]);
                        ancova_censor_right.push_back(outcome_censor_selected[j]);                      
                        cont_covariates_right(k_high, _) = cont_covariates_selected(j, _);
                        class_covariates_right(k_high, _) = class_covariates_selected(j, _);
                        k_high++;

                    }

                }
            
            }

            ancova_result_left = SurvANCOVA(ancova_treatment_left, ancova_outcome_left, ancova_censor_left, cont_covariates_left, class_covariates_left, direction, cont_covariates_flag, class_covariates_flag);

            ancova_result_right = SurvANCOVA(ancova_treatment_right, ancova_outcome_right, ancova_censor_right, cont_covariates_right, class_covariates_right, direction, cont_covariates_flag, class_covariates_flag);

            criterion = CriterionFunction(ancova_result_left.test_stat, ancova_result_right.test_stat, criterion_type);

            // Find the maximum value of the splitting criterion if no error are detected
            if (error_flag_left == 0 && error_flag_right == 0) { 
                if (max_found == 0) {
                    criterion_max = criterion;
                    // Sorted by biomarker values
                    value_max = unique[i];
                    max_found = 1;
                    if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                        sign = 1; // <=
                        test_statistic = ancova_result_left.test_stat;
                        pvalue = ancova_result_left.pvalue;
                        prom_estimate = ancova_result_left.estimate;
                        prom_sderr = ancova_result_left.sderr;
                    } else {
                        sign = 2; // >             
                        test_statistic = ancova_result_right.test_stat;
                        pvalue = ancova_result_right.pvalue;
                        prom_estimate = ancova_result_right.estimate;
                        prom_sderr = ancova_result_right.sderr;
                    }

                } else {
                    if (criterion_max < criterion) { 
                        criterion_max = criterion;      
                        // Sorted by biomarker values
                        value_max = unique[i];
                       if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                            sign = 1; // <=
                            test_statistic = ancova_result_left.test_stat;
                            pvalue = ancova_result_left.pvalue;
                            prom_estimate = ancova_result_left.estimate;
                            prom_sderr = ancova_result_left.sderr;
                        } else {
                            sign = 2; // >             
                            test_statistic = ancova_result_right.test_stat;
                            pvalue = ancova_result_right.pvalue;
                            prom_estimate = ancova_result_right.estimate;
                            prom_sderr = ancova_result_right.sderr;
                        }
                    }
                }
            }

        }

    }

    // Save the split which corresponds to the maximal spitting criterion
    vector<double> value;
    value.clear();
    value.push_back(value_max);

    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows(parent_rows);

    // Left subgroup (x <= c)
    if (sign == 1) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] > value_max) subgroup_rows[i] = 0;
        }
    }

    // Right subgroup (x > c)
    if (sign == 2) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] <= value_max) subgroup_rows[i] = 0;
        }    
    }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    res.prom_estimate = prom_estimate;
    res.prom_sderr = prom_sderr;
    res.prom_sd = -1.0;   
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = sign;
    // Size of the current patient subgroup
    res.size = size_control + size_treatment;
    res.size_control = size_control;
    res.size_treatment = size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change) 
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change) 
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change) 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors) (this temporary value will change) 
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change) 
    res.terminal_subgroup = 0;

    return res;

}

SingleSubgroup ContOutNomBioANCOVA(const std::vector<int> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter, const ModelCovariates &model_covariates, const int &analysis_method) {

    SingleSubgroup res;

    int i, j;

    double temp;

    // Total number of observations
    int n = biomarker.size();

    // Extract the covariates
    std::vector<double> ancova_cov1, ancova_cov2, ancova_cov3, ancova_cov4, ancova_center;

    ancova_cov1 = model_covariates.cov1; 
    ancova_cov2 = model_covariates.cov2; 
    ancova_cov3 = model_covariates.cov3; 
    ancova_cov4 = model_covariates.cov4; 
    ancova_center = model_covariates.cov_class; 

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<int> biomarker_selected;
    vector<double> biomarker_selected_double;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    vector<double> ancova_cov1_selected;
    vector<double> ancova_cov2_selected;
    vector<double> ancova_cov3_selected;
    vector<double> ancova_cov4_selected;
    vector<double> ancova_center_selected;
    for(i = 0; i < n; i++) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            biomarker_selected_double.push_back(biomarker[i] + 0.0);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            ancova_cov1_selected.push_back(ancova_cov1[i]);
            ancova_cov2_selected.push_back(ancova_cov2[i]);
            ancova_cov3_selected.push_back(ancova_cov3[i]);
            ancova_cov4_selected.push_back(ancova_cov4[i]);
            ancova_center_selected.push_back(ancova_center[i]);
            n_parent_rows++;    
        }
    }

    // List of all possible partitions

    // Compute the number of unique values 
    int nlevels = CountUniqueValues(biomarker_selected_double);
    int ncombinations = pow(2, nlevels) - 2;

    ublas::matrix<int> h(ncombinations, nlevels);

    for (i = 0; i < nlevels; i++) {
        for (j = 0; j < ncombinations; j++) {
            temp=floor((j+1.0)/(pow(2.0, nlevels-i-1)));
            if (temp/2.0==floor(temp/2.0)) h(j,i)=1; else h(j,i)=0;
        }
    }

    double criterion_max=numeric_limits<double>::quiet_NaN();
    int code = 1;//no partition found
    int max_found = 0;
    int error_flag_left, error_flag_right;
    int sign=-1;

    error_flag_left = 0; error_flag_right = 0;

    ModelCovariates model_covariates_local;

    double pvalue, test_statistic, criterion, prom_estimate, prom_sderr;
    int count_low, count_high;

    vector<double> ancova_treatment_left, ancova_outcome_left, ancova_cov1_left, ancova_cov2_left, ancova_cov3_left, ancova_cov4_left, ancova_center_left, ancova_treatment_right, ancova_outcome_right, ancova_cov1_right, ancova_cov2_right, ancova_cov3_right, ancova_cov4_right, ancova_center_right, value;

    ANCOVAResult ancova_result_left, ancova_result_right;

    // Vectors of left and right indices
    vector<int> left_index, right_index, value_int;

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)

    for(i = 0; i < ncombinations; i++) {

        left_index.clear(); right_index.clear();

        ancova_treatment_left.clear(); ancova_outcome_left.clear(); ancova_cov1_left.clear(); ancova_cov2_left.clear(); ancova_cov3_left.clear(); ancova_cov4_left.clear(); ancova_center_left.clear(); ancova_treatment_right.clear(); ancova_outcome_right.clear(); ancova_cov1_right.clear(); ancova_cov2_right.clear(); ancova_cov3_right.clear(); ancova_cov4_right.clear(); ancova_center_right.clear(); 

        // Left and right indices
        for(j = 0; j < nlevels; j++) {
            if (h(i, j) == 1) left_index.push_back(j + 1); else right_index.push_back(j + 1);
        }

        count_low = 0; count_high = 0;

        for(j = 0; j < biomarker_selected.size(); j++) {

            // Non-missing
            if (isnan(biomarker_selected[j]) < 1.0) { 

                if (Included(left_index, biomarker_selected[j]) == 1) {

                    ancova_treatment_left.push_back(treatment_selected[j]);
                    ancova_outcome_left.push_back(outcome_selected[j]);
                    ancova_cov1_left.push_back(ancova_cov1_selected[j]);
                    ancova_cov2_left.push_back(ancova_cov2_selected[j]);
                    ancova_cov3_left.push_back(ancova_cov3_selected[j]);
                    ancova_cov4_left.push_back(ancova_cov4_selected[j]);
                    ancova_center_left.push_back(ancova_center_selected[j]);
                    count_low++;

                }
                if (Included(right_index, biomarker_selected[j]) == 1) {

                    ancova_treatment_right.push_back(treatment_selected[j]);
                    ancova_outcome_right.push_back(outcome_selected[j]);
                    ancova_cov1_right.push_back(ancova_cov1_selected[j]);
                    ancova_cov2_right.push_back(ancova_cov2_selected[j]);
                    ancova_cov3_right.push_back(ancova_cov3_selected[j]);
                    ancova_cov4_right.push_back(ancova_cov4_selected[j]);
                    ancova_center_right.push_back(ancova_center_selected[j]);
                    count_high++;

                }

            }

        }


        if (count_low >= nmin && count_high >= nmin) {

            code = 0;

            model_covariates_local.cov1 = ancova_cov1_left; 
            model_covariates_local.cov2 = ancova_cov2_left; 
            model_covariates_local.cov3 = ancova_cov3_left; 
            model_covariates_local.cov4 = ancova_cov4_left; 
            model_covariates_local.cov_class = ancova_center_left; 

            ancova_result_left = ContANCOVA(ancova_treatment_left, ancova_outcome_left, model_covariates_local, analysis_method, direction);

            model_covariates_local.cov1 = ancova_cov1_right; 
            model_covariates_local.cov2 = ancova_cov2_right; 
            model_covariates_local.cov3 = ancova_cov3_right; 
            model_covariates_local.cov4 = ancova_cov4_right; 
            model_covariates_local.cov_class = ancova_center_right; 

            ancova_result_right = ContANCOVA(ancova_treatment_right, ancova_outcome_right, model_covariates_local, analysis_method, direction);

            criterion = CriterionFunction(ancova_result_left.test_stat, ancova_result_right.test_stat, criterion_type);

            // Find the maximum value of the splitting criterion if no error are detected
            if (error_flag_left == 0 && error_flag_right == 0) { 
                if (max_found == 0) {
                    criterion_max = criterion;
                    max_found = 1;
                    if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                        sign = 3; // ==
                        value_int = left_index;
                        test_statistic = ancova_result_left.test_stat;
                        pvalue = ancova_result_left.pvalue;
                        prom_estimate = ancova_result_left.estimate;
                        prom_sderr = ancova_result_left.sderr;
                    } else {
                        sign = 3; // ==           
                        value_int = right_index;  
                        test_statistic = ancova_result_right.test_stat;
                        pvalue = ancova_result_right.pvalue;
                        prom_estimate = ancova_result_right.estimate;
                        prom_sderr = ancova_result_right.sderr;
                    }

                } else {
                    if (criterion_max < criterion) { 
                        criterion_max = criterion;      
                       if (ancova_result_left.test_stat > ancova_result_right.test_stat) {
                            sign = 3; // ==
                            value_int = left_index;
                            test_statistic = ancova_result_left.test_stat;
                            pvalue = ancova_result_left.pvalue;
                            prom_estimate = ancova_result_left.estimate;
                            prom_sderr = ancova_result_left.sderr;
                        } else {
                            sign = 3; // ==             
                            value_int = right_index;
                            test_statistic = ancova_result_right.test_stat;
                            pvalue = ancova_result_right.pvalue;
                            prom_estimate = ancova_result_right.estimate;
                            prom_sderr = ancova_result_right.sderr;
                        }
                    }
                }
            }

        }

    }

    // Save the optimal partition
    value.clear();
    for (i = 0; i < value_int.size(); i++) value.push_back(value_int[i] + 0.0);

    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows(parent_rows);

    // Selected subgroup
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && Included(value_int, biomarker[i]) == 0) subgroup_rows[i] = 0;
        }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    res.prom_estimate = prom_estimate;
    res.prom_sderr = prom_sderr;   
    res.prom_sd = -1.0;   
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = 3;
    // Size of the current patient subgroup
    res.size_control = size_control;
    res.size_treatment = size_treatment;
    res.size = size_control + size_treatment;   
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change) 
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change) 
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change) 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors) (this temporary value will change) 
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change) 
    res.terminal_subgroup = 0;

    return res;

}


struct TestData{
    double treatment_count;
    double control_count;

    double treatment_sum;
    double control_sum;

    double treatment_sum_squares;
    double control_sum_squares;
    TestData():treatment_count(0.0),control_count(0.0),treatment_sum(0.0),control_sum(0.0),treatment_sum_squares(0.0),control_sum_squares(0.0){}
};

map<int,TestData> NomPar(const vector<int> &marker, const std::vector<double> &treatment, const std::vector<double> &outcome)
{
    map<int,TestData> dict;
    int x;
    double t,y;
    for (int i = 0; i < marker.size();++i)
    {
        x = marker[i];
        t = treatment[i];
        y = outcome[i];
        if(dict.find(x)==dict.end())
        {
            TestData d;
            if( t == 1){
                d.treatment_count = 1;
                d.treatment_sum = y;
                d.treatment_sum_squares = y*y;
            } else {
                d.control_count=1;
                d.control_sum = y;
                d.control_sum_squares = y*y;
            }
            dict.insert(pair<int,TestData>(x,d));
        } else {
            TestData *d = &dict[x];
            if( t == 1){
                d->treatment_count++;
                d->treatment_sum += y;
                d->treatment_sum_squares += y*y;
            } else {
                d->control_count++;
                d->control_sum += y;
                d->control_sum_squares += y*y;
            }
        }
    }
    return dict;
}

void NomSelectInc(vector<int> &sel)
{
    int n = sel.size();
    for (int i=0; i<n; ++i)
    {
        if (sel[i] == 1)
            sel[i] =0;
        else{
            sel[i]=1;
            break;
        }
    }
}

// Find all child subgroups for a nominal biomarker when the outcome variable is continuous 
SingleSubgroup ContOutNomBio(const std::vector<int> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter) {

    
    SingleSubgroup res;

    int i, j;

    // Total number of observations
    int n = biomarker.size();

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<int> biomarker_selected;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            n_parent_rows++;
        }
    }

    // Structures to sort the columns by biomarker values
    double criterion_max=numeric_limits<double>::quiet_NaN();
    int code = 1; // no partition was found
    int sign=-1;
    int max_found = 0;
    int error_flag_right, error_flag_left;

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)
    map<int,TestData> dat = NomPar(biomarker_selected,treatment_selected,outcome_selected);
    int tot = dat.size();
    vector<int> nom,sel(tot,0),sel_best;
    nom.reserve(tot);
    for(map<int,TestData>::iterator it = dat.begin(); it != dat.end(); ++it) {
      nom.push_back(it->first);
  }


    double treatment_mean_left,treatment_var_left, control_mean_left, control_var_left, treatment_mean_right, treatment_var_right, control_mean_right, control_var_right, test_left, test_right, test_statistic, criterion, prom_estimate, prom_sd;


    for(i = 0; i < pow(2,tot-1)-1; ++i) {
        double treatment_count_left = 0;
        double treatment_count_right = 0;
        double control_count_left = 0;
        double control_count_right = 0;

        double treatment_sum_left = 0;
        double treatment_sum_right = 0;
        double control_sum_left = 0;
        double control_sum_right = 0;

        double treatment_sum_squares_left = 0;
        double treatment_sum_squares_right = 0;
        double control_sum_squares_left = 0;
        double control_sum_squares_right = 0;


        NomSelectInc(sel);
        for(j = 0; j <tot;++j){

            TestData d = dat[nom[j]];
            if (sel[j]==1){
                treatment_count_right += d.treatment_count;
                control_count_right += d.control_count;

                treatment_sum_right += d.treatment_sum;
                control_sum_right += d.control_sum;

                treatment_sum_squares_right += d.treatment_sum_squares;
                control_sum_squares_right += d.control_sum_squares;
            } else {
                treatment_count_left += d.treatment_count;
                control_count_left += d.control_count;

                treatment_sum_left += d.treatment_sum;
                control_sum_left += d.control_sum;

                treatment_sum_squares_left += d.treatment_sum_squares;
                control_sum_squares_left += d.control_sum_squares;
            }

        }
        if(treatment_count_left+control_count_left < nmin || treatment_count_right+control_count_right < nmin)
            continue;
        code=0;

        treatment_mean_left = treatment_sum_left / treatment_count_left;
        treatment_var_left = (treatment_sum_squares_left / treatment_count_left - treatment_mean_left * treatment_mean_left) * treatment_count_left / (treatment_count_left - 1);
        treatment_mean_right = treatment_sum_right / treatment_count_right;
        treatment_var_right = (treatment_sum_squares_right / treatment_count_right - treatment_mean_right * treatment_mean_right) * treatment_count_right / (treatment_count_right - 1);
        control_mean_left = control_sum_left / control_count_left;
        control_var_left = (control_sum_squares_left / control_count_left - control_mean_left * control_mean_left) * control_count_left / (control_count_left - 1);
        control_mean_right = control_sum_right / control_count_right;
        control_var_right = (control_sum_squares_right / control_count_right - control_mean_right * control_mean_right) * control_count_right / (control_count_right - 1);

        error_flag_left = 0;
        error_flag_right = 0;

        test_left = TTestStatistic(treatment_mean_left,
                                            treatment_var_left,
                                            control_mean_left,
                                            control_var_left,
                                            treatment_count_left,
                                            control_count_left,
                                            direction, error_flag_left);
        test_right = TTestStatistic(treatment_mean_right,
                                            treatment_var_right,
                                            control_mean_right,
                                            control_var_right,
                                            treatment_count_right,
                                            control_count_right,
                                            direction, error_flag_right);


        criterion = CriterionFunction(test_left, test_right, criterion_type);

        // Find the maximum value of the splitting criterion if no error are detected 
        //if (error_flag_left == 0 && error_flag_right == 0) { 
            if (max_found == 0) {
                criterion_max = criterion;
                sel_best = sel;
                // Sorted by biomarker values
                max_found = 1;
                if (test_left > test_right) {
                    sign = 0; // <=
                    test_statistic = test_left;
                    prom_estimate =  (treatment_mean_left - control_mean_left);
                    prom_sd = sqrt((treatment_var_left * (treatment_count_left - 1) + control_var_left * (control_count_left - 1))/(treatment_count_left + control_count_left- 2));                    
                } else {
                    sign = 1; // >
                    test_statistic = test_right;
                    prom_estimate =  (treatment_mean_right - control_mean_right);
                    prom_sd = sqrt((treatment_var_right * (treatment_count_right - 1) + control_var_right * (control_count_right - 1))/(treatment_count_right + control_count_right- 2));
                }

            } else {
                if (criterion_max < criterion) {
                    criterion_max = criterion;
                    // Sorted by biomarker values
                    sel_best = sel;
                    if (test_left > test_right) {
                        sign = 0; // <=
                        test_statistic = test_left;
                        prom_estimate =  (treatment_mean_left - control_mean_left);
                        prom_sd = sqrt((treatment_var_left * (treatment_count_left - 1) + control_var_left * (control_count_left - 1))/(treatment_count_left + control_count_left- 2));                        
                    } else {
                        sign = 1; // >
                        test_statistic = test_right;
                        prom_estimate =  (treatment_mean_right - control_mean_right);
                        prom_sd = sqrt((treatment_var_right * (treatment_count_right - 1) + control_var_right * (control_count_right - 1))/(treatment_count_right + control_count_right- 2));
                    }
                }
        }
    }

    vector<int> subgroup_rows(parent_rows);
    vector<double> value;
    if (code == 0){
        // Save the split which corresponds to the maximal spitting criterion
        vector<int> value_int;
        value.clear();
        int tar = sign;
        for(i = 0; i < tot;++i){
            if (sel_best[i] == tar){
                value.push_back(nom[i]);
                value_int.push_back(nom[i]);
            }
        }
        int val_n = value_int.size();

        // Indexes of rows in the original data set that define the current patient subgroup

        // Left subgroup (x <= c)

        bool tt = false;
        for(i = 0; i < n; i++) {
            tt = false;
            for(j = 0; j<val_n ;++j)
                if (biomarker[i] == value_int[j]){
                    tt = true;
                    break;
                }
            if ( !tt) subgroup_rows[i] = 0;
        }
    }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    double pvalue = 1.0;
    if (test_statistic == test_statistic) pvalue = 1.0 - rcpp_pnorm(test_statistic);

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    
    res.prom_estimate = prom_estimate; 
    res.prom_sderr = -1.0;   
    res.prom_sd = prom_sd;   
    res.size_control = size_control;
    res.size_treatment = size_treatment;

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = 3;
    // Size of the current patient subgroup
    res.size = size_control + size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change)
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change)
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change)
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;
    // Error code for testing (0 if no errors) (this temporary value will change)
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change)
    res.terminal_subgroup = 0;
    return res;

}


// Find all child subgroups for a nominal biomarker when the outcome variable is binary. 
SingleSubgroup BinOutNomBio(const std::vector<int> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter) {
    
    SingleSubgroup res;

    int i, j;

    // Total number of observations
    int n = biomarker.size();

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<int> biomarker_selected;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            n_parent_rows++;
        }
    }

    // Structures to sort the columns by biomarker values
    double criterion_max=numeric_limits<double>::quiet_NaN();
    int code = 1;//no partition found
    int sign=-1;
    int max_found = 0;
    int error_flag_left, error_flag_right;

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)
    map<int,TestData> dat = NomPar(biomarker_selected,treatment_selected,outcome_selected);
    int tot = dat.size();
    vector<int> nom,sel(tot,0),sel_best;
    nom.reserve(tot);
    for(map<int,TestData>::iterator it = dat.begin(); it != dat.end(); ++it) {
      nom.push_back(it->first);
    }


    double test_left, test_right, test_statistic, criterion, prom_estimate;

    test_statistic = 0;
    criterion = 0;
    criterion_max = 0;

    for(i = 0; i < pow(2,tot-1)-1; ++i) {

        double treatment_count_left = 0;
        double treatment_count_right = 0;
        double control_count_left = 0;
        double control_count_right = 0;

        double treatment_sum_left = 0;
        double treatment_sum_right = 0;
        double control_sum_left = 0;
        double control_sum_right = 0;


        NomSelectInc(sel);

        for(j = 0; j <tot;++j){
            TestData d = dat[nom[j]];
            if (sel[j]==1){
                treatment_count_right += d.treatment_count;
                control_count_right += d.control_count;

                treatment_sum_right += d.treatment_sum;
                control_sum_right += d.control_sum;

            } else {
                treatment_count_left += d.treatment_count;
                control_count_left += d.control_count;

                treatment_sum_left += d.treatment_sum;
                control_sum_left += d.control_sum;

            }

        }

        if(treatment_count_left+control_count_left < nmin || treatment_count_right+control_count_right < nmin)
            continue;
        code=0;

        error_flag_left = 0;
        error_flag_right = 0;
        test_left = PropTestStatistic(treatment_sum_left, control_sum_left,
                                      treatment_count_left, control_count_left, direction, error_flag_left);
        test_right = PropTestStatistic(treatment_sum_right, control_sum_right,
                                      treatment_count_right, control_count_right, direction, error_flag_right);

        criterion = CriterionFunction(test_left, test_right, criterion_type);

        // Find the maximum value of the splitting criterion if no error are detected
        //if (error_flag_left == 0 && error_flag_right == 0) { 
            if (max_found == 0) {
                criterion_max = criterion;
                sel_best = sel;
                // Sorted by biomarker values
                max_found = 1;
                if (test_left > test_right) {
                    sign = 0; // <=
                    test_statistic = test_left;
                    prom_estimate =  (treatment_sum_left / treatment_count_left - control_sum_left / control_count_left);
                } else {
                    sign = 1; // >
                    test_statistic = test_right;
                    prom_estimate =  (treatment_sum_right / treatment_count_right - control_sum_right / control_count_right);
                }

            } else {
                if (criterion_max < criterion) {
                    criterion_max = criterion;
                    // Sorted by biomarker values
                    sel_best = sel;
                    if (test_left > test_right) {
                        sign = 0; // <=
                        test_statistic = test_left;
                        prom_estimate =  (treatment_sum_left / treatment_count_left - control_sum_left / control_count_left);
                    } else {
                        sign = 1; // >
                        test_statistic = test_right;
                        prom_estimate =  (treatment_sum_right / treatment_count_right - control_sum_right / control_count_right);
                    }
                }
            }
        //}
    }


    vector<int> subgroup_rows(parent_rows);
    // Save the split which corresponds to the maximal spitting criterion
    vector<double> value;

    if (code == 0){
        vector<int> value_int;
        value.clear();
        int tar = sign;
        for(i = 0; i < tot;++i){
            if (sel_best[i] == tar){
                value.push_back(nom[i]);
                value_int.push_back(nom[i]);
            }
        }
        int val_n = value_int.size();

        // Indexes of rows in the original data set that define the current patient subgroup

        // Left subgroup (x <= c)

        bool tt = false;
        for(i = 0; i < n; i++) {
            tt = false;
            for(j = 0; j<val_n ;++j)
                if (biomarker[i] == value_int[j]){
                    tt = true;
                    break;
                }
            if ( !tt) subgroup_rows[i] = 0;
        }
    }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    double pvalue = 1.0;
    if (test_statistic == test_statistic) pvalue = 1.0 - rcpp_pnorm(test_statistic);

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    
    res.prom_estimate = prom_estimate; 
    res.prom_sderr = -1.0;   
    res.prom_sd = -1.0;   
    res.size_control = size_control;
    res.size_treatment = size_treatment;

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = 3;
    // Size of the current patient subgroup
    res.size = size_control + size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change)
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change)
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change)
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;
    // Error code for testing (0 if no errors) (this temporary value will change)
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change)
    res.terminal_subgroup = 0;

    return res;

}

// Find all child subgroups for a numerical biomarker when the outcome variable is binary 
SingleSubgroup BinOutContBio(const std::vector<double> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter) {

    
    SingleSubgroup res;

    int i;

    // Total number of observations
    int n = biomarker.size();

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<double> biomarker_selected;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            n_parent_rows++;    
        }
    }

    // Structures to sort the columns by biomarker values
    vector<ddpair> treatment_sorted;
    vector<ddpair> outcome_sorted;
    vector<ddpair> outcome_censor_sorted;

    for(i = 0; i < n_parent_rows; i++) {
        treatment_sorted.push_back(ddpair(biomarker_selected[i], treatment_selected[i]));
        outcome_sorted.push_back(ddpair(biomarker_selected[i], outcome_selected[i]));
        outcome_censor_sorted.push_back(ddpair(biomarker_selected[i], outcome_censor_selected[i]));
    }    

    sort(treatment_sorted.begin(), treatment_sorted.end(), DDPairSortUp);
    sort(outcome_sorted.begin(), outcome_sorted.end(), DDPairSortUp);
    sort(outcome_censor_sorted.begin(), outcome_censor_sorted.end(), DDPairSortUp);

    double criterion_max=numeric_limits<double>::quiet_NaN(), value_max=numeric_limits<double>::quiet_NaN();
    int code = 1;//no partition found
    int sign=-1;
    int max_found = 0;

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)

    double treatment_count_left = 0; 
    double treatment_count_right = 0;
    double control_count_left = 0;
    double control_count_right = 0;

    double treatment_event_count_left = 0;  
    double treatment_event_count_right = 0;
    double control_event_count_left = 0;
    double control_event_count_right = 0;
    int error_flag_right, error_flag_left;

    double test_left, test_right, test_statistic, criterion, prom_estimate;

    for(i = 0; i < nmin-1; i++) {
        // Observations (rows) are sorted by biomarker values
        if (treatment_sorted[i].second == 0) {
            control_event_count_left += outcome_sorted[i].second;
            control_count_left ++;
        }
        if (treatment_sorted[i].second == 1) {
            treatment_event_count_left += outcome_sorted[i].second;
            treatment_count_left ++;
        }
    }

    for(i = nmin-1; i < n_parent_rows; i++) {
        if (treatment_sorted[i].second == 0) {
            control_event_count_right += outcome_sorted[i].second;
            control_count_right ++;
        }
        if (treatment_sorted[i].second == 1) {
            treatment_event_count_right += outcome_sorted[i].second;
            treatment_count_right ++;
        }
    }

    for(i = nmin-1; i < n_parent_rows - nmin; i++) {
        // If the current biomarker value is equal to the next value, recompute statistics
            if (treatment_sorted[i].second == 0) {
                control_event_count_left += outcome_sorted[i].second;
                control_event_count_right -= outcome_sorted[i].second;
                control_count_left ++;
                control_count_right --;
            }
            if (treatment_sorted[i].second == 1) {
                treatment_event_count_left += outcome_sorted[i].second;
                treatment_event_count_right -= outcome_sorted[i].second;
                treatment_count_left ++;
                treatment_count_right --;
            }

        // If the current biomarker value is not equal to the next value, recompute statistics, save the value of the splitting criterion and biomarker value
        if (treatment_sorted[i].first != treatment_sorted[i + 1].first) {
            code = 0;
            error_flag_left = 0;
            error_flag_right = 0;
            test_left = PropTestStatistic(treatment_event_count_left, control_event_count_left, 
                                          treatment_count_left, control_count_left, direction, error_flag_left);
            test_right = PropTestStatistic(treatment_event_count_right, control_event_count_right, 
                                          treatment_count_right, control_count_right, direction, error_flag_right);
            criterion = CriterionFunction(test_left, test_right, criterion_type);

            // Find the maximum value of the splitting criterion if no error are detected
            if (error_flag_left == 0 && error_flag_right == 0) { 
                if (max_found == 0) {
                    criterion_max = criterion;
                    // Sorted by biomarker values
                    value_max = treatment_sorted[i].first;
                    max_found = 1;
                    if (test_left > test_right) {
                        sign = 1; // <=
                        test_statistic = test_left;
                        prom_estimate =  (treatment_event_count_left / treatment_count_left - control_event_count_left / control_count_left); 
                    } else {
                        sign = 2; // >             
                        test_statistic = test_right;
                        prom_estimate =  (treatment_event_count_right / treatment_count_right - control_event_count_right / control_count_right);
                    }

                } else {
                    if (criterion_max < criterion) { 
                        criterion_max = criterion;      
                        // Sorted by biomarker values
                        value_max = treatment_sorted[i].first;
                        if (test_left > test_right) {
                            sign = 1; // <=
                            test_statistic = test_left;
                            prom_estimate =  (treatment_event_count_left / treatment_count_left - control_event_count_left / control_count_left);
                        } else {
                            sign = 2; // >             
                            test_statistic = test_right;
                            prom_estimate =  (treatment_event_count_right / treatment_count_right - control_event_count_right / control_count_right);
                        }
                    }
                }
            }
        }
    }

    // Save the split which corresponds to the maximal spitting criterion
    vector<double> value;
    value.clear();
    value.push_back(value_max);

    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows(parent_rows);

    // Left subgroup (x <= c)
    if (sign == 1) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] > value_max) subgroup_rows[i] = 0;
        }
    }

    // Right subgroup (x > c)
    if (sign == 2) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] <= value_max) subgroup_rows[i] = 0;
        }    
    }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    double pvalue = 1.0;
    if (test_statistic==test_statistic) pvalue = 1.0 - rcpp_pnorm(test_statistic);

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    
    res.prom_estimate = prom_estimate; 
    res.prom_sderr = -1.0;   
    res.prom_sd = -1.0;
    res.size_control = size_control;
    res.size_treatment = size_treatment;

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = sign;
    // Size of the current patient subgroup
    res.size = size_control + size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change) 
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change) 
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change) 
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;     
    // Error code for testing (0 if no errors) (this temporary value will change) 
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change) 
    res.terminal_subgroup = 0;

    return res;

}

// Find all child subgroups for a numerical biomarker when the outcome variable is time-to-event 
SingleSubgroup SurvOutContBio(const std::vector<double> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter) {

    
    SingleSubgroup res;

    int i;

    // Total number of observations
    int n = biomarker.size();

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<double> biomarker_selected;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            n_parent_rows++;
        }
    }

    // Structures to sort the columns by biomarker values
    vector<ddpair> treatment_sorted;
    vector<ddpair> outcome_sorted;
    vector<ddpair> outcome_censor_sorted;

    for(i = 0; i < n_parent_rows; i++) {
        treatment_sorted.push_back(ddpair(biomarker_selected[i], treatment_selected[i]));
        outcome_sorted.push_back(ddpair(biomarker_selected[i], outcome_selected[i]));
        outcome_censor_sorted.push_back(ddpair(biomarker_selected[i], outcome_censor_selected[i]));
    }

    sort(treatment_sorted.begin(), treatment_sorted.end(), DDPairSortUp);
    sort(outcome_sorted.begin(), outcome_sorted.end(), DDPairSortUp);
    sort(outcome_censor_sorted.begin(), outcome_censor_sorted.end(), DDPairSortUp);

    double criterion_max=numeric_limits<double>::quiet_NaN(), value_max=numeric_limits<double>::quiet_NaN(), prom_estimate;
    int code = 1;//no partition found
    int sign=-1;
    int max_found = 0;

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)

    vector<double> treatment_left, treatment_right,
            outcome_left, outcome_right,
            outcome_c_left, outcome_c_right;
    treatment_left.reserve(n);treatment_left.reserve(n);
    outcome_left.reserve(n);outcome_right.reserve(n);
    outcome_c_left.reserve(n);outcome_c_right.reserve(n);

    for(i = 0; i < nmin-1; i++) {
        // Observations (rows) are sorted by biomarker values
        treatment_left.push_back(treatment_sorted[i].second);
        outcome_left.push_back(outcome_sorted[i].second);
        outcome_c_left.push_back(outcome_censor_sorted[i].second);
    }

    for(i = nmin-1; i < n_parent_rows; i++) {
        treatment_right.push_back(treatment_sorted[i].second);
        outcome_right.push_back(outcome_sorted[i].second);
        outcome_c_right.push_back(outcome_censor_sorted[i].second);
    }

    double test_left, test_right, test_statistic, criterion;
    for(i = nmin-1; i < n_parent_rows - nmin; i++) {

                treatment_left.push_back(treatment_sorted[i].second);
                outcome_left.push_back(outcome_sorted[i].second);
                outcome_c_left.push_back(outcome_censor_sorted[i].second);
                treatment_right.erase(treatment_right.begin());
                outcome_right.erase(outcome_right.begin());
                outcome_c_right.erase(outcome_c_right.begin());

        // If the current biomarker value is not equal to the next value, recompute statistics, save the value of the splitting criterion and biomarker value
        if (treatment_sorted[i].first != treatment_sorted[i + 1].first) {
            code = 0;
 
            test_left = LRTest(outcome_left,outcome_c_left,treatment_left, direction);
            test_right = LRTest(outcome_right,outcome_c_right,treatment_right, direction);
            criterion = CriterionFunction(test_left, test_right, criterion_type);

            // Find the maximum value of splitting criterion qq
            if (max_found == 0) {
                criterion_max = criterion;
                // Sorted by biomarker values
                value_max = treatment_sorted[i].first;
                max_found = 1;
                if (test_left > test_right) {
                    sign = 1; // <=
                    test_statistic = test_left;
                    prom_estimate = HazardRatio(outcome_left, outcome_c_left, treatment_left, direction);
                } else {
                    sign = 2; // >
                    test_statistic = test_right;
                    prom_estimate = HazardRatio(outcome_right, outcome_c_right, treatment_right, direction);
                }

            } else {
                if (criterion_max < criterion) {
                    criterion_max = criterion;
                    // Sorted by biomarker values
                    value_max = treatment_sorted[i].first;
                    if (test_left > test_right) {
                        sign = 1; // <=
                        test_statistic = test_left;
                        prom_estimate = HazardRatio(outcome_left, outcome_c_left, treatment_left, direction);
                    } else {
                        sign = 2; // >
                        test_statistic = test_right;
                        prom_estimate = HazardRatio(outcome_right, outcome_c_right, treatment_right, direction);
                    }
                }
            }
        }
    }

    // Save the split which corresponds to the maximal spitting criterion
    vector<double> value;
    value.clear();
    value.push_back(value_max);

    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows(parent_rows);

    // Left subgroup (x <= c)
    if (sign == 1) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] > value_max) subgroup_rows[i] = 0;
        }
    }

    // Right subgroup (x > c)
    if (sign == 2) {
    for(i = 0; i < n; i++) {
            if (parent_rows[i] == 1 && biomarker[i] <= value_max) subgroup_rows[i] = 0;
        }
    }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    double pvalue = 1.0;
    if (test_statistic == test_statistic) pvalue = 1.0 - rcpp_pnorm(test_statistic);

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    
    res.prom_estimate = prom_estimate; 
    res.prom_sderr = -1.0;   
    res.prom_sd = -1.0;   
    res.size_control = size_control;
    res.size_treatment = size_treatment;

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = sign;
    // Size of the current patient subgroup
    res.size = size_control + size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change)
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change)
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change)
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;
    // Error code for testing (0 if no errors) (this temporary value will change)
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change)
    res.terminal_subgroup = 0;

    return res;

}

// Find all child subgroups for a nominal biomarker when the outcome variable is time-to-event 
SingleSubgroup SurvOutNomBio(const std::vector<int> &biomarker, const std::vector<double> &treatment, const std::vector<double> &outcome, const std::vector<double> &outcome_censor, const std::vector<int> &parent_rows, const int &nmin, const int &direction, const int &criterion_type, const double &adj_parameter) {

    
    SingleSubgroup res;

    int i, j;

    // Total number of observations
    int n = biomarker.size();

    // Select the observations (rows) included in the parent group
    int n_parent_rows = 0;
    vector<int> biomarker_selected;
    vector<double> treatment_selected;
    vector<double> outcome_selected;
    vector<double> outcome_censor_selected;
    for(i = 0; i < n; ++i) {
        if (parent_rows[i] == 1) {
            biomarker_selected.push_back(biomarker[i]);
            treatment_selected.push_back(treatment[i]);
            outcome_selected.push_back(outcome[i]);
            outcome_censor_selected.push_back(outcome_censor[i]);
            n_parent_rows++;
        }
    }

    // Structures to sort the columns by biomarker values


    double criterion_max=numeric_limits<double>::quiet_NaN(), prom_estimate;
    int code = 1;//no partition found
    int sign=-1;
    int max_found = 0;

    // Compute the splitting criterion for all possible splits (starting with the smallest left subgroup and ending with the smallest right subgroup)

    map<int,TestData> dat = NomPar(biomarker_selected,treatment_selected,outcome_selected);
    int tot = dat.size();
    vector<int> nom,sel(tot,0),sel_best;
    nom.reserve(tot);
    for(map<int,TestData>::iterator it = dat.begin(); it != dat.end(); ++it) {
      nom.push_back(it->first);
    }

    vector<int> biomarker_sorted;
    vector<double> treatment_sorted, outcome_sorted, outcome_censor_sorted;
    biomarker_sorted.reserve(n_parent_rows);treatment_sorted.reserve(n_parent_rows);
    outcome_sorted.reserve(n_parent_rows);outcome_censor_sorted.reserve(n_parent_rows);
    vector<vector<int>::iterator> bio_begin(tot),bio_end(tot);
    vector<vector<double>::iterator> treatment_begin(tot),treatment_end(tot),outcome_begin(tot),outcome_end(tot),outcome_c_begin(tot),outcome_c_end(tot);

    for(i = 0; i < tot;++i){
        bio_begin[i] = biomarker_sorted.end();
        treatment_begin[i] = treatment_sorted.end();
        outcome_begin[i] = outcome_sorted.end();
        outcome_c_begin[i] = outcome_censor_sorted.end();

        for(j = 0; j<n_parent_rows;++j){
            if(biomarker_selected[j]==nom[i]){
                biomarker_sorted.push_back(biomarker_selected[j]);
                treatment_sorted.push_back(treatment_selected[j]);
                outcome_sorted.push_back(outcome_selected[j]);
                outcome_censor_sorted.push_back(outcome_censor_selected[j]);
            }
        }

        bio_end[i] = biomarker_sorted.end();
        treatment_end[i] = treatment_sorted.end();
        outcome_end[i] = outcome_sorted.end();
        outcome_c_end[i] = outcome_censor_sorted.end();

    }

    double test_left, test_right, test_statistic, criterion;

    vector<double> outcome_left, outcome_c_left, outcome_right, outcome_c_right;
    vector<double> treatment_left, treatment_right;
    outcome_left.reserve(n_parent_rows); outcome_c_left.reserve(n_parent_rows); outcome_right.reserve(n_parent_rows); outcome_c_right.reserve(n_parent_rows);
    treatment_left.reserve(n_parent_rows); treatment_right.reserve(n_parent_rows);

    for(i = 0; i < pow(2,tot-1)-1; ++i) {

        outcome_left.erase(outcome_left.begin(),outcome_left.end()); outcome_c_left.erase(outcome_c_left.begin(),outcome_c_left.end());
        outcome_right.erase(outcome_right.begin(),outcome_right.end()); outcome_c_right.erase(outcome_c_right.begin(),outcome_c_right.end());
        treatment_left.erase(treatment_left.begin(),treatment_left.end()); treatment_right.erase(treatment_right.begin(),treatment_right.end());

        NomSelectInc(sel);

        for(j = 0; j <tot;++j){
            if (sel[j]==1){
                treatment_right.insert(treatment_right.end(),treatment_begin[j],treatment_end[j] );
                outcome_right.insert(outcome_right.end(),outcome_begin[j],outcome_end[j] );
                outcome_c_right.insert(outcome_c_right.end(),outcome_c_begin[j],outcome_c_end[j] );
            } else {
                treatment_left.insert(treatment_left.end(),treatment_begin[j],treatment_end[j] );
                outcome_left.insert(outcome_left.end(),outcome_begin[j],outcome_end[j] );
                outcome_c_left.insert(outcome_c_left.end(),outcome_c_begin[j],outcome_c_end[j] );
            }

        }
        if(treatment_left.empty() || treatment_right.empty())
            continue;
        if(treatment_left.size()<nmin || treatment_right.size()<nmin)
            continue;
        code=0;

        test_left = LRTest(outcome_left,outcome_c_left,treatment_left, direction);
        test_right = LRTest(outcome_right,outcome_c_right,treatment_right, direction);

        criterion = CriterionFunction(test_left, test_right, criterion_type);
        // Find the maximum value of splitting criterion
        if (max_found == 0) {
            criterion_max = criterion;
            sel_best = sel;
            // Sorted by biomarker values
            max_found = 1;
            if (test_left > test_right) {
                sign = 0; // <=
                test_statistic = test_left;
                prom_estimate = HazardRatio(outcome_left, outcome_c_left, treatment_left, direction);

            } else {
                sign = 1; // >
                test_statistic = test_right;
                prom_estimate = HazardRatio(outcome_right, outcome_c_right, treatment_right, direction);
            }

        } else {
            if (criterion_max < criterion) {
                criterion_max = criterion;
                // Sorted by biomarker values
                sel_best = sel;
                if (test_left > test_right) {
                    sign = 0; // <=
                    test_statistic = test_left;
                    prom_estimate = HazardRatio(outcome_left, outcome_c_left, treatment_left, direction);

                } else {
                    sign = 1; // >
                    test_statistic = test_right;
                    prom_estimate = HazardRatio(outcome_right, outcome_c_right, treatment_right, direction);
                }
            }
        }
    }


    // Save the split which corresponds to the maximal spitting criterion
    vector<double> value;
    // Indexes of rows in the original data set that define the current patient subgroup
    vector<int> subgroup_rows(parent_rows);
    if (code == 0){
        vector<int> value_int;
        value.clear();
        int tar = sign;
        for(i = 0; i < tot;++i){
            if (sel_best[i] == tar){
                value.push_back(nom[i]);
                value_int.push_back(nom[i]);
            }
        }
        int val_n = value_int.size();

        // Left subgroup (x <= c)

        bool tt = false;
        for(i = 0; i < n; i++) {
            tt = false;
            for(j = 0; j<val_n ;++j)
                if (biomarker[i] == value_int[j]){
                    tt = true;
                    break;
                }
            if ( !tt) subgroup_rows[i] = 0;
        }
    }

    // Size of the current patient subgroup
    int size_control = 0, size_treatment = 0;
    for(i = 0; i < n; ++i) {
        if (treatment[i] == 0) size_control += subgroup_rows[i]; else size_treatment += subgroup_rows[i];
    }

    double criterion_pvalue = 1.0;
    if (criterion_max==criterion_max)
        criterion_pvalue = CriterionPvalue(criterion_max, criterion_type); 

    double adjusted_criterion_pvalue;
    adjusted_criterion_pvalue = Log1minusPSupp(criterion_pvalue, criterion_max, adj_parameter, 2); 

    double pvalue = 1.0;
    if (test_statistic==test_statistic) pvalue = 1.0 - rcpp_pnorm(test_statistic);

    // Save the results

    // Splitting criterion
    res.criterion = criterion_max;
    // Splitting criterion on p-value scale
    res.criterion_pvalue = criterion_pvalue;
    // Splitting criterion on p-value scale with a local multiplicity adjustment
    res.adjusted_criterion_pvalue = adjusted_criterion_pvalue;
    res.test_statistic = test_statistic;
    res.pvalue = pvalue;
    // This temporary value will change
    res.adjusted_pvalue = -1.0;    
    res.prom_estimate = prom_estimate; 
    res.prom_sderr = -1.0;   
    res.prom_sd = -1.0;   
    res.size_control = size_control;
    res.size_treatment = size_treatment;

    // Vector of biomarker values to define the current patient subgroup
    res.value = value;
    // 1 if <=, 2 if >, 3 if =
    res.sign = 3;
    // Size of the current patient subgroup
    res.size = size_control + size_treatment;
    // Indexes of rows in the original data set that define the current patient subgroup
    res.subgroup_rows = subgroup_rows;
    // Index of biomarker used in the current patient subgroup (this temporary value will change)
    res.biomarker_index = 1;
    // Level of the current patient subgroup (this temporary value will change)
    res.level  = 1;
    // Parent index for the current patient subgroup (this temporary value will change)
    res.parent_index = 1;
    // Indexes of child subgroups for the current patient subgroup (this temporary value will change)
    vector<int> child_index;
    child_index.clear();
    child_index.push_back(1);
    res.child_index = child_index;
    // Error code for testing (0 if no errors) (this temporary value will change)
    res.code = code;
    // Is the current patient subgroup terminal (0/1) (this temporary value will change)
    res.terminal_subgroup = 0;

    return res;

}

// Recursive function for computing variable importance for all biomarkers
void ComputeVarImp(const vector<SingleSubgroup> &sub, vector<double> &imp, int &nterm, int &subnterm, vector<int> &signat){
    int l=0, totl=0;
    for(int i = 0; i<sub.size();++i){
        bool skip = false;
        for(int j =0; j < signat.size();++j){
            if (sub[i].signat == signat[j]){
                skip = true;
                break;
            }
        }
        signat.push_back(sub[i].signat);
        if(skip){
            continue;
        }
        if (!sub[i].subgroups.empty()){
            ComputeVarImp(sub[i].subgroups,imp,nterm,l,signat);
            totl += l;
        }
        if(sub[i].terminal_subgroup == 1){
            ++nterm;
            ++totl;
            l = 1;
        }
        //imp[sub[i].biomarker_index-1] -= log(sub[i].adjusted_criterion_pvalue)*l;
        imp[sub[i].biomarker_index-1] -= sub[i].adjusted_criterion_pvalue * l;
    }
    subnterm =totl;
}

// Find all patient subgroups 
std::vector<SingleSubgroup> FindSubgroups(const vector<double> df, const ModelCovariates &model_covariates,  
            const int ncol, const int nrow, const std::vector<int> &parent_rows, const int &nmin, const int &direction, 
            const int &criterion_type, const int &outcome_type, const int &analysis_method, 
            const vector<int> &biomarker_type, vector<double> &adj_parameter, const int &width, const int &depth, 
            const int &cur_depth, vector<int> depth_hist, const std::vector<double> gamma, double parent_pvalue, 
            double parent_test_statistic, const NumericMatrix &cont_covariates, const NumericMatrix &class_covariates, 
            const int &cont_covariates_flag, const int &class_covariates_flag
) {

    int i;
    int n_biomarkers = ncol-3;
    vector<double> current_biomarker;
   
    // List of all subgroups found at this level
    vector<SingleSubgroup> all_subgroups(n_biomarkers);

    // List of selected (promising) subgroups
    vector<SingleSubgroup> best_subgroups(width);

    // Structure for sorting subgroups by multiplicity-adjusted splitting criterion on a p-value scale (if no local multiplicity adjustment, unadjusted splitting criterion on a p-value scale is actually used)
    vector<dipair> adjusted_criterion_pvalue_sorted;

    adjusted_criterion_pvalue_sorted.clear();

        // Select a biomarker and find the best patient subgroup for this biomarker
        for(i = 0; i < n_biomarkers; i++) {

            if (biomarker_level[i] <= cur_depth) {

                current_biomarker = vector<double>(df.begin()+i*nrow,df.begin()+(i+1)*nrow);

                vector<int> bio_rows(parent_rows);
                na_index_rows(current_biomarker,bio_rows);

                vector<double> treatment(df.begin()+(ncol - 3)*nrow,df.begin()+(ncol - 2)*nrow);
                vector<double> outcome(df.begin()+(ncol - 2)*nrow,df.begin()+(ncol - 1)*nrow);
                vector<double> outcome_censor(df.begin()+(ncol - 1)*nrow,df.begin()+(ncol)*nrow);

                // Numeric biomarker
                if (biomarker_type[i] == 1) {

                    if (outcome_type == 1 && analysis_method == 1)
                        all_subgroups[i] = ContOutContBio(current_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i]);
                    if (outcome_type == 2 && analysis_method == 1)
                        all_subgroups[i] = BinOutContBio(current_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i]);
                    if (outcome_type == 3 && analysis_method == 1)
                        all_subgroups[i] = SurvOutContBio(current_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i]);    

                    if (outcome_type == 1 && analysis_method >= 2) all_subgroups[i] = ContOutContBioANCOVA(current_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i], model_covariates, analysis_method);

                    if (outcome_type == 2 && analysis_method >= 2) all_subgroups[i] = BinOutContBioANCOVA(current_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i], cont_covariates, class_covariates, cont_covariates_flag, class_covariates_flag);

                    if (outcome_type == 3 && analysis_method >= 2) all_subgroups[i] = SurvOutContBioANCOVA(current_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i], cont_covariates, class_covariates, cont_covariates_flag, class_covariates_flag);

                }
                // Nominal biomarker
                if (biomarker_type[i] == 2){
                    vector<int> int_biomarker(current_biomarker.begin(),current_biomarker.end());

                    if (outcome_type == 1 && analysis_method == 1)
                        all_subgroups[i] = ContOutNomBio(int_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i]);
                    if (outcome_type == 2 && analysis_method == 1)
                        all_subgroups[i] = BinOutNomBio(int_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i]);
                    if (outcome_type == 3 && analysis_method == 1)
                        all_subgroups[i] = SurvOutNomBio(int_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i]);

                    if (outcome_type == 1 && analysis_method >= 2) all_subgroups[i] = ContOutNomBioANCOVA(int_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i], model_covariates, analysis_method);

                    if (outcome_type == 2 && analysis_method >= 2) all_subgroups[i] = BinOutNomBioANCOVA(current_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i], cont_covariates, class_covariates, cont_covariates_flag, class_covariates_flag);

                    if (outcome_type == 3 && analysis_method >= 2) all_subgroups[i] = SurvOutNomBioANCOVA(current_biomarker, treatment, outcome, outcome_censor, bio_rows, nmin, direction, criterion_type, adj_parameter[i], cont_covariates, class_covariates, cont_covariates_flag, class_covariates_flag);

                }

                all_subgroups[i].biomarker_index = i + 1;
                if (cur_depth == depth)
                    all_subgroups[i].terminal_subgroup = 1;

                // Save the value of the splitting criterion on a p-value scale for the current patient subgroup
                if (all_subgroups[i].code == 0) {
                    adjusted_criterion_pvalue_sorted.push_back(dipair(all_subgroups[i].adjusted_criterion_pvalue, i));
                
                }
            }
        }

    // Sort the patient subgroups by multiplicity-adjusted splitting criterion on a p-value scale (from smallest to largest)
    sort(adjusted_criterion_pvalue_sorted.begin(), adjusted_criterion_pvalue_sorted.end(), DIPairSortUp); 

    // Select the subgroups with a high splitting criterion
    best_subgroups.clear();

    int m = width;
    if (width > adjusted_criterion_pvalue_sorted.size()) m = adjusted_criterion_pvalue_sorted.size();

    double log_parent, log_child;

    for(i = 0; i < m; i++) {
        int ind = adjusted_criterion_pvalue_sorted[i].second;

        // Apply complexity control if gamma is not equal to -1 (which means that complexity control is disabled)
        if (gamma[cur_depth-1] > 0) {
            // Compare the treatment effect p-value within the current subgroup to the parent p-value 
            // The complexity criterion is applied on a log scale
            log_parent = Log1minusPSupp(parent_pvalue, parent_test_statistic, 1, 1);
            log_child = Log1minusPSupp(all_subgroups[ind].pvalue, all_subgroups[ind].test_statistic, 1, 1);

            // if(parent_pvalue*gamma[cur_depth-1] >= all_subgroups[ind].pvalue){
            if (log_parent + log(gamma[cur_depth-1]) >= log_child) {
                all_subgroups[ind].signat = 10000000 * all_subgroups[ind].size + round(1000000 * all_subgroups[ind].test_statistic);
                best_subgroups.push_back(all_subgroups[ind]);
            }
        }
        else {        
                all_subgroups[ind].signat = 10000000 * all_subgroups[ind].size + round(1000000 * all_subgroups[ind].test_statistic);
                best_subgroups.push_back(all_subgroups[ind]);          
        }
    }

    if(cur_depth < depth){
        for(vector<SingleSubgroup>::iterator isubs =best_subgroups.begin(); isubs!=best_subgroups.end();++isubs){
            if (isubs->terminal_subgroup !=1 && isubs->size > 2*nmin){
                depth_hist[cur_depth-1]=isubs->biomarker_index;
                std::vector<SingleSubgroup> xxx = FindSubgroups(df, model_covariates, ncol, nrow,isubs->subgroup_rows,nmin,direction,criterion_type,outcome_type, analysis_method,biomarker_type,adj_parameter,width,depth,cur_depth+1,depth_hist,gamma,isubs->pvalue, isubs->test_statistic, cont_covariates, class_covariates, cont_covariates_flag, class_covariates_flag);
                isubs->subgroups = xxx;
            }
            if (isubs->subgroups.empty()){
                isubs->terminal_subgroup = 1;
            }

        }
    }

    return best_subgroups;

}

	
// [[Rcpp::export]]
double SIDES(const NumericVector &ancova_outcome_arg, const NumericVector &ancova_censor_arg, const NumericVector &ancova_treatment_arg, 
            const NumericMatrix &cont_covariates, const NumericMatrix &class_covariates, 
            const int &n_cont_covariates, const int &n_class_covariates,
            const std::string project_filename,
            const std::string output_filename
) {
    
	ofstream out(output_filename);

    out << fixed << showpoint;

    int i, j, k, n_subgroups;
			
	xml_document<> doc;	//create xml_document object (http://rostakagmfun.github.io/xml-file-parsing-cpp.html)
	file<> xmlFile(project_filename.c_str()); //open file
	doc.parse<0>(xmlFile.data()); //parse the contents of file
	xml_node<>* root = doc.first_node("head");//find our root node
	xml_node<>* nodeStructure = root->first_node("structure");
	xml_node<>* nodeParameters = root->first_node("parameters");

	int nrow = atoi(root->first_node("data")->first_attribute("nrow")->value());
	int ncol = atoi(root->first_node("data")->first_attribute("ncol")->value());
	int skip = atoi(root->first_node("data")->first_attribute("skipfirstrow")->value());
	string std_data_filename = root->first_node("data")->first_attribute("stddata")->value();

	// Read in standardized data
	std::ifstream std_data_file(std_data_filename.c_str());
	std::string data_set((std::istreambuf_iterator<char>(std_data_file)),
		                 (std::istreambuf_iterator<char>()));

	//[X]min_subgroup_size (number)
	int min_subgroup_size = atoi(nodeParameters->first_node("min_subgroup_size")->first_attribute("value")->value());
	//[X]outcome_variable_direction (number)
	string outcome_variable_direction_str = nodeStructure->first_node("outcome")->first_attribute("direction")->value();
	int outcome_variable_direction = 1;	// larger
	if (outcome_variable_direction_str == "smaller")
		outcome_variable_direction = -1;
	//[X]outcome_type (number)
	string outcome_type_str = nodeStructure->first_node("outcome")->first_attribute("type")->value();
	int outcome_type = 1; // continuous
	if (outcome_type_str == "binary")
		outcome_type = 2;
	if (outcome_type_str == "time")
		outcome_type = 3;

	//[x]covariate_types (vector)
	xml_node<>* nodeBiomarkers = nodeStructure->first_node("biomarkers");
	vector<int> covariate_types;
	covariate_types.clear();
    biomarker_level.clear();
	if (nodeBiomarkers) {
		for (xml_node<> *child = nodeBiomarkers->first_node(); child; child = child->next_sibling()) {
			covariate_types.push_back(atoi(child->first_attribute("numeric")->value()));
            biomarker_level.push_back(atoi(child->first_attribute("level")->value()));
        }
	}

	//[X]criterion_type (number)
	int criterion_type = atoi(nodeParameters->first_node("criterion_type")->first_attribute("value")->value());
	//[X]width (number)
	int width = atoi(nodeParameters->first_node("width")->first_attribute("value")->value());
	//[X]depth (number)
	int depth = atoi(nodeParameters->first_node("depth")->first_attribute("value")->value());
	//[X]gamma (vector)
	xml_node<>* nodeComplexityControl = nodeParameters->first_node("complexity_control");
	bool complexity_control = (nodeComplexityControl != NULL);
	vector<double> gamma;
	gamma.clear();
	if (complexity_control) {
		for (xml_node<> *child = nodeComplexityControl->first_node(); child; child = child->next_sibling())
			gamma.push_back(atof(child->first_attribute("value")->value()));
	} else {
		for (i = 0; i < depth; i++)
			gamma.push_back(-1.0f);
	}
	//[X]pvalue_max (number)
	// double pvalue_max = atof(nodeParameters->first_node("pvalue_max")->first_attribute("value")->value());
	//[X]local_mult_adj (number)
	double local_mult_adj = atof(nodeParameters->first_node("local_mult_adj")->first_attribute("value")->value());
	//[X]n_perms_mult_adjust (number)
	int n_perms_mult_adjust = atoi(nodeParameters->first_node("n_perms_mult_adjust")->first_attribute("value")->value());
    //[X] perm_type (number)
    int perm_type = atoi(nodeParameters->first_node("perm_type")->first_attribute("value")->value());
	//[X]subgroup_search_algorithm (number)
	int subgroup_search_algorithm = atoi(nodeParameters->first_node("subgroup_search_algorithm")->first_attribute("value")->value());
	//[X]n_top_biomarkers (number)
	int n_top_biomarkers = atoi(nodeParameters->first_node("n_top_biomarkers")->first_attribute("value")->value());
	//[X]multiplier (number)
	double multiplier = atof(nodeParameters->first_node("multiplier")->first_attribute("value")->value());
	//[X]n_perms_vi_score (number)
	int n_perms_vi_score = atoi(nodeParameters->first_node("n_perms_vi_score")->first_attribute("value")->value());
    int analysis_method = atoi(nodeParameters->first_node("analysis_method")->first_attribute("value")->value());
    int nperc = atoi(nodeParameters->first_node("nperc")->first_attribute("value")->value()) - 1;

    // Global parameters
    precision = atof(nodeParameters->first_node("precision")->first_attribute("value")->value()) - 1;
    max_iter = atoi(nodeParameters->first_node("max_iter")->first_attribute("value")->value()) - 1;

    //***************************************************************************************************************

	// Read in the main data set
	vector<string> CsvLines;
	boost::split(CsvLines, data_set, boost::is_any_of("\n"));
	// Remove empty lines
	for (i = CsvLines.size()-1; i >= 0; i--) {
		boost::trim(CsvLines[i]);
		if (CsvLines[i].size() == 0)
			CsvLines.erase(CsvLines.begin() + i);
	}
	// Remove header lines
	for (i = 0; i < skip; i++)
		CsvLines.erase(CsvLines.begin());

	vector<double> treatment;
	vector<double> outcome;
	vector<double> outcome_censor;

    vector<double> data_all_biomarkers;
    vector<double> data_all_biomarkers_permuted;
    vector<double> data_selected_biomarkers;
    vector<double> data_selected_biomarkers_permuted;

    // All covariates placed in the same vector
    // vector<double> ancova_all_covariates;
    // vector<double> ancova_all_covariates_permuted;

    ModelCovariates model_covariates;

	vector<double> yy, yy2, yy3;

    int n_biomarkers = ncol - 3;
    vector<double> adj_parameter(n_biomarkers);
    vector<int> biomarker_type(n_biomarkers);
    vector<double> current_biomarker, transformed_biomarker;
    int n_levels;

    data_all_biomarkers.clear();

	vector<string> row;

    // Type I error rate
    double error_rate = 0.0;

	for (i = 0; i < CsvLines.size(); i++) {
		boost::split(row, CsvLines[i], boost::is_any_of(","), boost::token_compress_on);	
		for (j = 0; j < ncol; j++) {
            if (row[j] != "." && row[j] != " ." && row[j] != "NA") {
                yy2.push_back(atof(row[j].c_str()));
            }
            else {                
                yy2.push_back(numeric_limits<double>::quiet_NaN());
            }
        }
		treatment.push_back( atof(row[ncol - 3].c_str()) );
		outcome.push_back( atof(row[ncol - 2].c_str()) );
		outcome_censor.push_back( atof(row[ncol - 1].c_str()) );
	}

	for (int i1 = 0; i1<ncol; i1++) {
        current_biomarker.clear();
        transformed_biomarker.clear();
        for (int j1 = 0; j1<nrow; j1++) {
            current_biomarker.push_back(yy2[i1 + j1*ncol]);

        }

        // Convert numeric biomarkers to percentiles
        if (nperc > 0 && CountUniqueValues(current_biomarker) > nperc && covariate_types[i1] == 1) {
            transformed_biomarker = QuantileTransform(current_biomarker, nperc); 
            data_all_biomarkers.insert(data_all_biomarkers.end(), transformed_biomarker.begin(), transformed_biomarker.end());
        } else {
            data_all_biomarkers.insert(data_all_biomarkers.end(), current_biomarker.begin(), current_biomarker.end());
        }

    }

    // Define the biomarker type and compute local multiplicity adjustment for each biomarker (if local multiplicity adjustment is enabled)
    for (i = 0; i < n_biomarkers; i++) {

        // Numerical biomarker
        if (covariate_types[i] == 1) {
            biomarker_type[i] = 1;
        }
        // Nominal biomarker
        if (covariate_types[i] == 0) {
            biomarker_type[i] = 2;
        }

        // Local multiplicity adjustment
        current_biomarker = vector<double>(data_all_biomarkers.begin()+i*nrow, data_all_biomarkers.begin()+(i+1)*nrow);
        n_levels = CountUniqueValues(current_biomarker);

        // Compute local multiplicity adjustment parameter
        adj_parameter[i] = (1.0 - local_mult_adj) + local_mult_adj * AdjParameter(biomarker_type[i], n_levels, criterion_type);


    }

    //***************************************************************************************************************

    // Read in the ANCOVA data set

    std::vector<double> ancova_outcome, ancova_treatment, ancova_censor, temp_vec(1);

    NumericVector temp_numeric_vec;

    // Continuous endpoint
    if (outcome_type == 1 && analysis_method >= 2) {

        ancova_outcome = as<vector<double>>(ancova_outcome_arg);
        ancova_treatment = as<vector<double>>(ancova_treatment_arg);

        // analysis_method == 2: cov1
        // analysis_method == 3: cov1, cov2
        // analysis_method == 4: cov1, cov2, factor
        // analysis_method == 5: cov1, cov2, cov3, factor
        // analysis_method == 6: cov1, cov2, cov3, cov4, factor

        if (analysis_method == 2) {
            model_covariates.cov2 = temp_vec;
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            model_covariates.cov_class = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 3) {
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            model_covariates.cov_class = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
        }
        
        if (analysis_method == 4) {
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 5) {
            model_covariates.cov4 = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 2);
            model_covariates.cov3 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 6) {
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 2);
            model_covariates.cov3 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 3);
            model_covariates.cov4 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

    }

    // Binary endpoint
    if (outcome_type == 2 && analysis_method >= 2) {

        ancova_outcome = as<vector<double>>(ancova_outcome_arg);
        ancova_treatment = as<vector<double>>(ancova_treatment_arg);


    }

    // Survival endpoint
    if (outcome_type == 3 && analysis_method >= 2) {

        ancova_outcome = as<vector<double>>(ancova_outcome_arg);
        ancova_censor = as<vector<double>>(ancova_censor_arg);
        ancova_treatment = as<vector<double>>(ancova_treatment_arg);

    }

    //***************************************************************************************************************

	// Start subgroup search

	vector<int> parent_rows(nrow, 1);
	vector<double> best_subgroup_pvalue(n_perms_mult_adjust);

    // List of found subgroups
    vector<SingleSubgroup> single_level, single_level_permuted;
    SingleSubgroup parent_group, parent_group_permuted;

    vector<double> vi_max (n_perms_vi_score);

    // Random seed
    int random_seed = 74228;
    //srand(random_seed);
    set_seed(random_seed);

    // Biomarker indices
    vector<int> biomarker_index;
    biomarker_index.clear();

    // Threshold for variable importance used in the adaptive two-stage procedure
    double vi_threshold = -1.0;
    double mean_vi_max = -1.0;
    double sd_vi_max = -1.0;

    // Second-stage parameters
    int width_second_stage(width);
    int depth_second_stage(depth);
    std::vector<double> gamma_second_stage(gamma);

    // Randomly shuffled indices
    int n_patients = nrow;
    vector<int> shuffle_index(n_patients);
    vector<double> treatment_permuted(n_patients), outcome_permuted(n_patients), outcome_censor_permuted(n_patients);

    // Analysis of overall group of patients

    // Basic test without covariates 
    if (analysis_method == 1) parent_group = OverallAnalysis(treatment, outcome, outcome_censor, outcome_type, outcome_variable_direction); 

    // Continuous outcomes
    if (outcome_type == 1 && analysis_method >= 2) parent_group = OverallAnalysisCont(ancova_treatment, ancova_outcome, model_covariates, analysis_method, outcome_variable_direction);

    // Binary outcomes
    if (outcome_type == 2 && analysis_method >= 2) parent_group = OverallAnalysisBin(ancova_treatment, ancova_outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);

    // Survival outcomes
    if (outcome_type == 3 && analysis_method >= 2) parent_group = OverallAnalysisSurv(ancova_treatment, ancova_outcome, ancova_censor, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);

    // XML file header
    out<<"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    out<<"<head>\n";
    out<<" <subgroups>\n";
    out<<"  <subgroup>\n";
    out<<"   <definition>\n";
    out<<"    <component description=\"Overall population\" biomarker=\"\" sign=\"\" value=\"\"/>\n";
    out<<"   </definition>\n";
    out<<"   <parameters size_control=\""<< parent_group.size_control<<"\" size_treatment=\""<<parent_group.size_treatment <<"\" splitting_criterion=\"\" splitting_criterion_log_p_value=\"\" p_value=\""<< parent_group.pvalue<<"\" prom_estimate=\""<< parent_group.prom_estimate<<"\" prom_sderr=\""<< parent_group.prom_sderr<<"\" prom_sd=\""<< parent_group.prom_sd<<"\" adjusted_p_value=\"-1\"/>\n";
    out<<"  </subgroup>\n";
 
    vector<int> depth_hist(depth, 0);

    // Find subgroups 
    single_level = FindSubgroups(data_all_biomarkers, model_covariates, n_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type, adj_parameter, width, depth, 1, depth_hist, gamma, parent_group.pvalue, parent_group.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);


    // Treatment effect p-values from original subgroups
	vector<double> pvalue;

    // Treatment effect p-values from permuted subgroups
	vector<double> pvalue_permuted;

 	int total_number_subgroups = 0, w = 0;
	TreeSize(single_level, total_number_subgroups, w);
	int iter = 0;
	vector<int> signat1(total_number_subgroups);

    // Compute variable importance for all biomarkers
	int nterm = 0, subnterm = 0;
	vector<double> variable_importance(n_biomarkers, 0);
	signat1.clear();
	ComputeVarImp(single_level, variable_importance, nterm, subnterm, signat1);

    // Stop everything if all VI scores are missing
    int stop_everything = 0;

    for(i=0; i < variable_importance.size(); ++i) {
        variable_importance[i] /= nterm; 
        if (nterm == 0) {
            variable_importance[i] = 0.0;
            stop_everything = 1;
        }
        if (std::isnan(variable_importance[i]) == 1) {
            variable_importance[i] = 0.0;
            stop_everything = 1;
        }
    }

    if (stop_everything == 1) {
        vi_threshold = 0.0;
        mean_vi_max = 0.0;
        sd_vi_max = 0.0;
    }



    string par_info;
    vector<double> adj_pvalue;

    // Regular one-stage SIDES procedure
    if (subgroup_search_algorithm == 1 && stop_everything == 0) {

    	total_number_subgroups = 0; 
    	w = 0;
		TreeSize(single_level, total_number_subgroups, w);
		iter = 0;
		vector<int> signat2(total_number_subgroups);

        // Extract the vector of treatment effect p-values from the set of subgroup
		par_info = "";
		pvalue.clear();
		ExtractPvalues(single_level, par_info, iter, 0, signat2, pvalue);

        // Run permutations to compute multiplicity-adjusted treatment effect p-value for each group (under the null distribution)
  	    n_subgroups = pvalue.size();
        vector<double> adj_pvalue(n_subgroups, 0);

        data_all_biomarkers_permuted = data_all_biomarkers;
        treatment_permuted = treatment;
        outcome_permuted = outcome;


        for (i = 0; i < n_perms_mult_adjust; i++) {

            // Permute the treatment column only
            if (perm_type == 1) {

                shuffle_vector(treatment_permuted);
                for(j = 0; j<n_patients; ++j){
                    data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
                }

                // Analysis of overall group of patients
                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

                // Continuous outcomes
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);


            }

            // Permute the treatment, outcome and outcome censor columns
            if (perm_type == 2) {

                for (j = 0; j < n_patients; j++) shuffle_index[j] = j;
                shuffle_vector(shuffle_index);

                for(j = 0; j<n_patients; j++){

                    treatment_permuted[j] = treatment[shuffle_index[j]];
                    outcome_permuted[j] = outcome[shuffle_index[j]];
                    outcome_censor_permuted[j] = outcome_censor[shuffle_index[j]];
                    data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
                    data_all_biomarkers_permuted[(n_biomarkers + 1) * n_patients + j] = outcome_permuted[j];
                    data_all_biomarkers_permuted[(n_biomarkers + 2) * n_patients + j] = outcome_censor_permuted[j];

                }

                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome_permuted, outcome_censor_permuted, outcome_type, outcome_variable_direction);

                // Continuous outcomes
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome_permuted, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome_permuted, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);
                // parent_group_permuted = OverallAnalysisContOld(treatment_permuted, outcome_permuted, ancova_all_covariates_permuted, analysis_method, outcome_variable_direction); 

            }

            // Permute the outcome column only
            if (perm_type == 3) {

                shuffle_vector(outcome_permuted);
                for(j = 0; j<n_patients; ++j){
                    data_all_biomarkers_permuted[(n_biomarkers + 1) * n_patients + j] = outcome_permuted[j];
                }

                // Analysis of overall group of patients
                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment, outcome_permuted, outcome_censor, outcome_type, outcome_variable_direction);

                // Continuous outcomes
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment, outcome_permuted, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment, outcome_permuted, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

            }

	    	vector<int> depth_hist(depth, 0);        

	        // List of all subgroups (under the null distribution)
	        single_level_permuted = FindSubgroups(data_all_biomarkers_permuted, model_covariates, n_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type, adj_parameter, width, depth, 1, depth_hist, gamma, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

	        // Extract the vector of treatment effect p-values from the set of subgroup
	        pvalue_permuted.clear();
	        par_info = "";
	        total_number_subgroups = 0;
	        w = 0;
			TreeSize(single_level_permuted, total_number_subgroups, w);
			iter = 0;
			vector<int> signat3(total_number_subgroups);

	        ExtractPvalues(single_level_permuted, par_info, iter, 0, signat3, pvalue_permuted);

            if (pvalue_permuted.size() > 0) {

    		    // Find the most significant p-value in the subgroups
                best_subgroup_pvalue[i] = *std::min_element(pvalue_permuted.begin(), pvalue_permuted.end());
                for (j = 0; j < n_subgroups; j++) {
                    if (best_subgroup_pvalue[i] <= pvalue[j]) adj_pvalue[j]++;
                }

                if (best_subgroup_pvalue[i] <= 0.025) error_rate++;

            } 
            else {

                // Increment the adjusted p-value if no subgroups were not found in the current permutation
                for (j = 0; j < n_subgroups; j++) {
                    adj_pvalue[j]++;
                }

                error_rate++;

            }


	    }

	    // Compute multiplicity-adjusted treatment effect p-value
        for (j = 0; j < n_subgroups; j++) {
            adj_pvalue[j] /= n_perms_mult_adjust;
        }

        error_rate /= n_perms_mult_adjust;

        // Save the multiplicity-adjusted treatment effect p-value

        // Indices of biomarkers used in the found subgroups
        for (i = 0; i < n_biomarkers; i++) {
            biomarker_index.push_back(i + 1);
        }

        par_info = "";
        total_number_subgroups = 0;
        w = 0;
        TreeSize(single_level, total_number_subgroups, w);
        iter = 0;
        vector<int> signat4(total_number_subgroups);

        IterateSubgroupSummaryCSV(single_level, out, par_info, iter, 0, signat4, biomarker_index, adj_pvalue);


    }

    // Sort biomarkers based on variable importance 
    vector<dipair> vi_sorted;
    vi_sorted.clear();

    // Fixed two-stage SIDES procedure (SIDEScreen procedure)
    if (subgroup_search_algorithm == 2 && stop_everything == 0) {

        // Prepare biomarkers for sorting
        for (i = 0; i < n_biomarkers; i++) {
            vi_sorted.push_back(dipair(variable_importance[i], i));
        }

        // Sort biomarkers based on variable importance
        sort(vi_sorted.begin(), vi_sorted.end(), DIPairDown);

        // Create a new data set which includes only the top biomarkers, treatment, outcome and outcome_censor columns
        data_selected_biomarkers.clear();

        // Create a vector of biomarker types for the top biomarkers
        vector<int> biomarker_type_second_stage; 

        // Create a vector of adjustment parameters types for the top biomarkers
        vector<double> adj_parameter_second_stage; 

        biomarker_type_second_stage.clear();
        adj_parameter_second_stage.clear();

        biomarker_index.clear();

        // Select the top biomarkers based on variable importance
        for (i = 0; i < n_top_biomarkers; i++) {
            
            j = vi_sorted[i].second;
            current_biomarker = vector<double>(data_all_biomarkers.begin()+j * nrow, data_all_biomarkers.begin() + (j + 1) * nrow);
            biomarker_index.push_back(j + 1);
            biomarker_type_second_stage.push_back(biomarker_type[j]);
            adj_parameter_second_stage.push_back(adj_parameter[j]);
            data_selected_biomarkers.insert(data_selected_biomarkers.end(), current_biomarker.begin(), current_biomarker.end());
        
        }

        // Add treatment column
        data_selected_biomarkers.insert(data_selected_biomarkers.end(), treatment.begin(), treatment.end());

        // Add outcome column
        data_selected_biomarkers.insert(data_selected_biomarkers.end(), outcome.begin(), outcome.end());

        // Add outcome_censor column
        data_selected_biomarkers.insert(data_selected_biomarkers.end(), outcome_censor.begin(), outcome_censor.end());

        vector<int> depth_hist_second_stage(depth_second_stage, 0);

        // Find all subgroups in the second stage of the subgroup search algorithm
        single_level = FindSubgroups(data_selected_biomarkers, model_covariates, n_top_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type_second_stage, adj_parameter_second_stage, width_second_stage, depth_second_stage, 1, depth_hist_second_stage, gamma_second_stage, parent_group.pvalue, parent_group.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);  

   		total_number_subgroups = 0; 
    	w = 0;
		TreeSize(single_level, total_number_subgroups, w);
		iter = 0;
		vector<int> signat5(total_number_subgroups);

        // Extract the vector of treatment effect p-values from the set of subgroup
		par_info = "";
		pvalue.clear();
		ExtractPvalues(single_level, par_info, iter, 0, signat5, pvalue);

        // Run permutations to compute multiplicity-adjusted treatment effect p-value for each group (under the null distribution)
  	    n_subgroups = pvalue.size();
        vector<double> adj_pvalue(n_subgroups, 0);

        data_all_biomarkers_permuted = data_all_biomarkers;
        treatment_permuted = treatment;
        outcome_permuted = outcome;

        for (i = 0; i < n_perms_mult_adjust; i++) {

            // Permute the treatment column only
            if (perm_type == 1) {

                shuffle_vector(treatment_permuted);
                for(j = 0; j<n_patients; ++j){
                    data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
                }

                // Analysis of overall group of patients
                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

                // Continuous outcome
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

            }

            // Permute the treatment, outcome and outcome censor columns
            if (perm_type == 2) {

                for (j = 0; j < n_patients; j++) shuffle_index[j] = j;
                shuffle_vector(shuffle_index);

                for(j = 0; j<n_patients; j++){

                    treatment_permuted[j] = treatment[shuffle_index[j]];
                    outcome_permuted[j] = outcome[shuffle_index[j]];
                    outcome_censor_permuted[j] = outcome_censor[shuffle_index[j]];
                    data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
                    data_all_biomarkers_permuted[(n_biomarkers + 1) * n_patients + j] = outcome_permuted[j];
                    data_all_biomarkers_permuted[(n_biomarkers + 2) * n_patients + j] = outcome_censor_permuted[j];

                }

                // Analysis of overall group of patients
                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome_permuted, outcome_censor_permuted, outcome_type, outcome_variable_direction);

                // Continuous outcome               
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome_permuted, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome_permuted, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

                // parent_group_permuted = OverallAnalysisContOld(treatment_permuted, outcome_permuted, ancova_all_covariates_permuted, analysis_method, outcome_variable_direction); 

            }

            vector<int> depth_hist(depth,0);

            // Find all subgroups
            single_level_permuted = FindSubgroups(data_all_biomarkers_permuted, model_covariates, n_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type, adj_parameter, width, depth, 1, depth_hist, gamma, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

            // Compute variable importance based on the permuted set
            int nterm = 0, subnterm = 0;
            vector<double> imp_permuted(n_biomarkers, 0);
            vector<int> signat;
            ComputeVarImp(single_level_permuted, imp_permuted, nterm, subnterm, signat);
            for(k=0; k < n_biomarkers; ++k)
                imp_permuted[k] /= nterm;

            vi_sorted.clear();

            // Prepare biomarkers for sorting
            for (k = 0; k < n_biomarkers; k++) {
                vi_sorted.push_back(dipair(imp_permuted[k], k));
            }

            // Sort biomarkers based on variable importance
            sort(vi_sorted.begin(), vi_sorted.end(), DIPairDown);

            // Create a new data set which includes only the top biomarkers, treatment, outcome and outcome_censor columns
            data_selected_biomarkers_permuted.clear();

            // Create a vector of biomarker types for the top biomarkers
            vector<int> biomarker_type_second_stage; 

            // Create a vector of adjustment parameters types for the top biomarkers
            vector<double> adj_parameter_second_stage; 

            biomarker_type_second_stage.clear();
            adj_parameter_second_stage.clear();

            // Select the top biomarkers based on variable importance
            for (k = 0; k < n_top_biomarkers; k++) {
                j = vi_sorted[k].second;
                current_biomarker = vector<double>(data_all_biomarkers.begin()+j * nrow, data_all_biomarkers.begin() + (j + 1) * nrow);
                // biomarker_index.push_back(j + 1);
                biomarker_type_second_stage.push_back(biomarker_type[j]);
                adj_parameter_second_stage.push_back(adj_parameter[j]);
                data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), current_biomarker.begin(), current_biomarker.end());
            }

            // Permute the treatment column only
            if (perm_type == 1) {

                // Add treatment column
                data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), treatment_permuted.begin(), treatment_permuted.end());

                // Add outcome column
                data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome.begin(), outcome.end());

                // Add outcome_censor column
                data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome_censor.begin(), outcome_censor.end());

                // Analysis of overall group of patients
                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

                // Continuous outcome
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);
            }

            // Permute the treatment, outcome and outcome censor columns
            if (perm_type == 2) {

                // Add treatment column
                data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), treatment_permuted.begin(), treatment_permuted.end());

                // Add outcome column
                data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome_permuted.begin(), outcome_permuted.end());

                // Add outcome_censor column
                data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome_censor_permuted.begin(), outcome_censor_permuted.end());

                // Analysis of overall group of patients
                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome_permuted, outcome_censor_permuted, outcome_type, outcome_variable_direction);

                // Continuous outcome
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome_permuted, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome_permuted, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

                // parent_group_permuted = OverallAnalysisContOld(treatment_permuted, outcome_permuted, ancova_all_covariates_permuted, analysis_method, outcome_variable_direction); 


            }

            vector<int> depth_hist_second_stage(depth_second_stage, 0);

            // List of all subgroups (under the null distribution)
            single_level_permuted = FindSubgroups(data_selected_biomarkers_permuted, model_covariates, n_top_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type_second_stage, adj_parameter_second_stage, width_second_stage, depth_second_stage, 1, depth_hist_second_stage, gamma_second_stage, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

            // Extract the vector of treatment effect p-values from the set of subgroup
            pvalue_permuted.clear();
            par_info = "";
            total_number_subgroups = 0;
            w = 0;
            TreeSize(single_level_permuted, total_number_subgroups, w);
            iter = 0;
            vector<int> signat13(total_number_subgroups);

            ExtractPvalues(single_level_permuted, par_info, iter, 0, signat13, pvalue_permuted);

            if (pvalue_permuted.size() > 0) {

                // Find the most significant p-value in the subgroups
                best_subgroup_pvalue[i] = *std::min_element(pvalue_permuted.begin(), pvalue_permuted.end());
                for (j = 0; j < n_subgroups; j++) {
                    if (best_subgroup_pvalue[i] <= pvalue[j]) adj_pvalue[j]++;
                }

                if (best_subgroup_pvalue[i] <= 0.025) error_rate++;
            
            }
            else {

                // Increment the adjusted p-value if no subgroups were not found in the current permutation
                for (j = 0; j < n_subgroups; j++) {
                    adj_pvalue[j]++;
                }

                error_rate++;

            }


	    }

	    // Compute multiplicity-adjusted treatment effect p-value
        for (j = 0; j < n_subgroups; j++) {
            adj_pvalue[j] /= n_perms_mult_adjust;
        }
 
        error_rate /= n_perms_mult_adjust;

        par_info = "";
        total_number_subgroups = 0;
        w = 0;
        TreeSize(single_level, total_number_subgroups, w);
        iter = 0;
        vector<int> signat7(total_number_subgroups);

        IterateSubgroupSummaryCSV(single_level, out, par_info, iter, 0, signat7, biomarker_index, adj_pvalue);

    }

    // Adaptive two-stage SIDES procedure (SIDEScreen procedure)
    if (subgroup_search_algorithm == 3 && stop_everything == 0) {

        data_all_biomarkers_permuted = data_all_biomarkers;
        treatment_permuted = treatment;
        outcome_permuted = outcome;

        // Run permutations to compute mean and standard deviation of variable importance under the null distribution
        for (i = 0; i < n_perms_vi_score; i++) {

            // Permute the treatment column only
            if (perm_type == 1) {

                shuffle_vector(treatment_permuted);
                for(j = 0; j<n_patients; ++j){
                    data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
                }

                // Analysis of overall group of patients
                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

                // Continuous outcome
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

            }

            // Permute the treatment, outcome and outcome censor columns
            if (perm_type == 2) {

                for (j = 0; j < n_patients; j++) shuffle_index[j] = j;
                shuffle_vector(shuffle_index);

                for(j = 0; j<n_patients; j++){

                    treatment_permuted[j] = treatment[shuffle_index[j]];
                    outcome_permuted[j] = outcome[shuffle_index[j]];
                    outcome_censor_permuted[j] = outcome_censor[shuffle_index[j]];
                    data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
                    data_all_biomarkers_permuted[(n_biomarkers + 1) * n_patients + j] = outcome_permuted[j];
                    data_all_biomarkers_permuted[(n_biomarkers + 2) * n_patients + j] = outcome_censor_permuted[j];

                }


                // Analysis of overall group of patients
                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome_permuted, outcome_censor_permuted, outcome_type, outcome_variable_direction);

                // Continuous outcome
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome_permuted, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome_permuted, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

                    // parent_group_permuted = OverallAnalysisContOld(treatment_permuted, outcome_permuted, ancova_all_covariates_permuted, analysis_method, outcome_variable_direction);

            }

            // Permute the outcome column only
            if (perm_type == 3) {

                shuffle_vector(outcome_permuted);
                for(j = 0; j<n_patients; ++j){
                    data_all_biomarkers_permuted[(n_biomarkers + 1) * n_patients + j] = outcome_permuted[j];
                }

                // Analysis of overall group of patients
                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment, outcome_permuted, outcome_censor, outcome_type, outcome_variable_direction);

                // Continuous outcome
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment, outcome_permuted, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment, outcome_permuted, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

            }            

            vector<int> depth_hist(depth,0);

            // Find all subgroups
            single_level_permuted = FindSubgroups(data_all_biomarkers_permuted, model_covariates, n_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type, adj_parameter, width, depth, 1, depth_hist, gamma, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

            // Compute variable importance based on the permuted set
            int nterm = 0, subnterm = 0;
            vector<double> imp_permuted(n_biomarkers, 0);
            vector<int> signat8;
            ComputeVarImp(single_level_permuted, imp_permuted, nterm, subnterm, signat8);
            for(k=0; k < n_biomarkers; ++k)
                imp_permuted[k] /= nterm;

            // Find the maximum variable importance 
            vi_max[i] = *std::max_element(imp_permuted.begin(), imp_permuted.end());

            // Set the maximum variable importance to 0 if it is not computed
            if (std::isnan(vi_max[i]) == 1) {
                vi_max[i] = 0.0;
            }

            // Set the maximum variable importance to 0 if it is not computed
            if (std::isinf(vi_max[i]) == 1) {
                vi_max[i] = 0.0;
            }

        }

        // Compute mean and standard deviation of maximum variable importance 
        mean_vi_max = std::accumulate(vi_max.begin(), vi_max.end(), 0.0) / n_perms_vi_score;
        double sq_sum = std::inner_product(vi_max.begin(), vi_max.end(), vi_max.begin(), 0.0);
        sd_vi_max = sq_sum / n_perms_vi_score - mean_vi_max * mean_vi_max;
		if (sd_vi_max < 0.0)
			sd_vi_max = 0.0;
		else
			sd_vi_max = std::sqrt(sd_vi_max);

        // Threshold for variable importance used in the adaptive two-stage procedure
        vi_threshold = mean_vi_max + multiplier * sd_vi_max;

        // Create a new data set which includes only the top biomarkers, treatment, outcome and outcome_censor columns
        data_selected_biomarkers.clear();

        // Number of selected biomarkers with variable importance above the threshold
        int n_selected_biomarkers = 0;

        // Create a vector of biomarker types for the top biomarkers
        vector<int> biomarker_type_second_stage; 

        // Create a vector of adjustment parameters for the top biomarkers
        vector<double> adj_parameter_second_stage; 

        biomarker_type_second_stage.clear();
        adj_parameter_second_stage.clear();

        // Select biomarkers with variable importance above the threshold
        for (i = 0; i < n_biomarkers; i++) {

            if (variable_importance[i] >= vi_threshold) {
                n_selected_biomarkers++;
                current_biomarker = vector<double>(data_all_biomarkers.begin()+i * nrow, data_all_biomarkers.begin() + (i + 1) * nrow); 
                biomarker_index.push_back(i + 1);
                biomarker_type_second_stage.push_back(biomarker_type[i]);
                adj_parameter_second_stage.push_back(adj_parameter[i]);
                data_selected_biomarkers.insert(data_selected_biomarkers.end(), current_biomarker.begin(), current_biomarker.end());               
            }

        }

        // Find subgroups in the second stage  of the subgroup search algorithm (if at least one biomarker is selected)
        if (n_selected_biomarkers > 0) {
         
            // Add treatment column
            data_selected_biomarkers.insert(data_selected_biomarkers.end(), treatment.begin(), treatment.end());

            // Add outcome column
            data_selected_biomarkers.insert(data_selected_biomarkers.end(), outcome.begin(), outcome.end());

            // Add outcome_censor column
            data_selected_biomarkers.insert(data_selected_biomarkers.end(), outcome_censor.begin(), outcome_censor.end());

            // Analysis of overall group of patients
            if (analysis_method == 1) parent_group = OverallAnalysis(treatment, outcome, outcome_censor, outcome_type, outcome_variable_direction);

            // Continuous outcomes
            if (outcome_type == 1 && analysis_method >= 2) parent_group = OverallAnalysisCont(treatment, outcome, model_covariates, analysis_method, outcome_variable_direction);

            // Binary outcomes
            if (outcome_type == 2 && analysis_method >= 2) parent_group = OverallAnalysisBin(treatment, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               


            vector<int> depth_hist_second_stage(depth_second_stage, 0);

            // Find all subgroups in the second stage of the subgroup search algorithm
            single_level = FindSubgroups(data_selected_biomarkers, model_covariates, n_selected_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type_second_stage, adj_parameter_second_stage, width_second_stage, depth_second_stage, 1, depth_hist_second_stage, gamma_second_stage, parent_group.pvalue, parent_group.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

            total_number_subgroups = 0; 
            w = 0;
            TreeSize(single_level, total_number_subgroups, w);
            iter = 0;
            vector<int> signat8(total_number_subgroups);

            // Extract the vector of treatment effect p-values from the set of subgroup
            par_info = "";
            pvalue.clear();
            ExtractPvalues(single_level, par_info, iter, 0, signat8, pvalue);

            // Run permutations to compute multiplicity-adjusted treatment effect p-value for each group (under the null distribution)
            n_subgroups = pvalue.size();
            vector<double> adj_pvalue(n_subgroups, 0);

            data_all_biomarkers_permuted = data_all_biomarkers;
            treatment_permuted = treatment;
            outcome_permuted = outcome;

            for (i = 0; i < n_perms_mult_adjust; i++) {

                // Permute the treatment column only
                if (perm_type == 1) {

                    shuffle_vector(treatment_permuted);
                    for(j = 0; j<n_patients; ++j){
                        data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
                    }

                    // Analysis of overall group of patients
                    if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

                    // Continuous outcome
                    if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

                    // Binary outcomes
                    if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

                }

                // Permute the treatment, outcome and outcome censor columns
                if (perm_type == 2) {

                    for (j = 0; j < n_patients; j++) shuffle_index[j] = j;
                    shuffle_vector(shuffle_index);

                    for(j = 0; j<n_patients; j++){

                        treatment_permuted[j] = treatment[shuffle_index[j]];
                        outcome_permuted[j] = outcome[shuffle_index[j]];
                        outcome_censor_permuted[j] = outcome_censor[shuffle_index[j]];
                        data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
                        data_all_biomarkers_permuted[(n_biomarkers + 1) * n_patients + j] = outcome_permuted[j];
                        data_all_biomarkers_permuted[(n_biomarkers + 2) * n_patients + j] = outcome_censor_permuted[j];

                    }

                    // Analysis of overall group of patients
                    if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome_permuted, outcome_censor_permuted, outcome_type, outcome_variable_direction);

                    // Continuous outcome
                    if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome_permuted, model_covariates, analysis_method, outcome_variable_direction);

                    // Binary outcomes
                    if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome_permuted, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

                    // parent_group_permuted = OverallAnalysisContOld(treatment_permuted, outcome_permuted, ancova_all_covariates_permuted, analysis_method, outcome_variable_direction);

                }

                // Permute the outcome column only
                if (perm_type == 3) {

                    shuffle_vector(outcome_permuted);
                    for(j = 0; j<n_patients; ++j){
                        data_all_biomarkers_permuted[(n_biomarkers + 1) * n_patients + j] = outcome_permuted[j];
                    }

                    // Analysis of overall group of patients
                    if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment, outcome_permuted, outcome_censor, outcome_type, outcome_variable_direction);

                    // Continuous outcome
                    if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment, outcome_permuted, model_covariates, analysis_method, outcome_variable_direction);

                    // Binary outcomes
                    if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment, outcome_permuted, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

                }

                vector<int> depth_hist(depth,0);

                // List of all subgroups (under the null distribution)
                single_level_permuted = FindSubgroups(data_all_biomarkers_permuted, model_covariates, n_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type, adj_parameter, width, depth, 1, depth_hist, gamma, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

                // Compute variable importance based on the permuted set
                int nterm = 0, subnterm = 0;
                vector<double> imp_permuted(n_biomarkers, 0);
                vector<int> signat18;
                ComputeVarImp(single_level_permuted, imp_permuted, nterm, subnterm, signat18);
                for(k=0; k < n_biomarkers; ++k)
                    imp_permuted[k] /= nterm;

                // Create a vector of biomarker types for the biomarkers above the threshold
                vector<int> biomarker_type_second_stage; 

                // Create a vector of adjustment parameters types for the biomarkers above the threshold
                vector<double> adj_parameter_second_stage; 

                data_selected_biomarkers_permuted.clear();

                biomarker_type_second_stage.clear();
                adj_parameter_second_stage.clear();

                int n_selected_biomarkers_permuted = 0;

                // Select biomarkers with variable importance above the threshold
                for (k = 0; k < n_biomarkers; k++) {

                    if (imp_permuted[k] >= vi_threshold) {
                        n_selected_biomarkers_permuted++; 
                        current_biomarker = vector<double>(data_all_biomarkers.begin()+k * nrow, data_all_biomarkers.begin() + (k + 1) * nrow);
                        // index.push_back(k + 1);
                        biomarker_type_second_stage.push_back(biomarker_type[k]);                
                        adj_parameter_second_stage.push_back(adj_parameter[k]);  
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), current_biomarker.begin(), current_biomarker.end());               
                    }
                }

                if (n_selected_biomarkers_permuted > 0) {

                // Permute the treatment column only
                if (perm_type == 1) {

                        // Add treatment column
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), treatment_permuted.begin(), treatment_permuted.end());

                        // Add outcome column
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome.begin(), outcome.end());

                        // Add outcome_censor column
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome_censor.begin(), outcome_censor.end());

                        // Analysis of overall group of patients
                        if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

                        // Continuous outcome
                        if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

                        // Binary outcomes
                        if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               


                    }

                // Permute the treatment, outcome and outcome censor columns
                if (perm_type == 2) {

                        // Add treatment column
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), treatment_permuted.begin(), treatment_permuted.end());

                        // Add outcome column
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome_permuted.begin(), outcome_permuted.end());

                        // Add outcome_censor column
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome_censor_permuted.begin(), outcome_censor_permuted.end());

                        // Analysis of overall group of patients
                        if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome_permuted, outcome_censor_permuted, outcome_type, outcome_variable_direction);

                        // Continuous outcome
                        if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome_permuted, model_covariates, analysis_method, outcome_variable_direction);

                        // Binary outcomes
                        if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome_permuted, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates); 

                        // parent_group_permuted = OverallAnalysisContOld(treatment_permuted, outcome_permuted, ancova_all_covariates_permuted, analysis_method, outcome_variable_direction);
                    }

                    // Permute the outcome column only
                    if (perm_type == 3) {

                        // Add treatment column
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), treatment.begin(), treatment.end());

                        // Add outcome column
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome_permuted.begin(), outcome_permuted.end());

                        // Add outcome_censor column
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome_censor.begin(), outcome_censor.end());

                        // Analysis of overall group of patients
                        if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment, outcome_permuted, outcome_censor, outcome_type, outcome_variable_direction);

                        // Continuous outcome
                        if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment, outcome_permuted, model_covariates, analysis_method, outcome_variable_direction);

                        // Binary outcomes
                        if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment, outcome_permuted, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates); 

                    }

                    vector<int> depth_hist_second_stage(depth_second_stage, 0);

                    // Find all subgroups in the second stage of the subgroup search algorithm
                    single_level_permuted = FindSubgroups(data_selected_biomarkers_permuted, model_covariates, n_selected_biomarkers_permuted + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type_second_stage, adj_parameter_second_stage, width_second_stage, depth_second_stage, 1, depth_hist_second_stage, gamma_second_stage, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

                    // Extract the vector of treatment effect p-values from the set of subgroup
                    pvalue_permuted.clear();
                    par_info = "";
                    total_number_subgroups = 0;
                    w = 0;
                    TreeSize(single_level_permuted, total_number_subgroups, w);
                    iter = 0;
                    vector<int> signat14(total_number_subgroups);

                    ExtractPvalues(single_level_permuted, par_info, iter, 0, signat14, pvalue_permuted);

                    if (pvalue_permuted.size() > 0) {

                        // Find the most significant p-value in the subgroups
                        best_subgroup_pvalue[i] = *std::min_element(pvalue_permuted.begin(), pvalue_permuted.end());
                        for (j = 0; j < n_subgroups; j++) {
                            if (best_subgroup_pvalue[i] <= pvalue[j]) adj_pvalue[j]++;
                        }

                        if (best_subgroup_pvalue[i] <= 0.025) error_rate++;

                    }
                    else {

                    // Increment the adjusted p-value if no subgroups were not found in the current permutation
                    for (j = 0; j < n_subgroups; j++) {
                        adj_pvalue[j]++;
                    }

                    error_rate++;

                    }


                }

            }

            // Compute multiplicity-adjusted treatment effect p-value
            for (j = 0; j < n_subgroups; j++) {
                adj_pvalue[j] /= n_perms_mult_adjust;
            }

            error_rate /= n_perms_mult_adjust;
     
            par_info = "";
            total_number_subgroups = 0;
            w = 0;
            TreeSize(single_level, total_number_subgroups, w);
            iter = 0;
            vector<int> signat10(total_number_subgroups);

            IterateSubgroupSummaryCSV(single_level, out, par_info, iter, 0, signat10, biomarker_index, adj_pvalue);

        }

    }

    out<<" </subgroups>\n";
    out<<" <vi_scores>\n";
	for(int i=0; i < n_biomarkers; ++i) {
		out <<"  <vi_score biomarker=\""<<i + 1<<"\" value=\""<<variable_importance[i]<< "\"/>\n";
	}
    out<<" </vi_scores>\n";

    // Save the VI threshold if adaptive two-stage SIDES procedure (SIDEScreen procedure)
    if (subgroup_search_algorithm == 3) {

        out<<" <vi_threshold mean=\""<<mean_vi_max<<"\" sd=\""<<sd_vi_max<<"\" threshold=\""<<vi_threshold<<"\" />\n";

    }

    // Type I error rate 
    out<<" <error_rate value=\""<<error_rate<<"\" />\n"; 
    out<<"</head>\n";
	out.close();

    return 37.0;
}

// [[Rcpp::export]]
NumericVector SIDESAdjP(
    const NumericVector &ancova_outcome_arg, const NumericVector &ancova_censor_arg, 
    const NumericVector &ancova_treatment_arg, const NumericMatrix &cont_covariates, 
    const NumericMatrix &class_covariates, 
    const int &n_cont_covariates, const int &n_class_covariates, const int &random_seed,
    const std::string project_filename
) {

    int i, j;
            
    xml_document<> doc; //create xml_document object (http://rostakagmfun.github.io/xml-file-parsing-cpp.html)
    file<> xmlFile(project_filename.c_str()); //open file
    doc.parse<0>(xmlFile.data()); //parse the contents of file
    xml_node<>* root = doc.first_node("head");//find our root node
    xml_node<>* nodeStructure = root->first_node("structure");
    xml_node<>* nodeParameters = root->first_node("parameters");

    int nrow = atoi(root->first_node("data")->first_attribute("nrow")->value());
    int ncol = atoi(root->first_node("data")->first_attribute("ncol")->value());
    int skip = atoi(root->first_node("data")->first_attribute("skipfirstrow")->value());
    string std_data_filename = root->first_node("data")->first_attribute("stddata")->value();

    // Read in standardized data
    std::ifstream std_data_file(std_data_filename.c_str());
    std::string data_set((std::istreambuf_iterator<char>(std_data_file)),
                         (std::istreambuf_iterator<char>()));

    //[X]min_subgroup_size (number)
    int min_subgroup_size = atoi(nodeParameters->first_node("min_subgroup_size")->first_attribute("value")->value());
    //[X]outcome_variable_direction (number)
    string outcome_variable_direction_str = nodeStructure->first_node("outcome")->first_attribute("direction")->value();
    int outcome_variable_direction = 1; // larger
    if (outcome_variable_direction_str == "smaller")
        outcome_variable_direction = -1;
    //[X]outcome_type (number)
    string outcome_type_str = nodeStructure->first_node("outcome")->first_attribute("type")->value();
    int outcome_type = 1; // continuous
    if (outcome_type_str == "binary")
        outcome_type = 2;
    if (outcome_type_str == "time")
        outcome_type = 3;
    // if (outcome_type_str == "ancova")
    //     outcome_type = 4;

    //[x]covariate_types (vector)
    xml_node<>* nodeBiomarkers = nodeStructure->first_node("biomarkers");
    vector<int> covariate_types;
    covariate_types.clear();
    biomarker_level.clear();
    if (nodeBiomarkers) {
        for (xml_node<> *child = nodeBiomarkers->first_node(); child; child = child->next_sibling()) {
            covariate_types.push_back(atoi(child->first_attribute("numeric")->value()));
            biomarker_level.push_back(atoi(child->first_attribute("level")->value()));
        }
    }

    //[X]criterion_type (number)
    int criterion_type = atoi(nodeParameters->first_node("criterion_type")->first_attribute("value")->value());
    //[X]width (number)
    int width = atoi(nodeParameters->first_node("width")->first_attribute("value")->value());
    //[X]depth (number)
    int depth = atoi(nodeParameters->first_node("depth")->first_attribute("value")->value());
    //[X]gamma (vector)
    xml_node<>* nodeComplexityControl = nodeParameters->first_node("complexity_control");
    bool complexity_control = (nodeComplexityControl != NULL);
    vector<double> gamma;
    gamma.clear();
    if (complexity_control) {
        for (xml_node<> *child = nodeComplexityControl->first_node(); child; child = child->next_sibling())
            gamma.push_back(atof(child->first_attribute("value")->value()));
    } else {
        for (i = 0; i < depth; i++)
            gamma.push_back(-1.0f);
    }
    //[X]pvalue_max (number)
    // double pvalue_max = atof(nodeParameters->first_node("pvalue_max")->first_attribute("value")->value());
    //[X]local_mult_adj (number)
    double local_mult_adj = atof(nodeParameters->first_node("local_mult_adj")->first_attribute("value")->value());
    //[X]n_perms_mult_adjust (number)
    int n_perms_mult_adjust = atoi(nodeParameters->first_node("n_perms_mult_adjust")->first_attribute("value")->value());
    int analysis_method = atoi(nodeParameters->first_node("analysis_method")->first_attribute("value")->value());
    int nperc = atoi(nodeParameters->first_node("nperc")->first_attribute("value")->value()) - 1;

    //***************************************************************************************************************

    // Read in the main data set
    vector<string> CsvLines;
    boost::split(CsvLines, data_set, boost::is_any_of("\n"));
    // Remove empty lines
    for (i = CsvLines.size()-1; i >= 0; i--) {
        boost::trim(CsvLines[i]);
        if (CsvLines[i].size() == 0)
            CsvLines.erase(CsvLines.begin() + i);
    }
    // Remove header lines
    for (i = 0; i < skip; i++)
        CsvLines.erase(CsvLines.begin());


    vector<double> treatment;
    vector<double> outcome;
    vector<double> outcome_censor;

    vector<double> data_all_biomarkers;
    vector<double> data_all_biomarkers_permuted;
    vector<double> data_selected_biomarkers;
    vector<double> data_selected_biomarkers_permuted;

    // All covariates placed in the same vector
    // vector<double> ancova_all_covariates;
    // vector<double> ancova_all_covariates_permuted;

    ModelCovariates model_covariates;

    vector<double> yy, yy2, yy3;

    int n_biomarkers = ncol - 3;
    vector<double> adj_parameter(n_biomarkers);
    vector<int> biomarker_type(n_biomarkers);
    vector<double> current_biomarker, transformed_biomarker;
    int n_levels;

    data_all_biomarkers.clear();

    vector<string> row;

    for (i = 0; i < CsvLines.size(); i++) {
        boost::split(row, CsvLines[i], boost::is_any_of(","), boost::token_compress_on);    
        for (j = 0; j < ncol; j++) {
            if (row[j] != "." && row[j] != " .") {
                yy2.push_back(atof(row[j].c_str()));
            }
            else {                
                yy2.push_back(numeric_limits<double>::quiet_NaN());
            }
        }
        treatment.push_back( atof(row[ncol - 3].c_str()) );
        outcome.push_back( atof(row[ncol - 2].c_str()) );
        outcome_censor.push_back( atof(row[ncol - 1].c_str()) );
    }

    for (int i1 = 0; i1<ncol; i1++) {
        current_biomarker.clear();
        transformed_biomarker.clear();
        for (int j1 = 0; j1<nrow; j1++) {
            current_biomarker.push_back(yy2[i1 + j1*ncol]);

        }

        // Convert numeric biomarkers to percentiles
        if (nperc > 0 && CountUniqueValues(current_biomarker) > nperc && covariate_types[i1] == 1) {
            transformed_biomarker = QuantileTransform(current_biomarker, nperc); 
            data_all_biomarkers.insert(data_all_biomarkers.end(), transformed_biomarker.begin(), transformed_biomarker.end());
        } else {
            data_all_biomarkers.insert(data_all_biomarkers.end(), current_biomarker.begin(), current_biomarker.end());
        }

    }

    // Define the biomarker type and compute local multiplicity adjustment for each biomarker (if local multiplicity adjustment is enabled)
    for (i = 0; i < n_biomarkers; i++) {

        // Numerical biomarker
        if (covariate_types[i] == 1) {
            biomarker_type[i] = 1;
        }
        // Nominal biomarker
        if (covariate_types[i] == 0) {
            biomarker_type[i] = 2;
        }

        // Local multiplicity adjustment
        current_biomarker = vector<double>(data_all_biomarkers.begin()+i*nrow, data_all_biomarkers.begin()+(i+1)*nrow);
        n_levels = CountUniqueValues(current_biomarker);

        // Compute local multiplicity adjustment parameter
        adj_parameter[i] = (1.0 - local_mult_adj) + local_mult_adj * AdjParameter(biomarker_type[i], n_levels, criterion_type);

    }

    //***************************************************************************************************************

    // Read in the ANCOVA data set

    std::vector<double> ancova_outcome, ancova_treatment, temp_vec(1);

    NumericVector temp_numeric_vec;

    // Continuous endpoint
    if (outcome_type == 1 && analysis_method >= 2) {

        ancova_outcome = as<vector<double>>(ancova_outcome_arg);
        ancova_treatment = as<vector<double>>(ancova_treatment_arg);

        // analysis_method == 2: cov1
        // analysis_method == 3: cov1, cov2
        // analysis_method == 4: cov1, cov2, factor
        // analysis_method == 5: cov1, cov2, cov3, factor
        // analysis_method == 6: cov1, cov2, cov3, cov4, factor

        if (analysis_method == 2) {
            model_covariates.cov2 = temp_vec;
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            model_covariates.cov_class = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 3) {
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            model_covariates.cov_class = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
        }
        
        if (analysis_method == 4) {
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 5) {
            model_covariates.cov4 = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 2);
            model_covariates.cov3 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 6) {
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 2);
            model_covariates.cov3 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 3);
            model_covariates.cov4 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }


    }

    // Binary endpoint
    if (outcome_type == 2 && analysis_method >= 2) {

        ancova_outcome = as<vector<double>>(ancova_outcome_arg);
        ancova_treatment = as<vector<double>>(ancova_treatment_arg);

    }

    // Survival endpoint
    if (outcome_type == 3 && analysis_method >= 2) {

        ancova_outcome = as<vector<double>>(ancova_outcome_arg);
        ancova_treatment = as<vector<double>>(ancova_treatment_arg);

    }

    //***************************************************************************************************************   

    vector<int> parent_rows(nrow, 1);
    NumericVector best_subgroup_pvalue(n_perms_mult_adjust);

    // List of found subgroups
    vector<SingleSubgroup> single_level, single_level_permuted;
    SingleSubgroup parent_group, parent_group_permuted;

    int total_number_subgroups, w, iter;

    // Random seed
    // int random_seed = 74228;
    //srand(random_seed);
    set_seed(random_seed);

    // Randomly shuffled indices
    int n_patients = nrow;
    vector<int> shuffle_index(n_patients);
    vector<double> treatment_permuted(n_patients), outcome_permuted(n_patients), outcome_censor_permuted(n_patients);

    string par_info;
    vector<double> adj_pvalue, pvalue_permuted;

    data_all_biomarkers_permuted = data_all_biomarkers;
    treatment_permuted = treatment;
    outcome_permuted = outcome;

    for (i = 0; i < n_perms_mult_adjust; i++) {

        shuffle_vector(treatment_permuted);

        for(j = 0; j<n_patients; ++j){
            data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
        }

        // Analysis of overall group of patients
        if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

        // Continuous outcome
        if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

        // Binary outcomes
        if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates); 

        vector<int> depth_hist(depth, 0);        

        // List of all subgroups (under the null distribution)
        single_level_permuted = FindSubgroups(data_all_biomarkers_permuted, model_covariates, n_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type, adj_parameter, width, depth, 1, depth_hist, gamma, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

        // Extract the vector of treatment effect p-values from the set of subgroup
        pvalue_permuted.clear();
        par_info = "";
        total_number_subgroups = 0;
        w = 0;
        TreeSize(single_level_permuted, total_number_subgroups, w);
        iter = 0;
        vector<int> signat(total_number_subgroups);

        ExtractPvalues(single_level_permuted, par_info, iter, 0, signat, pvalue_permuted);

        if (pvalue_permuted.size() > 0) {

            // Find the most significant p-value in the subgroups
            best_subgroup_pvalue[i] = *std::min_element(pvalue_permuted.begin(), pvalue_permuted.end());

        } else {

            best_subgroup_pvalue[i] = 0.0;

        }


    }

   // return result
   return(best_subgroup_pvalue);
}

// [[Rcpp::export]]
NumericVector FixedSIDEScreenAdjP(
    const NumericVector &ancova_outcome_arg, const NumericVector &ancova_censor_arg, 
    const NumericVector &ancova_treatment_arg, const NumericMatrix &cont_covariates, 
    const NumericMatrix &class_covariates, 
    const int &n_cont_covariates, const int &n_class_covariates, const int &random_seed,
    const std::string project_file
) {

    int i, j, k;
            
    xml_document<> doc; //create xml_document object (http://rostakagmfun.github.io/xml-file-parsing-cpp.html)
    file<> xmlFile(project_file.c_str()); //open file
    doc.parse<0>(xmlFile.data()); //parse the contents of file
    xml_node<>* root = doc.first_node("head");//find our root node
    xml_node<>* nodeStructure = root->first_node("structure");
    xml_node<>* nodeParameters = root->first_node("parameters");

    //string data_set = root->first_node("data")->value();
    int nrow = atoi(root->first_node("data")->first_attribute("nrow")->value());
    int ncol = atoi(root->first_node("data")->first_attribute("ncol")->value());
    int skip = atoi(root->first_node("data")->first_attribute("skipfirstrow")->value());
    string std_data_filename = root->first_node("data")->first_attribute("stddata")->value();

    // Read in standardized data
    std::ifstream std_data_file(std_data_filename.c_str());
    std::string data_set((std::istreambuf_iterator<char>(std_data_file)),
                         (std::istreambuf_iterator<char>()));

    //[X]min_subgroup_size (number)
    int min_subgroup_size = atoi(nodeParameters->first_node("min_subgroup_size")->first_attribute("value")->value());
    //[X]outcome_variable_direction (number)
    string outcome_variable_direction_str = nodeStructure->first_node("outcome")->first_attribute("direction")->value();
    int outcome_variable_direction = 1; // larger
    if (outcome_variable_direction_str == "smaller")
        outcome_variable_direction = -1;
    //[X]outcome_type (number)
    string outcome_type_str = nodeStructure->first_node("outcome")->first_attribute("type")->value();
    int outcome_type = 1; // continuous
    if (outcome_type_str == "binary")
        outcome_type = 2;
    if (outcome_type_str == "time")
        outcome_type = 3;
    // if (outcome_type_str == "ancova")
    //     outcome_type = 4;

    //[x]covariate_types (vector)
    xml_node<>* nodeBiomarkers = nodeStructure->first_node("biomarkers");
    vector<int> covariate_types;
    covariate_types.clear();
    biomarker_level.clear();
    if (nodeBiomarkers) {
        for (xml_node<> *child = nodeBiomarkers->first_node(); child; child = child->next_sibling()) {
            covariate_types.push_back(atoi(child->first_attribute("numeric")->value()));
            biomarker_level.push_back(atoi(child->first_attribute("level")->value()));
        }
    }

    //[X]criterion_type (number)
    int criterion_type = atoi(nodeParameters->first_node("criterion_type")->first_attribute("value")->value());
    //[X]width (number)
    int width = atoi(nodeParameters->first_node("width")->first_attribute("value")->value());
    //[X]depth (number)
    int depth = atoi(nodeParameters->first_node("depth")->first_attribute("value")->value());
    //[X]gamma (vector)
    xml_node<>* nodeComplexityControl = nodeParameters->first_node("complexity_control");
    bool complexity_control = (nodeComplexityControl != NULL);
    vector<double> gamma;
    gamma.clear();
    if (complexity_control) {
        for (xml_node<> *child = nodeComplexityControl->first_node(); child; child = child->next_sibling())
            gamma.push_back(atof(child->first_attribute("value")->value()));
    } else {
        for (i = 0; i < depth; i++)
            gamma.push_back(-1.0f);
    }
    //[X]pvalue_max (number)
    // double pvalue_max = atof(nodeParameters->first_node("pvalue_max")->first_attribute("value")->value());
    //[X]local_mult_adj (number)
    double local_mult_adj = atof(nodeParameters->first_node("local_mult_adj")->first_attribute("value")->value());
    //[X]n_perms_mult_adjust (number)
    int n_perms_mult_adjust = atoi(nodeParameters->first_node("n_perms_mult_adjust")->first_attribute("value")->value());
    //[X]n_top_biomarkers (number)
    int n_top_biomarkers = atoi(nodeParameters->first_node("n_top_biomarkers")->first_attribute("value")->value());
    //[X]n_perms_vi_score (number)
    int n_perms_vi_score = atoi(nodeParameters->first_node("n_perms_vi_score")->first_attribute("value")->value());
    int analysis_method = atoi(nodeParameters->first_node("analysis_method")->first_attribute("value")->value());
    int nperc = atoi(nodeParameters->first_node("nperc")->first_attribute("value")->value()) - 1;

    //***************************************************************************************************************

    // Read in the main data set
    vector<string> CsvLines;
    boost::split(CsvLines, data_set, boost::is_any_of("\n"));
    // Remove empty lines
    for (i = CsvLines.size()-1; i >= 0; i--) {
        boost::trim(CsvLines[i]);
        if (CsvLines[i].size() == 0)
            CsvLines.erase(CsvLines.begin() + i);
    }
    // Remove header lines
    for (i = 0; i < skip; i++)
        CsvLines.erase(CsvLines.begin());


    vector<double> treatment;
    vector<double> outcome;
    vector<double> outcome_censor;

    vector<double> data_all_biomarkers;
    vector<double> data_all_biomarkers_permuted;
    vector<double> data_selected_biomarkers;
    vector<double> data_selected_biomarkers_permuted;

    ModelCovariates model_covariates;

    vector<double> yy, yy2, yy3;

    int n_biomarkers = ncol - 3;
    vector<double> adj_parameter(n_biomarkers);
    vector<int> biomarker_type(n_biomarkers);
    vector<double> current_biomarker, transformed_biomarker;
    int n_levels;

    data_all_biomarkers.clear();

    vector<string> row;

    for (i = 0; i < CsvLines.size(); i++) {
        boost::split(row, CsvLines[i], boost::is_any_of(","), boost::token_compress_on);    
        for (j = 0; j < ncol; j++) {
            if (row[j] != "." && row[j] != " .") {
                yy2.push_back(atof(row[j].c_str()));
            }
            else {                
                yy2.push_back(numeric_limits<double>::quiet_NaN());
            }
        }
        treatment.push_back( atof(row[ncol - 3].c_str()) );
        outcome.push_back( atof(row[ncol - 2].c_str()) );
        outcome_censor.push_back( atof(row[ncol - 1].c_str()) );
    }

    for (int i1 = 0; i1<ncol; i1++) {
        current_biomarker.clear();
        transformed_biomarker.clear();
        for (int j1 = 0; j1<nrow; j1++) {
            current_biomarker.push_back(yy2[i1 + j1*ncol]);

        }

        // Convert numeric biomarkers to percentiles
        if (nperc > 0 && CountUniqueValues(current_biomarker) > nperc && covariate_types[i1] == 1) {
            transformed_biomarker = QuantileTransform(current_biomarker, nperc); 
            data_all_biomarkers.insert(data_all_biomarkers.end(), transformed_biomarker.begin(), transformed_biomarker.end());
        } else {
            data_all_biomarkers.insert(data_all_biomarkers.end(), current_biomarker.begin(), current_biomarker.end());
        }

    }

    // Define the biomarker type and compute local multiplicity adjustment for each biomarker (if local multiplicity adjustment is enabled)
    for (i = 0; i < n_biomarkers; i++) {

        // Numerical biomarker
        if (covariate_types[i] == 1) {
            biomarker_type[i] = 1;
        }
        // Nominal biomarker
        if (covariate_types[i] == 0) {
            biomarker_type[i] = 2;
        }

        // Local multiplicity adjustment
        current_biomarker = vector<double>(data_all_biomarkers.begin()+i*nrow, data_all_biomarkers.begin()+(i+1)*nrow);
        n_levels = CountUniqueValues(current_biomarker);

        // Compute local multiplicity adjustment parameter
        adj_parameter[i] = (1.0 - local_mult_adj) + local_mult_adj * AdjParameter(biomarker_type[i], n_levels, criterion_type);

    }

    //***************************************************************************************************************

    // Read in the ANCOVA data set

    std::vector<double> ancova_outcome, ancova_treatment, temp_vec(1);

    NumericVector temp_numeric_vec;

    // Continuous endpoint
    if (outcome_type == 1 && analysis_method >= 2) {

        ancova_outcome = as<vector<double>>(ancova_outcome_arg);
        ancova_treatment = as<vector<double>>(ancova_treatment_arg);

        if (analysis_method == 2) {
            model_covariates.cov2 = temp_vec;
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            model_covariates.cov_class = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 3) {
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            model_covariates.cov_class = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
        }
        
        if (analysis_method == 4) {
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 5) {
            model_covariates.cov4 = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 2);
            model_covariates.cov3 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 6) {
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 2);
            model_covariates.cov3 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 3);
            model_covariates.cov4 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

    }

    // Binary endpoint
    if (outcome_type == 2 && analysis_method >= 2) {

    }

    // Survival endpoint
    if (outcome_type == 3 && analysis_method >= 2) {

    }

    //***************************************************************************************************************   

    vector<int> parent_rows(nrow, 1);
    NumericVector best_subgroup_pvalue(n_perms_mult_adjust);

    // List of found subgroups
    vector<SingleSubgroup> single_level, single_level_permuted;
    SingleSubgroup parent_group, parent_group_permuted;

    int total_number_subgroups, w, iter;

    // Random seed
    set_seed(random_seed);

    // Treatment effect p-values from original subgroups
    vector<double> pvalue;

    vector<double> variable_importance(n_biomarkers, 0);

    // Biomarker indices
    vector<int> biomarker_index;
    biomarker_index.clear();

    // Second-stage parameters
    int width_second_stage(width);
    int depth_second_stage(depth);
    std::vector<double> gamma_second_stage(gamma);

    vector<double> vi_max (n_perms_vi_score);

    // Randomly shuffled indices
    int n_patients = nrow;
    vector<int> shuffle_index(n_patients);
    vector<double> treatment_permuted(n_patients), outcome_permuted(n_patients), outcome_censor_permuted(n_patients);

    string par_info;
    vector<double> pvalue_permuted;

    int n_subgroups;

    data_all_biomarkers_permuted = data_all_biomarkers;
    treatment_permuted = treatment;
    outcome_permuted = outcome;

    // Sort biomarkers based on variable importance 
    vector<dipair> vi_sorted;
    vi_sorted.clear();

    // Prepare biomarkers for sorting
    for (i = 0; i < n_biomarkers; i++) {
        vi_sorted.push_back(dipair(variable_importance[i], i));
    }

    // Sort biomarkers based on variable importance
    sort(vi_sorted.begin(), vi_sorted.end(), DIPairDown);

    // Create a new data set which includes only the top biomarkers, treatment, outcome and outcome_censor columns
    data_selected_biomarkers.clear();

    // Create a vector of biomarker types for the top biomarkers
    vector<int> biomarker_type_second_stage; 

    // Create a vector of adjustment parameters types for the top biomarkers
    vector<double> adj_parameter_second_stage; 

    biomarker_type_second_stage.clear();
    adj_parameter_second_stage.clear();

    biomarker_index.clear();

    // Select the top biomarkers based on variable importance
    for (i = 0; i < n_top_biomarkers; i++) {
        
        j = vi_sorted[i].second;
        current_biomarker = vector<double>(data_all_biomarkers.begin()+j * nrow, data_all_biomarkers.begin() + (j + 1) * nrow);
        biomarker_index.push_back(j + 1);
        biomarker_type_second_stage.push_back(biomarker_type[j]);
        adj_parameter_second_stage.push_back(adj_parameter[j]);
        data_selected_biomarkers.insert(data_selected_biomarkers.end(), current_biomarker.begin(), current_biomarker.end());
    
    }

    // Add treatment column
    data_selected_biomarkers.insert(data_selected_biomarkers.end(), treatment.begin(), treatment.end());

    // Add outcome column
    data_selected_biomarkers.insert(data_selected_biomarkers.end(), outcome.begin(), outcome.end());

    // Add outcome_censor column
    data_selected_biomarkers.insert(data_selected_biomarkers.end(), outcome_censor.begin(), outcome_censor.end());

    vector<int> depth_hist_second_stage(depth_second_stage, 0);

    // Find all subgroups in the second stage of the subgroup search algorithm
    single_level = FindSubgroups(data_selected_biomarkers, model_covariates, n_top_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type_second_stage, adj_parameter_second_stage, width_second_stage, depth_second_stage, 1, depth_hist_second_stage, gamma_second_stage, parent_group.pvalue, parent_group.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);  

    total_number_subgroups = 0; 
    w = 0;
    TreeSize(single_level, total_number_subgroups, w);
    iter = 0;
    vector<int> signat5(total_number_subgroups);

    // Extract the vector of treatment effect p-values from the set of subgroup
    par_info = "";
    pvalue.clear();
    ExtractPvalues(single_level, par_info, iter, 0, signat5, pvalue);

    // Run permutations to compute multiplicity-adjusted treatment effect p-value for each group (under the null distribution)
    n_subgroups = pvalue.size();

    data_all_biomarkers_permuted = data_all_biomarkers;
    treatment_permuted = treatment;
    outcome_permuted = outcome;

    for (i = 0; i < n_perms_mult_adjust; i++) {

        shuffle_vector(treatment_permuted);
        for(j = 0; j<n_patients; ++j){
            data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
        }

        // Analysis of overall group of patients
        if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

        // Continuous outcome
        if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

        // Binary outcomes
        if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates); 

        vector<int> depth_hist(depth,0);

        // Find all subgroups
        single_level_permuted = FindSubgroups(data_all_biomarkers_permuted, model_covariates, n_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type, adj_parameter, width, depth, 1, depth_hist, gamma, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

        // Compute variable importance based on the permuted set
        int nterm = 0, subnterm = 0;
        vector<double> imp_permuted(n_biomarkers, 0);
        vector<int> signat;
        ComputeVarImp(single_level_permuted, imp_permuted, nterm, subnterm, signat);
        for(k=0; k < n_biomarkers; ++k)
            imp_permuted[k] /= nterm;

        vi_sorted.clear();

        // Prepare biomarkers for sorting
        for (k = 0; k < n_biomarkers; k++) {
            vi_sorted.push_back(dipair(imp_permuted[k], k));
        }

        // Sort biomarkers based on variable importance
        sort(vi_sorted.begin(), vi_sorted.end(), DIPairDown);

        // Create a new data set which includes only the top biomarkers, treatment, outcome and outcome_censor columns
        data_selected_biomarkers_permuted.clear();

        // Create a vector of biomarker types for the top biomarkers
        vector<int> biomarker_type_second_stage; 

        // Create a vector of adjustment parameters types for the top biomarkers
        vector<double> adj_parameter_second_stage; 

        biomarker_type_second_stage.clear();
        adj_parameter_second_stage.clear();

        // Select the top biomarkers based on variable importance
        for (k = 0; k < n_top_biomarkers; k++) {
            j = vi_sorted[k].second;
            current_biomarker = vector<double>(data_all_biomarkers.begin()+j * nrow, data_all_biomarkers.begin() + (j + 1) * nrow);
            // biomarker_index.push_back(j + 1);
            biomarker_type_second_stage.push_back(biomarker_type[j]);
            adj_parameter_second_stage.push_back(adj_parameter[j]);
            data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), current_biomarker.begin(), current_biomarker.end());
        }

        // Add treatment column
        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), treatment_permuted.begin(), treatment_permuted.end());

        // Add outcome column
        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome.begin(), outcome.end());

        // Add outcome_censor column
        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome_censor.begin(), outcome_censor.end());

        // Analysis of overall group of patients
        if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

        // Continuous outcome
        if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

        // Binary outcomes
        if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates); 

        vector<int> depth_hist_second_stage(depth_second_stage, 0);

        // List of all subgroups (under the null distribution)
        single_level_permuted = FindSubgroups(data_selected_biomarkers_permuted, model_covariates, n_top_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type_second_stage, adj_parameter_second_stage, width_second_stage, depth_second_stage, 1, depth_hist_second_stage, gamma_second_stage, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

        // Extract the vector of treatment effect p-values from the set of subgroup
        pvalue_permuted.clear();
        par_info = "";
        total_number_subgroups = 0;
        w = 0;
        TreeSize(single_level_permuted, total_number_subgroups, w);
        iter = 0;
        vector<int> signat13(total_number_subgroups);

        ExtractPvalues(single_level_permuted, par_info, iter, 0, signat13, pvalue_permuted);

        if (pvalue_permuted.size() > 0) {

            // Find the most significant p-value in the subgroups
            best_subgroup_pvalue[i] = *std::min_element(pvalue_permuted.begin(), pvalue_permuted.end());

        } 
        else {

        best_subgroup_pvalue[i] = 0.0;

        }


    }


   // return result
   return(best_subgroup_pvalue);    

}

// [[Rcpp::export]]
NumericVector AdaptiveSIDEScreenAdjP(
    const NumericVector &ancova_outcome_arg, const NumericVector &ancova_censor_arg, 
    const NumericVector &ancova_treatment_arg, const NumericMatrix &cont_covariates, 
    const NumericMatrix &class_covariates, 
    const int &n_cont_covariates, const int &n_class_covariates, const int &random_seed,
    const std::string project_file
) {

    int i, j, k, n_subgroups;
            
    xml_document<> doc; //create xml_document object (http://rostakagmfun.github.io/xml-file-parsing-cpp.html)
    file<> xmlFile(project_file.c_str()); //open file
    doc.parse<0>(xmlFile.data()); //parse the contents of file
    xml_node<>* root = doc.first_node("head");//find our root node
    xml_node<>* nodeStructure = root->first_node("structure");
    xml_node<>* nodeParameters = root->first_node("parameters");

    //string data_set = root->first_node("data")->value();
    int nrow = atoi(root->first_node("data")->first_attribute("nrow")->value());
    int ncol = atoi(root->first_node("data")->first_attribute("ncol")->value());
    int skip = atoi(root->first_node("data")->first_attribute("skipfirstrow")->value());
    string std_data_filename = root->first_node("data")->first_attribute("stddata")->value();

    // Read in standardized data
    std::ifstream std_data_file(std_data_filename.c_str());
    std::string data_set((std::istreambuf_iterator<char>(std_data_file)),
                         (std::istreambuf_iterator<char>()));

    //[X]min_subgroup_size (number)
    int min_subgroup_size = atoi(nodeParameters->first_node("min_subgroup_size")->first_attribute("value")->value());
    //[X]outcome_variable_direction (number)
    string outcome_variable_direction_str = nodeStructure->first_node("outcome")->first_attribute("direction")->value();
    int outcome_variable_direction = 1; // larger
    if (outcome_variable_direction_str == "smaller")
        outcome_variable_direction = -1;
    //[X]outcome_type (number)
    string outcome_type_str = nodeStructure->first_node("outcome")->first_attribute("type")->value();
    int outcome_type = 1; // continuous
    if (outcome_type_str == "binary")
        outcome_type = 2;
    if (outcome_type_str == "time")
        outcome_type = 3;
    
    //[x]covariate_types (vector)
    xml_node<>* nodeBiomarkers = nodeStructure->first_node("biomarkers");
    vector<int> covariate_types;
    covariate_types.clear();
    biomarker_level.clear();
    if (nodeBiomarkers) {
        for (xml_node<> *child = nodeBiomarkers->first_node(); child; child = child->next_sibling()) {
            covariate_types.push_back(atoi(child->first_attribute("numeric")->value()));
            biomarker_level.push_back(atoi(child->first_attribute("level")->value()));
        }
    }

    //[X]criterion_type (number)
    int criterion_type = atoi(nodeParameters->first_node("criterion_type")->first_attribute("value")->value());
    //[X]width (number)
    int width = atoi(nodeParameters->first_node("width")->first_attribute("value")->value());
    //[X]depth (number)
    int depth = atoi(nodeParameters->first_node("depth")->first_attribute("value")->value());
    //[X]gamma (vector)
    xml_node<>* nodeComplexityControl = nodeParameters->first_node("complexity_control");
    bool complexity_control = (nodeComplexityControl != NULL);
    vector<double> gamma;
    gamma.clear();
    if (complexity_control) {
        for (xml_node<> *child = nodeComplexityControl->first_node(); child; child = child->next_sibling())
            gamma.push_back(atof(child->first_attribute("value")->value()));
    } else {
        for (i = 0; i < depth; i++)
            gamma.push_back(-1.0f);
    }
    //[X]pvalue_max (number)
    // double pvalue_max = atof(nodeParameters->first_node("pvalue_max")->first_attribute("value")->value());
    //[X]local_mult_adj (number)
    double local_mult_adj = atof(nodeParameters->first_node("local_mult_adj")->first_attribute("value")->value());
    //[X]n_perms_mult_adjust (number)
    int n_perms_mult_adjust = atoi(nodeParameters->first_node("n_perms_mult_adjust")->first_attribute("value")->value());
    //[X]subgroup_search_algorithm (number)
    int subgroup_search_algorithm = atoi(nodeParameters->first_node("subgroup_search_algorithm")->first_attribute("value")->value());
    //[X]n_perms_vi_score (number)
    int n_perms_vi_score = atoi(nodeParameters->first_node("n_perms_vi_score")->first_attribute("value")->value());
    int analysis_method = atoi(nodeParameters->first_node("analysis_method")->first_attribute("value")->value());
    int nperc = atoi(nodeParameters->first_node("nperc")->first_attribute("value")->value()) - 1;

    //***************************************************************************************************************

    // Open the SIDES output file            
    xml_document<> doc_output; //create xml_document object (http://rostakagmfun.github.io/xml-file-parsing-cpp.html)
    file<> xmlFile_output("sides_output.xml"); //open file
    doc_output.parse<0>(xmlFile_output.data()); //parse the contents of file
    xml_node<>* root_output = doc_output.first_node("head");//find our root node
    double vi_threshold = atof(root_output->first_node("vi_threshold")->first_attribute("threshold")->value()); 

    // VI scores
    vector<double> variable_importance;
    variable_importance.clear();

    xml_node<>* vi_scores = root_output->first_node("vi_scores");
    if (vi_scores) {
        for (xml_node<> *child = vi_scores->first_node(); child; child = child->next_sibling()) {
            variable_importance.push_back(atof(child->first_attribute("value")->value()));
        }
    }

    //***************************************************************************************************************

    // Read in the main data set
    vector<string> CsvLines;
    boost::split(CsvLines, data_set, boost::is_any_of("\n"));
    // Remove empty lines
    for (i = CsvLines.size()-1; i >= 0; i--) {
        boost::trim(CsvLines[i]);
        if (CsvLines[i].size() == 0)
            CsvLines.erase(CsvLines.begin() + i);
    }
    // Remove header lines
    for (i = 0; i < skip; i++)
        CsvLines.erase(CsvLines.begin());


    vector<double> treatment;
    vector<double> outcome;
    vector<double> outcome_censor;

    vector<double> data_all_biomarkers;
    vector<double> data_all_biomarkers_permuted;
    vector<double> data_selected_biomarkers;
    vector<double> data_selected_biomarkers_permuted;

    ModelCovariates model_covariates;

    vector<double> yy, yy2, yy3;

    int n_biomarkers = ncol - 3;
    vector<double> adj_parameter(n_biomarkers);
    vector<int> biomarker_type(n_biomarkers);
    vector<double> current_biomarker, transformed_biomarker;
    int n_levels;

    data_all_biomarkers.clear();

    vector<string> row;

    // Type I error rate
    double error_rate = 0.0;

    for (i = 0; i < CsvLines.size(); i++) {
        boost::split(row, CsvLines[i], boost::is_any_of(","), boost::token_compress_on);    
        for (j = 0; j < ncol; j++) {
            if (row[j] != "." && row[j] != " .") {
                yy2.push_back(atof(row[j].c_str()));
            }
            else {                
                yy2.push_back(numeric_limits<double>::quiet_NaN());
            }
        }
        treatment.push_back( atof(row[ncol - 3].c_str()) );
        outcome.push_back( atof(row[ncol - 2].c_str()) );
        outcome_censor.push_back( atof(row[ncol - 1].c_str()) );
    }

    for (int i1 = 0; i1<ncol; i1++) {
        current_biomarker.clear();
        transformed_biomarker.clear();
        for (int j1 = 0; j1<nrow; j1++) {
            current_biomarker.push_back(yy2[i1 + j1*ncol]);

        }

        // Convert numeric biomarkers to percentiles
        if (nperc > 0 && CountUniqueValues(current_biomarker) > nperc && covariate_types[i1] == 1) {
            transformed_biomarker = QuantileTransform(current_biomarker, nperc); 
            data_all_biomarkers.insert(data_all_biomarkers.end(), transformed_biomarker.begin(), transformed_biomarker.end());
        } else {
            data_all_biomarkers.insert(data_all_biomarkers.end(), current_biomarker.begin(), current_biomarker.end());
        }

    }

    // Define the biomarker type and compute local multiplicity adjustment for each biomarker (if local multiplicity adjustment is enabled)
    for (i = 0; i < n_biomarkers; i++) {

        // Numerical biomarker
        if (covariate_types[i] == 1) {
            biomarker_type[i] = 1;
        }
        // Nominal biomarker
        if (covariate_types[i] == 0) {
            biomarker_type[i] = 2;
        }

        // Local multiplicity adjustment
        current_biomarker = vector<double>(data_all_biomarkers.begin()+i*nrow, data_all_biomarkers.begin()+(i+1)*nrow);
        n_levels = CountUniqueValues(current_biomarker);

        // Compute local multiplicity adjustment parameter
        adj_parameter[i] = (1.0 - local_mult_adj) + local_mult_adj * AdjParameter(biomarker_type[i], n_levels, criterion_type);


    }

    //***************************************************************************************************************

    // Read in the ANCOVA data set

    std::vector<double> ancova_outcome, ancova_treatment, temp_vec(1);

    NumericVector temp_numeric_vec;

    // Continuous endpoint
    if (outcome_type == 1 && analysis_method >= 2) {

        ancova_outcome = as<vector<double>>(ancova_outcome_arg);
        ancova_treatment = as<vector<double>>(ancova_treatment_arg);

        if (analysis_method == 2) {
            model_covariates.cov2 = temp_vec;
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            model_covariates.cov_class = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 3) {
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            model_covariates.cov_class = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
        }
        
        if (analysis_method == 4) {
            model_covariates.cov3 = temp_vec;
            model_covariates.cov4 = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 5) {
            model_covariates.cov4 = temp_vec;
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 2);
            model_covariates.cov3 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

        if (analysis_method == 6) {
            temp_numeric_vec = cont_covariates(_, 0);
            model_covariates.cov1 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 1);
            model_covariates.cov2 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 2);
            model_covariates.cov3 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = cont_covariates(_, 3);
            model_covariates.cov4 = as<vector<double>>(temp_numeric_vec);
            temp_numeric_vec = class_covariates(_, 0);
            model_covariates.cov_class = as<vector<double>>(temp_numeric_vec);
        }

    }

    // Binary endpoint
    if (outcome_type == 2 && analysis_method >= 2) {

        ancova_outcome = as<vector<double>>(ancova_outcome_arg);
        ancova_treatment = as<vector<double>>(ancova_treatment_arg);

    }

    // Survival endpoint
    if (outcome_type == 3 && analysis_method >= 2) {

        ancova_outcome = as<vector<double>>(ancova_outcome_arg);
        ancova_treatment = as<vector<double>>(ancova_treatment_arg);

    }

    //***************************************************************************************************************

    // Start subgroup search

    vector<int> parent_rows(nrow, 1);
    // vector<double> best_subgroup_pvalue(n_perms_mult_adjust, 0);
    NumericVector best_subgroup_pvalue(n_perms_mult_adjust);

    // List of found subgroups
    vector<SingleSubgroup> single_level, single_level_permuted;
    SingleSubgroup parent_group, parent_group_permuted;

    vector<double> vi_max (n_perms_vi_score);

    // Random seed
    //srand(random_seed);
    set_seed(random_seed);

    // Biomarker indices
    vector<int> biomarker_index;
    biomarker_index.clear();

    // Second-stage parameters
    int width_second_stage(width);
    int depth_second_stage(depth);
    std::vector<double> gamma_second_stage(gamma);

    // Randomly shuffled indices
    int n_patients = nrow;
    vector<int> shuffle_index(n_patients);
    vector<double> treatment_permuted(n_patients), outcome_permuted(n_patients), outcome_censor_permuted(n_patients);

    // Analysis of overall group of patients
    if (analysis_method == 1) parent_group = OverallAnalysis(treatment, outcome, outcome_censor, outcome_type, outcome_variable_direction);
    // if (outcome_type == 1 && analysis_method >= 2) parent_group = OverallAnalysisContOld(ancova_treatment, ancova_outcome, ancova_all_covariates, analysis_method);

    // Continuous outcome
    if (outcome_type == 1 && analysis_method >= 2) parent_group = OverallAnalysisCont(ancova_treatment, ancova_outcome, model_covariates, analysis_method, outcome_variable_direction);

    // Binary outcomes
    if (outcome_type == 2 && analysis_method >= 2) parent_group = OverallAnalysisBin(ancova_treatment, ancova_outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);  

    vector<int> depth_hist(depth, 0);

    // Find subgroups 
    single_level = FindSubgroups(data_all_biomarkers, model_covariates, n_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type, adj_parameter, width, depth, 1, depth_hist, gamma, parent_group.pvalue, parent_group.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);


    // Treatment effect p-values from original subgroups
    vector<double> pvalue;

    // Treatment effect p-values from permuted subgroups
    vector<double> pvalue_permuted;

     int total_number_subgroups = 0, w = 0, iter = 0;

    string par_info;
    vector<double> adj_pvalue;

   // Adaptive two-stage SIDES procedure (SIDEScreen procedure)
    if (subgroup_search_algorithm == 3) {

        data_all_biomarkers_permuted = data_all_biomarkers;
        treatment_permuted = treatment;
        outcome_permuted = outcome;

        // Create a new data set which includes only the top biomarkers, treatment, outcome and outcome_censor columns
        data_selected_biomarkers.clear();

        // Number of selected biomarkers with variable importance above the threshold
        int n_selected_biomarkers = 0;

        // Create a vector of biomarker types for the top biomarkers
        vector<int> biomarker_type_second_stage; 

        // Create a vector of adjustment parameters for the top biomarkers
        vector<double> adj_parameter_second_stage; 

        biomarker_type_second_stage.clear();
        adj_parameter_second_stage.clear();

        // Select biomarkers with variable importance above the threshold
        for (i = 0; i < n_biomarkers; i++) {

            if (variable_importance[i] >= vi_threshold) {
                n_selected_biomarkers++;
                current_biomarker = vector<double>(data_all_biomarkers.begin()+i * nrow, data_all_biomarkers.begin() + (i + 1) * nrow); 
                biomarker_index.push_back(i + 1);
                biomarker_type_second_stage.push_back(biomarker_type[i]);
                adj_parameter_second_stage.push_back(adj_parameter[i]);
                data_selected_biomarkers.insert(data_selected_biomarkers.end(), current_biomarker.begin(), current_biomarker.end());               
            }

        }

        // Find subgroups in the second stage  of the subgroup search algorithm (if at least one biomarker is selected)
        if (n_selected_biomarkers > 0) {
         
            // Add treatment column
            data_selected_biomarkers.insert(data_selected_biomarkers.end(), treatment.begin(), treatment.end());

            // Add outcome column
            data_selected_biomarkers.insert(data_selected_biomarkers.end(), outcome.begin(), outcome.end());

            // Add outcome_censor column
            data_selected_biomarkers.insert(data_selected_biomarkers.end(), outcome_censor.begin(), outcome_censor.end());

            // Analysis of overall group of patients
            if (analysis_method == 1) parent_group = OverallAnalysis(treatment, outcome, outcome_censor, outcome_type, outcome_variable_direction);

            // Continuous outcome
            if (outcome_type == 1 && analysis_method >= 2) parent_group = OverallAnalysisCont(treatment, outcome, model_covariates, analysis_method, outcome_variable_direction);

            // Binary outcomes
            if (outcome_type == 2 && analysis_method >= 2) parent_group = OverallAnalysisBin(treatment, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

            vector<int> depth_hist_second_stage(depth_second_stage, 0);

            // Find all subgroups in the second stage of the subgroup search algorithm
            single_level = FindSubgroups(data_selected_biomarkers, model_covariates, n_selected_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type_second_stage, adj_parameter_second_stage, width_second_stage, depth_second_stage, 1, depth_hist_second_stage, gamma_second_stage, parent_group.pvalue, parent_group.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

            total_number_subgroups = 0; 
            w = 0;
            TreeSize(single_level, total_number_subgroups, w);
            iter = 0;
            vector<int> signat8(total_number_subgroups);

            // Extract the vector of treatment effect p-values from the set of subgroup
            par_info = "";
            pvalue.clear();
            ExtractPvalues(single_level, par_info, iter, 0, signat8, pvalue);

            // Run permutations to compute multiplicity-adjusted treatment effect p-value for each group (under the null distribution)
            n_subgroups = pvalue.size();
            vector<double> adj_pvalue(n_subgroups, 0);

            data_all_biomarkers_permuted = data_all_biomarkers;
            treatment_permuted = treatment;
            outcome_permuted = outcome;

            for (i = 0; i < n_perms_mult_adjust; i++) {

                best_subgroup_pvalue[i] = 1.0;

                shuffle_vector(treatment_permuted);
                for(j = 0; j<n_patients; ++j){
                    data_all_biomarkers_permuted[n_biomarkers * n_patients + j] = treatment_permuted[j];
                }

                // Analysis of overall group of patients
                if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

                // Continuous outcome
                if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

                // Binary outcomes
                if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               

                vector<int> depth_hist(depth,0);

                // List of all subgroups (under the null distribution)
                single_level_permuted = FindSubgroups(data_all_biomarkers_permuted, model_covariates, n_biomarkers + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type, adj_parameter, width, depth, 1, depth_hist, gamma, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

                // Compute variable importance based on the permuted set
                int nterm = 0, subnterm = 0;
                vector<double> imp_permuted(n_biomarkers, 0);
                vector<int> signat18;
                ComputeVarImp(single_level_permuted, imp_permuted, nterm, subnterm, signat18);
                for(k=0; k < n_biomarkers; ++k)
                    imp_permuted[k] /= nterm;

                // Create a vector of biomarker types for the biomarkers above the threshold
                vector<int> biomarker_type_second_stage; 

                // Create a vector of adjustment parameters types for the biomarkers above the threshold
                vector<double> adj_parameter_second_stage; 

                data_selected_biomarkers_permuted.clear();

                biomarker_type_second_stage.clear();
                adj_parameter_second_stage.clear();

                int n_selected_biomarkers_permuted = 0;

                // Select biomarkers with variable importance above the threshold
                for (k = 0; k < n_biomarkers; k++) {

                    if (imp_permuted[k] >= vi_threshold) {
                        n_selected_biomarkers_permuted++; 
                        current_biomarker = vector<double>(data_all_biomarkers.begin()+k * nrow, data_all_biomarkers.begin() + (k + 1) * nrow);
                        // index.push_back(k + 1);
                        biomarker_type_second_stage.push_back(biomarker_type[k]);                
                        adj_parameter_second_stage.push_back(adj_parameter[k]);  
                        data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), current_biomarker.begin(), current_biomarker.end());               
                    }
                }

                if (n_selected_biomarkers_permuted > 0) {

                    // Add treatment column
                    data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), treatment_permuted.begin(), treatment_permuted.end());

                    // Add outcome column
                    data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome.begin(), outcome.end());

                    // Add outcome_censor column
                    data_selected_biomarkers_permuted.insert(data_selected_biomarkers_permuted.end(), outcome_censor.begin(), outcome_censor.end());

                    // Analysis of overall group of patients
                    if (analysis_method == 1) parent_group_permuted = OverallAnalysis(treatment_permuted, outcome, outcome_censor, outcome_type, outcome_variable_direction);

                    // Continuous outcome
                    if (outcome_type == 1 && analysis_method >= 2) parent_group_permuted = OverallAnalysisCont(treatment_permuted, outcome, model_covariates, analysis_method, outcome_variable_direction);

                    // Binary outcomes
                    if (outcome_type == 2 && analysis_method >= 2) parent_group_permuted = OverallAnalysisBin(treatment_permuted, outcome, cont_covariates, class_covariates, outcome_variable_direction, n_cont_covariates, n_class_covariates);               


                    vector<int> depth_hist_second_stage(depth_second_stage, 0);

                    // Find all subgroups in the second stage of the subgroup search algorithm
                    single_level_permuted = FindSubgroups(data_selected_biomarkers_permuted, model_covariates, n_selected_biomarkers_permuted + 3, n_patients, parent_rows, min_subgroup_size, outcome_variable_direction, criterion_type, outcome_type, analysis_method, biomarker_type_second_stage, adj_parameter_second_stage, width_second_stage, depth_second_stage, 1, depth_hist_second_stage, gamma_second_stage, parent_group_permuted.pvalue, parent_group_permuted.test_statistic, cont_covariates, class_covariates, n_cont_covariates, n_class_covariates);

                    // Extract the vector of treatment effect p-values from the set of subgroup
                    pvalue_permuted.clear();
                    par_info = "";
                    total_number_subgroups = 0;
                    w = 0;
                    TreeSize(single_level_permuted, total_number_subgroups, w);
                    iter = 0;
                    vector<int> signat14(total_number_subgroups);

                    ExtractPvalues(single_level_permuted, par_info, iter, 0, signat14, pvalue_permuted);

                    if (pvalue_permuted.size() > 0) {

                        // Find the most significant p-value in the subgroups
                        best_subgroup_pvalue[i] = *std::min_element(pvalue_permuted.begin(), pvalue_permuted.end());
                        for (j = 0; j < n_subgroups; j++) {
                            if (best_subgroup_pvalue[i] <= pvalue[j]) adj_pvalue[j]++;
                        }

                        if (best_subgroup_pvalue[i] <= 0.025) error_rate++;

                    }
                    else {

                    // Increment the adjusted p-value if no subgroups were not found in the current permutation
                    for (j = 0; j < n_subgroups; j++) {
                        adj_pvalue[j]++;
                    }

                    error_rate++;
                    best_subgroup_pvalue[i] = 1.0;

                    }


                } else {
                    best_subgroup_pvalue[i] = 1.0;
                }

            }

            // Compute multiplicity-adjusted treatment effect p-value
            for (j = 0; j < n_subgroups; j++) {
                adj_pvalue[j] /= n_perms_mult_adjust;
            }

            error_rate /= n_perms_mult_adjust;
     
        }

    }

    // return result
    return(best_subgroup_pvalue);


}

// [[Rcpp::export]]
List Quant(const NumericVector &vec_arg, const int &nperc) {

    vector<double> vec = as<vector<double>>(vec_arg);

    int i, j, m = vec.size();
    vector<double> complete_vec, res, transform(m);

    double prob;

    for (i = 0; i < m; i++) {
        if (!::isnan(vec[i])) complete_vec.push_back(vec[i]);
    }

    for (i = 0; i < nperc - 1; i++) {
        prob = (i + 1.0) / (nperc + 0.0);
        res.push_back(Quantile(vec, prob));
    }

    for (i = 0; i < m; i++) {
        if (!::isnan(vec[i])) {    
            if (vec[i] <= res[0]) transform[i] = res[0];
            if (vec[i] > res[nperc - 2]) transform[i] = res[nperc - 2];
            for (j = 0; j < nperc - 2; j++) {
                if (vec[i] > res[j] && vec[i] <= res[j + 1]) transform[i] = res[j + 1];
            }
        } else {
            transform[i] = vec[i];
        }

    }

    return List::create(Named("quant") = res,
                        Named("transform") = transform);
}
