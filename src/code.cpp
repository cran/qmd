#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <list>

using namespace Rcpp;

//' Returns a vector 
//' 
//' @param x A usually sorted vector
//' @return A sequence along x. If consecutive values in x are equal the maximal value is used.
//' @export
// [[Rcpp::export(seq_until_changes)]]
IntegerVector seq_until_changes(const NumericVector& x) {
    R_xlen_t n = x.length();
    IntegerVector result(n);

    // As sortdata is sorted by value the position of (x, orig) is the rank of x
    // disregarding ties.
    double last_x = 0;
    R_xlen_t current_ties = 0; // How many times did we have the last x?
    for (R_xlen_t i = 0; i < n; i++) {
        if (x[i] == last_x) // Keep counting
            current_ties++;
        else {
            for (R_xlen_t j = 0; j < current_ties; j++)
                result[i-1-j] = i; // max here, could also do min and avg
            last_x = x[i];
            current_ties = 1;
        }
    }

    // The last streak should not be written yet, so we need to do it here
    if (current_ties > 0)
        for (R_xlen_t j = 0; j < current_ties; j++)
            result[n-1-j] = n;

    return result;        
}


// Returns a map between the values of x and the numbers of ties of each element of x
// Assumes that x contains elements within 1:n where n is the length of x
IntegerVector range(const IntegerVector& x) {
    IntegerVector result(x.length());

    for (R_xlen_t i = 0; i < x.length(); i++) {
        if(i%100000 == 0) Rcpp::checkUserInterrupt();
        result[(R_xlen_t) x[i] - 1] = result[(R_xlen_t) x[i] - 1] + 1;
    }

    return result;
}

// Does the same as range, but for a two dimensional array of such vectors
IntegerMatrix array_range(const IntegerVector& x, R_xlen_t cols, R_xlen_t rows) {
    IntegerMatrix result(cols, rows);

    for (R_xlen_t i = 0; i < cols; i++) {
        Rcpp::checkUserInterrupt();
        for (R_xlen_t j = 0; j < rows; j++) {
            if(j%100000 == 0) Rcpp::checkUserInterrupt();
            result(i, x(j, i) - 1) = result(i, x(j, i) - 1) + 1;
        }
    }

    return result;
}

// Iterates throw the catesian product lowers:uppers. current_index returns the position
// within the final array, current_val is an array, representing the current vector,
// aggr_lambda contains a cumulatively aggregated product of the lambdas, lowers, uppers
// are the limits and rho, resolution are the number of variables and the dimension
// of the final result.
R_xlen_t update_values(R_xlen_t current_index, R_xlen_t *current_val, double *aggr_lambda, const R_xlen_t *lowers, const R_xlen_t *uppers, const std::vector<double> *lambda, R_xlen_t rho, R_xlen_t resolution) {
    R_xlen_t k = 0;
    do {
        current_val[k]++;
        current_index = current_index / resolution;
        if (current_val[k] > uppers[k]) {
            current_val[k] = lowers[k];
            k++;
        } else {
            if (k == rho-1) aggr_lambda[k] = lambda[k][current_val[k] - lowers[k]];
            else aggr_lambda[k] = aggr_lambda[k+1] * lambda[k][current_val[k] - lowers[k]];
            current_index = current_index*resolution + current_val[k]-1;
            for (R_xlen_t l = k-1; l >= 0; l--) {
                aggr_lambda[l] = aggr_lambda[l+1] * lambda[l][current_val[l] - lowers[l]];
                current_index = current_index*resolution + current_val[l]-1;
            }
            return current_index;
        }
    } while (k < rho);
    return -1; // This should only be reached if we are done
}

//' Calculates the empirical checkerboard approximation to some data.
//' 
//' @param X A nxrho matrix of n samples of rho variables
//' @param resolution The resolution of the CB approximation
//' @return A matrix of dimension resolution^rho
//' @export
// [[Rcpp::export(.ECBC)]]
NumericVector ECBC (const IntegerMatrix& X, R_xlen_t resolution) {
    IntegerVector dim = X.attr("dim");
    R_xlen_t sample_size = dim[0];
    R_xlen_t rho = dim[1];
    IntegerVector outdim = rep((int) resolution, rho);
    NumericVector result(std::pow(resolution, rho));
    IntegerMatrix ranges = array_range(X, rho, sample_size);

    std::vector<double> *lambda = new std::vector<double>[rho];
    R_xlen_t rZ;
    R_xlen_t upper;
    R_xlen_t lower;
    R_xlen_t *uppers = new R_xlen_t[rho];
    R_xlen_t *lowers = new R_xlen_t[rho];
    R_xlen_t *current_val = new R_xlen_t[rho];
    double *aggr_lambda = new double[rho];
    R_xlen_t current_index;
    for (R_xlen_t i = 0; i < sample_size; i++) {
        if(i%100000 == 0) Rcpp::checkUserInterrupt();

        for (R_xlen_t j = 0; j < rho; j++) {
            lambda[j].clear();
            rZ = ranges(j, static_cast<R_xlen_t>(X(i, j))-1);
            upper = static_cast<R_xlen_t>(std::ceil(X(i, j) / static_cast<double>(sample_size) * resolution));
            uppers[j] = upper;
            lower = std::max(static_cast<R_xlen_t>(1 + std::floor((X(i, j) - rZ) / static_cast<double>(sample_size) * resolution)), (R_xlen_t) 1); // Use floor+1 instead of ceil, as if ceil(x)=0 then lambda=0.
            lowers[j] = lower;
            for (R_xlen_t z = lower; z <= upper; z++)
                lambda[j].push_back((std::min(static_cast<double>(X(i, j)), z / static_cast<double>(resolution) * sample_size) - std::max(static_cast<double>(X(i, j)) - rZ, (z - 1) / static_cast<double>(resolution) * sample_size)) / static_cast<double>(rZ));
        }
        std::copy(lowers, lowers+rho, current_val);
        current_index = 0;
        for (R_xlen_t k = rho-1; k >= 0; k--) current_index = current_index*resolution + current_val[k]-1;
        aggr_lambda[rho-1] = lambda[rho-1][0];
        for (R_xlen_t k = rho-2; k >= 0; k--) aggr_lambda[k] = aggr_lambda[k+1] * lambda[k][0];
        do {
            result[current_index] = result[current_index] + aggr_lambda[0] / sample_size;
            current_index = update_values(current_index, current_val, aggr_lambda, lowers, uppers, lambda, rho, resolution);
        } while (current_index >= 0);
    }

    delete[] lambda;
    delete[] uppers;
    delete[] lowers;
    delete[] current_val;
    delete[] aggr_lambda;
    result.attr("dim") = outdim;

    return result;
}

//' Calculates an empirical CB approximation with adaptive bin sizes. This will be faster on data with many ties.
//' 
//' @param X A nxrho matrix of n samples of rho variables
//' @param resolution The resolution of the CB approximation
//' @return A matrix of dimension resolution^rho
//' @export
// [[Rcpp::export(.EACBC)]]
NumericVector EACBC (const IntegerMatrix& X, R_xlen_t resolution) {
    IntegerVector dim = X.attr("dim");
    R_xlen_t sample_size = dim[0];
    R_xlen_t rho = dim[1];
    IntegerVector outdim = rep((int) resolution, rho);
    NumericVector result(std::pow(resolution, rho));
    IntegerMatrix ranges = array_range(X, rho, sample_size);

    R_xlen_t rZ;
    double lower_block;
    double upper_block;
    R_xlen_t mean_block;
    R_xlen_t current_idx;
    for (R_xlen_t i = 0; i < sample_size; i++) {
        if(i%100000 == 0) Rcpp::checkUserInterrupt();

        current_idx = 0;
        for (R_xlen_t j = rho-1; j >= 0; j--) {
            // round X(i,j) to the nearest CB-Block. This has the same effect as rounding the CB-Block to the 1/n-Grid, but doesn’t require a mapping
            // Considering each CB-Block is given by it's maximum edges this means ceiling the mean of start end end point in CB-coordinates.
            // In case of ties we round the mean of the theoretical extremal CB-Blocks to the CB-Block that would theoretically contain the longest part of the streak, effectively enlargening it to fit the whole range
            rZ = ranges(j, static_cast<R_xlen_t>(X(i, j))-1);
            lower_block = (X(i, j) - rZ) / static_cast<double>(sample_size) * resolution;
            upper_block = X(i, j) / static_cast<double>(sample_size) * resolution;
            mean_block = std::max(std::ceil((lower_block + upper_block) / 2), 1.);
            current_idx = current_idx*resolution + mean_block - 1;
        }
        result[current_idx] = result[current_idx] + 1./sample_size;
    }

    result.attr("dim") = outdim;

    return result;
}

//' Returns the sizes of the adaptive bins used for the adaptive ECBC for one vector.
//' 
//' @param X A vector, representing one sample of one variable
//' @param resolution The resolution of the CB approximation
//' @return A numeric vector of bin sizes
//' @export
// [[Rcpp::export(.adaptive_masses)]]
NumericVector adaptive_masses (const IntegerVector& X, R_xlen_t resolution) {
    R_xlen_t sample_size = X.length();
    NumericVector result(resolution);
    IntegerVector r = range(X);

    R_xlen_t rZ;
    double lower_block;
    double upper_block;
    R_xlen_t mean_block;
    for (R_xlen_t i = 0; i < sample_size; i++) {
        // round X(i) to the nearest 1/res-Block. This has the same effect as rounding the CB-Block to the 1/n-Grid, but doesn’t require a mapping
        // Considering each CB-Block is given by it's maximum edges this means ceiling the mean of start end end point in CB-coordinates.
        // In case of ties we round the mean of the theoretical extremal CB-Blocks to the CB-Block that would theoretically contain the longest part of the streak, effectively enlargening it to fit the whole range
        rZ = r[X[i]-1];
        lower_block = (X[i] - rZ) / static_cast<double>(sample_size) * resolution;
        upper_block = X[i] / static_cast<double>(sample_size) * resolution;
        mean_block = std::max(std::ceil((lower_block + upper_block) / 2.), 1.);
        result[mean_block-1] = result[mean_block-1] + 1./sample_size;
    }

    return result;
}

//' Returns a list of reverse cumulative margins of a CB copula. The nth entry is thus the copula of X1,...,Xn
//' 
//' @param CB A matrix of CB weights.
//' @return A list of CB weight matrixes of ascending dimension
//' @export
// [[Rcpp::export(.CB_make_cumulative_df)]]
List CB_make_cumulative_df (NumericVector& CB) {
    IntegerVector dim = CB.attr("dim");
    R_xlen_t res = dim[1];
    R_xlen_t rho = dim.length();
    std::list<NumericVector> outT;
    R_xlen_t currentpow = std::pow(res, rho-1);
    NumericVector *last_frame;
    for (R_xlen_t i = rho-1; i >= 0; i--) {
        // each index is of the form (xrho-1, ..., xi, ..., x0)res
        // thus the possible coordinates for i are given by A + xires^i + B
        // where A < res^i and res^i | B, thus B = res^i+1 b where b < res^(rho-i-1)
        // Also the upper margin of b has b = res^(rho-i-1)-1
        NumericVector this_frame(currentpow*res);
        if (i < rho-1) {
            last_frame = &(outT.front());
        }
        for (R_xlen_t A = 0; A < currentpow; A++) {
            for (R_xlen_t xi = 0; xi < res; xi++) {
                if (i < rho-1) {
                    this_frame[A+currentpow*xi] = (*last_frame)[A+currentpow*xi];
                    for (R_xlen_t xj = 1; xj < res; xj++)
                        this_frame[A+currentpow*xi] = (*last_frame)[A+currentpow*xi+res*currentpow*xj] + this_frame[A+currentpow*xi];
                } else {
                    this_frame[A+currentpow*xi] = CB[A+currentpow*xi];
                }
            }
        }
        if (i>0) currentpow = currentpow / res;
        IntegerVector outdim = rep((int) res, i+1);
        this_frame.attr("dim") = outdim;
        outT.push_front(this_frame);
    }

    return wrap(outT);
}

NumericVector indep_CB(R_xlen_t rho, R_xlen_t steps) {
    R_xlen_t res = 1 << steps;
    R_xlen_t N = 1 << (steps*rho);
    NumericVector RES = rep(1. / N, N);
    RES.attr("dim") = rep((int) res, rho);

    return RES;
}

NumericVector random_CB_dep(R_xlen_t rho, R_xlen_t steps) {
    R_xlen_t res = 1 << steps;
    IntegerMatrix S(res, rho);
    IntegerVector s(res);
    for (R_xlen_t j = 0; j < rho; j++) {
        s = sample(res, res);
        for (R_xlen_t i = 0; i < res; i++) {
            S(i, j) = s[i];
        }
    }

    return ECBC(S, res);
}


//' Creates a random CB copula of resolution 2^steps
//' 
//' @param rho The number of variables
//' @param steps Number of iteration steps, the final resolution will be 2^steps
//' @param de Exponent to increase dependence
//' @param ie Exponent to increase independence
//' @return A matrix of dimension (2^steps)^rho
//' @export
// [[Rcpp::export(.random_CB)]]
NumericVector random_CB(R_xlen_t rho, R_xlen_t steps, double de, double ie) {
    Rcpp::checkUserInterrupt();
    double my = std::pow(runif(1)[0], de);
    double lambda = std::pow(runif(1)[0], ie);
    NumericVector depres = random_CB_dep(rho, 1);
    NumericVector result = runif(1 << rho);
    double S = sum(result);
    for (R_xlen_t i = 0; i < result.length(); i++) {
        result[i] = result[i] / S;
    }
    R_xlen_t pow = 1;
    R_xlen_t cpow = 1 << (rho-1);
    for (R_xlen_t j = 0; j < rho; j++) {
        double margin = 0;
        for (R_xlen_t A = 0; A < pow; A++) {
            for (R_xlen_t B = 0; B < cpow; B++) {
                margin = margin + result[A+B*2*pow];
            }
        }
        if (margin > 0.5) {
            for (R_xlen_t A = 0; A < pow; A++) {
                for (R_xlen_t B = 0; B < cpow; B++) {
                    double tmp = result[A+B*2*pow];
                    result[A+B*2*pow] = tmp*(0.5/margin);
                    result[A+B*2*pow+pow] = result[A+B*2*pow+pow] + tmp*(1-0.5/margin);
                }
            }
        } else {
            for (R_xlen_t A = 0; A < pow; A++) {
                for (R_xlen_t B = 0; B < cpow; B++) {
                    double tmp = result[A+B*2*pow+pow];
                    result[A+B*2*pow+pow] = tmp*(0.5/(1-margin));
                    result[A+B*2*pow] = result[A+B*2*pow] + tmp*(1-0.5/(1-margin));
                }
            }
        }
        pow <<= 1;
        cpow >>= 1;
    }
    result = (my*result + (1-my)*depres)*lambda + (1-lambda)*indep_CB(rho, 1);
    result.attr("dim") = rep(2, rho);
    if (steps <= 1) return result;

    NumericVector RESULT(1 << (steps*rho));
    NumericVector sub_CB;
    for (R_xlen_t i = 0; i < result.length(); i++) {
        sub_CB = random_CB(rho, steps-1, de, ie);
        for (R_xlen_t j = 0; j < sub_CB.length(); j++) {
            R_xlen_t J = 0;
            R_xlen_t _j = j;
            R_xlen_t _i = i;
            R_xlen_t p = 1 << (steps-1);
            R_xlen_t q = 2;
            for (R_xlen_t l = 0; l < rho; l++) {
                R_xlen_t m1 = _j%p;
                R_xlen_t m2 = _i%q;
                R_xlen_t d1 = _j/p;
                R_xlen_t d2 = _i/q;
                J = J + ((m1 << (l*steps)) + (m2 << (l*steps + (steps-1))));
                _j = d1;
                _i = d2;
            }
            RESULT[J] = result[i] * sub_CB[j];
        }
    }
    RESULT.attr("dim") = rep((int) (1 << steps), rho);
    return RESULT;
}

//' Generate a sample of some CB copula-
//' 
//' @param CB A weight matrix of a CB copula
//' @param n The number of samples to be generated
//' @return Matrix of dimension nxm where m is the dimension of CB
//' @export
// [[Rcpp::export(.sample_CB)]]
NumericVector sample_CB(NumericVector &CB, R_xlen_t n) {
    IntegerVector D = CB.attr("dim");
    R_xlen_t rho = D.length();
    R_xlen_t res = D[0];
    List CDF = CB_make_cumulative_df(CB);
    NumericVector result(n*rho);
    IntegerVector indexes(n*rho);
    IntegerVector idx;
    NumericVector x;
    idx = sample(res, n, true);
    x = runif(n);
    for (R_xlen_t j = 0; j<n; j++) {
        result[j] = (idx[j]-1 + x[j])/res;
        indexes[j] = idx[j];
    }
    for (R_xlen_t i = 1; i < rho; i++) {
        NumericVector L = CDF[i];
        x = runif(n);
        for (R_xlen_t j = 0; j < n; j++) {
            R_xlen_t IDX = 0;
            for (R_xlen_t k = i-1; k >=0; k--) {
                IDX = IDX*res + indexes[j + n*k]-1;
            }
            NumericVector probs(res);
            for (R_xlen_t k = 0; k < res; k++) {
                probs[k] = L[IDX + k*std::pow(res, i)];
            }
            IntegerVector this_idx = sample(res, 1, true, probs);
            idx[j] = this_idx[0];
        }
        for (R_xlen_t j = 0; j<n; j++) {
            result[j+i*n] = (idx[j]-1 + x[j])/res;
            indexes[j+i*n] = idx[j];
        }
    }
    IntegerVector odim(2);
    odim[0] = n;
    odim[1] = rho;
    result.attr("dim") = odim;
    return(result);
}

//' Computes the D1-difference of two CB matrizes on a local CB dimension
//' 
//' @param k1 Vector of local CB weights of first matrix
//' @param k2 Vector of local CB weights of second matrix
//' @param y Vector indicating the bin sizes of the local dimension
//' @return number indicating the difference between k1 and k2
//' @export
// [[Rcpp::export(.local_kernel_integral)]]
double local_kernel_integral(const NumericVector& k1, const NumericVector& k2, const NumericVector& y) {
    R_xlen_t N = std::min(std::min(k1.length(), k2.length()), y.length());
    double k1sum = 0;
    double pk1sum = 0;
    double k2sum = 0;
    double pk2sum = 0;
    double ysum = 0;
    double x0;
    double xr;
    double I1;
    double I1b;
    double I2;
    double I2b;
    double I;

    for (R_xlen_t i = 0; i < N; i++) {
        k1sum = k1sum + k1[i];
        k2sum = k2sum + k2[i];
        // case 1: no intersection
        if ((k1sum - k2sum) * (pk1sum - pk2sum) >= 0) {
            I1 = k1[i] / (2 * N) + pk1sum / N;
            I2 = k2[i] / (2 * N) + pk2sum / N;
            I = std::abs(I1 - I2);
        } else {
            x0 = (pk1sum - pk2sum) / ((k2[i] - k1[i]) * N);
            xr = ((double) 1)/N - x0;

            I1 = (N * k1[i] * x0 * x0) / 2 + pk1sum * x0;
            I1b = (N * k2[i] * x0 * x0) / 2 + pk2sum * x0;
            I2 = (N * k1[i] * xr * xr) / 2;
            I2b = ((N * k2[i] * xr * xr)) / 2;
            I = std::abs(I1 - I1b) + std::abs(I2 - I2b);
        }
        ysum = ysum + I * y[i] * N;
        pk1sum = k1sum;
        pk2sum = k2sum;
    }

    return ysum;
}


void lexsort(std::vector<std::vector<R_xlen_t>*>::iterator begin, std::vector<std::vector<R_xlen_t>*>::iterator end, R_xlen_t cbegin, R_xlen_t cend) {
    R_xlen_t min = (*begin)->at(cbegin);
    R_xlen_t max = (*begin)->at(cbegin);
    for (std::vector<std::vector<R_xlen_t>*>::iterator it = begin; it < end; it++) {
        R_xlen_t val = (*it)->at(cbegin);
        if (val < min) min = val;
        if (val > max) max = val;
    }
    std::vector<std::vector<R_xlen_t>*>* D = new std::vector<std::vector<R_xlen_t>*>[max-min+1];
    for(std::vector<std::vector<R_xlen_t>*>::iterator it = begin; it < end; it++) {
        D[(*it)->at(cbegin)-min].push_back(*it);
    }
    if (cbegin < cend)
        for(R_xlen_t i = 0; i <= max-min; i++)
            if(D[i].size() > 1)
                lexsort(D[i].begin(), D[i].end(), cbegin+1, cend);
    R_xlen_t current = 0;
    std::vector<std::vector<R_xlen_t>*>::iterator cb = D[current].begin();
    std::vector<std::vector<R_xlen_t>*>::iterator ce = D[current].end();
    for(std::vector<std::vector<R_xlen_t>*>::iterator it = begin; it < end; it++) {
        while(cb==ce) {
            current++;
            cb = D[current].begin();
            ce = D[current].end();
        }
        *it = *cb;
        cb++;
    }
}

//' Returns non 0 entries of the EACBC
//' 
//' @param X A nxrho matrix of n samples of rho variables
//' @param resolution The resolution of the CB approximation
//' @return A list of local kernel masses
//' @export
// [[Rcpp::export(.EACBC_nonzero)]]
List EACBC_nonzero (const IntegerMatrix& X, R_xlen_t resolution) {
    IntegerVector dim = X.attr("dim");
    R_xlen_t sample_size = dim[0];
    R_xlen_t rho = dim[1];
    std::vector<NumericVector> *result = new std::vector<NumericVector>;
    std::vector<std::vector<R_xlen_t>*> blocks(sample_size);
    IntegerMatrix ranges = array_range(X, rho, sample_size);

    R_xlen_t rZ;
    double lower_block;
    double upper_block;
    R_xlen_t mean_block;
    for (R_xlen_t i = 0; i < sample_size; i++) {
        blocks[i] = new std::vector<R_xlen_t>(rho);
        if(i%100000 == 0) Rcpp::checkUserInterrupt();
        for (R_xlen_t j = 0; j < rho; j++) {
            // round X(i,j) to the nearest CB-Block. This has the same effect as rounding the CB-Block to the 1/n-Grid, but doesn’t require a mapping
            // Considering each CB-Block is given by it's maximum edges this means ceiling the mean of start end end point in CB-coordinates.
            // In case of ties we round the mean of the theoretical extremal CB-Blocks to the CB-Block that would theoretically contain the longest part of the streak, effectively enlargening it to fit the whole range
            rZ = ranges(j, static_cast<R_xlen_t>(X(i, j))-1);
            lower_block = (X(i, j) - rZ) / static_cast<double>(sample_size) * resolution;
            upper_block = X(i, j) / static_cast<double>(sample_size) * resolution;
            mean_block = std::max(std::ceil((lower_block + upper_block) / 2), 1.);
            (*blocks[i])[j] = mean_block;
        }
    }

    lexsort(blocks.begin(), blocks.end(), 0, rho-1);

    R_xlen_t *prev_val = new R_xlen_t[rho-1];
    for (R_xlen_t j = 0; j < rho-1; j++) prev_val[j] = (*blocks[0])[j];

    NumericVector this_vector(resolution);
    bool ready_for_write = false;
    for (R_xlen_t i = 0; i < sample_size; i++) {
        if(i%100000 == 0) Rcpp::checkUserInterrupt();
        for (R_xlen_t j = 0; j < rho; j++) {
            if (j == rho-1) {
                this_vector[(*blocks[i])[j]-1] = this_vector[(*blocks[i])[j]-1] + ((double) 1) / sample_size;
                ready_for_write = true;
            } else if((*blocks[i])[j] != prev_val[j]) {
                prev_val[j] = (*blocks[i])[j];
                if (ready_for_write) {
                    result->push_back(clone(this_vector));
                    this_vector.fill(0);
                    ready_for_write = false; // don’t write on further changes in this sample, wait for vector to populate
                }
            }
        }
    }

    // push out last value:
    result->push_back(this_vector);

    for (std::vector<std::vector<R_xlen_t>*>::iterator it = blocks.begin(); it < blocks.end(); it++) delete (*it);

    delete[] prev_val;

    List ret = wrap(*result);
    return ret;
}
