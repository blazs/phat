#include <phat/compute_persistence_pairs.h>

// main data structure (choice affects performance)
#include <phat/representations/vector_vector.h>

// algorithm (choice affects performance)
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/twist_reduction.h>

#include <fstream>
#include <iostream>
#include <string>

class Simplex {
public:
    Simplex(const int& i, const double& t, const int& d, const bool& L)
        : idx_(i), time_(t), dim_(d), inL_(L)
    { }

    void add_bd(const int& idx) { boundary_.push_back(idx); }

    int idx() const { return idx_; }
    double time() const { return time_; }
    int dim() const { return dim_; }
    bool inL() const { return inL_; }
    std::vector<int> bd() const { return boundary_; }

    bool operator<(const Simplex& rhs) const {
        bool returnP = (time() < rhs.time()) ||
            (time() == rhs.time() && dim() < rhs.dim()) ||
            (time() == rhs.time() && dim() == rhs.dim() && !inL() && rhs.inL());
        return returnP;
    }
    bool operator<=(const Simplex& rhs) const {
        return *this == rhs || *this < rhs;
    }
    bool operator==(const Simplex& rhs) const {
        // It should suffice to only check whether idx()==rhs.idx()
        return idx() == rhs.idx() && time() == rhs.time() && dim() == rhs.dim() && inL() == rhs.inL() && bd() == rhs.bd();
    }
    bool operator!=(const Simplex& rhs) const {
        return !(*this == rhs);
    }
    std::string str() const {
        std::string s = std::to_string(idx_) + " " + std::to_string(time_) + " " + std::to_string(dim_) + " " + std::to_string(inL_) + " ";
        for (int i = 0; i < boundary_.size(); ++i) { s += std::to_string(boundary_[i]) + " "; }
        return s;
    }
private:
    int idx_; // index representing the simplex
    double time_; // the time the simplex enters the filtration
    int dim_; // the dimension of the simplex
    bool inL_; // whether the simplex is in L; if false, it is in K
    std::vector<int> boundary_; // boundary of the simplex
};

// Given a filtration in the vector simplices, construct the boundary matrix and compute relative persistence
void compute_relative_persistence(std::vector<Simplex>& simplices, const std::map<int, int>& L, std::vector<phat::persistence_pairs>& pairs);
void example();
void load(std::string prefix, std::vector<Simplex>& v, std::map<int, int>& L, int& mx_dim);
void experiments(std::string prefix);
void build_bdm(std::vector<Simplex>& v, std::map<int, int>& L, phat::boundary_matrix< phat::vector_vector >& bd_m);
void test();

int main(int argc, char** argv) {
    // experiments("hi_3_0p9");
    if (argc <= 1) { return 1; }
    experiments(argv[1]);
    // experiments("hi_0p94_0p9_0p5_3");
    // test();
    // experiments("test-2");
    return 1;
    example();
    return 1;
    // first define a boundary matrix with the chosen internal representation
    phat::boundary_matrix< phat::vector_vector > bd_m;

    // set the number of columns, which equals the number of simplices in the
    // simplicial complex
    bd_m.set_num_cols(8);

    // for each column, i, set the dimension of the simplex it represents
    bd_m.set_dim(0, 0); // a_K
    bd_m.set_dim(1, 0); // b_K
    bd_m.set_dim(2, 0); // c_K
    bd_m.set_dim(3, 1); // ab_K
    bd_m.set_dim(4, 1); // ac_K
    bd_m.set_dim(5, 1); // bc_K
    bd_m.set_dim(6, 0); // a_L
    bd_m.set_dim(7, 0); // b_L

    std::vector< phat::index > tmp_col;

    // Define boundary of the 0-dimensional simplices, labeled with integers 0 through 2
    bd_m.set_col(0, tmp_col);
    bd_m.set_col(1, tmp_col);
    bd_m.set_col(2, tmp_col);

    // Now add the 2-dimensional simplex 3 whose boundary consists of 0-dimensional simplices 0 and 1
    tmp_col.push_back(0);
    tmp_col.push_back(1);
    bd_m.set_col(3, tmp_col);

    // Empty the vector 
    tmp_col.clear();

    // Repeat for 1-dimensional simpliex 4
    tmp_col.push_back(0);
    tmp_col.push_back(2);
    bd_m.set_col(4, tmp_col);

    tmp_col.clear();

    tmp_col.push_back(1);
    tmp_col.push_back(2);
    bd_m.set_col(5, tmp_col);

    tmp_col.clear();

    tmp_col.push_back(0);
    bd_m.set_col(6, tmp_col);

    tmp_col.clear();

    tmp_col.push_back(1);
    bd_m.set_col(7, tmp_col);

    // Now compute the homology via matrix reduction
    std::cout << "The boundary matrix has " << bd_m.get_num_cols() << " columns." << std::endl;
    for (phat::index col_idx = 0; col_idx < bd_m.get_num_cols(); ++col_idx) {
        std::cout << "Processing column corresponding to " << static_cast<int>(bd_m.get_dim(col_idx)) << "-dimensional simplex." << std::endl;
        if (!bd_m.is_empty(col_idx)) { // If column at the given index is not empty 
            std::vector<phat::index> tmp_col;
            bd_m.get_col(col_idx, tmp_col);
            std::cout << "Its boundary:" << std::endl;
            for (phat::index idx = 0; idx < static_cast<phat::index>(tmp_col.size()); ++idx) {
                std::cout << tmp_col[idx] << " ";
            }
            std::cout << std::endl;
        } else {
            std::cout << "Its boudary is the empty chain." << std::endl;
        }
    }
    std::cout << "Overall, the matrix has " << bd_m.get_num_entries() << " entries." << std::endl;

    // Compute persistence
    std::cout << "Persistence" << std::endl;
    std::map<int, int> L; L[6] = 0; L[7] = 1;
    std::vector<phat::persistence_pairs> pp_v(2);
    phat::compute_relative_persistence_pairs<phat::standard_reduction>(pp_v, bd_m, L);
    for (int d = 0; d < pp_v.size(); ++d) {
        std::cout << "Dimension " << d << std::endl;
        for (int idx = 0; idx < pp_v[d].get_num_pairs(); ++idx) {
            phat::index birth = pp_v[d].get_pair(idx).first;
            phat::index death = pp_v[d].get_pair(idx).second;
            std::cout << "(" << birth << ", " << death << ")" << std::endl;
        }
    }

    return 0;
}

// simplices is a vector of simplices, which, when ordered, represents the filtration
// L maps the label of each L-simplex to the label of the corresponding K-simplex
// pairs is a vector of persistence_pairs, whose d-th component corresponds to d-dimensional persistence pairs
void compute_relative_persistence(std::vector<Simplex>& simplices, const std::map<int, int>& L, std::vector<phat::persistence_pairs>& pairs) {
    std::map<int, int> label_to_idx; // Maps each simplex to the column that represents that simplex
    std::sort(simplices.begin(), simplices.end());
    const int N = simplices.size();
    int mx_dim = 0; // Dimension of the final complex
    phat::boundary_matrix< phat::vector_vector > bd_m;
    bd_m.set_num_cols(N);
    // Set dimension of each of the simplices
    for (int idx = 0; idx < simplices.size(); ++idx) {
        bd_m.set_dim(idx, simplices[idx].dim());
        label_to_idx[simplices[idx].idx()] = idx;
        mx_dim = std::max(mx_dim, simplices[idx].dim());
    }
    // Set boundary for each simplex
    std::vector< phat::index > tmp_col;
    for (int idx = 0; idx < simplices.size(); ++idx) {
        tmp_col.clear();
        if (simplices[idx].dim() > 0) {
            // The boudary consists of column indices---use label_to_idx for that
            std::vector<int> tmp_bd = simplices[idx].bd();
            for (int jdx = 0; jdx < tmp_bd.size(); ++jdx) {
                tmp_col.push_back(label_to_idx[tmp_bd[jdx]]);
            }
        }
        // If an L-simplex, add the corresponding K-simplex to its boundary
        if (simplices[idx].inL()) {
            // Need to find the corresponding K-simplex
            const int idx_K = L.find(simplices[idx].idx())->second;
            tmp_col.push_back(label_to_idx[idx_K]);
        }
        sort(tmp_col.begin(), tmp_col.end());
        bd_m.set_col(idx, tmp_col);
    }
    // Compute persistence
    //std::vector<phat::persistence_pairs> pp_v(mx_dim+1);
    std::map<int, int> L2;
    for (std::map<int, int>::const_iterator it = L.begin(); it != L.end(); ++it) { L2[label_to_idx[it->first]] = label_to_idx[it->second]; }

    // Now compute the homology via matrix reduction
/*
    std::cout << "The boundary matrix has " << bd_m.get_num_cols() << " columns." << std::endl;
    for (phat::index col_idx = 0; col_idx < bd_m.get_num_cols(); ++col_idx) {
        std::cout << "Processing column corresponding to " << static_cast<int>(bd_m.get_dim(col_idx)) << "-dimensional simplex." << std::endl;
        if (!bd_m.is_empty(col_idx)) { // If column at the given index is not empty 
            std::vector<phat::index> tmp_col;
            bd_m.get_col(col_idx, tmp_col);
            std::cout << "Its boundary:" << std::endl;
            for (phat::index idx = 0; idx < static_cast<phat::index>(tmp_col.size()); ++idx) {
                std::cout << tmp_col[idx] << " ";
            }
            std::cout << std::endl;
        } else {
            std::cout << "Its boudary is the empty chain." << std::endl;
        }
    }
    std::cout << "Overall, the matrix has " << bd_m.get_num_entries() << " entries." << std::endl;
*/
    phat::compute_relative_persistence_pairs<phat::standard_reduction>(pairs, bd_m, L);
}

void example() {
    std::vector<Simplex> v;
    v.push_back(Simplex(0, 0.1, 0, true)); // a_K
    v.push_back(Simplex(1, 0.1, 0, false)); // a _L
    v.push_back(Simplex(2, 0.4, 0, false)); // c_K
    v.push_back(Simplex(3, 0.3, 0, true)); // b_L
    v.push_back(Simplex(5, 0.3, 0, false)); // b_K
    Simplex s(4, 0.3, 1, false); s.add_bd(1); s.add_bd(5); // ab_K
    v.push_back(s);
    std::sort(v.begin(), v.end());
    for (std::vector<Simplex>::iterator it = v.begin(); it != v.end(); ++it) { std::cout << "(" << it->idx() << ", " << it->time() << ", " << it->dim() << ") " << it->inL() << std::endl; }

    std::map<int, int> ml; ml[0]=1; ml[3]=5;
    std::vector<phat::persistence_pairs> pp_v(2);

    compute_relative_persistence(v, ml, pp_v);

    for (int d = 0; d < pp_v.size(); ++d) {
        std::cout << "Dimension " << d << std::endl;
        for (int idx = 0; idx < pp_v[d].get_num_pairs(); ++idx) {
            phat::index birth = pp_v[d].get_pair(idx).first;
            phat::index death = pp_v[d].get_pair(idx).second;
            std::cout << "(" << birth << ", " << death << ")" << std::endl;
        }
    }
}

// XXX: Crashes if the file is not in the right format
void load(std::string prefix, std::vector<Simplex>& v, std::map<int, int>& L, int& mx_dim) {
    std::ifstream fin(std::string(std::string("data/") + prefix + std::string(".cx")).c_str());
    int nr_lines;
    int idx, dim, isL;
    int iK, iL;
    double time;

    mx_dim = 0;

    v.clear();
    fin >> nr_lines;
    while (nr_lines-- > 0) {
        fin >> idx >> time >> dim >> isL;
        Simplex s(idx, time, dim, isL);
        mx_dim = std::max(mx_dim, dim); // Keep track of the highest-dimensional simplex
        if (dim != 0 ) {
            while (dim+1 > 0) { fin >> idx; --dim; s.add_bd(idx); }
        }
        v.push_back(s);
    }
    fin.close();

    L.clear();
    fin.open(std::string(std::string("data/") + prefix + std::string(".L")).c_str());
    fin >> nr_lines;
    while (nr_lines-- > 0) {
        fin >> iL >> iK;
        L[iL] = iK;
    }
    fin.close();
}

void experiments(std::string prefix) {
    std::vector<Simplex> v;
    std::map<int, int> L;
    int mx_dim;
    std::cout << "Loading..." << std::endl;
    load(prefix, v, L, mx_dim);

    std::vector<phat::persistence_pairs> pp_v(mx_dim+2);
    std::cout << "Compuiting..." << std::endl;
    compute_relative_persistence(v, L, pp_v);
    std::cout << "Writing out..." << std::endl;
    for (int d = 0; d < pp_v.size(); ++d) {
        // std::cout << "Dimension " << d << std::endl;
        std::string fname = std::string("data/")+prefix+std::string("_pp_")+std::to_string(d)+std::string(".pp");
        std::ofstream fout(fname.c_str());
        fout << pp_v[d].get_num_pairs() << std::endl;
        for (int idx = 0; idx < pp_v[d].get_num_pairs(); ++idx) {
            phat::index birth = pp_v[d].get_pair(idx).first;
            phat::index death = pp_v[d].get_pair(idx).second;
            double b_tm = v[birth].time();
            double d_tm = (death >= 0 ? v[death].time() : death);
            fout << b_tm << " " << d_tm << std::endl;
            // std::cout << "(" << b_tm << ", " << d_tm << ")" << std::endl;
        }
        fout.close();
    }
}

void build_bdm(std::vector<Simplex>& v, std::map<int, int>& L, phat::boundary_matrix< phat::vector_vector >& bd_m) {
    const int N = v.size();

    bd_m.set_num_cols(N);

    for (int idx = 0; idx < v.size(); ++idx) {
        bd_m.set_dim(idx, v[idx].dim());
    }

    for (int idx = 0; idx < v.size(); ++idx) {
        std::vector< phat::index > tmp_col; // tmp_col.clear();
        if (v[idx].dim() > 0) {
            // The boudary consists of column indices---use label_to_idx for that
            std::vector<int> tmp_bd = v[idx].bd();
            for (int jdx = 0; jdx < tmp_bd.size(); ++jdx) {
                tmp_col.push_back(tmp_bd[jdx]);
            }
        }
        // If an L-simplex, add the corresponding K-simplex to its boundary
        if (v[idx].inL()) {
            // Need to find the corresponding K-simplex
            const int idx_K = L.find(v[idx].idx())->second;
            tmp_col.push_back(idx_K);
        }
        sort(tmp_col.begin(), tmp_col.end());
        bd_m.set_col(idx, tmp_col);
    }
}

void test() {
    // first define a boundary matrix with the chosen internal representation
    phat::boundary_matrix< phat::vector_vector > bd_m;

    // set the number of columns, which equals the number of simplices in the
    // simplicial complex
    bd_m.set_num_cols(19);

    // for each column, i, set the dimension of the simplex it represents
    bd_m.set_dim(0, 0); // a_K
    bd_m.set_dim(1, 0); // b_K
    bd_m.set_dim(2, 0); // c_K
    bd_m.set_dim(3, 0); // a_K
    bd_m.set_dim(4, 0); // b_K
    bd_m.set_dim(5, 0); // c_K
    bd_m.set_dim(6, 0); // a_K
    bd_m.set_dim(7, 0); // b_K
    bd_m.set_dim(8, 0); // c_K
    bd_m.set_dim(9, 0); // a_K
    bd_m.set_dim(10, 0); // b_K
    bd_m.set_dim(11, 0); // c_K
    bd_m.set_dim(12, 0); // a_K
    bd_m.set_dim(13, 0); // b_K
    bd_m.set_dim(14, 0); // b_K
    bd_m.set_dim(15, 0); // b_K
    bd_m.set_dim(16, 0); // b_K
    bd_m.set_dim(17, 1); // b_K
    bd_m.set_dim(18, 1); // b_K

    std::vector< phat::index > tmp_col;

    // Define boundary of the 0-dimensional simplices, labeled with integers 0 through 2
    bd_m.set_col(0, tmp_col);
    bd_m.set_col(1, tmp_col);
    bd_m.set_col(2, tmp_col);
    bd_m.set_col(3, tmp_col);
    bd_m.set_col(4, tmp_col);
    bd_m.set_col(5, tmp_col);
    bd_m.set_col(6, tmp_col);
    bd_m.set_col(7, tmp_col);
    bd_m.set_col(8, tmp_col);

    // Now add the 2-dimensional simplex 3 whose boundary consists of 0-dimensional simplices 0 and 1
    tmp_col.push_back(6);
    bd_m.set_col(9, tmp_col);
    tmp_col.clear();

    // Repeat for 1-dimensional simpliex 4
    tmp_col.push_back(7);
    bd_m.set_col(10, tmp_col);
    tmp_col.clear();

    tmp_col.push_back(2);
    bd_m.set_col(11, tmp_col);
    tmp_col.clear();

    tmp_col.push_back(0);
    bd_m.set_col(12, tmp_col);
    tmp_col.clear();

    tmp_col.push_back(1);
    bd_m.set_col(13, tmp_col);
    tmp_col.clear();

    tmp_col.push_back(3);
    bd_m.set_col(14, tmp_col);
    tmp_col.clear();

    tmp_col.push_back(4);
    bd_m.set_col(15, tmp_col);
    tmp_col.clear();

    tmp_col.push_back(8);
    bd_m.set_col(16, tmp_col);
    tmp_col.clear();

    tmp_col.push_back(1);
    tmp_col.push_back(6);
    bd_m.set_col(17, tmp_col);
    // tmp_col.clear();

    tmp_col.push_back(17);
    bd_m.set_col(18, tmp_col);
    tmp_col.clear();

    // Now compute the homology via matrix reduction
    std::cout << "The boundary matrix has " << bd_m.get_num_cols() << " columns." << std::endl;
    for (phat::index col_idx = 0; col_idx < bd_m.get_num_cols(); ++col_idx) {
        std::cout << "Processing column corresponding to " << static_cast<int>(bd_m.get_dim(col_idx)) << "-dimensional simplex." << std::endl;
        if (!bd_m.is_empty(col_idx)) { // If column at the given index is not empty 
            std::vector<phat::index> tmp_col;
            bd_m.get_col(col_idx, tmp_col);
            std::cout << "Its boundary:" << std::endl;
            for (phat::index idx = 0; idx < static_cast<phat::index>(tmp_col.size()); ++idx) {
                std::cout << tmp_col[idx] << " ";
            }
            std::cout << std::endl;
        } else {
            std::cout << "Its boudary is the empty chain." << std::endl;
        }
    }
    std::cout << "Overall, the matrix has " << bd_m.get_num_entries() << " entries." << std::endl;

    // Compute persistence
    std::cout << "Persistence" << std::endl;
    std::map<int, int> L; L[13] = 1; L[12] = 0; L[11] = 2; L[10] = 7; L[9] = 6; L[14]=3; L[15]=4; L[16]=8; L[18]=17;
    std::vector<phat::persistence_pairs> pp_v(2);
    phat::compute_relative_persistence_pairs<phat::standard_reduction>(pp_v, bd_m, L);
    for (int d = 0; d < pp_v.size(); ++d) {
        std::cout << "Dimension " << d << std::endl;
        for (int idx = 0; idx < pp_v[d].get_num_pairs(); ++idx) {
            phat::index birth = pp_v[d].get_pair(idx).first;
            phat::index death = pp_v[d].get_pair(idx).second;
            std::cout << "(" << birth << ", " << death << ")" << std::endl;
        }
    }
}

