#include <phat/compute_persistence_pairs.h>

// main data structure (choice affects performance)
#include <phat/representations/vector_vector.h>

// algorithm (choice affects performance)
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/twist_reduction.h>

int main(int argc, char** argv) {
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
            std::cout << "Its boudary is the empty chain.";
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

