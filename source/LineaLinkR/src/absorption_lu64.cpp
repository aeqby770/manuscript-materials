#include <Eigen/SparseLU>
#include <cstdint>
#include <cstdio>
#include <exception>

// CSC input and dense RHS/output use R's zero-based, column-major layout.
// Factorization uses 64-bit indices while preserving the supplied equation.
extern "C" int llr_absorption_lu64(
    int size, const int* outer, const int* inner, const double* values,
    const double* rhs, int n_rhs, double* output, char* error, int error_size) {
    using Sparse32 = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
    using Sparse64 = Eigen::SparseMatrix<double, Eigen::ColMajor, std::int64_t>;
    static_assert(sizeof(Sparse64::StorageIndex) == 8, "64-bit factor indices required");
    try {
        Eigen::Map<const Sparse32> mapped(size, size, outer[size], outer, inner, values);
        Sparse64 matrix = mapped;
        Eigen::Map<const Eigen::MatrixXd> right_hand_side(rhs, size, n_rhs);
        Eigen::SparseLU<Sparse64, Eigen::COLAMDOrdering<std::int64_t>> factor;
        factor.analyzePattern(matrix);
        factor.factorize(matrix);
        if (factor.info() != Eigen::Success) {
            std::snprintf(error, error_size, "%s", factor.lastErrorMessage().c_str());
            return 1;
        }
        Eigen::Map<Eigen::MatrixXd> solution(output, size, n_rhs);
        solution = factor.solve(right_hand_side);
        if (factor.info() != Eigen::Success) {
            std::snprintf(error, error_size, "64-bit SparseLU solve failed");
            return 1;
        }
        return 0;
    } catch (const std::exception& failure) {
        std::snprintf(error, error_size, "%s", failure.what());
        return 1;
    } catch (...) {
        std::snprintf(error, error_size, "64-bit SparseLU failed with a native exception");
        return 1;
    }
}
