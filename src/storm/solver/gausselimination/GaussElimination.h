#pragma once

#include "storm/storage/SparseAccessForDenseMatrix.h"
#include "storm/storage/FlexibleSparseMatrix.h"
#include "storm/storage/HTMLMatrixExport.h"
#include "storm/storage/SCCInfo.h"
#include <memory>
#include <string>


namespace storm {
namespace solver {
namespace gausselimination {

template <typename Matrix, typename Worker>
class GaussElimination {
public:
    static void eliminate(Matrix& matrix,
                          boost::optional<const storm::storage::SCCInfo&> sccInfo,
                          Worker& worker,
                          bool debug = false,
                          std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport = nullptr) {
        if (sccInfo && storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().isIndividualSCCsSet()) {
            STORM_LOG_INFO("Performing Gauss elimination for each SCC separately...");
            eliminateWithSCCs(matrix, sccInfo.get(), worker, debug, htmlExport);
        } else {
            eliminateAll(matrix, worker, debug, htmlExport);
        }
    }

    static void eliminateAll(Matrix& matrix,
                             Worker& worker,
                             bool debug,
                             std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport) {
        std::size_t rows = matrix.getRowCount();
        std::size_t cols = matrix.getColumnCount();

        assert(rows <= cols);

        if (htmlExport) {
            htmlExport->slideBegin();
            htmlExport->print("<h2>Initial</h2>");
            htmlExport->printMatrix(matrix);
            htmlExport->print("<hr>");
            htmlExport->slideEnd();
        }

        typedef typename Matrix::value_type value_type;

        worker.startTriangulation(matrix);
        for (std::size_t row = 0; row < rows; row++) {
            const value_type& a_kk = matrix(row, row);
            worker.startTriangulation(matrix, row);
            for (std::size_t i = row + 1; i < rows; i++) {
                const value_type& a_ik = matrix(i,row);
                worker.startTriangulation(matrix, row, i, a_kk, a_ik);
                for (std::size_t j = row + 1; j < cols; j++) {
                    const value_type& a_ij = matrix(i,j);
                    const value_type& a_kj = matrix(row,j);
                    matrix(i,j) = worker.doTriangulation(matrix, row, i, j, a_kk, a_ij, a_ik, a_kj);
                }
                matrix(i, row) = std::move(value_type());
                // worker.endTriangulation(matrix, row, i);
            }

            worker.endTriangulation(matrix, row, a_kk);
            if (debug) std::cout << "After row " << row << ":\n" << matrix << "\n";

            if (htmlExport) {
                htmlExport->slideBegin();
                htmlExport->print("<h2>After triangulation for row " + std::to_string(row) + "</h2>");
                htmlExport->printMatrix(matrix, row);
                htmlExport->print("<hr>");
                htmlExport->slideEnd();
            }
        }

        worker.endTriangulation(matrix);
    }

    static void eliminateWithSCCs(Matrix& matrix,
                                  const storm::storage::SCCInfo& sccInfo,
                                  Worker& worker,
                                  bool debug,
                                  std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport) {
        std::size_t rows = matrix.getRowCount();
        std::size_t cols = matrix.getColumnCount();

        assert(rows <= cols);

        if (htmlExport) {
            htmlExport->slideBegin();
            htmlExport->print("<h2>Initial</h2>");
            htmlExport->printMatrix(matrix);
            htmlExport->print("<hr>");
            htmlExport->slideEnd();
        }

        worker.startTriangulation(matrix);
        for (std::size_t row = 0; row < rows; row++) {
            worker.startTriangulation(matrix, row);
            const typename Matrix::value_type& a_kk = matrix(row, row);
            std::size_t cur_scc = sccInfo.getSCC(row);
            for (std::size_t i = row + 1; i < rows; i++) {
                if (sccInfo.getSCC(i) != cur_scc) {
                    // state i is already in another SCC
                    break;
                }
                const typename Matrix::value_type& a_ik = matrix(i,row);
                worker.startTriangulation(matrix, row, i, a_kk, a_ik);
                for (std::size_t j = row + 1; j < cols; j++) {
                    const typename Matrix::value_type& a_ij = matrix(i,j);
                    const typename Matrix::value_type& a_kj = matrix(row,j);
                    matrix(i,j) = worker.doTriangulation(matrix, row, i, j, a_kk, a_ij, a_ik, a_kj);
                }
            }

            worker.endTriangulation(matrix, row, a_kk);
            if (debug) std::cout << "After row " << row << ":\n" << matrix << "\n";

            if (htmlExport) {
                htmlExport->slideBegin();
                htmlExport->print("<h2>After triangulation for row " + std::to_string(row) + "</h2>");
                htmlExport->printMatrix(matrix, row);
                htmlExport->print("<hr>");
                htmlExport->slideEnd();
            }
        }

        worker.endTriangulation(matrix);
    }

    template <typename OtherValueType>
    static void solve(Matrix& matrix, Worker& worker, std::vector<OtherValueType>& x, bool debug = false) {
        storm::storage::SparseAccessForDenseMatrix<Matrix> safdm(matrix);

        worker.startSolve_sparse(safdm, x);
        for (std::size_t row = safdm.getRowCount(); row > 0; row--) {
            if (debug) std::cout << "Solving row " << (row-1) << ":\n";
            x.at(row - 1) = std::move(worker.solve_sparse(safdm, row - 1, x));
            if (debug) std::cout << "  " << x.at(row-1) << "\n";
        }
        worker.postprocessSolution(safdm, x);
    }

};

/// ----------------- specialization for FlexibleSparseMatrix ----------------------------------------

template <typename ValueType, typename Worker>
class GaussElimination<storm::storage::FlexibleSparseMatrix<ValueType>, Worker> {
public:
    typedef storm::storage::FlexibleSparseMatrix<ValueType> FlexibleMatrix;

    static void eliminate(FlexibleMatrix& matrix,
                          boost::optional<const storm::storage::SCCInfo&> sccInfo,
                          Worker& worker,
                          bool debug = false,
                          std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport = nullptr) {
        if (sccInfo && storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().isIndividualSCCsSet()) {
            STORM_LOG_INFO("Performing Gauss elimination for each SCC separately...");
            eliminateWithSCCs(matrix, sccInfo.get(), worker, debug, htmlExport);
        } else {
            eliminateAll(matrix, worker, debug, htmlExport);
        }
    }

    static void eliminateAll(FlexibleMatrix& matrix,
                             Worker& worker,
                             bool debug,
                             std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport) {
        std::size_t rows = matrix.getRowCount();
        std::size_t cols = matrix.getColumnCount();

        assert(rows <= cols);

        if (htmlExport) {
            htmlExport->slideBegin();
            htmlExport->print("<h2>Initial</h2>");
            htmlExport->printMatrix(matrix);
            htmlExport->print("<hr>");
            htmlExport->slideEnd();
        }

        typedef typename FlexibleMatrix::value_type value_type;

        worker.startTriangulation(matrix);
        for (std::size_t row = 0; row < rows; row++) {
            const value_type& a_kk = matrix.getValueLeft(row, row);
            worker.startTriangulation(matrix, row);
            for (std::size_t i = row + 1; i < rows; i++) {
                const value_type& a_ik = matrix.getValueLeft(i,row);
                worker.startTriangulation(matrix, row, i, a_kk, a_ik);

                auto it_row_i  = matrix.row_begin(i, row+1);
                auto end_row_i = matrix.row_end(i);

                auto it_row_k  = matrix.row_begin(row, row+1);
                auto end_row_k = matrix.row_end(row);

                typename FlexibleMatrix::row_type tmpRow;

                while (it_row_i != end_row_i || it_row_k != end_row_k) {
                    std::size_t col_i = it_row_i != end_row_i ? it_row_i->getColumn() : cols;
                    std::size_t col_k = it_row_k != end_row_k ? it_row_k->getColumn() : cols;
                    std::size_t j = std::min(col_i, col_k);
                    const value_type& a_ij = (j == col_i ? it_row_i->getValue() : matrix.getZeroValue());
                    const value_type& a_kj = (j == col_k ? it_row_k->getValue() : matrix.getZeroValue());

                    tmpRow.emplace_back(j, worker.doTriangulation(matrix, row, i, j, a_kk, a_ij, a_ik, a_kj));

                    if (j == col_i) ++it_row_i;
                    if (j == col_k) ++it_row_k;
                }

                matrix.replaceRow(i, std::move(tmpRow));
            }

            worker.endTriangulation(matrix, row, a_kk);
            if (debug) std::cout << "After row " << row << ":\n" << matrix << "\n";

            if (htmlExport) {
                htmlExport->slideBegin();
                htmlExport->print("<h2>After triangulation for row " + std::to_string(row) + "</h2>");
                htmlExport->printMatrix(matrix, row);
                htmlExport->print("<hr>");
                htmlExport->slideEnd();
            }
        }

        worker.endTriangulation(matrix);
    }

    static void eliminateWithSCCs(FlexibleMatrix& matrix,
                                  const storm::storage::SCCInfo& sccInfo,
                                  Worker& worker,
                                  bool debug,
                                  std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport) {
        std::size_t rows = matrix.getRowCount();
        std::size_t cols = matrix.getColumnCount();

        assert(rows <= cols);

        if (htmlExport) {
            htmlExport->slideBegin();
            htmlExport->print("<h2>Initial</h2>");
            htmlExport->printMatrix(matrix);
            htmlExport->print("<hr>");
            htmlExport->slideEnd();
        }

        typedef typename FlexibleMatrix::value_type value_type;

        worker.startTriangulation(matrix);
        for (std::size_t row = 0; row < rows; row++) {
            worker.startTriangulation(matrix, row);
            const value_type& a_kk = matrix.getValueLeft(row, row);
            std::size_t cur_scc = sccInfo.getSCC(row);
            for (std::size_t i = row + 1; i < rows; i++) {
                if (sccInfo.getSCC(i) != cur_scc) {
                    // state i is already in another SCC
                    break;
                }

                const value_type& a_ik = matrix.getValueLeft(i,row);
                worker.startTriangulation(matrix, row, i, a_kk, a_ik);
                auto it_row_i  = matrix.row_begin(i, row+1);
                auto end_row_i = matrix.row_end(i);

                auto it_row_k  = matrix.row_begin(row, row+1);
                auto end_row_k = matrix.row_end(row);

                typename FlexibleMatrix::row_type tmpRow;

                while (it_row_i != end_row_i || it_row_k != end_row_k) {
                    std::size_t col_i = it_row_i != end_row_i ? it_row_i->getColumn() : cols;
                    std::size_t col_k = it_row_k != end_row_k ? it_row_k->getColumn() : cols;
                    std::size_t j = std::min(col_i, col_k);
                    const value_type& a_ij = (j == col_i ? it_row_i->getValue() : matrix.getZeroValue());
                    const value_type& a_kj = (j == col_k ? it_row_k->getValue() : matrix.getZeroValue());

                    tmpRow.emplace_back(j, worker.doTriangulation(matrix, row, i, j, a_kk, a_ij, a_ik, a_kj));

                    if (j == col_i) ++it_row_i;
                    if (j == col_k) ++it_row_k;
                }

                matrix.replaceRow(i, std::move(tmpRow));
            }

            worker.endTriangulation(matrix, row, a_kk);
            if (debug) std::cout << "After row " << row << ":\n" << matrix << "\n";

            if (htmlExport) {
                htmlExport->slideBegin();
                htmlExport->print("<h2>After triangulation for row " + std::to_string(row) + "</h2>");
                htmlExport->printMatrix(matrix, row);
                htmlExport->print("<hr>");
                htmlExport->slideEnd();
            }
        }

        worker.endTriangulation(matrix);
    }

    template <typename OtherValueType>
    static void solve(FlexibleMatrix& matrix, Worker& worker, std::vector<OtherValueType>& x, bool debug = false) {
        worker.startSolve_sparse(matrix, x);
        for (std::size_t row = matrix.getRowCount(); row > 0; row--) {
            if (debug) std::cout << "Solving row " << (row-1) << ":\n";
            x.at(row - 1) = std::move(worker.solve_sparse(matrix, row - 1, x));
            if (debug) std::cout << "  " << x.at(row-1) << "\n";
        }
        worker.postprocessSolution(matrix, x);
    }

};

}
}
}
