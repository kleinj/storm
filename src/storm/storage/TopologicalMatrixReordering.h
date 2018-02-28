#pragma once

#include "storm/storage/SparseMatrix.h"
#include "storm/storage/FlexibleSparseMatrix.h"
#include "storm/storage/PermutatedDenseMatrix.h"
#include "storm/storage/Permutation.h"

#include <tuple>

namespace storm {
namespace storage {

class TopologicalMatrixReordering {
public:
    struct TopologicalInfo {
        TopologicalInfo(std::size_t rowCount) :
            sccInfo(new storm::storage::SCCInfo),
            permutation(new storm::storage::Permutation(rowCount)) {
        }

        storm::storage::SCCInfo::ptr sccInfo;
        storm::storage::Permutation::ptr permutation;
    };

    template <typename ValueType, typename SCCDecomposition>
    static std::tuple<PermutatedDenseMatrix<ValueType>, TopologicalInfo>
    topologicallySorted(DenseMatrix<ValueType>& matrix, const SCCDecomposition& sccs) {
        TopologicalInfo info(matrix.getRowCount());
        constructPermutation(sccs, info);

        PermutatedDenseMatrix<ValueType> sortedMatrix(matrix, info.permutation);

        return std::make_tuple(std::move(sortedMatrix), std::move(info));
    }

    template <typename ValueType, typename SCCDecomposition>
    static std::tuple<FlexibleSparseMatrix<ValueType>&, TopologicalInfo>
    topologicallySorted(FlexibleSparseMatrix<ValueType>& matrix, const SCCDecomposition& sccs) {
        TopologicalInfo info(matrix.getRowCount());

        constructPermutation(sccs, info);
        matrix.permute(*info.permutation);

        return std::tuple<FlexibleSparseMatrix<ValueType>&, TopologicalInfo>(matrix, std::move(info));
    }

    template <typename ValueType, typename SCCDecomposition>
    static std::tuple<SparseMatrix<ValueType>, TopologicalInfo>
    topologicallySorted(const SparseMatrix<ValueType>& matrix, const SCCDecomposition& sccs) {
        TopologicalInfo info(matrix.getRowCount());
        constructPermutation(sccs, info);

        SparseMatrix<ValueType> sortedMatrix = matrix.permute(*info.permutation);

        return std::make_tuple(std::move(sortedMatrix), std::move(info));
    }


private:
    template <typename SCCDecomposition>
    static void constructPermutation(const SCCDecomposition& sccs, TopologicalInfo& info) {
        storm::storage::SCCInfo& sccInfo = *info.sccInfo;
        storm::storage::Permutation& permutation = *info.permutation;

        std::size_t index = 0;
        for (auto scc_it = sccs.rbegin(); scc_it != sccs.rend(); scc_it++) {
            auto& scc = *scc_it;
            sccInfo.startSCC(scc.isTrivial());

            if (debug) std::cout << " -----------------\n";
            for (auto state_it = scc.rbegin(); state_it != scc.rend(); state_it++) {
                std::size_t state = *state_it;
                if (debug) std::cout << " " << state << " => index " << index << "\n";
                permutation.set(state, index);
                sccInfo.addState(index);
                index++;
            }

            sccInfo.endSCC();
        }

        if (debug) {
            std::cout << " -----------------\n\n";
            std::cout << permutation;
        }
    }


private:
    static const bool debug = false;
};

}
}

