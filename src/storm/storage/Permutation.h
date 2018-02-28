#pragma once

#include <vector>
#include <memory>

namespace storm {
namespace storage {

class Permutation {
public:
    typedef std::shared_ptr<Permutation> ptr;

    Permutation(std::size_t size) : perm(size), inv(size) {
    }

    Permutation(const Permutation& other) = default;
    Permutation(Permutation&& other) = default;
    Permutation& operator=(const Permutation& other) = default;
    Permutation& operator=(Permutation&& other) = default;

    std::size_t size() const {
        return perm.size();
    }

    void set(std::size_t from, std::size_t to) {
        perm.at(from) = to;
        inv.at(to) = from;
    }

    std::size_t get(std::size_t from) const {
        if (from >= size()) {
            // element out of range of permutation -> identity
            return from;
        }
        return perm[from];
    }

    std::size_t getInv(std::size_t to) const {
        if (to >= size()) {
            // element out of range of permutation -> identity
            return to;
        }
        return inv[to];
    }

    template <typename T>
    void permute(std::vector<T>& v) const {
        std::vector<T> tmp;
        tmp.reserve(v.size());
        for (std::size_t i = 0; i < v.size(); i++) {
            tmp.push_back(std::move(v.at(getInv(i))));
        }
        v.swap(tmp);
    }

    template <typename T>
    void permuteInv(std::vector<T>& v) const {
        std::vector<T> tmp;
        tmp.reserve(v.size());
        for (std::size_t i = 0; i < v.size(); i++) {
            tmp.push_back(std::move(v.at(get(i))));
        }
        v.swap(tmp);
    }

    friend
    std::ostream& operator<<(std::ostream& out, const Permutation& permutation);

private:
    std::vector<std::size_t> perm;
    std::vector<std::size_t> inv;
};

inline std::ostream& operator<<(std::ostream& out, const Permutation& permutation) {
    out << "Permutation:\n";
    for (std::size_t i = 0; i < permutation.size(); i++) {
        out << i << " -> " << permutation.get(i) << "\n";
    }
    return out;
}

}
}

