#ifndef SHASTA_KMER_CHECKER_HPP
#define SHASTA_KMER_CHECKER_HPP

// As we transition to longer markers, we will no longer be able to store a k-mer
// table. Instead, the KmerChecker will be used to find out if a given k-mer
// is a marker. The initial implementation of the KmerChecker is table based,
// but later we will switch to hashing.

#include "Kmer.hpp"
#include "MemoryMappedVector.hpp"

namespace shasta {
    class KmerChecker;
}

class shasta::KmerChecker {
public:

    bool isMarker(KmerId kmerId) const
    {
        return data[kmerId];
    }

    void create(
        const MemoryMapped::Vector<KmerInfo>& kmerTable,
        const string& name,
        uint64_t pageSize);
    void access(const string& name);

private:
    MemoryMapped::Vector<bool> data;
};


#endif
