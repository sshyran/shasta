#include "KmerChecker.hpp"
using namespace shasta;


void KmerChecker::create(
    const MemoryMapped::Vector<KmerInfo>& kmerTable,
    const string& name,
    uint64_t pageSize)
{
    data.createNew(name, pageSize);
    data.resize(kmerTable.size());
    for(uint64_t i=0; i<data.size(); i++) {
        data[i] = kmerTable[i].isMarker;
    }
}

void KmerChecker::access(const string& name)
{
    data.accessExistingReadOnly(name);
}
