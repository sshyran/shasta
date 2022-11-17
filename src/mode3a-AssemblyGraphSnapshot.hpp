#ifndef SHASTA_MODE3A_ASSEMBLY_GRAPH_SNAPSHOT_HPP
#define SHASTA_MODE3A_ASSEMBLY_GRAPH_SNAPSHOT_HPP

// Class mode3a::AssemblyGraphSnapshot creates a persistent
// snapshot of the mode3a::Assembly graph.
// It is stored using MemoryMapped classes.

// Shasta
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVector.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

// Standard library.
#include "cstdint.hpp"
#include "string.hpp"

namespace shasta {
    namespace mode3a {
        class AssemblyGraphSnapshot;

        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;
    }
}



class shasta::mode3a::AssemblyGraphSnapshot :
    public MappedMemoryOwner {
public:
    AssemblyGraphSnapshot(
        const AssemblyGraph&,
        const string& name,
        const MappedMemoryOwner& mappedMemoryOwner);

    const string name;

    // Class SegmentInfo stores segmentId and id from AssemblyGraphVertex.
    // For the purpose of the snapshot, the index info this vector
    // is used as a vertex id.
    class SegmentInfo{
    public:
        uint64_t segmentId;
        uint64_t id;
        SegmentInfo(const AssemblyGraphVertex&);
        SegmentInfo() {}
    };
    MemoryMapped::Vector<SegmentInfo> segments;

    // Class LinkInfo stores vertex ids (indexes into the vertices vector).
    class LinkInfo {
    public:
        uint64_t vertexId0;
        uint64_t vertexId1;
        LinkInfo(uint64_t vertexId0, uint64_t vertexId1);
        LinkInfo() {}
    };
    MemoryMapped::Vector<LinkInfo> links;
};

#endif
