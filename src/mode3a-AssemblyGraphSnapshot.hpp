#ifndef SHASTA_MODE3A_ASSEMBLY_GRAPH_SNAPSHOT_HPP
#define SHASTA_MODE3A_ASSEMBLY_GRAPH_SNAPSHOT_HPP

// Class mode3a::AssemblyGraphSnapshot creates a persistent
// snapshot of the mode3a::Assembly graph.
// It is stored using MemoryMapped classes.

// Shasta
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVector.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

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

    // Class Segment stores segmentId and id from AssemblyGraphVertex.
    // For the purpose of the snapshot, the index into this vector
    // is used as a vertex id.
    class Segment{
    public:
        uint64_t segmentId;
        uint64_t id;
        Segment(const AssemblyGraphVertex&);
        Segment() {}
    };
    MemoryMapped::Vector<Segment> segments;

    // Class Link stores vertex ids (indexes into the vertices vector).
    class Link {
    public:
        uint64_t vertexId0;
        uint64_t vertexId1;
    };
    MemoryMapped::Vector<Link> links;

    // The path of each oriented read, stored as a sequence of vertex ids
    // (that is,indexes into the segments vector above).
    // Indexed by OrientedRead::getValue().
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> paths;

    // The path entries for each vertex.
    class PathEntry {
    public:
        OrientedReadId orientedReadId;
        uint64_t position;  // The position i the path for this oriented read
    };
    MemoryMapped::VectorOfVectors<PathEntry, uint64_t> vertexPathEntries;

    // The transitions of each link don't need to be stored.
    // They can be computed from the vertexPathEntries.
    class Transition {
    public:
        OrientedReadId orientedReadId;
        uint64_t position;  // The transition if between position and position+1
    };
    void getTransitions(uint64_t linkId, vector<Transition>&) const;
};

#endif
