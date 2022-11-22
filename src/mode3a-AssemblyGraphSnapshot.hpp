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

        class PackedMarkerGraph;
    }
}



class shasta::mode3a::AssemblyGraphSnapshot :
    public MappedMemoryOwner {
public:

    // This creates a snapshot of the AssemblyGraph in the current state.
    AssemblyGraphSnapshot(
        const AssemblyGraph&,
        const string& name,
        const MappedMemoryOwner& mappedMemoryOwner);

    // This accesses an existing snapshot.
    AssemblyGraphSnapshot(
        const string& name,
        const MappedMemoryOwner& mappedMemoryOwner,
        const PackedMarkerGraph& packedMarkerGraph);

    const string name;

    // The MarkerGraph and PackedMarkerGraph that this AssemblyGraphSnapshot refers to.
    const PackedMarkerGraph& packedMarkerGraph;

    // Class Vertex stores segmentId and id from AssemblyGraphVertex.
    // For the purpose of the snapshot, the index into this vector
    // is used as a vertex id.
    class Vertex {
    public:
        uint64_t segmentId;
        uint64_t segmentCopyIndex;
        Vertex(const AssemblyGraphVertex&);
        Vertex() {}
        string stringId() const
        {
            string s = to_string(segmentId);
            if(segmentCopyIndex != 0) {
                s += "." + to_string(segmentCopyIndex);
            }
            return s;
        }
    };
    MemoryMapped::Vector<Vertex> vertexVector;  // Can't call it vertices due to boost graph macros.

    // Class Edge stores vertex ids (indexes into the vertices vector).
    class Edge {
    public:
        uint64_t vertexId0;
        uint64_t vertexId1;
    };
    MemoryMapped::Vector<Edge> edgeVector;  // Can't call it edges due to boost graph macros.

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

    // The transitions of each edge (link) don't need to be stored.
    // They can be computed from the vertexPathEntries.
    class Transition {
    public:
        OrientedReadId orientedReadId;
        uint64_t position;  // The transition if between position and position+1
    };
    void getEdgeTransitions(uint64_t edgeId, vector<Transition>&) const;
    uint64_t getEdgeCoverage(uint64_t edgeId) const;

    // Find out if the segments of an edge are adjacent in the marker graph.
    bool segmentsAreAdjacent(uint64_t edgeId) const;

    // Connectivity.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> edgesBySource;
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> edgesByTarget;

    // A data structure that allows to get a vertexId (index in vertexVector)
    // given a segmentId and segmentCopyIndex.
    // Indexed by segmentId.
    // vertexTable[segmentId][segmentCopyIndex] contains the vertexId
    // with the given segmentId and segmentCopyIndex, or invalid<uint64_t>
    // if no such vertex.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> vertexTable;
    void createVertexTable(const PackedMarkerGraph&);

    // Get the length of assembled sequence for a vertex.
    uint64_t getVertexAssembledSequenceLength(uint64_t vertexId) const;

    void write() const;
    void writeGfa(uint64_t minLinkCoverage) const;
    void writePaths() const;
};

#endif
