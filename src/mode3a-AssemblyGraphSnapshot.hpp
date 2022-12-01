#ifndef SHASTA_MODE3A_ASSEMBLY_GRAPH_SNAPSHOT_HPP
#define SHASTA_MODE3A_ASSEMBLY_GRAPH_SNAPSHOT_HPP

// Class mode3a::AssemblyGraphSnapshot creates a persistent
// snapshot of the mode3a::Assembly graph.
// It is stored using MemoryMapped classes.

// Shasta.
#include "mode3a-AssemblyGraph.hpp"
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

    // Class Vertex stores the segmentId and segmentReplicaIndex from AssemblyGraphVertex.
    // For the purpose of the snapshot, the index into this vector
    // is used as a vertex id.
    class Vertex {
    public:
        uint64_t segmentId;
        uint64_t segmentReplicaIndex;
        Vertex(const AssemblyGraphVertex&);
        Vertex() {}
        string stringId() const
        {
            string s = to_string(segmentId);
            if(segmentReplicaIndex != 0) {
                s += "." + to_string(segmentReplicaIndex);
            }
            return s;
        }
    };
    MemoryMapped::Vector<Vertex> vertexVector;  // Can't call it vertices due to boost graph macros.

    // Get the stringId for a given vertexId, or "None" if vertexId is invalid<uint64_t>.
    string vertexStringId(uint64_t vertexId) const;

    // Class Edge stores vertex ids (indexes into the vertices vector).
    class Edge {
    public:
        uint64_t vertexId0;
        uint64_t vertexId1;
    };
    MemoryMapped::Vector<Edge> edgeVector;  // Can't call it edges due to boost graph macros.

    // The journey of each oriented read, stored as a sequence of vertex ids
    // (that is, indexes into the segments vector above).
    // Indexed by OrientedRead::getValue().
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> journeys;

    // The path entries for each vertex.
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



    // Analyze the simple tangle at a given vertex vertexId1 by
    // "following the reads" one step forward or backward.
    // On return:
    // - The previousVertices and nextVertices vectors
    //   have size equal to the number of path entries for vertexId1,
    //   and contain the vertices vertexId0 and vertexId2 each of
    //   the oriented reads visits immediately before/after visiting vertexId1.
    //   This can be invalid<uint64_t> if vertexId1 is the
    //   first or last vertex of a path.
    // - inCoverage[vertexId0] gives the number of oriented reads that
    //   visit vertexId0 immediately before visiting vertexId1.
    // - outCoverage[vertexId2] gives the number of oriented reads that
    //   visit vertexId0 immediately after visiting vertexId1.
    // - tangle_matrix[make_pair(vertexId0, vertexId2)] is
    //   the number of oriented reads that visit vertexId0 immediately before vertexId1, and
    //   visit vertexId2 immediately after  vertexId1.
    // The three maps can include entries for which vertexId0 and/or vertexId2
    // are invalid<uint64_t>, if vertexId1 is the first or last
    // vertex of one or more oriented read paths.
    void analyzeSimpleTangleAtVertex(
        uint64_t vertexId1,
        vector<uint64_t>& inVertices,
        vector<uint64_t>& outVertices,
        std::map<uint64_t, uint64_t>& inCoverage,
        std::map<uint64_t, uint64_t>& outCoverage,
        std::map< pair<uint64_t, uint64_t>, uint64_t>& tangleMatrix) const;



    void write() const;
    void writeGfa(uint64_t minLinkCoverage) const;
    void writeJourneys() const;
    void writePathEntries() const;
    void writeTransitions() const;
};

#endif
