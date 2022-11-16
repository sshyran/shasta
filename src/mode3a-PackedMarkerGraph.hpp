#ifndef SHASTA_MODE3A_PACKED_MARKER_GRAPH_HPP
#define SHASTA_MODE3A_PACKED_MARKER_GRAPH_HPP

// See comments in mode3a.hpp.

// Shasta.
#include "Base.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"

// Standard library.
#include "string.hpp"

namespace shasta {

    namespace mode3a {
        class PackedMarkerGraph;
    }
    extern template class MultithreadedObject<mode3a::PackedMarkerGraph>;

    class CompressedMarker;
    class MarkerGraph;
}



// In the PackedMarkerGraph, each segment corresponds exactly to a
// path in the marker graph.
class shasta::mode3a::PackedMarkerGraph :
    public MappedMemoryOwner,
    public MultithreadedObject<PackedMarkerGraph> {
public:
    PackedMarkerGraph(
        const MappedMemoryOwner&,
        const string& name,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&);

    // Objects stored by the constructor.
    string name;
    uint64_t k;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;

    // For each segment, we store its path in the marker graph.
    // The path is a sequence of marker graph edge ids.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> segments;

    void createSegmentsFromMarkerGraph();

    // Get the first or last marker graph vertex of a segment.
    uint64_t getFirstSegmentVertex(uint64_t segmentId) const;
    uint64_t getLastSegmentVertex (uint64_t segmentId) const;

    // Assembled sequence of each segment.
    MemoryMapped::VectorOfVectors<Base, uint64_t> segmentSequences;
    void assembleSegmentSequences();

    // The segmentSequences stores, for each segment,
    // the entire sequence from the AssembledSegment.
    // This includes the entire sequence of the first
    // and last vertex of each segment.
    // This returns the clipped sequence of each segment,
    // which excludes the first and last k/2 bases.
    span<const Base> clippedSequence(uint64_t segmentId) const;

    // A link between segments s0 and s1 is created if
    // the last marker graph vertex of s0 coincides with the
    // first marker graph vertex of s1.
    class Link {
    public:
        uint64_t segmentId0;
        uint64_t segmentId1;
        Link() {}
        Link(
            uint64_t segmentId0,
            uint64_t segmentId1) :
            segmentId0(segmentId0),
            segmentId1(segmentId1) {}

    };
    MemoryMapped::Vector<Link> links;

    // This creates the links once the segments are available.
    void createLinks();

    // Given the links, create the connectivity.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> linksBySource;
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> linksByTarget;
    void createConnectivity();

    void writeGfa() const;

    void remove();

    // Keep track of the segment and position each marker graph edge corresponds to.
    // For each marker graph edge, store in the marker graph edge table
    // the corresponding segment id and position in the path, if any.
    // Indexed by the edge id in the marker graph.
    // This is needed when computing markerGraphJourneys below.
    MemoryMapped::Vector< pair<uint64_t, uint64_t> > markerGraphEdgeTable;
    void createMarkerGraphEdgeTable(uint64_t threadCount);
    void createMarkerGraphEdgeTableThreadFunction(uint64_t threadId);
};

#endif

