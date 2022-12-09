#ifndef SHASTA_MODE3A_PACKED_MARKER_GRAPH_HPP
#define SHASTA_MODE3A_PACKED_MARKER_GRAPH_HPP

// See comments in mode3a.hpp.

// Shasta.
#include "Base.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"

namespace shasta {

    namespace mode3a {
        class PackedMarkerGraph;
    }
    extern template class MultithreadedObject<mode3a::PackedMarkerGraph>;

    class CompressedMarker;
    class MarkerGraph;
    class OrientedReadId;
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
        const MarkerGraph&,
        bool accessExisting);

    // Objects stored by the constructor.
    string name;
    uint64_t k;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;

    // For each segment, we store its path in the marker graph.
    // The path is a sequence of marker graph edge ids.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> segments;
    void writeSegments();

    void createSegmentsFromMarkerGraph();

    // Get the first or last marker graph vertex of a segment.
    uint64_t getFirstSegmentVertex(uint64_t segmentId) const;
    uint64_t getLastSegmentVertex (uint64_t segmentId) const;

    // Assembled sequence of each segment.
    // This stores, for each segment, the sequence from the AssembledSegment
    // with the first and last k/2 bases removed.
    MemoryMapped::VectorOfVectors<Base, uint64_t> segmentSequences;
    void assembleSegmentSequences();

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

    void accessSegments();
    void accessSegmentSequences();
    void accessLinks();
    void accessJourneys();

    // Given the links, create the connectivity.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> linksBySource;
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> linksByTarget;
    void createConnectivity();

    void writeGfa() const;

    void remove();



    // The journey of an oriented read is the sequence of segments it encounters.
    // Indexed by OrientedReadId::getValue().
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> journeys;
    void computeJourneys(uint64_t threadCount);
    class ComputeJourneysData {
    public:

        // Keep track of the segment each marker graph edge corresponds to.
        // For each marker graph edge, store in the marker graph edge table
        // the corresponding segment id, if any.
        // Indexed by the edge id in the marker graph.
        // This is needed when computing oriented read journeys below.
        MemoryMapped::Vector<uint64_t> markerGraphEdgeTable;

        // In passes 1 and 2, we create pairs(ordinal, segmentId)
        // for each oriented read.
        // Indexed by OrientedReadId::getValue().
        // We do this by looping over the marker graph edges of each segment.
        // Pass 1 is for counting and pass 2 is for storing.
        // These pairs can contain duplicates.
        MemoryMapped::VectorOfVectors<pair<uint32_t, uint64_t>, uint64_t> journeyPairs;

        // Temporary storage of the journeys of each oriented read.
        // Indexed by OrientedReadId::getValue().
        // This is created in pass 3.
        vector< vector<uint64_t> > journeys;
    };
    ComputeJourneysData computeJourneysData;
    void createMarkerGraphEdgeTable(uint64_t threadCount);
    void createMarkerGraphEdgeTableThreadFunction(uint64_t threadId);
    void computeJourneysPass1ThreadFunction(uint64_t threadId);
    void computeJourneysPass2ThreadFunction(uint64_t threadId);
    void computeJourneysPass12ThreadFunction(uint64_t pass);
    void computeJourneysPass3ThreadFunction(uint64_t threadId);


    void writeJourneys() const;



};

#endif

