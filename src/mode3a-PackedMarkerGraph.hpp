#ifndef SHASTA_MODE3A_PACKED_MARKER_GRAPH_HPP
#define SHASTA_MODE3A_PACKED_MARKER_GRAPH_HPP

// See comments in mode3a.hpp.

// Shasta.
#include "Base.hpp"
#include "invalid.hpp"
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
    class Reads;
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
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&,
        bool accessExisting);

    // Objects stored by the constructor.
    string name;
    uint64_t k;
    const Reads& reads;
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

    // Compute average marker graph edge coverage for a segment.
    double averageMarkerGraphEdgeCoverage(uint64_t segmentId) const;

    // Assembled sequence of each segment.
    // This stores, for each segment, the entire sequence from the AssembledSegment.
    MemoryMapped::VectorOfVectors<Base, uint64_t> segmentSequences;
    MemoryMapped::VectorOfVectors<uint32_t, uint64_t> segmentVertexOffsets; // Filled in by assembleSegment.
    void assembleSegmentSequences();

    // Accessors for the segment sequences.
    // The first one returns the complete sequence as stored and as obtained from
    // the AssembledSegment. It includes the entire sequences of the first
    // and last vertex of the segment.
    span<const Base> segmentCompleteSequence(uint64_t segmentId) const
    {
        return segmentSequences[segmentId];
    }

    // This returns the sequence with the first and last k/2 bases removed.
    span<const Base> segmentClippedSequence(uint64_t segmentId) const
    {
        return span<const Base>(
            segmentSequences.begin(segmentId) + k/2,
            segmentSequences.end  (segmentId) - k/2);
    }

    uint64_t totalSegmentLength() const
    {
        return segmentSequences.totalSize();
    }



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



    // Class used to store a detailed representation of the journeys.
    class JourneyStep {
    public:
        uint64_t segmentId = invalid<uint64_t>;

        // The oriented read appears in a number of marker graph edges in this segment.
        // We store information on the source vertex of the first of those edges
        // and the target vertex of the last of those edges.
        // For each of those two vertices, we store the position of the vertex
        // in the segment and the ordinal of the oriented read in the vertex.
        array<uint64_t, 2> positions;
        array<uint32_t, 2> ordinals;

    };



    // The journey of an oriented read is the sequence of segments it encounters.
    // For each segment we also store some additional information (see class JourneyStep above).
    // Indexed by OrientedReadId::getValue().
    MemoryMapped::VectorOfVectors<JourneyStep, uint64_t> journeys;
    void computeJourneys(uint64_t threadCount);
    class ComputeJourneysData {
    public:

        // In passes 1 and 2, we compute the journey of
        // each oriented read in the marker graph - that is,
        // the sequence of marker graph edges it encounters.
        // Indexed by OrientedReadId::getValue().
        class MarkerGraphJourneyStep {
        public:
            uint64_t segmentId;
            uint64_t positionInSegment;
            uint64_t markerGraphEdgeId;
            uint32_t ordinal0;
            uint32_t ordinal1;

            // Order by position in the oriented read.
            bool operator<(const MarkerGraphJourneyStep& that) const
            {
                return ordinal0 < that.ordinal0;
            }
        };
        MemoryMapped::VectorOfVectors<MarkerGraphJourneyStep, uint64_t> markerGraphJourneys;

        // More detailed representation of the journeys.
        vector< vector<JourneyStep> > journeys;

    };
    ComputeJourneysData computeJourneysData;
    void computeJourneysPass1ThreadFunction(uint64_t threadId);
    void computeJourneysPass2ThreadFunction(uint64_t threadId);
    void computeJourneysPass12ThreadFunction(uint64_t pass);
    void computeJourneysPass3ThreadFunction(uint64_t threadId);


    void writeJourneys() const;



};

#endif

