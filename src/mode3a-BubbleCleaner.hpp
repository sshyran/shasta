#ifndef SHASTA_MODE3A_BUBBLE_CLEANER_HPP
#define SHASTA_MODE3A_BUBBLE_CLEANER_HPP

/*******************************************************************************

Class BubbleCleaner takes as input the initial PackedMarkerGraph
and flags marker graph edges of branch bubbles that are
likely to be caused by errors.

*******************************************************************************/

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3a {

        class BubbleCleaner;
        class BubbleCleanerVertex;
        class BubbleCleanerEdge;
        using BubbleCleanerBaseClass = boost::adjacency_list<
            boost::listS, boost::listS, boost::bidirectionalS,
            BubbleCleanerVertex, BubbleCleanerEdge>;

        class PackedMarkerGraph;
    }

    class Base;
    class MarkerGraph;
}



// Each vertex of the BubbleCleaner corresponds to a marker graph vertex.
class shasta::mode3a::BubbleCleanerVertex {
public:
    uint64_t markerGraphVertexId;
    BubbleCleanerVertex(uint64_t markerGraphVertexId) :
        markerGraphVertexId(markerGraphVertexId) {}
};



// Each edge of the BubbleCleaner corresponds to a path in the marker graph.
// The path is described as a sequence of segments in the PackedMarkerGraph.
class shasta::mode3a::BubbleCleanerEdge {
public:
    vector<uint64_t> segments;
    BubbleCleanerEdge(uint64_t segmentId) : segments(1, segmentId) {}
    BubbleCleanerEdge() {}

    string representation() const;
    void assembledSequence(const PackedMarkerGraph&, vector<Base>&) const;
};



class shasta::mode3a::BubbleCleaner : public BubbleCleanerBaseClass {
public:

    // Construct the BubbleCleaner from the PackedMarkerGraph.
    BubbleCleaner(const PackedMarkerGraph&);

    // Clean up the bubbles causes by errors.
    // This flags marker graph edges of bubble branches
    // likely to be errors.
    void cleanup(MarkerGraph&);

private:

    // A reference to the PackedMarkerGraph before bubble cleanup.
    const PackedMarkerGraph& packedMarkerGraph;

    // Get the vertex corresponding to a given marker graph vertex,
    // creating if necessary.
    vertex_descriptor getVertex(uint64_t markerGraphVertexId);

    // Map marker graph vertices to vertex descriptors.
    std::map<uint64_t, vertex_descriptor> vertexMap;

    // Hide BubbleCleanerBaseClass::Base.
    using Base = shasta::Base;

    // A bubble is a set of parallel edges between two vertices (v0, v1).
    using Bubble = vector<edge_descriptor>;
    using VertexPair = pair<vertex_descriptor, vertex_descriptor>;
    using Bubbles = std::map<VertexPair, Bubble>;
    void findBubbles(Bubbles&) const;

    // Get assembled sequences of the branches of a bubble.
    void getBubbleSequences(
        const Bubble&,
        vector< vector<Base> >&) const;

    // Compute average edge coverage for the branches of a bubble.
    void computeBubbleCoverage(
        const Bubble&,
        vector<double>&) const;

    // Given assembled sequences of the branches of a bubble,
    // figure out if this is a bubble caused by copy number
    // differences in repeats of period up to maxPeriod.
    // If this is the case, returns the shortest period for which this is true.
    // Otherwise, returns 0.
    static uint64_t computeCopyNumberDifferencePeriod(
        const vector< vector<Base> >& sequences,
        uint64_t maxPeriod);

    // Merge an edge with its only previous/next edge, if possible.
    // If the merge  was done this returns true and eNew is the
    // newly created edge.
    // Otherwise, this returns false and eNew is set equal to e.
    bool mergeWithPreviousIfPossible(
        edge_descriptor e,
        edge_descriptor& eNew
        );
    bool mergeWithNextIfPossible(
        edge_descriptor e,
        edge_descriptor& eNew
        );

    // Compute average marker graph edge coverage for an edge.
    double averageEdgeCoverage(edge_descriptor e) const;
};



#endif
