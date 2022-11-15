// Shasta.
#include "mode3a-BubbleCleaner.hpp"
#include "copyNumber.hpp"
#include "deduplicate.hpp"
#include "MarkerGraph.hpp"
#include "invalid.hpp"
#include "mode3a.hpp"
using namespace shasta;using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "utility.hpp"



BubbleCleaner::BubbleCleaner(
    const PackedMarkerGraph& packedMarkerGraph) :
    packedMarkerGraph(packedMarkerGraph)
{
    BubbleCleaner& bubbleCleaner = *this;

    // Construct an edge for each segment of the PackedMarkerGraph.
    // Add vertices as we go.
    for(uint64_t segmentId=0; segmentId<packedMarkerGraph.segments.size(); segmentId++) {
        const uint64_t vertexId0 = packedMarkerGraph.getFirstSegmentVertex(segmentId);
        const uint64_t vertexId1 = packedMarkerGraph.getLastSegmentVertex(segmentId);
        const vertex_descriptor v0 = getVertex(vertexId0);
        const vertex_descriptor v1 = getVertex(vertexId1);
        add_edge(v0, v1, BubbleCleanerEdge(segmentId), bubbleCleaner);
    }

}



string BubbleCleanerEdge::representation() const
{
    string s;
    for(uint64_t i=0; i<segments.size(); i++) {
        if(i != 0) {
            s += "-";
        }
        s += to_string(segments[i]);
    }
    return s;
}



// Compute average marker graph edge coverage for an edge.
double BubbleCleaner::averageEdgeCoverage(edge_descriptor e) const
{
    const BubbleCleanerEdge& edge = (*this)[e];

    uint64_t sum = 0;
    uint64_t n = 0;
    for(const uint64_t segmentId: edge.segments) {
        const auto segment = packedMarkerGraph.segments[segmentId];
        for(const uint64_t edgeId: segment) {
            const MarkerGraph::Edge& edge = packedMarkerGraph.markerGraph.edges[edgeId];
            sum += edge.coverage;
            n++;
        }
    }
    return double(sum) / double(n);
}



void BubbleCleanerEdge::assembledSequence(
    const PackedMarkerGraph& packedMarkerGraph,
    vector<Base>& sequence) const
{
    sequence.clear();
    for(const uint64_t segmentId: segments) {
        const span<const Base> segmentSequence =
            packedMarkerGraph.clippedSequence(segmentId);
        copy(segmentSequence.begin(), segmentSequence.end(), back_inserter(sequence));

    }
}





// Get the vertex corresponding to a given marker graph vertex,
// creating if necessary.
BubbleCleaner::vertex_descriptor BubbleCleaner::getVertex(uint64_t markerGraphVertexId)
{
    BubbleCleaner& bubbleCleaner = *this;

    auto it = vertexMap.find(markerGraphVertexId);
    if(it == vertexMap.end()) {
        const vertex_descriptor v = add_vertex(BubbleCleanerVertex(markerGraphVertexId), bubbleCleaner);
        vertexMap.insert(make_pair(markerGraphVertexId, v));
        return v;
    } else {
        return it->second;
    }
}



// Given the initial PackedMarkerGraph, cleanup the bubbles
// due to errors and store the result in a new PackedMarkerGraph.
void BubbleCleaner::cleanup(MarkerGraph& markerGraph)
{
    // The maximum period length that this will cleanup.
    // EXPOSE AFTER CODE STABILIZES?   ************
    const uint64_t maxPeriod = 4;

    const bool debug = false;
    BubbleCleaner& bubbleCleaner = *this;

    // A bubble is a set of parallel edges. Find them all.
    std::map<pair<vertex_descriptor, vertex_descriptor>, vector<edge_descriptor> > bubbles;
    BGL_FORALL_EDGES(e, bubbleCleaner, BubbleCleaner) {
        const vertex_descriptor v0 = source(e, bubbleCleaner);
        const vertex_descriptor v1 = target(e, bubbleCleaner);
        bubbles[make_pair(v0, v1)].push_back(e);
    }
    for(auto it=bubbles.begin(); it!=bubbles.end(); /* Increment later */) {
        auto itNext = it;
        ++itNext;
        if(it->second.size() < 2) {
            bubbles.erase(it);
        }
        it = itNext;
    }

    if(debug) {
        cout << "Found " << bubbles.size() << " bubbles." << endl;
    }



    // Cleanup the bubbles we found.
    // For now, this is not recursive.
    // This means it will only deal with one level of "bubble within a bubble"
    // situations.
    for(const auto& p: bubbles) {
        vertex_descriptor v0;
        vertex_descriptor v1;
        tie(v0, v1) = p.first;
        const vector<edge_descriptor>& bubble = p.second;

        // Get the sequence of the branches.
        vector< vector<Base> > sequences(bubble.size());
        for(uint64_t i=0; i<bubble.size(); i++) {
            const edge_descriptor e = bubble[i];
            bubbleCleaner[e].assembledSequence(packedMarkerGraph, sequences[i]);
        }


        if(debug) {
            cout << "Working on a bubble with " << bubble.size() << " branches.\n";
            cout << "MarkerGraph vertices " <<
                bubbleCleaner[v0].markerGraphVertexId << " " <<
                bubbleCleaner[v1].markerGraphVertexId << "\n";

            for(uint64_t i=0; i<bubble.size(); i++) {
                const edge_descriptor e = bubble[i];
                copy(sequences[i].begin(), sequences[i].end(), ostream_iterator<Base>(cout));
                cout << " " << bubbleCleaner[e].representation() << "\n";
            }
        }


        const uint64_t period = computeCopyNumberDifferencePeriod(sequences, maxPeriod);
        if(period == 0) {
            continue;
        }

        if(period > maxPeriod) {
            continue;
        }

        if(debug) {
            cout << "This bubble describes copy number changes in a repeat of period " << period << "\n";
        }

        // Compute average edge coverage for the branches of this bubble.
        vector<double> coverage(bubbles.size());
        for(uint64_t i=0; i<bubble.size(); i++) {
            const edge_descriptor e = bubble[i];
            coverage[i] = averageEdgeCoverage(e);
            if(debug) {
                cout << bubbleCleaner[e].representation() <<
                    " has length " << sequences[i].size() << " and coverage " << coverage[i] << "\n";
            }
        }

        // Compute a coverage-weighted average of the lengths.
        double sum = 0.;
        double sumCoverage = 0.;
        for(uint64_t i=0; i<bubble.size(); i++) {
            sum += coverage[i] * double(sequences[i].size());
            sumCoverage += coverage[i];
        }
        const double weightedAverageLength = double(sum) / double(sumCoverage);
        if(debug) {
            cout << "Coverage-weighted average length " << weightedAverageLength << "\n";
        }

        // Find the index of the branch that is closest to the weighted average.
        uint64_t iBest = invalid<uint64_t>;
        double bestDelta = 0.;
        for(uint64_t i=0; i<bubble.size(); i++) {
            const double delta = fabs(double(sequences[i].size()) - weightedAverageLength);
            if(i == 0 or delta < bestDelta) {
                iBest = i;
                bestDelta = delta;
            }
        }
        if(debug) {
            cout << "Best approximation branch is " << bubbleCleaner[bubble[iBest]].representation() << "\n";
        }


        // Flag as removed the marker graph edges of the other branches.
        for(uint64_t i=0; i<bubble.size(); i++) {
            if(i == iBest) {
                continue;
            }
            const edge_descriptor e = bubble[i];
            const BubbleCleanerEdge& edge = bubbleCleaner[e];
            for(const uint64_t segmentId: edge.segments) {
                const auto segment = packedMarkerGraph.segments[segmentId];
                for(const uint64_t edgeId: segment) {
                    markerGraph.edges[edgeId].isSuperBubbleEdge = 1;
                }
            }
        }


    }

}



// Given assembled sequences of the branches of a bubble,
// figure out if this is a bubble caused by copy number
// differences in repeats of period up to maxPeriod.
// If this is the case, returns the shortest period for which this is true.
// Otherwise, returns 0.
uint64_t BubbleCleaner::computeCopyNumberDifferencePeriod(
    const vector< vector<Base> >& sequences,
    uint64_t maxPeriod)
{

    // Check all pairs of branches.
    vector<uint64_t> periods;
    for(uint64_t i=0; i<sequences.size()-1; i++) {
        for(uint64_t j=i+1; j<sequences.size(); j++) {
            const uint64_t pairPeriod =
                shasta::isCopyNumberDifference(sequences[i], sequences[j], maxPeriod);
            if(pairPeriod == 0) {
                return 0;
            }
            periods.push_back(pairPeriod);
        }
    }
    deduplicate(periods);


    if(periods.size() == 1) {
        return periods.front();
    } else {
        return 0;
    }
}
