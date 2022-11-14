// Shasta.
#include "mode3a.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "AssembledSegment.hpp"
#include "copyNumber.hpp"
#include "deduplicate.hpp"
#include "invalid.hpp"
#include "MarkerGraph.hpp"
#include "Reads.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <map>


Assembler::Assembler(
    uint64_t threadCount,
    uint64_t k, // Marker length
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    MarkerGraph& markerGraph) :
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    markers(markers),
    markerGraph(markerGraph)
{

    // This requires the marker length k to be even.
    SHASTA_ASSERT((k % 2) == 0);

    // This does not work with RLE.
    SHASTA_ASSERT(reads.representation == 0);

    // Clear all the superbubble fags in marker graph edges -
    // just in case we alreayd ran this before.
    for( MarkerGraph::Edge& edge: markerGraph.edges) {
        edge.isSuperBubbleEdge = 0;
    }

    // Create the initial PackedMarkerGraph.
    const string name0 = "Mode3a-PackedMarkerGraph-0";
    PackedMarkerGraph packedMarkerGraph0(name0, k, MappedMemoryOwner(*this), markers, markerGraph);
    packedMarkerGraph0.assembleSegmentSequences(name0);
    cout << "The initial PackedMarkerGraph has " <<
        packedMarkerGraph0.segments.size() << " segments, " <<
        packedMarkerGraph0.links.size() << " links, and " <<
        packedMarkerGraph0.segmentSequences.totalSize() <<
        " bases of assembled sequence." << endl;
    packedMarkerGraph0.writeGfa(name0);

    // Clean up the bubbles causes by errors.
    // This keeps one branch of each bubble.
    // The marker graph edges of the remaining branches are flagged as removed.
    BubbleCleaner cleaner(packedMarkerGraph0);
    cleaner.cleanup(markerGraph);

    // Create the cleaned up PackedMarkerGraph.
    const string name1 = "Mode3a-PackedMarkerGraph-1";
    PackedMarkerGraph packedMarkerGraph1(name1, k, MappedMemoryOwner(*this), markers, markerGraph);
    packedMarkerGraph1.assembleSegmentSequences(name1);
    cout << "After bubble cleanup, the PackedMarkerGraph has " <<
        packedMarkerGraph1.segments.size() << " segments, " <<
        packedMarkerGraph1.links.size() << " links, and " <<
        packedMarkerGraph1.segmentSequences.totalSize() <<
        " bases of assembled sequence." << endl;
    packedMarkerGraph1.writeGfa(name1);
}



// Initial creation from the marker graph.
PackedMarkerGraph::PackedMarkerGraph(
    const string& name,
    uint64_t k,
    const MappedMemoryOwner& mappedMemoryOwner,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph) :
    MappedMemoryOwner(mappedMemoryOwner),
    k(k),
    markers(markers),
    markerGraph(markerGraph)
{
    createSegmentsFromMarkerGraph(name);
    createLinks(name);
}



void PackedMarkerGraph::createSegmentsFromMarkerGraph(const string& name)
{
    createNew(segments, name + "-segments");

    const bool debug = false;

    const uint64_t edgeCount = markerGraph.edges.size();
    vector<bool> wasFound(edgeCount, false);

    using MarkerGraphPath = vector<uint64_t>;
    MarkerGraphPath nextEdges;
    MarkerGraphPath previousEdges;
    MarkerGraphPath path;
    MarkerGraphPath reverseComplementedPath;

    // Main loop over all edges of the marker graph.
    // At each iteration we find a new linear path of edges.
    for(uint64_t startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {

        if(markerGraph.edges[startEdgeId].wasRemoved()) {
            continue;
        }

        // If we already found this edge, skip it.
        // It is part of a path we already found.
        if(wasFound[startEdgeId]) {
            continue;
        }

        if(debug) {
            cout << "Starting a new path at edge " << startEdgeId << endl;
        }

        // Follow the path forward.
        nextEdges.clear();
        MarkerGraph::EdgeId edgeId = startEdgeId;
        bool isCircular = false;
        while(true) {
            const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
            SHASTA_ASSERT(not edge.wasRemoved());
            const MarkerGraph::VertexId v1 = edge.target;
            if(markerGraph.outDegree(v1) != 1) {
                break;
            }
            if(markerGraph.inDegree(v1) != 1) {
                break;
            }
            edgeId = markerGraph.getFirstNonRemovedOutEdge(v1);
            if(edgeId == startEdgeId) {
                isCircular = true;
                break;
            }
            nextEdges.push_back(edgeId);
            SHASTA_ASSERT(not wasFound[edgeId]);
            if(debug) {
                cout << "Moving forward: added " << edgeId << endl;
            }
        }

        // Follow the path backward.
        previousEdges.clear();
        if(!isCircular) {
            edgeId = startEdgeId;
            while(true) {
                const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
                SHASTA_ASSERT(not edge.wasRemoved());
                const MarkerGraph::VertexId v0 = edge.source;
                if(markerGraph.outDegree(v0) != 1) {
                    break;
                }
                if(markerGraph.inDegree(v0) != 1) {
                    break;
                }
                edgeId = markerGraph.getFirstNonRemovedInEdge(v0);
                previousEdges.push_back(edgeId);
                SHASTA_ASSERT(not wasFound[edgeId]);
                if(debug) {
                    cout << "Moving backward: added " << edgeId << endl;
                }
            }
        }

        // Gather the path.
        path.clear();
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter(path));
        path.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter(path));

        // Mark all the edges in the path as found.
        for(const MarkerGraph::EdgeId edgeId: path) {
            if(wasFound[edgeId]) {
                cout << "Assertion failed at " << edgeId << endl;
                SHASTA_ASSERT(0);
            }
            wasFound[edgeId] = true;
        }

        // Store this path as a new segment.
        segments.appendVector(path);
    }



    // Check that all non-removed edges of the marker graph were found.
    for(uint64_t edgeId=0; edgeId<markerGraph.edges.size(); edgeId++) {
        const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
        SHASTA_ASSERT(edge.wasRemoved() or wasFound[edgeId]);
    }
}



// Get the first or last marker graph vertex of a PackedMarkerGraph segment.
uint64_t PackedMarkerGraph::getFirstSegmentVertex(uint64_t segmentId) const
{
    const auto segment = segments[segmentId];
    SHASTA_ASSERT(not segment.empty());
    const uint64_t firstEdgeId = segment.front();
    const MarkerGraph::Edge& edge = markerGraph.edges[firstEdgeId];
    return edge.source;
}
uint64_t PackedMarkerGraph::getLastSegmentVertex (uint64_t segmentId) const
{
    const auto segment = segments[segmentId];
    SHASTA_ASSERT(not segment.empty());
    const uint64_t lastEdgeId = segment.back();
    const MarkerGraph::Edge& edge = markerGraph.edges[lastEdgeId];
    return edge.target;
}




void PackedMarkerGraph::createLinks(const string& name)
{
    // Index the segments by their first vertex.
    // This can be made faster.
    std::map<uint64_t, vector<uint64_t> > segmentMap;    // Key=first vertex, Value=segments
    for(uint64_t segmentId=0; segmentId!=segments.size(); segmentId++) {
        segmentMap[getFirstSegmentVertex(segmentId)].push_back(segmentId);
    }


    // Create the links.
    createNew(links, name + "-links");
    for(uint64_t segmentId0=0; segmentId0!=segments.size(); segmentId0++) {
        const auto& v = segmentMap[getLastSegmentVertex(segmentId0)];
        for(const uint64_t segmentId1: v) {
            links.push_back(Link(segmentId0, segmentId1));
        }
    }
}



void PackedMarkerGraph::assembleSegmentSequences(const string& name)
{
    createNew(segmentSequences, name + "-sequences");

    for(uint64_t segmentId=0; segmentId!=segments.size(); segmentId++) {

        AssembledSegment assembledSegment;
        assembleMarkerGraphPath(
            0,
            k,
            markers,
            markerGraph,
            segments[segmentId],
            false,
            assembledSegment);

        segmentSequences.appendVector(assembledSegment.rawSequence);
    }
}



// The PackedMarkerGraph::segmentSequences stores, for each segment,
// the entire sequence from the AssembledSegment.
// This includes the entire sequence of the first
// and last vertex of each segment.
// This returns the clipped sequence of each segment,
// which excludes the first and last k/2 bases.
span<const Base> PackedMarkerGraph::clippedSequence(uint64_t segmentId) const
{
    const Base* begin = segmentSequences.begin(segmentId);
    const Base* end = segmentSequences.end(segmentId);
    return span<const Base>(begin + k / 2, end - k / 2);
}



void PackedMarkerGraph::writeGfa(const string& name) const
{
    ofstream gfa(name + ".gfa");
    ofstream csv(name + ".csv");

    // Write the headers.
    gfa << "H\tVN:Z:1.0\n";
    csv << "Segment,Sequence Length,Path Length,First marker graph vertex,Last marker graph vertex\n";

    // Write the segments.
    for(uint64_t segmentId=0; segmentId<segmentSequences.size(); segmentId++) {

        const auto sequence = segmentSequences[segmentId];
        gfa <<"S\t" << segmentId << "\t";
        copy(sequence.begin()+k/2, sequence.end()-k/2, ostream_iterator<Base>(gfa));
        gfa << "\n";

        const auto path = segments[segmentId];
        csv << segmentId << ",";
        csv << sequence.size() << ",";
        csv << path.size() << ",";
        csv << getFirstSegmentVertex(segmentId) << ",";
        csv << getLastSegmentVertex(segmentId) << "\n";
    }

    // Write the links.
    for(const Link& link: links) {
        gfa << "L\t" <<
            link.segmentId0 << "\t+\t" <<
            link.segmentId1 << "\t+\t0M\n";
    }


}



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
