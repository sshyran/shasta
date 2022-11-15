// Shasta.
#include "mode3a.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "AssembledSegment.hpp"
#include "MarkerGraph.hpp"
#include "mode3a-BubbleCleaner.hpp"
#include "Reads.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3a;

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

