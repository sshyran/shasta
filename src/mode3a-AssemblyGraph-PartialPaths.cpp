// Shasta
#include "mode3a-AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"



void AssemblyGraph::computePartialPaths(
    uint64_t threadCount,
    uint64_t minLinkCoverage)
{
    SHASTA_ASSERT(threadCount > 0);

    // Store the arguments so the thread function can see them.
    computePartialPathsData.minLinkCoverage = minLinkCoverage;

    const uint64_t segmentCount = verticesBySegment.size();
    const uint64_t batchSize = 10;
    setupLoadBalancing(segmentCount, batchSize);
    runThreads(&AssemblyGraph::computePartialPathsThreadFunction, threadCount);
}



void AssemblyGraph::computePartialPathsThreadFunction(uint64_t threadId)
{
    ofstream debugOut;
    debugOut.open("ComputePartialPathsDebug-Thread-" + to_string(threadId));

    const uint64_t minLinkCoverage = computePartialPathsData.minLinkCoverage;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; segmentId++) {

            // Loop over all vertices (replicas) of this segment.
            for(const vertex_descriptor v: verticesBySegment[segmentId]) {
                if(v != null_vertex()) {
                    computePartialPath(v, minLinkCoverage, debugOut);
                }
            }

        }
    }
}



// Compute partial paths (forward and backward) starting from a given vertex.
void AssemblyGraph::computePartialPath(
    vertex_descriptor vStart,
    uint64_t minLinkCoverage,
    ostream& debugOut
    )
{
    AssemblyGraph& assemblyGraph = *this;

    if(debugOut) {
        debugOut << "Following reads at " << vertexStringId(vStart) << "\n";
    }

    // Access the start vertex.
    SHASTA_ASSERT(vStart != null_vertex());
    AssemblyGraphVertex& startVertex = assemblyGraph[vStart];

    // Clear the partial paths.
    vector<vertex_descriptor>& forwardPartialPath = startVertex.forwardPartialPath;
    vector<vertex_descriptor>& backwardPartialPath = startVertex.backwardPartialPath;
    forwardPartialPath.clear();
    backwardPartialPath.clear();

    // The vertices we encounter when following the reads.
    vector<vertex_descriptor> verticesEncountered;

    // The transitions we encounter when following the reads.
    vector< pair<vertex_descriptor, vertex_descriptor> > transitionsEncountered;

    // Loop over JourneyEntry's of the start vertex.
    for(const JourneyEntry& journeyEntry: startVertex.journeyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;

        // Store the vertices encountered in the journey of this read.
        const auto journey = journeys[orientedReadId.getValue()];
        for(uint64_t position=0; position<journey.size(); position++) {
            const vertex_descriptor v = journey[position];
            if(v != null_vertex()) {
                verticesEncountered.push_back(v);
            }
        }

        // Also store the transitions.
        for(uint64_t position=1; position<journey.size(); position++) {
            const vertex_descriptor v0 = journey[position-1];
            const vertex_descriptor v1 = journey[position];
            if(v0 != null_vertex() and v1 != null_vertex()) {
                transitionsEncountered.push_back(make_pair(v0, v1));
            }
        }
    }

    // Count how many times we encountered each vertex.
    vector<uint64_t> vertexFrequency;
    deduplicateAndCount(verticesEncountered, vertexFrequency);

    // Count how many times we encountered each transition.
    // Keep only the ones that appear at least minLinkCoverage times.
    vector<uint64_t> transitionFrequency;
    deduplicateAndCountWithThreshold(
        transitionsEncountered, transitionFrequency, minLinkCoverage);

    if(false) {
        debugOut << "Segments:\n";
        for(uint64_t i=0; i<verticesEncountered.size(); i++) {
            const vertex_descriptor v = verticesEncountered[i];
            debugOut << vertexStringId(v) << " " << vertexFrequency[i] << "\n";
        }

        debugOut << "Transitions:\n";
        for(uint64_t i=0; i<transitionsEncountered.size(); i++) {
            const auto& p = transitionsEncountered[i];
            const vertex_descriptor v0 = p.first;
            const vertex_descriptor v1 = p.second;
            debugOut << vertexStringId(v0) << "->" << vertexStringId(v1) <<
                " " << transitionFrequency[i] << "\n";
        }
    }



    // The transitions we kept define a graph.

    // Starting at the start vertex, follow the linear portion of the graph forward.
    // Stop when we encounter a branch or a vertex seen less than minSegmentCoverage times.
    sort(transitionsEncountered.begin(), transitionsEncountered.end(),  // Not strictly necessary.
        OrderPairsByFirstOnly<vertex_descriptor, vertex_descriptor>());
    vertex_descriptor v = vStart;
    while(true) {
        vector<pair<vertex_descriptor, vertex_descriptor> >::iterator it0, it1;
        tie(it0, it1) = std::equal_range(transitionsEncountered.begin(), transitionsEncountered.end(),
            make_pair(v, null_vertex()),
            OrderPairsByFirstOnly<vertex_descriptor, vertex_descriptor>());
        if( it0 == transitionsEncountered.end() or
            (it1 - it0) != 1) {
            break;
        }
        v = it0->second;
        if(v == vStart) {
            break;
        }
        forwardPartialPath.push_back(v);
    }

    if(debugOut) {
        debugOut << "Forward partial path:";
        for(const vertex_descriptor v: forwardPartialPath) {
            debugOut << " " << vertexStringId(v);
        }
        debugOut << "\n";
    }



    // Starting at the start vertex, follow the linear portion of the graph backward.
    // Stop when we encounter a branch or a vertex seen less than minSegmentCoverage times.
    sort(transitionsEncountered.begin(), transitionsEncountered.end(),
        OrderPairsBySecondOnly<vertex_descriptor, vertex_descriptor>());
    v = vStart;
    while(true) {
        vector<pair<vertex_descriptor, vertex_descriptor> >::iterator it0, it1;
        tie(it0, it1) = std::equal_range(transitionsEncountered.begin(), transitionsEncountered.end(),
            make_pair(null_vertex(), v),
            OrderPairsBySecondOnly<vertex_descriptor, vertex_descriptor>());
        if( it0 == transitionsEncountered.end() or
            (it1 - it0) != 1) {
            break;
        }
        v = it0->first;
        if(v == vStart) {
            break;
        }
        backwardPartialPath.push_back(v);
    }

    if(debugOut) {
        debugOut << "Backward partial path:";
        for(const vertex_descriptor v: backwardPartialPath) {
            debugOut << " " << vertexStringId(v);
        }
        debugOut << "\n";
    }

}



void AssemblyGraph::writePartialPaths() const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv("PartialPaths.csv");
    csv << "Start,Direction,Position,Vertex\n";

    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];

        for(uint64_t position=0; position<vertex0.forwardPartialPath.size(); position++) {
            const vertex_descriptor v1 = vertex0.forwardPartialPath[position];
            csv << vertexStringId(v0) << ",Forward,";
            csv << position << ",";
            csv << vertexStringId(v1) << "\n";
        }

        for(uint64_t position=0; position<vertex0.backwardPartialPath.size(); position++) {
            const vertex_descriptor v1 = vertex0.backwardPartialPath[position];
            csv << vertexStringId(v0) << ",Backward,";
            csv << position << ",";
            csv << vertexStringId(v1) << "\n";
        }
    }
}



void AssemblyGraph::analyzePartialPaths() const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;
    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {

        for(const vertex_descriptor v1: assemblyGraph[v0].forwardPartialPath) {
            forwardPairs.push_back(make_pair(v0, v1));
        }

        for(const vertex_descriptor v1: assemblyGraph[v0].backwardPartialPath) {
            backwardPairs.push_back(make_pair(v1, v0));
        }
    }
    sort(forwardPairs.begin(), forwardPairs.end());
    sort(backwardPairs.begin(), backwardPairs.end());

    // The pairs we found in both directions are the ones we want to keep.
    vector< pair<vertex_descriptor, vertex_descriptor> > bidirectionalPairs;
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs));

    cout << "Partial paths analysis:" << endl;
    cout << "Number of vertices " << num_vertices(assemblyGraph) << endl;
    cout << "Number of forward pairs " << forwardPairs.size() << endl;
    cout << "Number of backward pairs " << backwardPairs.size() << endl;
    cout << "Number of bidirectional pairs " << bidirectionalPairs.size() << endl;



    // Write the bidirectional pairs as a Graphviz graph.
    ofstream dot("PartialPaths.dot");
    dot << "digraph PartialPaths {\n";
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        dot << "\"" << vertexStringId(v) << "\";\n";
    }
    for(const auto& p: bidirectionalPairs) {
        const vertex_descriptor v0 = p.first;
        const vertex_descriptor v1 = p.second;
        dot << "\"" << vertexStringId(v0) << "\"->";
        dot << "\"" << vertexStringId(v1) << "\";\n";
    }


    dot << "}\n";
}
