// Shasta
#include "mode3a-AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "invalid.hpp"
#include "orderPairs.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/reverse_graph.hpp>
#include "dominatorTree.hpp"

// Standard library.
#include "array.hpp"
#include "fstream.hpp"



void AssemblyGraph::computePartialPaths(
    uint64_t threadCount,
    uint64_t minSegmentCoverage,
    uint64_t minLinkCoverage)
{
    SHASTA_ASSERT(threadCount > 0);

    // Store the arguments so the thread function can see them.
    computePartialPathsData.minSegmentCoverage = minSegmentCoverage;
    computePartialPathsData.minLinkCoverage = minLinkCoverage;

    const uint64_t segmentCount = verticesBySegment.size();
    const uint64_t batchSize = 10;
    setupLoadBalancing(segmentCount, batchSize);
    runThreads(&AssemblyGraph::computePartialPathsThreadFunction, threadCount);
}



void AssemblyGraph::computePartialPathsThreadFunction(uint64_t threadId)
{
    ofstream debugOut;
    // debugOut.open("ComputePartialPathsDebug-Thread-" + to_string(threadId));

    const uint64_t minSegmentCoverage = computePartialPathsData.minSegmentCoverage;
    const uint64_t minLinkCoverage = computePartialPathsData.minLinkCoverage;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; segmentId++) {

            // Loop over all vertices (replicas) of this segment.
            for(const vertex_descriptor v: verticesBySegment[segmentId]) {
                if(v != null_vertex()) {
                    computePartialPath2(v, minSegmentCoverage, minLinkCoverage, debugOut);
                }
            }

        }
    }
}



// Compute partial paths (forward and backward) starting from a given vertex.
// Simple version.
void AssemblyGraph::computePartialPath1(
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




// Compute partial paths (forward and backward) starting from a given vertex.
// More robust version that uses dominator trees.
void AssemblyGraph::computePartialPath2(
    vertex_descriptor vStart,
    uint64_t minSegmentCoverage,
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

    // Write the graph.
    if(debugOut) {
        debugOut << "digraph Graph_" << startVertex.segmentId << "_" << startVertex.segmentReplicaIndex << " {\n";
        for(uint64_t i=0; i<verticesEncountered.size(); i++) {
            const vertex_descriptor v = verticesEncountered[i];
            debugOut << "\"" << vertexStringId(v) << "\" [label=\"" <<
                vertexStringId(v) << "\\n" << vertexFrequency[i] <<
                "\"];\n";
        }

        for(uint64_t i=0; i<transitionsEncountered.size(); i++) {
            const auto& p = transitionsEncountered[i];
            const vertex_descriptor v0 = p.first;
            const vertex_descriptor v1 = p.second;
            debugOut <<
                "\"" << vertexStringId(v0) << "\"->\"" <<
                vertexStringId(v1) << "\";\n";
        }
        debugOut << "}\n";
    }



    // The transitions we kept define a graph.
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;
    Graph graph(verticesEncountered.size());
    for(const auto& p: transitionsEncountered) {
        array<vertex_descriptor, 2> v = {p.first, p.second};
        array<uint64_t, 2> iv;
        for(uint64_t k=0; k<2; k++) {
            const auto q = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v[k]);
            SHASTA_ASSERT(q.first != verticesEncountered.end());
            SHASTA_ASSERT(q.second - q.first == 1);
            iv[k] = q.first - verticesEncountered.begin();
        }
        add_edge(iv[0], iv[1], graph);
    }


    // To compute the forward partial path, compute the dominator tree of the graph,
    // with the start vertex as the entrance.
    const auto q = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), vStart);
    SHASTA_ASSERT(q.first != verticesEncountered.end());
    SHASTA_ASSERT(q.second - q.first == 1);
    const uint64_t ivStart = q.first - verticesEncountered.begin();
    std::map<uint64_t, uint64_t> predecessorMap;
    shasta::lengauer_tarjan_dominator_tree(
        graph,
        ivStart,
        boost::make_assoc_property_map(predecessorMap));

    // The allSuccessorMap stores all possible successors in the dominator tree.
    std::map<uint64_t, vector<uint64_t> > allSuccessorMap;
    for(const auto& p: predecessorMap) {
        // In the predecessor map, the key is the target vertex and the value is the source vertex.
        const uint64_t iv0 = p.second;
        const uint64_t iv1 = p.first;
        allSuccessorMap[iv0].push_back(iv1);
    }
    std::map<uint64_t, uint64_t> bestSuccessorMap;
    for(const auto& p: allSuccessorMap) {
        const uint64_t iv0 = p.first;
        uint64_t iv1Best = invalid<uint64_t>;
        uint64_t frequencyBest = 0;
        for(const uint64_t iv1: p.second) {
            const uint64_t frequency = vertexFrequency[iv1];
            if(frequency > frequencyBest) {
                frequencyBest = frequency;
                iv1Best = iv1;
            }
        }
        bestSuccessorMap.insert(make_pair(iv0, iv1Best));
    }

    // Follow the bestSuccessorMap to construct the forward path.
    uint64_t iv = ivStart;
    while(true) {
        auto it = bestSuccessorMap.find(iv);
        if(it == bestSuccessorMap.end()) {
            break;
        }
        iv = it->second;
        if(vertexFrequency[iv] < minSegmentCoverage) {
            break;
        }
        forwardPartialPath.push_back(verticesEncountered[iv]);
    }




    // Write the forward dominator tree.
    if(debugOut) {
        debugOut << "digraph Forward_Tree_" << startVertex.segmentId << "_" << startVertex.segmentReplicaIndex << " {\n";

        // Gather the vertices of the dominator tree.
        vector<vertex_descriptor> dominatorTreeVertices;
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the target vertex and the value is the source vertex.
            const uint64_t iv0 = p.second;
            const uint64_t iv1 = p.first;
            const vertex_descriptor v0 = verticesEncountered[iv0];
            const vertex_descriptor v1 = verticesEncountered[iv1];
            dominatorTreeVertices.push_back(v0);
            dominatorTreeVertices.push_back(v1);
        }
        deduplicate(dominatorTreeVertices);

        for(const vertex_descriptor v: dominatorTreeVertices) {
            const uint64_t i = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v).first -
                verticesEncountered.begin();
            debugOut << "\"" << vertexStringId(v) << "\" [label=\"" <<
                vertexStringId(v) << "\\n" << vertexFrequency[i] <<
                "\"];\n";

        }
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the target vertex and the value is the source vertex.
            const uint64_t iv0 = p.second;
            const uint64_t iv1 = p.first;
            const vertex_descriptor v0 = verticesEncountered[iv0];
            const vertex_descriptor v1 = verticesEncountered[iv1];
            debugOut <<
                "\"" << vertexStringId(v0) << "\"->\"" <<
                vertexStringId(v1) << "\";\n";
        }
        debugOut << "}\n";
    }



    // To compute the backward partial path, compute the backward dominator tree of the graph,
    // with the start vertex as the entrance.
    predecessorMap.clear();
    shasta::lengauer_tarjan_dominator_tree(
        boost::make_reverse_graph(graph),
        ivStart,
        boost::make_assoc_property_map(predecessorMap));



    // The allSuccessorMap stores all possible successors in the dominator tree.
    allSuccessorMap.clear();
    for(const auto& p: predecessorMap) {
        const uint64_t iv0 = p.second;
        const uint64_t iv1 = p.first;
        allSuccessorMap[iv0].push_back(iv1);
    }
    bestSuccessorMap.clear();
    for(const auto& p: allSuccessorMap) {
        const uint64_t iv0 = p.first;
        uint64_t iv1Best = invalid<uint64_t>;
        uint64_t frequencyBest = 0;
        for(const uint64_t iv1: p.second) {
            const uint64_t frequency = vertexFrequency[iv1];
            if(frequency > frequencyBest) {
                frequencyBest = frequency;
                iv1Best = iv1;
            }
        }
        bestSuccessorMap.insert(make_pair(iv0, iv1Best));
    }

    // Follow the bestSuccessorMap to construct the backward path.
    iv = ivStart;
    while(true) {
        auto it = bestSuccessorMap.find(iv);
        if(it == bestSuccessorMap.end()) {
            break;
        }
        iv = it->second;
        if(vertexFrequency[iv] < minSegmentCoverage) {
            break;
        }
        backwardPartialPath.push_back(verticesEncountered[iv]);
    }



    // Write the backward dominator tree.
    if(debugOut) {
        debugOut << "digraph Backward_Tree_" << startVertex.segmentId << "_" << startVertex.segmentReplicaIndex << " {\n";

        // Gather the vertices of the dominator tree.
        vector<vertex_descriptor> dominatorTreeVertices;
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the source vertex and the value is the target vertex.
            const uint64_t iv0 = p.first;
            const uint64_t iv1 = p.second;
            const vertex_descriptor v0 = verticesEncountered[iv0];
            const vertex_descriptor v1 = verticesEncountered[iv1];
            dominatorTreeVertices.push_back(v0);
            dominatorTreeVertices.push_back(v1);
        }
        deduplicate(dominatorTreeVertices);

        for(const vertex_descriptor v: dominatorTreeVertices) {
            const uint64_t i = std::equal_range(verticesEncountered.begin(), verticesEncountered.end(), v).first -
                verticesEncountered.begin();
            debugOut << "\"" << vertexStringId(v) << "\" [label=\"" <<
                vertexStringId(v) << "\\n" << vertexFrequency[i] <<
                "\"];\n";

        }
        for(const auto& p: predecessorMap) {
            // In the predecessor map, the key is the target vertex and the value is the source vertex.
            const uint64_t iv0 = p.first;
            const uint64_t iv1 = p.second;
            const vertex_descriptor v0 = verticesEncountered[iv0];
            const vertex_descriptor v1 = verticesEncountered[iv1];
            debugOut <<
                "\"" << vertexStringId(v0) << "\"->\"" <<
                vertexStringId(v1) << "\";\n";
        }
        debugOut << "}\n";
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
