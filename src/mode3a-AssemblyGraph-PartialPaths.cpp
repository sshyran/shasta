// Shasta
#include "mode3a-AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "invalid.hpp"
#include "orderPairs.hpp"
#include "transitiveReduction.hpp"
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
    uint64_t segmentCoverageThreshold1,
    uint64_t segmentCoverageThreshold2,
    uint64_t minLinkCoverage)
{
    SHASTA_ASSERT(threadCount > 0);

    // Store the arguments so the thread function can see them.
    computePartialPathsData.segmentCoverageThreshold1 = segmentCoverageThreshold1;
    computePartialPathsData.segmentCoverageThreshold2 = segmentCoverageThreshold2;
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

    const uint64_t segmentCoverageThreshold1 = computePartialPathsData.segmentCoverageThreshold1;
    const uint64_t segmentCoverageThreshold2 = computePartialPathsData.segmentCoverageThreshold2;
    const uint64_t minLinkCoverage = computePartialPathsData.minLinkCoverage;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; segmentId++) {

            // Loop over all vertices (replicas) of this segment.
            for(const vertex_descriptor v: verticesBySegment[segmentId]) {
                if(v != null_vertex()) {
                    computePartialPath2(v,
                        segmentCoverageThreshold1, segmentCoverageThreshold2, minLinkCoverage, debugOut);
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
    uint64_t segmentCoverageThreshold1,
    uint64_t segmentCoverageThreshold2,
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

    // Explicitly construct the forward dominator tree.
    Graph forwardTree(verticesEncountered.size());
    for(const auto& p: predecessorMap) {
        const uint64_t iv0 = p.second;
        const uint64_t iv1 = p.first;
        add_edge(iv0, iv1, forwardTree);
    }



    // To compute the forward partial path, follow the forward dominator tree.
    uint64_t iv = ivStart;
    while(true) {

        // Find the out-vertices and sort them by decreasing vertex frequency.
        vector< pair<uint64_t, uint64_t> > outVertices;
        BGL_FORALL_OUTEDGES(iv, e, forwardTree, Graph) {
            const uint64_t iv1 = target(e, forwardTree);
            outVertices.push_back(make_pair(iv1, vertexFrequency[iv1]));
        }
        sort(outVertices.begin(), outVertices.end(), OrderPairsBySecondOnlyGreater<uint64_t, uint64_t >());

        // If there are no out-vertices, the forward path ends here.
        if(outVertices.empty()) {
            break;
        }

        // If the strongest out-vertex is too weak, the forward path ends here.
        if(outVertices.front().second < segmentCoverageThreshold1) {
            break;
        }

        // If the strongest in-vertex loses too much coverage compared to iv, the backward path ends here.
        const uint64_t coverageLoss =
            (outVertices.front().second >= vertexFrequency[iv]) ? 0 :
            (vertexFrequency[iv] - outVertices.front().second);
        if(coverageLoss > segmentCoverageThreshold2) {
            break;
        }

        // In all other cases, we add the strongest out-vertex to the forward path.
        iv = outVertices.front().first;
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

    // Explicitly construct the backward dominator tree.
    Graph backwardTree(verticesEncountered.size());
    for(const auto& p: predecessorMap) {
        const uint64_t iv0 = p.first;
        const uint64_t iv1 = p.second;
        add_edge(iv0, iv1, backwardTree);
    }



    // To compute the backward partial path, follow the backward dominator tree.
    iv = ivStart;
    while(true) {

        // Find the in-vertices and sort them by decreasing vertex frequency.
        vector< pair<uint64_t, uint64_t> > inVertices;
        BGL_FORALL_INEDGES(iv, e, backwardTree, Graph) {
            const uint64_t iv1 = source(e, forwardTree);
            inVertices.push_back(make_pair(iv1, vertexFrequency[iv1]));
        }
        sort(inVertices.begin(), inVertices.end(), OrderPairsBySecondOnlyGreater<uint64_t, uint64_t >());

        // If there are no in-vertices, the backward path ends here.
        if(inVertices.empty()) {
            break;
        }

        // If the strongest in-vertex is too weak, the backward path ends here.
        if(inVertices.front().second < segmentCoverageThreshold1) {
            break;
        }

        // If the strongest in-vertex loses too much coverage compared to iv, the backward path ends here.
        const uint64_t coverageLoss =
            (inVertices.front().second >= vertexFrequency[iv]) ? 0 :
            (vertexFrequency[iv] - inVertices.front().second);
        if(coverageLoss > segmentCoverageThreshold2) {
            break;
        }

        // In all other cases, we add the strongest in-vertex to the backward path.
        iv = inVertices.front().first;
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

    // One line for each partial path entry.
    ofstream csv1("PartialPaths1.csv");
    csv1 << "Start,Direction,Position,Vertex\n";

    // One line for each partial path.
    ofstream csv2("PartialPaths2.csv");
    csv1 << "Start,Direction,Vertices\n";

    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];

        csv2 << vertexStringId(v0) << ",Forward,";
        for(uint64_t position=0; position<vertex0.forwardPartialPath.size(); position++) {
            const vertex_descriptor v1 = vertex0.forwardPartialPath[position];
            csv1 << vertexStringId(v0) << ",Forward,";
            csv1 << position << ",";
            csv1 << vertexStringId(v1) << "\n";
            csv2 << vertexStringId(v1) << ",";
        }
        csv2 << "\n";

        csv2 << vertexStringId(v0) << ",Backward,";
        for(uint64_t position=0; position<vertex0.backwardPartialPath.size(); position++) {
            const vertex_descriptor v1 = vertex0.backwardPartialPath[position];
            csv1 << vertexStringId(v0) << ",Backward,";
            csv1 << position << ",";
            csv1 << vertexStringId(v1) << "\n";
            csv2 << vertexStringId(v1) << ",";
        }
        csv2 << "\n";
    }
}



void AssemblyGraph::analyzePartialPaths() const
{
    // EXPOSE WHEN CODE STABILIZES
    const uint64_t n = 10;

    const AssemblyGraph& assemblyGraph = *this;

    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;
    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {

        for(uint64_t i=0; i<min(n, assemblyGraph[v0].forwardPartialPath.size()); i++) {
            forwardPairs.push_back(make_pair(v0, assemblyGraph[v0].forwardPartialPath[i]));
        }

        for(uint64_t i=0; i<min(n, assemblyGraph[v0].backwardPartialPath.size()); i++) {
            backwardPairs.push_back(make_pair(assemblyGraph[v0].backwardPartialPath[i], v0));
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


#if 0
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
#endif


    // The bidirectional pairs define a graph.
    // Create it and compute its transitive reduction.
    // This will fail if the graph has cycles.
    vector<vertex_descriptor> vertexTable;
    std::map<vertex_descriptor, uint64_t> vertexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexTable.push_back(v);
        vertexMap.insert(make_pair(v, vertexIndex++));
    }
    using Graph = boost::adjacency_list<boost::listS, boost::vecS, boost::directedS>;
    Graph graph(vertexTable.size());
    for(const auto& p: bidirectionalPairs) {
        const vertex_descriptor v0 = p.first;
        const vertex_descriptor v1 = p.second;
        add_edge(vertexMap[v0], vertexMap[v1], graph);
    }
    // The transitive reduction could be computed in parallel over connected components.
    transitiveReduction(graph);
    cout << "Number of bidirectional pairs after transitive reduction " <<
        num_edges(graph) << endl;

    // Write out the transitive reduction as a Graphviz graph.
    ofstream dot("PartialPaths.dot");
    dot << "digraph PartialPaths {\n";
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = vertexTable[source(e, graph)];
        const vertex_descriptor v1 = vertexTable[target(e, graph)];
        dot << "\"" << vertexStringId(v0) << "\"->";
        dot << "\"" << vertexStringId(v1) << "\";\n";
    }
    dot << "}\n";
}
