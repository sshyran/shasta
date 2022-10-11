#include "mode3-JaccardGraph.hpp"
#include "deduplicate.hpp"
#include "mode3.hpp"
#include "orderPairs.hpp"
#include "orderVectors.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/topological_sort.hpp>

// Standard library.
#include "fstream.hpp"



// Create a JaccardGraph with the given number of vertices
// (one for each segment) and no edges.
JaccardGraph::JaccardGraph(uint64_t segmentCount)
{
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        vertexTable.push_back(add_vertex(JaccardGraphVertex(segmentId), *this));
    }
}


void AssemblyGraph::createJaccardGraph(
    size_t threadCount
    )
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minComponentSize = 100; //Likely needs to be decreased. Keep high for debugging.

    cout << timestamp << "createJaccardGraph begins." << endl;

    // Create the JaccardGraph and its vertices.
    const uint64_t segmentCount = markerGraphPaths.size();
    cout << "The total number of segments in the assembly graph is " << segmentCount << endl;
    jaccardGraphPointer = make_shared<JaccardGraph>(segmentCount);
    JaccardGraph& jaccardGraph = *jaccardGraphPointer;

    // Compute edges, in parallel.
    jaccardGraph.threadEdges.resize(threadCount);
    const uint64_t batchSize = 100;
    setupLoadBalancing(segmentCount, batchSize);
    runThreads(&AssemblyGraph::createJaccardGraphThreadFunction, threadCount);
    jaccardGraph.storeEdges();
    jaccardGraph.writeGraphviz("JaccardGraph0.dot", false, false);
    jaccardGraph.writeGraphviz("JaccardGraph0-Labeled.dot", false, true);
    jaccardGraph.writeEdgesCsv("JaccardGraph0Edges.csv");
    cout << "The initial Jaccard graph has " << num_vertices(jaccardGraph) <<
        " vertices (segments) and " << num_edges(jaccardGraph) << " edges." << endl;

    // Clear all weak vertices.
    jaccardGraph.clearWeakVertices();
    cout << "After clearing weak vertices, the Jaccard graph has " << num_vertices(jaccardGraph) <<
        " vertices (segments) and " << num_edges(jaccardGraph) << " edges." << endl;
    jaccardGraph.writeGraphviz("JaccardGraph1.dot", false, false);
    jaccardGraph.writeGraphviz("JaccardGraph1-Labeled.dot", false, true);
    jaccardGraph.writeEdgesCsv("JaccardGraph1Edges.csv");

    // Compute all connected components of size at least minComponentSize.
    jaccardGraph.computeConnectedComponents(minComponentSize);

    // Store the cluster id of each segment.
    // Each connected component of the Jaccard graph with sufficient size
    // generates a cluster.
    createNew(clusterIds, "Mode3-ClusterIds");
    jaccardGraph.findClusters(clusterIds);

    // Compute assembly paths.
    jaccardGraph.computeAssemblyPaths();

    // Create the ExpandedJaccardGraph.
    ExpandedJaccardGraph expandedJaccardGraph(jaccardGraph);
    expandedJaccardGraph.writeGraphviz("ExpandedJaccardGraph0.dot");
    expandedJaccardGraph.merge();
    expandedJaccardGraph.writeGraphviz("ExpandedJaccardGraph1.dot");

    cout << timestamp << "createJaccardGraph ends." << endl;
}



void AssemblyGraph::createJaccardGraphThreadFunction(size_t threadId)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; ++segmentId) {
            createJaccardGraphEdges(segmentId, jaccardGraphPointer->threadEdges[threadId]);
        }
    }
}



void AssemblyGraph::createJaccardGraphEdges(
    uint64_t segmentId,
    vector<JaccardGraphEdgeInfo>& edges)
{
    for(uint64_t direction=0; direction<2; direction++) {
        createJaccardGraphEdges(segmentId, direction, edges);
    }
}



// This follows an algorithm similar to the one used by createAssemblyPath3.
void AssemblyGraph::createJaccardGraphEdges(
    uint64_t primarySegmentId,
    uint64_t direction,
    vector<JaccardGraphEdgeInfo>& edges)
{
    // EXPOSE WHEN CODE STABILIZES.
    // FOR NOW THESE SHOULD BE THE SAME AS IN AssemblyGraph::createAssemblyPath3.
    const uint64_t minCommonForLink = 3;
    const uint64_t minCommonForPrimary = 3;
    const double minJaccard = 0.75;
    const int32_t minLinkSeparation = -20;

    // We start from primarySegmentId
    // and move in the specified direction until we find segmentId1 with
    // sufficiently high Jaccard similarity and number of
    // common oriented reads with primarySegmentId.
    // At each step, we choose the link that has the most common oriented
    // reads with the primarySegmentId.
    SegmentOrientedReadInformation infoPrimary;
    getOrientedReadsOnSegment(primarySegmentId, infoPrimary);
    JaccardGraphEdgeInfo edge;
    edge.direction = direction;
    uint64_t segmentId0 = primarySegmentId;
    std::set<uint64_t> previousSegments;
    while(true) {

        // Loop over outgoing or incoming links of segmentId0.
        // Find the link with the most common reads with the primarySegmentId.
        const auto linkIds = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
        if(linkIds.empty()) {
            return;
        }
        uint64_t linkIdBest = invalid<uint64_t>;
        uint64_t commonOrientedReadCountBest = 0;
        for(const uint64_t linkId: linkIds) {

            // If link separation is too negative, skip it.
            // The goal here is to avoid cycles in paths.
            const Link& link = links[linkId];
            if(link.separation < minLinkSeparation) {
                continue;
            }

            // Count the number of common oriented reads between the reference segment and this link.
            uint64_t commonOrientedReadCount;
            analyzeSegmentLinkPair(primarySegmentId, linkId, commonOrientedReadCount);

            // If better than the one we have it, record it.
            if(commonOrientedReadCount > commonOrientedReadCountBest) {
                linkIdBest = linkId;
                commonOrientedReadCountBest = commonOrientedReadCount;
            }
        }
        if(commonOrientedReadCountBest < minCommonForLink) {
            return;
        }
        const uint64_t linkId = linkIdBest;

        // Get the segment at the other side of this link.
        const Link& link = links[linkId];
        const uint64_t segmentId1 = (direction==0) ? link.segmentId1 : link.segmentId0;

        // Check that we haven't been here before.
        if(previousSegments.contains(segmentId1)) {
            break;
        }
        previousSegments.insert(segmentId1);

        // Check segmentId1 against the primary segment.
        SegmentOrientedReadInformation info1;
        getOrientedReadsOnSegment(segmentId1, info1);
        if(direction == 0) {
            analyzeSegmentPair(
                    primarySegmentId, segmentId1,
                    infoPrimary, info1,
                    markers, edge.segmentPairInformation);
        } else {
            analyzeSegmentPair(
                segmentId1, primarySegmentId,
                info1, infoPrimary,
                markers, edge.segmentPairInformation);
        }

        // If the Jaccard similarity is high, we found the Jaccard graph edge
        // we were looking for.
        if( edge.segmentPairInformation.commonCount >= minCommonForPrimary and
            edge.segmentPairInformation.rawJaccard() >= minJaccard) {   // ****** USING RAWJACCARD INSTEAD OF JACCARD
            if(direction == 0) {
                edge.segmentId0 = primarySegmentId;
                edge.segmentId1 = segmentId1;
            } else {
                edge.segmentId0 = segmentId1;
                edge.segmentId1 = primarySegmentId;
                reverse(edge.segmentIds.begin(), edge.segmentIds.end());
            }
            edges.push_back(edge);
            return;
        }

        edge.segmentIds.push_back(segmentId1);
        segmentId0 = segmentId1;
    }
}



// This storesin the Jaccard graph the edges found by all threads.
void JaccardGraph::storeEdges()
{
    JaccardGraph& jaccardGraph = *this;

    for(const auto& threadEdges: threadEdges) {
        for(const JaccardGraphEdgeInfo& info: threadEdges) {

            const uint64_t segmentId0 = info.segmentId0;
            const uint64_t segmentId1 = info.segmentId1;
            const JaccardGraph::vertex_descriptor v0 = vertexTable[segmentId0];
            const JaccardGraph::vertex_descriptor v1 = vertexTable[segmentId1];

            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(v0, v1, jaccardGraph);
            if(not edgeExists) {
                boost::add_edge(v0, v1,
                    JaccardGraphEdge(info.segmentPairInformation, info.direction, info.segmentIds),
                    jaccardGraph);
            } else {
                jaccardGraph[e].wasFoundInDirection[info.direction] = true;
            }
        }
    }
    threadEdges.clear();
}



// A strong vertex is one that is incident to at least one strong edge.
bool JaccardGraph::isStrongVertex(vertex_descriptor v) const
{
    const JaccardGraph& jaccardGraph = *this;

    // Check the out-edges.
    BGL_FORALL_OUTEDGES(v, e, jaccardGraph, JaccardGraph) {
        if(jaccardGraph[e].isStrong()) {
            return true;
        }
    }

    // Check the in-edges.
    BGL_FORALL_INEDGES(v, e, jaccardGraph, JaccardGraph) {
        if(jaccardGraph[e].isStrong()) {
            return true;
        }
    }

    // We did not find any strong edges.
    return false;
}




// Remove all weak vertices.
void JaccardGraph::removeWeakVertices()
{
    JaccardGraph& jaccardGraph = *this;

    // Find the vertices we are going to remove.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        if(not isStrongVertex(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }

    // Remove the vertices we flagged.
    for(const vertex_descriptor v: verticesToBeRemoved) {
        removeVertex(v);
    }

}



// Remove all edges to/from weak vertices.
void JaccardGraph::clearWeakVertices()
{
    JaccardGraph& jaccardGraph = *this;

    vector<vertex_descriptor> verticesToBeCleared;
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        if(not isStrongVertex(v)) {
            verticesToBeCleared.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeCleared) {
        clear_vertex(v, jaccardGraph);
    }

}



// Remove a vertex, making sure to update the vertexTable.
void JaccardGraph::removeVertex(vertex_descriptor v)
{
    JaccardGraph& jaccardGraph = *this;
    const uint64_t segmentId = jaccardGraph[v].segmentId;
    vertexTable[segmentId] = null_vertex();
    clear_vertex(v, jaccardGraph);
    remove_vertex(v, jaccardGraph);
}



void JaccardGraph::writeGraphviz(
    const string& fileName,
    bool includeIsolatedVertices,
    bool writeLabels) const
{
    ofstream file(fileName);
    writeGraphviz(file, includeIsolatedVertices, writeLabels);
}



void JaccardGraph::writeGraphviz(
    ostream& graphOut,
    bool includeIsolatedVertices,
    bool writeLabels) const
{
    const JaccardGraph& jaccardGraph = *this;

    graphOut << "digraph JaccardGraph {" << endl;

    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        if( includeIsolatedVertices or
            in_degree(v, jaccardGraph) or
            out_degree(v, jaccardGraph)) {
            graphOut << jaccardGraph[v].segmentId;
            if(writeLabels) {
                graphOut << " [label=" << jaccardGraph[v].segmentId << "]";
            }
            graphOut << ";\n";
        }
    }

    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const JaccardGraphEdge& edge = jaccardGraph[e];
        const JaccardGraph::vertex_descriptor v0 = source(e, jaccardGraph);
        const JaccardGraph::vertex_descriptor v1 = target(e, jaccardGraph);
        const uint64_t segmentId0 = jaccardGraph[v0].segmentId;
        const uint64_t segmentId1 = jaccardGraph[v1].segmentId;

        graphOut << segmentId0 << "->" << segmentId1 << "[";

        // Color the edge based on the direction flags.
        if(edge.wasFoundInDirection[0]) {
            if(edge.wasFoundInDirection[1]) {
                // Found in both directions.
                graphOut << " color=black";
            } else {
                // Only found in the forward direction.
                graphOut << " color=red";
            }
        } else {
            if(edge.wasFoundInDirection[1]) {
                // Only found in the backward direction.
                graphOut << " color=green";
            } else {
                SHASTA_ASSERT(0);
            }
        }

        if(writeLabels) {
            graphOut << " label=\"";
            for(const uint64_t segmentId: edge.segmentIds) {
                graphOut << segmentId << "\\n";
            }
            graphOut << "\"";
        }
        graphOut << "];\n";
    }

    graphOut << "}" << endl;

}



// Write edges in csv format.
void JaccardGraph::writeEdgesCsv(const string& fileName) const
{
    ofstream file(fileName);
    writeEdgesCsv(file);
}
void JaccardGraph::writeEdgesCsv(ostream& csv) const
{
    const JaccardGraph& jaccardGraph = *this;

    csv << "SegmentId0,SegmentId1,FoundForward,FoundBackward,SegmentId\n";
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const JaccardGraphEdge& edge = jaccardGraph[e];
        const JaccardGraph::vertex_descriptor v0 = source(e, jaccardGraph);
        const JaccardGraph::vertex_descriptor v1 = target(e, jaccardGraph);
        const uint64_t segmentId0 = jaccardGraph[v0].segmentId;
        const uint64_t segmentId1 = jaccardGraph[v1].segmentId;

        for(const uint64_t segmentId: edge.segmentIds) {
            csv << segmentId0 << ",";
            csv << segmentId1 << ",";
            csv << int(edge.wasFoundInDirection[0]) << ",";
            csv << int(edge.wasFoundInDirection[1]) << ",";
            csv << segmentId << "\n";
        }
    }
}



// Compute all connected components of size at least minComponentSize.
// They are stored in order of decreasing size.
void JaccardGraph::computeConnectedComponents(uint64_t minComponentSize)
{
    const JaccardGraph& jaccardGraph = *this;

    // This must be called without removing any vertices.
    const uint64_t segmentCount = num_vertices(jaccardGraph);

    // Compute connected components.
    vector<uint64_t> rank(segmentCount);
    vector<uint64_t> parent(segmentCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        disjointSets.make_set(segmentId);
    }
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const JaccardGraph::vertex_descriptor v0 = source(e, jaccardGraph);
        const JaccardGraph::vertex_descriptor v1 = target(e, jaccardGraph);
        const uint64_t segmentId0 = jaccardGraph[v0].segmentId;
        const uint64_t segmentId1 = jaccardGraph[v1].segmentId;
        disjointSets.union_set(segmentId0, segmentId1);
    }

    // Gather the segments in each connected component.
    vector< vector<uint64_t> > allComponents(segmentCount);
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        const uint64_t componentId = disjointSets.find_set(segmentId);
        allComponents[componentId].push_back(segmentId);
    }

    // Create a table of the components of size at least minComponentSize,
    // sorted by decreasing size.
    vector< pair<uint64_t, uint64_t> > componentTable; // pair(componentId, componentSize)
    for(uint64_t componentId=0; componentId<segmentCount; componentId++) {
        const uint64_t componentSize = allComponents[componentId].size();
        if(componentSize >= minComponentSize) {
            componentTable.push_back(make_pair(componentId, componentSize));
        }
    }
    sort(componentTable.begin(), componentTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Store the connected components of size at least minComponentSize.
    components.clear();
    for(uint64_t newComponentId=0; newComponentId<componentTable.size(); newComponentId++) {
        const auto& p = componentTable[newComponentId];
        const uint64_t oldComponentId = p.first;
        const uint64_t componentSize = p.second;
        const vector<uint64_t>& component = allComponents[oldComponentId];
        SHASTA_ASSERT(component.size() == componentSize);
        components.push_back(component);
    }


    // Write a histogram of component sizes.
    vector<uint64_t> histogram;
    for(const auto& p: componentTable) {
        const uint64_t componentSize = p.second;
        if(componentSize >= histogram.size()) {
            histogram.resize(componentSize + 1, 0);
        }
        ++histogram[componentSize];
    }
    ofstream csv("JaccardGraphComponentSizeHistogram.csv");
    csv << "Size,Frequency,Vertices,\n";
    for(uint64_t componentSize=1; componentSize<histogram.size(); componentSize++) {
        const uint64_t frequency = histogram[componentSize];
        if(frequency > 0) {
            csv << componentSize << ",";
            csv << frequency << ",";
            csv << frequency * componentSize << ",";
            csv << "\n";
        }
    }

}



// Compute connected component and store the component
// (define as a cluster) that each segment belongs to.
void JaccardGraph::findClusters(
    MemoryMapped::Vector<uint64_t>& clusterIds)
{
    const JaccardGraph& jaccardGraph = *this;

    // This must be called without removing any vertices.
    const uint64_t segmentCount = num_vertices(jaccardGraph);

    clusterIds.resize(segmentCount);
    fill(clusterIds.begin(), clusterIds.end(), invalid<uint64_t>);
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<uint64_t>& component = components[componentId];
        for(const uint64_t segmentId: component) {
            clusterIds[segmentId] = componentId;
        }
    }

}



// Construction of the ExpandedJaccardGraph.
// Each vertex of the JaccardGraph generates a vertex in the ExpandedJaccardGraph.
// Each edge of the JaccardGraph generates a linear chain of vertices
// in the ExpandedJaccardGraph.
ExpandedJaccardGraph::ExpandedJaccardGraph(const JaccardGraph& jaccardGraph)
{
    using Graph = ExpandedJaccardGraph;
    Graph& graph = *this;

    // Generate the vertices.
    std::map<JaccardGraph::vertex_descriptor, Graph::vertex_descriptor> vertexMap;
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        const Graph::vertex_descriptor u = add_vertex(
            ExpandedJaccardGraphVertex(jaccardGraph[v].segmentId, true), graph);
        vertexMap.insert(make_pair(v, u));
    }



    // Each edge of the JaccardGraph generates a linear chain of vertices
    // in the ExpandedJaccardGraph.
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const JaccardGraph::vertex_descriptor v0 = source(e, jaccardGraph);
        const JaccardGraph::vertex_descriptor v1 = target(e, jaccardGraph);
        const Graph::vertex_descriptor u0 = vertexMap[v0];
        const Graph::vertex_descriptor u1 = vertexMap[v1];
        const vector<uint64_t>& segmentIds = jaccardGraph[e].segmentIds;

        Graph::vertex_descriptor u = u0;
        for(const uint64_t segmentId: segmentIds) {
            const Graph::vertex_descriptor w = add_vertex(
                ExpandedJaccardGraphVertex(segmentId, false), graph);
            add_edge(u, w, graph);
            u = w;
        }
        add_edge(u, u1, graph);
    }
}



void ExpandedJaccardGraph::writeGraphviz(const string& fileName) const
{
    ofstream s(fileName);
    writeGraphviz(s);
}
void ExpandedJaccardGraph::writeGraphviz(ostream& s) const
{
    using Graph = ExpandedJaccardGraph;
    const Graph& graph = *this;

    const bool debug = false;

    s << "digraph ExpandedJaccardGraph {" << endl;

    // We can't use the segment ids to identify vertices
    // because each segment id can appear multiple times.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const ExpandedJaccardGraphVertex& vertex = graph[v];
        const double primaryFraction = vertex.primaryFraction();
        s << "\"" << v << "\" [label=\"" << vertex.segmentId;
        if(debug) {
            s << "\\n" << v;
        }
        s << "\\n" << vertex.primaryCount << "/" << vertex.totalCount << "\"";
        const double H = primaryFraction / 3.;
        const double S = 0.5;
        const double V = 1.;
        s << " style=filled fillcolor=\"" << H << " " << " " << S << " "<< V << "\"";
        s << "];\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const Graph::vertex_descriptor v0 = source(e, graph);
        const Graph::vertex_descriptor v1 = target(e, graph);

        s << "\"" << v0 << "\"->\"" << v1 << "\";\n";
    }

    s << "}" << endl;

}



// Recursively merge pairs of vertices that have a common parent or child
// and that refer to the same segmentId.
void ExpandedJaccardGraph::merge()
{
    using Graph = ExpandedJaccardGraph;
    Graph& graph = *this;

    const bool debug = false;
    if(debug) {
        cout << "ExpandedJaccardGraph::merge begins." << endl;
    }

    std::set<Branch> branches;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(out_degree(v, graph) > 1) {
            branches.insert(make_pair(v, 0));
        }
        if(in_degree(v, graph) > 1) {
            branches.insert(make_pair(v, 1));
        }
    }



    // Recursive merge.
    vector<vertex_descriptor> neighbors;
    while(not branches.empty()) {
        const auto it = branches.begin();
        const vertex_descriptor v0 = it->first;
        const uint64_t direction = it->second;
        branches.erase(it);

        if(debug) {
            cout << "Working on branch " << v0 << " " << direction << endl;
        }

        // Gather the children or parents.
        neighbors.clear();
        if(direction == 0) {
            BGL_FORALL_OUTEDGES(v0, e, graph, Graph) {
                neighbors.push_back(target(e, graph));
            }
        } else if(direction == 1) {
            BGL_FORALL_INEDGES(v0, e, graph, Graph) {
                neighbors.push_back(source(e, graph));
            }

        } else {
            SHASTA_ASSERT(0);
        }
        if(debug) {
            cout << neighbors.size() << " neighbors:";
            for(const vertex_descriptor v: neighbors) {
                cout << " " << v;
            }
            cout << endl;
        }
        SHASTA_ASSERT(neighbors.size() > 1);

        // Find a pair of neighbors with the same segmentId.
        vertex_descriptor v1, v2;
        bool found = false;
        for(uint64_t i1=0; i1<neighbors.size()-1; i1++) {
            v1 = neighbors[i1];
            for(uint64_t i2=i1+1; i2<neighbors.size(); i2++) {
                v2 = neighbors[i2];
                if(graph[v1].segmentId == graph[v2].segmentId) {
                    found = true;
                    break;
                }
            }
            if(found) {
                break;
            }
        }

        // If we did not find a pair of neighbors with the same segmentId,
        // there is nothing to do. We already removed this branch, so we
        // are done.
        if(not found) {
            if(debug) {
                cout << "No pair can be merged for this branch." << endl;
            }
            continue;
        }
        if(debug) {
            cout << "Merging " << v1 << " " << v2 << endl;
        }

        // Merge v1 and v2, and update the branches.
        merge(v1, v2, branches, debug);

    }

    if(debug) {
        cout << "ExpandedJaccardGraph::merge ends." << endl;
    }
}



// Merge v1 and v2 while updating the set of branches.
void ExpandedJaccardGraph::merge(
    vertex_descriptor v1,
    vertex_descriptor v2,
    std::set<Branch>& branches,
    bool debug)
{
    using Graph = ExpandedJaccardGraph;
    Graph& graph = *this;

    const ExpandedJaccardGraphVertex& vertex1 = graph[v1];
    const ExpandedJaccardGraphVertex& vertex2 = graph[v2];

    // Check the segmentId.
    const uint64_t segmentId = vertex1.segmentId;
    SHASTA_ASSERT(segmentId == vertex2.segmentId);

    // Find the children of v1 and v2.
    // These will be the children of the merged vertex v3.
    vector<vertex_descriptor> children;
    BGL_FORALL_OUTEDGES(v1, e, graph, Graph) {
        children.push_back(target(e, graph));
    }
    BGL_FORALL_OUTEDGES(v2, e, graph, Graph) {
        children.push_back(target(e, graph));
    }
    deduplicate(children);

    // Find the parents of v1 and v2.
    // These will be the parents of the merged vertex v3.
    vector<vertex_descriptor> parents;
    BGL_FORALL_INEDGES(v1, e, graph, Graph) {
        parents.push_back(source(e, graph));
    }
    BGL_FORALL_INEDGES(v2, e, graph, Graph) {
        parents.push_back(source(e, graph));
    }
    deduplicate(parents);

    if(debug) {
        cout << "Merging " << v1 << " " << v2 << endl;
        cout << "Children:";
        for(const vertex_descriptor v: children) {
            cout << " " << v;
        }
        cout << endl;
        cout << "Parents:";
        for(const vertex_descriptor v: parents) {
            cout << " " << v;
        }
        cout << endl;
    }

    // Remove the branches that will be affected by the merge.
    // We will add branches back as necessary.
    for(const vertex_descriptor v: children) {
        branches.erase(make_pair(v, 1));
    }
    for(const vertex_descriptor v: parents) {
        branches.erase(make_pair(v, 0));
    }
    branches.erase(make_pair(v1, 0));
    branches.erase(make_pair(v1, 1));
    branches.erase(make_pair(v2, 0));
    branches.erase(make_pair(v2, 1));

    // Create the merged vertex.
    ExpandedJaccardGraphVertex vertex3;
    vertex3.segmentId = segmentId;
    vertex3.totalCount = vertex1.totalCount + vertex2.totalCount;
    vertex3.primaryCount = vertex1.primaryCount + vertex2.primaryCount;
    const vertex_descriptor v3 = add_vertex(vertex3, graph);
    if(debug) {
        cout << "Created merged vertex " << v3 << endl;
    }

    // Remove the vertices that were merged, v1 and v2.
    clear_vertex(v1, graph);
    clear_vertex(v2, graph);
    remove_vertex(v1, graph);
    remove_vertex(v2, graph);
    if(debug) {
        cout << "Removed the merged vertices " << v1 << " " << v2 << endl;
    }

    // Add the edges to/from the merged vertex.
    for(const vertex_descriptor v: children) {
        add_edge(v3, v, graph);
        if(debug) {
            cout << "Added edge " << v3 << " " << v << endl;
        }
    }
    for(const vertex_descriptor v: parents) {
        add_edge(v, v3, graph);
        if(debug) {
            cout << "Added edge " << v << " " << v3 << endl;
        }
    }

    // Add back any necessary branches.
    if(out_degree(v3, graph) > 1) {
        branches.insert(make_pair(v3, 0));
        if(debug) {
            cout << "Added branch " << v3 << " " << 0 << endl;
        }
    }
    if(in_degree(v3, graph) > 1) {
        branches.insert(make_pair(v3, 1));
        if(debug) {
            cout << "Added branch " << v3 << " " << 1 << endl;
        }
    }
    for(const vertex_descriptor v: children) {
        if(in_degree(v, graph) > 1) {
            branches.insert(make_pair(v, 1));
            if(debug) {
                cout << "Added branch " << v << " " << 1 << endl;
            }
        }
    }
    for(const vertex_descriptor v: parents) {
        if(out_degree(v, graph) > 1) {
            branches.insert(make_pair(v, 0));
            if(debug) {
                cout << "Added branch " << v << " " << 0 << endl;
            }
        }
    }
}



// Compute assembly paths.
void JaccardGraph::computeAssemblyPaths()
{
    assemblyPaths.clear();
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        computeAssemblyPaths(componentId);
    }
}
void JaccardGraph::computeAssemblyPaths(uint64_t componentId)
{
    const JaccardGraph& jaccardGraph = *this;

    const bool debug = true;
    const vector<uint64_t>& component = components[componentId];
    if(debug) {
        cout << "Computing assembly paths for component " << componentId <<
            " of size " << component.size() << endl;
    }

    // Create a Graph to represent just this component.
    // Each vertex of the Graph stores the corresponding
    // vertex descriptor in the JaccardGraph.
    using Graph = boost::adjacency_list<
            boost::listS, boost::vecS, boost::bidirectionalS,
            JaccardGraph::vertex_descriptor>;
    Graph graph;
    std::map<JaccardGraph::vertex_descriptor, Graph::vertex_descriptor> vertexMap;
    for(uint64_t segmentId: component) {
        const JaccardGraph::vertex_descriptor jv = vertexTable[segmentId];
        const Graph::vertex_descriptor gv = add_vertex(jv, graph);
        vertexMap.insert(make_pair(jv, gv));
    }
    BGL_FORALL_VERTICES(gv0, graph, Graph) {
        const JaccardGraph::vertex_descriptor jv0 = graph[gv0];
        BGL_FORALL_OUTEDGES(jv0, e, jaccardGraph, JaccardGraph) {
            const JaccardGraph::vertex_descriptor jv1 = target(e, jaccardGraph);
            add_edge(vertexMap[jv0], vertexMap[jv1], graph);
        }
    }
    if(debug) {
        cout << "This component has " << num_vertices(graph) <<
            " vertices and " << num_edges(graph) << " edges." << endl;
    }

    // Topological sort of this connected component.
    vector<Graph::vertex_descriptor> reverseTopologicalSort;
    try {
        boost::topological_sort(graph, back_inserter(reverseTopologicalSort));
    } catch (boost::not_a_dag&) {
        if(debug) {
            cout << "Topological sort for this connected component failed." << endl;
            cout << "Computation of assembly path will skip this connected component." << endl;
        }
        return;
    }



    // Find the longest path in this component.
    // See https://en.wikipedia.org/wiki/Longest_path_problem#Acyclic_graphs
    vector<uint64_t> pathLength(component.size(), 0);
    vector<Graph::vertex_descriptor> successor(component.size(), Graph::null_vertex());

    // Process vertices in reverse topological order.
    for(const Graph::vertex_descriptor gv0: reverseTopologicalSort) {
        BGL_FORALL_OUTEDGES(gv0, e, graph, Graph) {
            const Graph::vertex_descriptor gv1 = target(e, graph);
            if(pathLength[gv1] + 1 > pathLength[gv0]) {
                pathLength[gv0] = pathLength[gv1] + 1;
                successor[gv0] = gv1;
            }
        }
    }

    // Find the vertex with the longest pathLength.
    // This will be the first vertex of the longest path.
    Graph::vertex_descriptor gv0 =
            std::max_element(pathLength.begin(), pathLength.end()) - pathLength.begin();

    // Find the longest path by following the successors.
    vector<uint64_t> longestPath;
    longestPath.push_back(jaccardGraph[graph[gv0]].segmentId);
    while(true) {
        const Graph::vertex_descriptor gv1 = successor[gv0];
        if(gv1 == Graph::null_vertex()) {
            break;
        }
        longestPath.push_back(jaccardGraph[graph[gv1]].segmentId);
        gv0 = gv1;
    }

    // Store the longest path.
    assemblyPaths.push_back(longestPath);

    if(debug) {
        cout << "Longest path has " << longestPath.size() << " segments:" << endl;
        for(const uint64_t segmentId: longestPath) {
            cout << segmentId << " ";
        }
        cout << endl;
    }
}
