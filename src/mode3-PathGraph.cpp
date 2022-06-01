// Shasta.
#include "mode3-PathGraph.hpp"
#include "orderPairs.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "iostream.hpp"
#include <queue>



// Create the PathGraph from the AssemblyGraph.
// Start with a single segment for each vertex
// (that is, paths of length 1).
PathGraph::PathGraph(const AssemblyGraph& assemblyGraph) :
    MultithreadedObject<PathGraph>(*this),
    assemblyGraph(assemblyGraph)
{
    // HARDWIRED CONSTANTS TO BE EXPOSED WHEN CODE STABILIZES.
    const uint64_t minCoverage = 3;
    const uint64_t partitionMaxDistance = 10;

    PathGraph& pathGraph = *this;

    createVertices();
    createEdges(minCoverage);
    cout << "The initial path graph has " << num_vertices(pathGraph) <<
        " vertices and " << num_edges(pathGraph) << " edges." << endl;

    // Partition the PathGraph into subgraphs.
    partition(partitionMaxDistance);
}



// Initial creation of the vertices.
// Start with a single segment for each vertex
// (that is, paths of length 1).
void PathGraph::createVertices() {

    PathGraph& pathGraph = *this;


    // Create a vertex for each segment in the AssemblyGraph.
    for(uint64_t segmentId=0; segmentId<assemblyGraph.paths.size(); segmentId++) {

        // Create the vertex.
        const vertex_descriptor v = add_vertex(pathGraph);
        PathGraphVertex& vertex = pathGraph[v];
        vertex.id = nextVertexId++;

        // Store the path.
        vertex.path.push_back(segmentId);

        // Store the AssemblyGraphJourneyInterval's.
        const span<const pair<OrientedReadId, uint64_t> > journeyInfos =
            assemblyGraph.assemblyGraphJourneyInfos[segmentId];
        for(const pair<OrientedReadId, uint64_t>& p: journeyInfos) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t position = p.second;
            AssemblyGraphJourneyInterval interval;
            interval.orientedReadId = orientedReadId;
            interval.first = position;
            interval.last = position;
            vertex.journeyIntervals.push_back(interval);
        }
    }

}



// Recreate all edges from scratch, using only the
// information stored in the vertices.
void PathGraph::createEdges(uint64_t minCoverage)
{
    PathGraph& pathGraph = *this;

    // Gather AssemblyGraphJourneyInterval's for all oriented reads.
    vector< vector<pair<AssemblyGraphJourneyInterval, vertex_descriptor> > >
        journeyIntervals(2 * assemblyGraph.readCount());
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        for(const AssemblyGraphJourneyInterval& interval: pathGraph[v].journeyIntervals) {
            journeyIntervals[interval.orientedReadId.getValue()].push_back(
                make_pair(interval, v));
        }
    }
    for(auto& v: journeyIntervals) {
        sort(v.begin(), v.end(),
            OrderPairsByFirstOnly<AssemblyGraphJourneyInterval, vertex_descriptor>());
    }


    // Create the edges.
    for(const auto& orientedReadJourneyIntervals: journeyIntervals) {

        for(uint64_t i=1; i<orientedReadJourneyIntervals.size(); i++) {
            const vertex_descriptor v0 = orientedReadJourneyIntervals[i-1].second;
            const vertex_descriptor v1 = orientedReadJourneyIntervals[i  ].second;

            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(v0, v1, pathGraph);
            if(not edgeExists) {
                tie(e, edgeExists) = add_edge(v0, v1, pathGraph);
                SHASTA_ASSERT(edgeExists);
            }
            ++pathGraph[e].coverage;
        }
    }



    // Remove the low coverage edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        if(pathGraph[e].coverage < minCoverage) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, pathGraph);
    }
}



// Partition the PathGraph into subgraphs.
void PathGraph::partition(uint64_t maxDistance)
{
    PathGraph& pathGraph = *this;

    // Mark all vertices as not assigned to any partition.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        pathGraph[v].subgraphId = noSubgraph;
    }

    // Simple iteration to begin, but we can do better.
    uint64_t subgraphId = 0;
    vector<vertex_descriptor> boundaryVertices;
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        if(pathGraph[v].subgraphId == noSubgraph) {
            partitionIteration(v, maxDistance, subgraphId++, boundaryVertices);
        }
    }

    cout << "Partitioned the path graph into " << subgraphId << " subgraphs." << endl;
}



// A partition iteration does a single BFS starting at v.
// It moves forward from v, avoiding vertices already
// assigned to a subgraph, and up to maxDistance from v.
// It also returns the boundaryVertices, that is the
// vertices found in the process that are at distance maxDistance+1
// from v and are not yet assigned to a subgraph.
// These can then used as starting points new partition iterations.
void PathGraph::partitionIteration(
    vertex_descriptor v,
    uint64_t maxDistance,
    uint64_t subgraphId,
    vector<vertex_descriptor>& boundaryVertices)
{
    PathGraph& pathGraph = *this;

    boundaryVertices.clear();

    // Initialize the BFS.
    std::queue<vertex_descriptor> q;
    q.push(v);
    PathGraphVertex& vertex = pathGraph[v];
    SHASTA_ASSERT(vertex.subgraphId == noSubgraph);
    vertex.subgraphId = subgraphId;
    vertex.distance = 0;

    // BFS loop.
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        const uint64_t distance0 = pathGraph[v0].distance;
        const uint64_t distance1 = distance0 + 1;
        SHASTA_ASSERT(distance0 <= maxDistance);

        // Loop over edges starting at v0.
        BGL_FORALL_OUTEDGES(v0, e01, pathGraph, PathGraph) {
            const vertex_descriptor v1 = target(e01, pathGraph);
            PathGraphVertex& vertex1 = pathGraph[v1];

            // If v1 is already in a subgraph, skip it.
            if(vertex1.subgraphId != noSubgraph) {
                continue;
            }

            // Assign v1 to this subgraph, if it is within maxDistance.
            if(distance1 <= maxDistance) {
                vertex1.subgraphId = subgraphId;
                vertex1.distance = distance1;
            }

            // Queue it or add it to the boundary vertices.
            if(distance1 <= maxDistance) {
                q.push(v1);
            } else {
                SHASTA_ASSERT(distance1 == maxDistance + 1);
                boundaryVertices.push_back(v1);
            }

        }

    }
}


