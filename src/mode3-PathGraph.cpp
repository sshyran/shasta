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
#include <stack>



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
    const uint64_t minSubgraphSize = 8;

    PathGraph& pathGraph = *this;

    createVertices();
    createEdges(minCoverage);
    cout << "The initial path graph has " << num_vertices(pathGraph) <<
        " vertices and " << num_edges(pathGraph) << " edges." << endl;

    // Partition the PathGraph into subgraphs.
    partition(partitionMaxDistance, minSubgraphSize);
    writeGfa("PathGraph");
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
void PathGraph::partition(
    uint64_t maxDistance,
    uint64_t minSubgraphSize)
{
    PathGraph& pathGraph = *this;

    // Mark all vertices as not assigned to any partition.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        pathGraph[v].subgraphId = noSubgraph;
    }

    // Start at all vertices with zero in-degree,
    // plus the boundary vertices we find that way.
    vector<vertex_descriptor> boundaryVertices;
    std::stack<vertex_descriptor> s;
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        if(in_degree(v, pathGraph) == 0) {
            s.push(v);
        }
    }
    uint64_t subgraphId = 0;
    while(not s.empty()) {
        const vertex_descriptor v = s.top();
        s.pop();

        if(pathGraph[v].subgraphId == noSubgraph) {
            partitionIteration(v, maxDistance, subgraphId++, boundaryVertices);
            for(const vertex_descriptor v: boundaryVertices) {
                s.push(v);
            }
        }
    }



    // In exceptional cases, the above procedure might not assign all
    // vertices to a subgraph.
    // This code takes care of that.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        if(pathGraph[v].subgraphId == noSubgraph) {
            partitionIteration(v, maxDistance, subgraphId++, boundaryVertices);
        }
    }



    // Combine small subgraphs with adjacent subgraphs, if possible.
    // This can leave subgraphs with size 0, but we don't worry about that.
    while(true) {

        // Gather the subgraphs based on the current settings of
        // the vertices subgraphId.
        gatherSubgraphs();

        // Find the small subgraphs.
        std::set<uint64_t> smallSubgraphs;
        for(uint64_t subgraphId=0; subgraphId<subgraphs.size(); subgraphId++) {
            const vector<vertex_descriptor>& subgraph = subgraphs[subgraphId];
            const uint64_t subgraphSize = subgraph.size();
            if((subgraphSize != 0) and (subgraph.size() < minSubgraphSize)) {
                smallSubgraphs.insert(subgraphId);
            }
        }



        // Try and merge small subgraphs with adjacent subgraphs.

        // Loop over small subgraphs.
        bool changesWereMade = false;
        for(uint64_t subgraphId0: smallSubgraphs) {
            const vector<vertex_descriptor>& subgraph0 = subgraphs[subgraphId0];
            const uint64_t subgraph0Size = subgraph0.size();
            SHASTA_ASSERT(subgraph0Size < minSubgraphSize);

            // Find adjacent subgraphs and their sizes.
            vector< pair<uint64_t, uint64_t> > adjacentSubgraphsTable; // (size, subgraphId) of adjacent.
            for(const vertex_descriptor v0: subgraph0) {
                BGL_FORALL_OUTEDGES(v0, e, pathGraph, PathGraph) {
                    const vertex_descriptor v1 = target(e, pathGraph);
                    const uint64_t subgraphId1 = pathGraph[v1].subgraphId;
                    if(subgraphId1 != subgraphId0){
                        adjacentSubgraphsTable.push_back(make_pair(subgraphs[subgraphId1].size(), subgraphId1));
                    }
                }
                BGL_FORALL_INEDGES(v0, e, pathGraph, PathGraph) {
                    const vertex_descriptor v1 = source(e, pathGraph);
                    const uint64_t subgraphId1 = pathGraph[v1].subgraphId;
                    if(subgraphId1 != subgraphId0){
                        adjacentSubgraphsTable.push_back(make_pair(subgraphs[subgraphId1].size(), subgraphId1));
                    }
                }
            }
            sort(adjacentSubgraphsTable.begin(), adjacentSubgraphsTable.end());

            // Merge it with the smallest adjacent subgraph.
            const uint64_t subgraphId1 = adjacentSubgraphsTable.front().second;
            smallSubgraphs.erase(subgraphId1);
            for(const vertex_descriptor v0: subgraph0) {
                pathGraph[v0].subgraphId = subgraphId1;
            }
            changesWereMade = true;
        }

        if(not changesWereMade) {
            break;
        }
    }



    // Subgraph statistics.
    cout << "Partitioned the path graph into " << subgraphs.size() << " subgraphs." << endl;
    histogramSubgraphs();

    // Count the edges across subgraphs.
    uint64_t crossEdgeCount = 0;
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        const vertex_descriptor v0 = source(e, pathGraph);
        const vertex_descriptor v1 = target(e, pathGraph);
        if(pathGraph[v0].subgraphId != pathGraph[v1].subgraphId) {
            ++crossEdgeCount;
        }
    }
    cout << "Number of edges across subgraphs is " << crossEdgeCount << endl;
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



// Gather subgraphs using the subgraphId stored in each vertex.
void PathGraph::gatherSubgraphs()
{
    PathGraph& pathGraph = *this;

    subgraphs.clear();
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        const uint64_t subgraphId = pathGraph[v].subgraphId;
        SHASTA_ASSERT(subgraphId != noSubgraph);

        if(subgraphId >= subgraphs.size()) {
            subgraphs.resize(subgraphId + 1);
        }

        subgraphs[subgraphId].push_back(v);
    }
}



void PathGraph::histogramSubgraphs()
{
    vector<uint64_t> histogram;
    for(const vector<vertex_descriptor>& subgraph: subgraphs) {
        const uint64_t subgraphSize = subgraph.size();
        if(subgraphSize >= histogram.size()) {
            histogram.resize(subgraphSize + 1, 0);
        }
        ++histogram[subgraphSize];
    }

    ofstream csv("PathGraphSubgraphHistogram.csv");
    csv << "Size,Frequency,Vertices\n";
    for(uint64_t subgraphSize=0; subgraphSize<histogram.size(); subgraphSize++) {
        const uint64_t frequency = histogram[subgraphSize];
        csv << subgraphSize << ",";
        csv << frequency << ",";
        csv << subgraphSize*frequency << "\n";
    }
}




void PathGraph::writeGfa(const string& baseName) const
{
    const PathGraph& pathGraph = *this;

    // Open the gfa and write the header.
    ofstream gfa(baseName + ".gfa");
    gfa << "H\tVN:Z:1.0\n";

    // Open the csv and write the header.
    ofstream csv(baseName + ".csv");
    csv << "Segment,Color,SubgraphId\n";

    // Write each vertex as a segment.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        gfa <<
            "S\t" <<
            pathGraph[v].id << "\t" // Segment name
            "*"                     // Segment length
            "\n";


        // Color based on the subgraphId.
        const uint64_t subgraphId = pathGraph[v].subgraphId;
        string color = "LightGrey";
        if(subgraphId != noSubgraph) {
            const uint64_t r = MurmurHash2(&subgraphId,  sizeof(subgraphId),  231) &255;
            const uint64_t g = MurmurHash2(&subgraphId,  sizeof(subgraphId),  233) &255;
            const uint64_t b = MurmurHash2(&subgraphId,  sizeof(subgraphId),  235) &255;

            std::ostringstream s;
            s.fill('0');
            s << "#";
            s << hex << std::setw(2) << r;
            s << hex << std::setw(2) << g;
            s << hex << std::setw(2) << b;
            color = s.str();
        }

        csv << pathGraph[v].id << "," << color << "," << subgraphId << "\n";

    }

    // Write each edge as a link.
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        const vertex_descriptor v0 = source(e, pathGraph);
        const vertex_descriptor v1 = target(e, pathGraph);
        gfa <<
            "L\t" <<
            pathGraph[v0].id << "\t+\t" <<
            pathGraph[v1].id << "\t+\t0M\n";
    }

}