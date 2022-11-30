// Shasta
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "deduplicate.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"



AssemblyGraph::AssemblyGraph(
    const PackedMarkerGraph& packedMarkerGraph) :
    packedMarkerGraph(packedMarkerGraph)
{
    createSegmentsAndPaths();
    createLinks();
}



void AssemblyGraph::createSegmentsAndPaths()
{
    AssemblyGraph& assemblyGraph = *this;

    // Initially, create a vertex for each segment in the packedMarkerGraph.
    vector<vertex_descriptor> vertexTable(packedMarkerGraph.segments.size());
    for(uint64_t segmentId=0; segmentId<packedMarkerGraph.segments.size(); segmentId++) {
        vertexTable[segmentId] = add_vertex(AssemblyGraphVertex(segmentId), assemblyGraph);
    }
    vertexCountBySegment.resize(vertexTable.size(), 1);

    // The sequence of segments visited by each oriented read
    // is a path. Initially, we construct it from the corresponding journey
    // in the PackedMarkerGraph.
    // While constructing the paths, we also store path entries in the vertices.
    paths.resize(packedMarkerGraph.journeys.size());
    for(uint64_t i=0; i<paths.size(); i++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        const auto journey = packedMarkerGraph.journeys[i];
        auto& path = paths[i];
        for(uint64_t position=0; position<journey.size(); position++) {
            const uint64_t segmentId = journey[position];
            const vertex_descriptor v = vertexTable[segmentId];
            path.push_back(v);
            assemblyGraph[v].pathEntries.push_back({orientedReadId, position});
        }
    }
}



void AssemblyGraph::createLinks()
{
    AssemblyGraph& assemblyGraph = *this;

    // Gather transitions for all oriented reads.
    vector< pair<vertex_descriptor, vertex_descriptor> > transitions;
    for(uint64_t i=0; i<paths.size(); i++) {

        // Loop over the path for this oriented read.
        const auto& path = paths[i];
        for(uint64_t position1=1; position1<path.size(); position1++) {
            const uint64_t position0 = position1 - 1;
            transitions.push_back({path[position0], path[position1]});
        }
    }
    deduplicate(transitions);

    for(const pair<vertex_descriptor, vertex_descriptor>& transition: transitions) {
        const vertex_descriptor v0 = transition.first;
        const vertex_descriptor v1 = transition.second;
        add_edge(v0, v1, assemblyGraph);
    }
}



// Get the transitions for an edge.
void AssemblyGraph::getEdgeTransitions(
    edge_descriptor e,
    vector<Transition>& transitions) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Access the vertices of this edge.
    const vertex_descriptor v0 = source(e, assemblyGraph);
    const vertex_descriptor v1 = target(e, assemblyGraph);
    const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];



    // Loop over path entries of vertex1.
    transitions.clear();
    for(const PathEntry& pathEntry: vertex1.pathEntries) {
        const uint64_t position1 = pathEntry.position;
        if(position1 == 0) {
            // v1 is at the beginning of the path.
            // There is no previous vertex in the path.
            continue;
        }
        const uint64_t position0 = position1 - 1;

        // Access the path for this oriented read.
        const OrientedReadId orientedReadId = pathEntry.orientedReadId;
        const vector<vertex_descriptor>& path = paths[orientedReadId.getValue()];

        // If the previous entry is not on v0, this does not
        // correspond to a transiton for the edge we are working on.
        if(path[position0] != v0) {
            continue;
        }

        // Store this transition.
        transitions.push_back({position0, position1, orientedReadId});
    }
}



// This is similar to getTransitions above, but it
// just counts the transitions instead of storing them.
uint64_t AssemblyGraph::edgeCoverage(edge_descriptor e) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Access the vertices of this edge.
    const vertex_descriptor v0 = source(e, assemblyGraph);
    const vertex_descriptor v1 = target(e, assemblyGraph);
    const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];



    // Loop over path entries of vertex1.
    uint64_t coverage = 0;
    for(const PathEntry& pathEntry: vertex1.pathEntries) {
        const uint64_t position1 = pathEntry.position;
        if(position1 == 0) {
            // v1 is at the beginning of the path.
            // There is no previous vertex in the path.
            continue;
        }
        const uint64_t position0 = position1 - 1;

        // Access the path for this oriented read.
        const OrientedReadId orientedReadId = pathEntry.orientedReadId;
        const vector<vertex_descriptor>& path = paths[orientedReadId.getValue()];

        // If the previous entry is not on v0, this does not
        // correspond to a transiton for the edge we are working on.
        if(path[position0] != v0) {
            continue;
        }

        ++coverage;
    }

    return coverage;
}



void AssemblyGraph::write(const string& name) const
{
    for(uint64_t minLinkCoverage=2; minLinkCoverage<=6; minLinkCoverage++) {
        writeGfa(name + "-minLinkCoverage-" + to_string(minLinkCoverage) + ".gfa", minLinkCoverage);
    }
    writeLinkCoverageHistogram(name + "-LinkCoverageHistogram.csv");
    writePaths(name + "-paths.csv");
}



void AssemblyGraph::writeLinkCoverageHistogram(const string& name) const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector<uint64_t> histogram;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const uint64_t coverage = edgeCoverage(e);
        if(histogram.size() <= coverage) {
            histogram.resize(coverage+1, 0);
        }
        ++histogram[coverage];
    }

    ofstream csv(name);
    csv << "Coverage,Frequency\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        csv << coverage << "," << histogram[coverage] << "\n";
    }
}



void AssemblyGraph::writePaths(const string& name) const
{
    const AssemblyGraph& assemblyGraph = *this;
    ofstream csv(name);

    for(uint64_t i=0; i<paths.size(); i++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        const auto path = paths[orientedReadId.getValue()];

        csv << orientedReadId << ",";
        for(const vertex_descriptor v: path) {
            csv << assemblyGraph[v].stringId() << ",";
        }
        csv << "\n";
    }

}



void AssemblyGraph::writeGfa(const string& name, uint64_t minLinkCoverage) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream gfa(name);

    // Write the headers.
    gfa << "H\tVN:Z:1.0\n";

    // Write the segments.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        gfa <<"S\t" << assemblyGraph[v].stringId() << "\t*\n";
    }

    // Write the links.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(edgeCoverage(e) < minLinkCoverage) {
            continue;
        }
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        gfa << "L\t" <<
            assemblyGraph[v0].stringId() << "\t+\t" <<
            assemblyGraph[v1].stringId() << "\t+\t0M\n";
    }
}



void AssemblyGraph::simpleDetangle(
    uint64_t minLinkCoverage,
    uint64_t minTangleCoverage)
{
    AssemblyGraph& assemblyGraph = *this;

    vector<vertex_descriptor> allVertices;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        allVertices.push_back(v);
    }
    for(const vertex_descriptor v: allVertices) {
        simpleDetangle(v, minLinkCoverage, minTangleCoverage);
    }
}



void AssemblyGraph::simpleDetangle(
    vertex_descriptor v1,
    uint64_t minLinkCoverage,
    uint64_t minTangleCoverage)
{
    const bool debug = true;

    AssemblyGraph& assemblyGraph = *this;
    const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];

    // Find adjacent vertices by following the reads.
    vector< pair<vertex_descriptor, vertex_descriptor> > adjacentVertices;
    findAdjacentVertices(v1, adjacentVertices);

    // Group them.
    // Store path entry indexes, that is indexes into the pathEntries vector
    // of v1 and into the adjacentVertices vector.
    std::map<vertex_descriptor, vector<uint64_t> > map0;  // By previous vertex
    std::map<vertex_descriptor, vector<uint64_t> > map2;  // By next vertex.
    std::map< pair<vertex_descriptor, vertex_descriptor>, vector<uint64_t> > map02;    // By previous and next vertex.
    for(uint64_t i=0; i<adjacentVertices.size(); i++) {
        const pair<vertex_descriptor, vertex_descriptor>& v02 = adjacentVertices[i];
        const vertex_descriptor v0 = v02.first;
        const vertex_descriptor v2 = v02.second;
        map0[v0].push_back(i);
        map2[v2].push_back(i);
        map02[v02].push_back(i);
    }

    // For detangling, we only consider parent/children
    // that are not null_vertex() and with at least minLinkCoverage
    // oriented reads.
    vector<vertex_descriptor> parents;
    for(const auto& p: map0) {
        const vertex_descriptor v0 = p.first;
        if(v0 != null_vertex() and p.second.size() >= minLinkCoverage) {
            parents.push_back(v0);
        }
    }
    vector<vertex_descriptor> children;
    for(const auto& p: map2) {
        const vertex_descriptor v2 = p.first;
        if(v2 != null_vertex() and p.second.size() >= minLinkCoverage) {
            children.push_back(v2);
        }
    }

    // For now we only attempt detangling if there are at least two
    // parents and two children.
    if(parents.size() < 2 or children.size() < 2){
        return;
    }

    if(debug) {
        cout << "Detangling " << vertex1.stringId() << " with " <<
            parents.size() << " parents and " <<
            children.size() << " children.\n";
        cout << "Parents: ";
        for(const vertex_descriptor parent: parents) {
            cout << " " << assemblyGraph[parent].stringId();
        }
        cout << "\n";
        cout << "Children: ";
        for(const vertex_descriptor child: children) {
            cout << " " << assemblyGraph[child].stringId();
        }
        cout << "\n";

        if(false) {
            cout << "Details of path entries:\n";
            for(uint64_t i=0; i<vertex1.pathEntries.size(); i++) {
                cout << i << " " << vertex1.pathEntries[i].orientedReadId << " " <<
                    assemblyGraph[adjacentVertices[i].first].stringId() << " " <<
                    assemblyGraph[adjacentVertices[i].second].stringId() << "\n";
            }
        }
    }


    // Find the pairs that will generate a new vertex.
    // These are called the "active pairs" here.
    // They are the ones for which map02 contains at least
    // minTangleCoverage entries.
    vector< pair<vertex_descriptor, vertex_descriptor> > activePairs;
    for(const auto& p: map02) {
        const auto& v02 = p.first;
        const vertex_descriptor v0 = v02.first;
        const vertex_descriptor v2 = v02.second;
        if(v0 != null_vertex() and v2 != null_vertex() and p.second.size() >= minTangleCoverage) {
            activePairs.push_back(v02);
        }
    }


    // Each active pair generates a new vertex with the same segmentId as v1.
    vector<vertex_descriptor> newVertices;
    const uint64_t segmentId1 = vertex1.segmentId;
    for(uint64_t i=0; i<activePairs.size(); i++) {
        const vertex_descriptor vNew = add_vertex(
            AssemblyGraphVertex(segmentId1, vertexCountBySegment[segmentId1]++),
            assemblyGraph);
        newVertices.push_back(vNew);
    }
    if(false) {
        cout << "Active pairs:\n";
        for(uint64_t i=0; i<activePairs.size(); i++) {
            const auto& v02 = activePairs[i];
            const vertex_descriptor v0 = v02.first;
            const vertex_descriptor v2 = v02.second;
            cout << assemblyGraph[v0].stringId() << " " <<
                assemblyGraph[v2].stringId() << " new vertex " <<
                assemblyGraph[newVertices[i]].stringId() << "\n";
        }
    }

    // In addition, we create a vertex that will receive path entries
    // that are not in active pairs.
    const vertex_descriptor vNew = add_vertex(
        AssemblyGraphVertex(segmentId1, vertexCountBySegment[segmentId1]++),
        assemblyGraph);
    if(false) {
        cout << "New vertex not associated with an active pair: " << assemblyGraph[vNew].stringId() << endl;
    }



    // Assign each path entry of v1 to one of the new vertices we just created.
    for(uint64_t i=0; i<vertex1.pathEntries.size(); i++) {
        const PathEntry& pathEntry = vertex1.pathEntries[i];
        const auto& v02 = adjacentVertices[i];
        const vertex_descriptor v0 = v02.first;
        const vertex_descriptor v2 = v02.second;
        if(false) {
            cout << "Assigning path entry " << i << " " << flush <<
            assemblyGraph[v0].stringId() << " " << flush <<
            assemblyGraph[v2].stringId() << " to a new vertex." << endl;
        }

        // Look it up in the active pairs.
        auto it = find(activePairs.begin(), activePairs.end(), v02);

        // Find the vertex we are going to add this PathEntry to.
        const vertex_descriptor v = (it == activePairs.end()) ? vNew : newVertices[it - activePairs.begin()];
        if(false) {
            cout << "This path entry will be assigned to vertex " << assemblyGraph[v].stringId() << endl;
        }

        // Add the Path entry to this vertex.
        assemblyGraph[v].pathEntries.push_back(pathEntry);

        // Update the path to reflect this change.
        paths[pathEntry.orientedReadId.getValue()][pathEntry.position] = v;

        // Make sure the edges v0->v and v->v2 exist.
        if(v0 != null_vertex()) {
            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(v0, v, assemblyGraph);
            if(not edgeExists) {
                add_edge(v0, v, assemblyGraph);
            }
        }
        if(v2 != null_vertex()) {
            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(v, v2, assemblyGraph);
            if(not edgeExists) {
                add_edge(v, v2, assemblyGraph);
            }
        }
    }

    // Now we can remove v1.
    clear_vertex(v1, assemblyGraph);
    remove_vertex(v1, assemblyGraph);

}



// Find the previous and next vertex for each PathEntry in a given vertex.
// On return, adjacentVertices contains a pair of vertex descriptors for
// each PathEntry in vertex v, in the same order.
// Those vertex descriptors are the previous and next vertex visited
// by the oriented read for that PathEntry, and can be null_vertex()
// if v is at the beginning or end of the path of an oriented read.
void AssemblyGraph::findAdjacentVertices(
    vertex_descriptor v,
    vector< pair<vertex_descriptor, vertex_descriptor> >& adjacentVertices
) const
{
    const AssemblyGraph& assemblyGraph = *this;
    const AssemblyGraphVertex& vertex = assemblyGraph[v];

    adjacentVertices.clear();
    for(const PathEntry& pathEntry: vertex.pathEntries) {
        const OrientedReadId orientedReadId = pathEntry.orientedReadId;
        const vector<vertex_descriptor>& path = paths[orientedReadId.getValue()];

        const uint64_t position1 = pathEntry.position;

        // Find the previous vertex visited by this path, if any.
        vertex_descriptor v0 = null_vertex();
        if(position1 > 0) {
            const uint64_t position0 = position1 - 1;
            v0 = path[position0];
        }

        // Find the next vertex visited by this path, if any.
        vertex_descriptor v2 = null_vertex();
        const uint64_t position2 = position1 + 1;
        if(position2 < path.size()) {
            v2 = path[position2];
        }

        // Store this pair of adjacent vertices.
        adjacentVertices.push_back({v0, v2});
    }
}

