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

