// Shasta
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
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
        const auto& journey = packedMarkerGraph.journeys[i];
        auto& path = paths[i];
        for(const uint64_t segmentId: journey) {
            const vertex_descriptor v = vertexTable[segmentId];
            auto it = path.insert(path.end(), v);
            assemblyGraph[v].pathEntries.push_back(AssemblyGraphVertex::PathEntry(orientedReadId, it));
        }
    }
}



void AssemblyGraph::createLinks()
{
    AssemblyGraph& assemblyGraph = *this;

    vector<AssemblyGraphEdge::Transition> transitions;
    AssemblyGraphEdge::Transition transition;
    for(uint64_t i=0; i<paths.size(); i++) {
        transition.orientedReadId = OrientedReadId::fromValue(ReadId(i));
        auto& path = paths[i];  // Can't be const because we are storing non-const iterators to it.

        for(auto it=path.begin(); /* Check later*/; ++it) {
            auto itNext = it;
            ++itNext;
            if(itNext == path.end()) {
                break;
            }
            transition.it0 = it;
            transition.it1 = itNext;
            transitions.push_back(transition);
        }
    }
    sort(transitions.begin(), transitions.end());



    // Each streak of transitions with the same vertices generates a link.
    for(auto it=transitions.begin(); it!=transitions.end(); /* Increment later */) {
        const vertex_descriptor v0 = *(it->it0);
        const vertex_descriptor v1 = *(it->it1);

        // Find the streak for these two vertices.
        auto streakBegin = it;
        auto streakEnd = it;
        ++streakEnd;
        while(streakEnd != transitions.end()) {
            if(*(streakEnd->it0) != v0 or *(streakEnd->it1) != v1) {
                break;
            }
            ++streakEnd;
        }

        // The transitions of this streak generate an edge.
        edge_descriptor e;
        tie(e, ignore) = add_edge(v0, v1, assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];
        for(auto jt=streakBegin; jt!=streakEnd; ++jt) {
            edge.transitions.push_back(*jt);
        }

        it = streakEnd;
    }

}



void AssemblyGraph::write(const string& name) const
{
    for(uint64_t minLinkCoverage=2; minLinkCoverage<=6; minLinkCoverage++) {
        writeGfa(name + "-minLinkCoverage-" + to_string(minLinkCoverage) + ".gfa", minLinkCoverage);
    }
    writeLinkCoverageHistogram(name + "-LinkCoverageHistogram.csv");
}



void AssemblyGraph::writeLinkCoverageHistogram(const string& name) const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector<uint64_t> histogram;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const uint64_t coverage = assemblyGraph[e].transitions.size();
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
        if(assemblyGraph[e].coverage() < minLinkCoverage) {
            continue;
        }
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        gfa << "L\t" <<
            assemblyGraph[v0].stringId() << "\t+\t" <<
            assemblyGraph[v1].stringId() << "\t+\t0M\n";
    }
}

