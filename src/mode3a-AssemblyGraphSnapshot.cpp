#include "mode3a-AssemblyGraphSnapshot.hpp"
#include "mode3a-AssemblyGraph.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>



AssemblyGraphSnapshot::AssemblyGraphSnapshot(
    const AssemblyGraph& assemblyGraph,
    const string& name,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner),
    name(name)
{
    // Store the segments.
    createNew(segments, name + "-segments");
    std::map<AssemblyGraph::vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexMap.insert(make_pair(v, segments.size()));
        segments.push_back(Segment(assemblyGraph[v]));
    }

    // Store the links.
    createNew(links, name + "-links");
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraph::vertex_descriptor v0 = source(e, assemblyGraph);
        const AssemblyGraph::vertex_descriptor v1 = target(e, assemblyGraph);
        links.push_back({vertexMap[v0], vertexMap[v1]});
    }

    // Create the paths and the path entries for each vertex.
    createNew(paths, name + "-paths");
    vector< vector<PathEntry> > pathEntries(segments.size());
    for(uint64_t i=0; i<paths.size(); i++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        const auto assemblyGraphPath = assemblyGraph.paths[i];
        paths.appendVector();
        uint64_t position = 0;
        for(const AssemblyGraph::vertex_descriptor v: assemblyGraphPath) {
            const uint64_t vertexId = vertexMap[v];
            paths.append(vertexId);
            pathEntries[orientedReadId.getValue()].push_back({orientedReadId, position++});
        }
    }
    createNew(vertexPathEntries, name + "-vertexPathEntries");
    for(const vector<PathEntry>& v: pathEntries) {
        vertexPathEntries.appendVector(v);
    }
}



AssemblyGraphSnapshot::Segment::Segment(const AssemblyGraphVertex& segment) :
    segmentId(segment.segmentId),
    id(segment.id)
{}



void AssemblyGraphSnapshot::getTransitions(
    uint64_t linkId,
    vector<Transition>& transitions) const
{
    transitions.clear();
    const Link& link = links[linkId];
    const uint64_t vertexId0 = link.vertexId0;
    const uint64_t vertexId1 = link.vertexId1;

    // Loop over path entries of vertex0.
    for(const PathEntry& pathEntry: vertexPathEntries[vertexId0]) {
        const OrientedReadId orientedReadId = pathEntry.orientedReadId;
        const uint64_t position0 = pathEntry.position;
        const uint64_t position1 = position0 + 1;
        const auto path = paths[orientedReadId.getValue()];
        if(position1 < path.size()) {
            if(path[position1] == vertexId1) {
                transitions.push_back({orientedReadId, position0});
            }
        }
    }
}


