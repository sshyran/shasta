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
    // Create the segments.
    createNew(segments, name + "-segments");
    std::map<AssemblyGraph::vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexMap.insert(make_pair(v, segments.size()));
        segments.push_back(SegmentInfo(assemblyGraph[v]));
    }

    // Create the links.
    createNew(links, name + "-links");
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraph::vertex_descriptor v0 = source(e, assemblyGraph);
        const AssemblyGraph::vertex_descriptor v1 = target(e, assemblyGraph);
        links.push_back(LinkInfo(vertexMap[v0], vertexMap[v1]));
    }
}



AssemblyGraphSnapshot::SegmentInfo::SegmentInfo(const AssemblyGraphVertex& segment) :
    segmentId(segment.segmentId),
    id(segment.id)
{}



AssemblyGraphSnapshot::LinkInfo::LinkInfo(uint64_t vertexId0, uint64_t vertexId1) :
    vertexId0(vertexId0),
    vertexId1(vertexId1)
{}


