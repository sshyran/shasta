// Shasta
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
using namespace shasta;
using namespace mode3a;



AssemblyGraph::AssemblyGraph(
    const PackedMarkerGraph& packedMarkerGraph) :
    packedMarkerGraph(packedMarkerGraph)
{
    createSegmentsAndPaths();
}



void AssemblyGraph::createSegmentsAndPaths()
{
    AssemblyGraph& assemblyGraph = *this;

    // Initially, create a vertex for each segment in the packedMarkerGraph.
    vector<vertex_descriptor> vertexTable(packedMarkerGraph.segments.size());
    for(uint64_t segmentId=0; segmentId<packedMarkerGraph.segments.size(); segmentId++) {
        vertexTable[segmentId] = add_vertex(AssemblyGraphVertex(segmentId), assemblyGraph);
    }

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
