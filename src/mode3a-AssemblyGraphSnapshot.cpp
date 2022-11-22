#include "mode3a-AssemblyGraphSnapshot.hpp"
#include "invalid.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"



// This creates a snapshot of the AssemblyGraph in the current state.
AssemblyGraphSnapshot::AssemblyGraphSnapshot(
    const AssemblyGraph& assemblyGraph,
    const string& name,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner),
    name(name),
    packedMarkerGraph(assemblyGraph.packedMarkerGraph)
{
    // Store the segments.
    createNew(vertexVector, name + "-vertices");
    std::map<AssemblyGraph::vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexMap.insert(make_pair(v, vertexVector.size()));
        vertexVector.push_back(Vertex(assemblyGraph[v]));
    }

    // Store the edges.
    createNew(edgeVector, name + "-edges");
    std::map<AssemblyGraph::edge_descriptor, uint64_t> edgeMap;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraph::vertex_descriptor v0 = source(e, assemblyGraph);
        const AssemblyGraph::vertex_descriptor v1 = target(e, assemblyGraph);
        edgeMap.insert(make_pair(e, edgeVector.size()));
        edgeVector.push_back({vertexMap[v0], vertexMap[v1]});
    }

    // Create the paths and the path entries for each vertex.
    createNew(paths, name + "-paths");
    vector< vector<PathEntry> > pathEntries(vertexVector.size());
    for(uint64_t i=0; i<assemblyGraph.paths.size(); i++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        const auto assemblyGraphPath = assemblyGraph.paths[i];
        paths.appendVector();
        uint64_t position = 0;
        for(const AssemblyGraph::vertex_descriptor v: assemblyGraphPath) {
            const uint64_t vertexId = vertexMap[v];
            paths.append(vertexId);
            pathEntries[vertexId].push_back({orientedReadId, position++});
        }
    }
    createNew(vertexPathEntries, name + "-vertexPathEntries");
    for(const vector<PathEntry>& v: pathEntries) {
        vertexPathEntries.appendVector(v);
    }


    // Compute connectivity.
    createNew(edgesBySource, name + "-edgesBySource");
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        edgesBySource.appendVector();
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
            edgesBySource.append(edgeMap[e]);
        }
    }
    createNew(edgesByTarget, name + "-edgesByTarget");
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        edgesByTarget.appendVector();
        BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
            edgesByTarget.append(edgeMap[e]);
        }
    }

    createVertexTable(assemblyGraph.packedMarkerGraph);
}



// The vertex table is a data structure that allows to get a vertexId (index in vertexVector)
// given a segmentId and segmentCopyIndex.
// Indexed by segmentId.
// vertexTable[segmentId][segmentCopyIndex] contains the vertexId
// with the given segmentId and segmentCopyIndex, or invalid<uint64_t>
// if no such vertex.
void AssemblyGraphSnapshot::createVertexTable(
    const PackedMarkerGraph& packedMarkerGraph)
{
    const uint64_t segmentCount = packedMarkerGraph.segments.size();

    // Create a temporary vertex table stored in a way that is easy to manipulate.
    vector< vector<uint64_t> > tmpVertexTable(segmentCount);
    for(uint64_t vertexId=0; vertexId<vertexVector.size(); vertexId++) {
        const Vertex& vertex = vertexVector[vertexId];
        auto& v = tmpVertexTable[vertex.segmentId];
        if(v.size() <= vertex.segmentCopyIndex) {
            v.resize(vertex.segmentCopyIndex + 1, invalid<uint64_t>);
        }
        v[vertex.segmentCopyIndex] = vertexId;
    }

    // Now copy it to its permanent location.
    createNew(vertexTable, name + "-vertexTable");
    for(const auto& v: tmpVertexTable) {
        vertexTable.appendVector(v);
    }
}



// This accesses an existing snapshot.
AssemblyGraphSnapshot::AssemblyGraphSnapshot(
    const string& name,
    const MappedMemoryOwner& mappedMemoryOwner,
    const PackedMarkerGraph& packedMarkerGraph) :
    MappedMemoryOwner(mappedMemoryOwner),
    name(name),
    packedMarkerGraph(packedMarkerGraph)
{
    accessExistingReadOnly(vertexVector, name + "-vertices");
    accessExistingReadOnly(edgeVector, name + "-edges");
    accessExistingReadOnly(paths, name + "-paths");
    accessExistingReadOnly(vertexPathEntries, name + "-vertexPathEntries");
    accessExistingReadOnly(edgesBySource, name + "-edgesBySource");
    accessExistingReadOnly(edgesByTarget, name + "-edgesByTarget");
    accessExistingReadOnly(vertexTable, name + "-vertexTable");
}



AssemblyGraphSnapshot::Vertex::Vertex(const AssemblyGraphVertex& vertex) :
    segmentId(vertex.segmentId),
    segmentCopyIndex(vertex.segmentCopyIndex)
{}



// Get the length of assembled sequence for a vertex.
uint64_t AssemblyGraphSnapshot::getVertexAssembledSequenceLength(uint64_t vertexId) const
{
    const Vertex& vertex = vertexVector[vertexId];
    return packedMarkerGraph.segmentSequences[vertex.segmentId].size();
}



void AssemblyGraphSnapshot::getTransitions(
    uint64_t linkId,
    vector<Transition>& transitions) const
{
    transitions.clear();
    const Edge& edge = edgeVector[linkId];
    const uint64_t vertexId0 = edge.vertexId0;
    const uint64_t vertexId1 = edge.vertexId1;

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

uint64_t AssemblyGraphSnapshot::linkCoverage(uint64_t linkId) const
{
    uint64_t coverage = 0;

    const Edge& edge = edgeVector[linkId];
    const uint64_t vertexId0 = edge.vertexId0;
    const uint64_t vertexId1 = edge.vertexId1;

    // Loop over path entries of vertex0.
    for(const PathEntry& pathEntry: vertexPathEntries[vertexId0]) {
        const OrientedReadId orientedReadId = pathEntry.orientedReadId;
        const uint64_t position0 = pathEntry.position;
        const uint64_t position1 = position0 + 1;
        const auto path = paths[orientedReadId.getValue()];
        if(position1 < path.size()) {
            if(path[position1] == vertexId1) {
                ++coverage;
            }
        }
    }

    // cout << "*** " << linkId << " " << coverage << "\n";
    return coverage;
}




void AssemblyGraphSnapshot::write() const
{
    writePaths();
    for(uint64_t minLinkCoverage=2; minLinkCoverage<=6; minLinkCoverage++) {
        writeGfa(minLinkCoverage);
    }
}



void AssemblyGraphSnapshot::writeGfa(uint64_t minLinkCoverage) const
{

    ofstream gfa(name + "-minLinkCoverage-" + to_string(minLinkCoverage) + ".gfa");

    // Write the header.
    gfa << "H\tVN:Z:1.0\n";

    // Write the segments.
    for(const Vertex& vertex: vertexVector) {
        gfa <<"S\t" << vertex.stringId() << "\t*\n";
    }

    // Write the links.
    for(uint64_t linkId=0; linkId<edgeVector.size(); linkId++) {
        if(linkCoverage(linkId) < minLinkCoverage) {
            continue;
        }
        const Edge& edge = edgeVector[linkId];
        gfa << "L\t" <<
            vertexVector[edge.vertexId0].stringId() << "\t+\t" <<
            vertexVector[edge.vertexId1].stringId() << "\t+\t0M\n";
    }

}



void AssemblyGraphSnapshot::writePaths() const
{
    ofstream csv(name + "-paths.csv");

    for(uint64_t i=0; i<paths.size(); i++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        const auto path = paths[orientedReadId.getValue()];

        csv << orientedReadId << ",";
        for(const uint64_t vertexId: path) {
            csv << vertexVector[vertexId].stringId() << ",";
        }
        csv << "\n";
    }

}

