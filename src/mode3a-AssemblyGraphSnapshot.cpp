#include "mode3a-AssemblyGraphSnapshot.hpp"
#include "invalid.hpp"
#include "MarkerGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
using namespace shasta;
using namespace mode3a;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "algorithm.hpp"
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
    createNew(vertexJourneyEntries, name + "-vertexJourneyEntries");
    std::map<AssemblyGraph::vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& assemblyGraphVertex = assemblyGraph[v];
        vertexMap.insert(make_pair(v, vertexVector.size()));
        vertexVector.push_back(Vertex(assemblyGraphVertex));
        vertexJourneyEntries.appendVector(assemblyGraphVertex.journeyEntries);
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



    // Copy the journeys to the snapshot.
    createNew(journeys, name + "-journeys");
    for(uint64_t i=0; i<assemblyGraph.journeys.size(); i++) {
        const vector<AssemblyGraph::vertex_descriptor>& assemblyGraphJourney = assemblyGraph.journeys[i];
        journeys.appendVector();
        for(const AssemblyGraph::vertex_descriptor v: assemblyGraphJourney) {
            uint64_t vertexId = invalid<uint64_t>;
            auto it = vertexMap.find(v);
            if(it != vertexMap.end()) {
                vertexId = it->second;
            }
            journeys.append(vertexId);
        }
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
        if(v.size() <= vertex.segmentReplicaIndex) {
            v.resize(vertex.segmentReplicaIndex + 1, invalid<uint64_t>);
        }
        v[vertex.segmentReplicaIndex] = vertexId;
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
    accessExistingReadOnly(journeys, name + "-journeys");
    accessExistingReadOnly(vertexJourneyEntries, name + "-vertexJourneyEntries");
    accessExistingReadOnly(edgesBySource, name + "-edgesBySource");
    accessExistingReadOnly(edgesByTarget, name + "-edgesByTarget");
    accessExistingReadOnly(vertexTable, name + "-vertexTable");
}



AssemblyGraphSnapshot::Vertex::Vertex(const AssemblyGraphVertex& vertex) :
    segmentId(vertex.segmentId),
    segmentReplicaIndex(vertex.segmentReplicaIndex)
{}

// Get the stringId for a given vertexId, or "None" if vertexId is invalid<uint64_t>.
string AssemblyGraphSnapshot::vertexStringId(uint64_t vertexId) const
{
    if(vertexId == invalid<uint64_t>) {
        return "None";
    } else {
        SHASTA_ASSERT(vertexId < vertexVector.size());
        const Vertex& vertex = vertexVector[vertexId];
        return vertex.stringId();
    }

}



// Get the length of assembled sequence for a vertex.
uint64_t AssemblyGraphSnapshot::getVertexAssembledSequenceLength(uint64_t vertexId) const
{
    const Vertex& vertex = vertexVector[vertexId];
    return packedMarkerGraph.segmentSequences[vertex.segmentId].size();
}



void AssemblyGraphSnapshot::getEdgeTransitions(
    uint64_t edgeId,
    vector<Transition>& transitions) const
{
    transitions.clear();
    const Edge& edge = edgeVector[edgeId];
    const uint64_t vertexId0 = edge.vertexId0;
    const uint64_t vertexId1 = edge.vertexId1;

    // Loop over journey entries of vertex0.
    for(const JourneyEntry& journeyEntry: vertexJourneyEntries[vertexId0]) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const uint64_t position0 = journeyEntry.position;
        const uint64_t position1 = position0 + 1;
        const auto journey = journeys[orientedReadId.getValue()];
        if(position1 < journey.size()) {
            if(journey[position1] == vertexId1) {
                transitions.push_back({orientedReadId, position0});
            }
        }
    }
}



uint64_t AssemblyGraphSnapshot::getEdgeCoverage(uint64_t edgeId) const
{
    uint64_t coverage = 0;

    const Edge& edge = edgeVector[edgeId];
    const uint64_t vertexId0 = edge.vertexId0;
    const uint64_t vertexId1 = edge.vertexId1;

    // Loop over journey entries of vertex0.
    for(const JourneyEntry& journeyEntry: vertexJourneyEntries[vertexId0]) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const uint64_t position0 = journeyEntry.position;
        const uint64_t position1 = position0 + 1;
        const auto journey = journeys[orientedReadId.getValue()];
        if(position1 < journey.size()) {
            if(journey[position1] == vertexId1) {
                ++coverage;
            }
        }
    }

    return coverage;
}



// Find out if the segments of an edge are adjacent in the marker graph.
bool AssemblyGraphSnapshot::segmentsAreAdjacent(uint64_t edgeId) const
{
    // Get this edge.
    const Edge& edge = edgeVector[edgeId];

    // Get the corresponding vertices.
    const uint64_t vertexId0 = edge.vertexId0;
    const uint64_t vertexId1 = edge.vertexId1;
    const Vertex& vertex0 = vertexVector[vertexId0];
    const Vertex& vertex1 = vertexVector[vertexId1];

    // Get the corresponding segments.
    const uint64_t segmentId0 = vertex0.segmentId;
    const uint64_t segmentId1 = vertex1.segmentId;

    // Get the marker graph paths of these segments.
    const auto path0 = packedMarkerGraph.segments[segmentId0];
    const auto path1 = packedMarkerGraph.segments[segmentId1];

    // The the last marker graph edge of path0
    // and the first marker graph edge of path1.
    const uint64_t markerGraphEdgeId0 = path0.back();
    const uint64_t markerGraphEdgeId1 = path1.front();
    const MarkerGraph::Edge& markerGraphEdge0 = packedMarkerGraph.markerGraph.edges[markerGraphEdgeId0];
    const MarkerGraph::Edge& markerGraphEdge1 = packedMarkerGraph.markerGraph.edges[markerGraphEdgeId1];

    return markerGraphEdge0.target == markerGraphEdge1.source;
}



void AssemblyGraphSnapshot::write() const
{
    writeJourneys();
    writeJourneyEntries();
    writeTransitions();

    for(uint64_t minLinkCoverage=2; minLinkCoverage<=10; minLinkCoverage++) {
        writeGfa(minLinkCoverage);
    }
}



void AssemblyGraphSnapshot::writeGfa(uint64_t minLinkCoverage) const
{

    ofstream gfa(name + "-minLinkCoverage-" + to_string(minLinkCoverage) + ".gfa");

    // Write the header.
    gfa << "H\tVN:Z:1.0\n";

    // Write the vertices.
    for(const Vertex& vertex: vertexVector) {
        const span<const Base> sequence = packedMarkerGraph.clippedSequence(vertex.segmentId);
        gfa <<"S\t" << vertex.stringId() << "\t";
        copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(gfa));
        gfa << "\n";
    }

    // Write the edges.
    for(uint64_t edgeId=0; edgeId<edgeVector.size(); edgeId++) {
        if(getEdgeCoverage(edgeId) < minLinkCoverage) {
            continue;
        }
        const Edge& edge = edgeVector[edgeId];
        gfa << "L\t" <<
            vertexVector[edge.vertexId0].stringId() << "\t+\t" <<
            vertexVector[edge.vertexId1].stringId() << "\t+\t0M\n";
    }

}



void AssemblyGraphSnapshot::writeJourneys() const
{
    ofstream csv(name + "-journeys.csv");

    for(uint64_t i=0; i<journeys.size(); i++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        const auto journey = journeys[orientedReadId.getValue()];

        csv << orientedReadId << ",";
        for(const uint64_t vertexId: journey) {
            csv << vertexVector[vertexId].stringId() << ",";
        }
        csv << "\n";
    }

}



void AssemblyGraphSnapshot::writeJourneyEntries() const
{
    ofstream csv(name + "-journeyEntries.csv");
    csv << "VertexId,SegmentId,Replica,OrientedReadId,Position\n";

    for(uint64_t vertexId=0; vertexId<vertexJourneyEntries.size(); vertexId++)
    {
        const Vertex& vertex = vertexVector[vertexId];
        const auto journeyEntries = vertexJourneyEntries[vertexId];

        for(const JourneyEntry& journeyEntry: journeyEntries) {
            csv << vertexId << ",";
            csv << vertex.segmentId << ",";
            csv << vertex.segmentReplicaIndex << ",";
            csv << journeyEntry.orientedReadId << ",";
            csv << journeyEntry.position << "\n";
        }
    }

}



void AssemblyGraphSnapshot::writeTransitions() const
{
    ofstream csv(name + "-transitions.csv");
    csv << "LinkId,SegmentReplica0,SegmentReplica1,OrientedReadId,Position0,Position1\n";

    vector<Transition> transitions;
    for(uint64_t edgeId=0; edgeId<edgeVector.size(); edgeId++) {
        const Edge& edge = edgeVector[edgeId];
        getEdgeTransitions(edgeId, transitions);
        for(const Transition& transition: transitions) {
            csv << edgeId << ",";
            csv << vertexVector[edge.vertexId0].stringId() << ",";
            csv << vertexVector[edge.vertexId1].stringId() << ",";
            csv << transition.orientedReadId << ",";
            csv << transition.position << ",";
            csv << transition.position + 1 << "\n";
        }

    }
}



// Analyze the simple tangle at a given vertex vertexId1 by
// "following the reads" one step forward or backward.
// On return:
// - The previousVertices and nextVertices vectors
//   have size equal to the number of journey entries for vertexId1,
//   and contain the vertices vertexId0 and vertexId2 each of
//   the oriented reads visits immediately before/after visiting vertexId1.
//   This can be invalid<uint64_t> if vertexId1 is the
//   first or last vertex of a journey.
// - inCoverage[vertexId0] gives the number of oriented reads that
//   visit vertexId0 immediately before visiting vertexId1.
// - outCoverage[vertexId2] gives the number of oriented reads that
//   visit vertexId0 immediately after visiting vertexId1.
// - tangle_matrix[make_pair(vertexId0, vertexId2)] is
//   the number of oriented reads that visit vertexId0 immediately before vertexId1, and
//   visit vertexId2 immediately after  vertexId1.
// The three maps can include entries for which vertexId0 and/or vertexId2
// are invalid<uint64_t>, if vertexId1 is the first or last
// vertex of one or more oriented read journeys.
void AssemblyGraphSnapshot::analyzeSimpleTangleAtVertex(
    uint64_t vertexId1,
    vector<uint64_t>& inVertices,
    vector<uint64_t>& outVertices,
    std::map<uint64_t, uint64_t>& inCoverage,
    std::map<uint64_t, uint64_t>& outCoverage,
    std::map< pair<uint64_t, uint64_t>, uint64_t>& tangleMatrix) const
{
    // Access the journey entries for vertexId1.
    auto journeyEntries = vertexJourneyEntries[vertexId1];

    // Initialize output arguments.
    inVertices.clear();
    inVertices.reserve(journeyEntries.size());
    outVertices.clear();
    outVertices.reserve(journeyEntries.size());
    inCoverage.clear();
    outCoverage.clear();
    tangleMatrix.clear();

    // Loop over the journey entries for this vertex.
    for(const JourneyEntry& journeyEntry: journeyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        const uint64_t position1 = journeyEntry.position;
        const auto journey = journeys[orientedReadId.getValue()];

        // Find the vertex this oriented read visits before vertexId1.
        uint64_t vertexId0 = invalid<uint64_t>;
        if(position1 > 0) {
            const uint64_t position0 = position1 - 1;
            vertexId0 = journey[position0];
        }

        // Find the vertex this oriented read visits after vertexId1.
        uint64_t vertexId2 = invalid<uint64_t>;
        if(position1 < journey.size() - 1) {
            const uint64_t position2 = position1 + 1;
            vertexId2 = journey[position2];
        }

        // Store vertexId0 and vertexId2.
        inVertices.push_back(vertexId0);
        outVertices.push_back(vertexId2);

        // Update the output maps.
        ++inCoverage[vertexId0];
        ++outCoverage[vertexId2];
        ++tangleMatrix[make_pair(vertexId0, vertexId2)];
    }
}



// Get deduplicated oriented read ids for a vertex.
// The journey entries are sorted by OrientedReadId,
// so this can be done quickly.
void AssemblyGraphSnapshot::getDeduplicatedOrientedReads(
    uint64_t vertexId,
    vector<OrientedReadId>& orientedReadIds
    ) const
{
    orientedReadIds.clear();
    for(const JourneyEntry& journeyEntry: vertexJourneyEntries) {
        const OrientedReadId orientedReadId = journeyEntry.orientedReadId;
        if(orientedReadIds.empty() or orientedReadId != orientedReadIds.back()) {
            orientedReadIds.push_back(orientedReadId);
        }
    }
}



// Compute the Jaccard similarity of the oriented read composition of two vertices.
// Duplicate oriented reads in the path entries for the vertices are ignored.
// This also computes the number of oriented reads in the union and intersection
// of the two read compositions, as well as vectors containing the deduplicated
// oriented reads for each of the two vertices.
double AssemblyGraphSnapshot::jaccard(
    uint64_t vertexId0,
    uint64_t vertexId1,
    uint64_t& unionCount,
    uint64_t& intersectionCount,
    vector<OrientedReadId>& orientedReadIds0,
    vector<OrientedReadId>& orientedReadIds1,
    vector<OrientedReadId>& unionOrientedReads,
    vector<OrientedReadId>& intersectionOrientedReads
) const
{
    // Get deduplicated OrientedReadId's.
    getDeduplicatedOrientedReads(vertexId0, orientedReadIds0);
    getDeduplicatedOrientedReads(vertexId1, orientedReadIds1);

    // Compute the union.
    unionOrientedReads.clear();
    std::set_union(
        orientedReadIds0.begin(), orientedReadIds0.end(),
        orientedReadIds1.begin(), orientedReadIds1.end(),
        back_inserter(unionOrientedReads));

    // Compute the intersection.
    intersectionOrientedReads.clear();
    std::set_intersection(
        orientedReadIds0.begin(), orientedReadIds0.end(),
        orientedReadIds1.begin(), orientedReadIds1.end(),
        back_inserter(intersectionOrientedReads));

    // Return the Jaccard similarity.
    return
        double(intersectionOrientedReads.size()) /
        double(unionOrientedReads.size());
}
