// Shasta.
#include "Assembler.hpp"
#include "invalid.hpp"
#include "mode3a-AssemblyGraphSnapshot.hpp"
#include "mode3a-LocalAssemblyGraph.hpp"
using namespace shasta;
using namespace mode3a;



void Assembler::exploreMode3aAssemblyGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters for the request.
    mode3a::LocalAssemblyGraph::SvgOptions options(request);

    uint64_t snapshotIndex = 0;
    getParameterValue(request, "snapshotIndex", snapshotIndex);

    uint64_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    uint64_t minLinkCoverage = 4;
    getParameterValue(request, "minLinkCoverage", minLinkCoverage);

    uint64_t startSegmentId;
    const bool startSegmentIdIsPresent = getParameterValue(request, "startSegmentId", startSegmentId);

    uint64_t startSegmentReplicaIndex = 0;
    getParameterValue(request, "startSegmentReplicaIndex", startSegmentReplicaIndex);

    double timeout = 30.;
    getParameterValue(request, "timeout", timeout);



    // Write the form.
    html <<
        "<h2>Local assembly graph</h2>"
        "<form>"
        "<table>"

        "<tr>"
        "<td>Snapshot index"
        "<td class=centered><input type=text required name=snapshotIndex size=8 style='text-align:center'"
        " value='" << snapshotIndex <<
        "'>"

        "<tr>"
        "<td>Start segment"
        "<td class=centered>"
        "<table>"
        "<tr><td class=left>Segment id<td>"
        "<input type=text required name=startSegmentId size=8 style='text-align:center'"
        " value='" << (startSegmentIdIsPresent ? to_string(startSegmentId) : "") <<
        "'>"
        "<tr><td class=left>Start segment replica index<td>"
        "<input type=text required name=startSegmentReplicaIndex size=8 style='text-align:center'"
        " value='" << startSegmentReplicaIndex <<
        "'>"
        "</table>"

        "<tr>"
        "<td>Maximum distance"
        "<td class=centered><input type=text name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance <<
        "'>"

        "<tr>"
        "<td>Minimum link coverage"
        "<td class=centered><input type=text name=minLinkCoverage size=8 style='text-align:center'"
        " value='" << minLinkCoverage <<
        "'>"

        "<tr>"
        "<td>Timeout for layout (seconds)"
        "<td class=centered><input type=text name=timeout size=8 style='text-align:center'"
        " value='" << timeout <<
        "'>";

    options.addFormRows(html);

    html <<
        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";



    if(not startSegmentIdIsPresent) {
        return;
    }


    // Access the requested snapshot.
    const uint64_t snapshotCount = mode3aAssemblyData.assemblyGraphSnapshots.size();
    if(snapshotIndex >= snapshotCount) {
        html << "<br>Invalid snapshot index. The number of available snapshots is " << snapshotCount <<
            ". Valid snapshot indexes are 0 through " << snapshotCount - 1 << ".";
        return;
    }
    const AssemblyGraphSnapshot& snapshot = *mode3aAssemblyData.assemblyGraphSnapshots[snapshotIndex];



    // Locate the vertexId for this segmentId and segmentCopyIndex.
    if(startSegmentId >= snapshot.vertexTable.size()) {
        html << "<br>Invalid segment id. Valid segment ids are 0 through " << snapshot.vertexTable.size() - 1 << ".";
        return;
    }
    auto v = snapshot.vertexTable[startSegmentId];
    if(startSegmentReplicaIndex >= v.size() or v[startSegmentReplicaIndex] == invalid<uint64_t> ) {
        html << "<br>Invalid segment replica index.<br>Valid segment replica indexes for this segment are:";
        for(uint64_t segmentReplicaIndex=0; segmentReplicaIndex<v.size(); segmentReplicaIndex++) {
            if(v[segmentReplicaIndex] != invalid<uint64_t>) {
                html << " " << segmentReplicaIndex;
            }
        }
        return;
    }
    const uint64_t startVertexId = v[startSegmentReplicaIndex];
    SHASTA_ASSERT(startVertexId < snapshot.vertexVector.size());
    const AssemblyGraphSnapshot::Vertex& startVertex = snapshot.vertexVector[startVertexId];
    SHASTA_ASSERT(startVertex.segmentId == startSegmentId);
    SHASTA_ASSERT(startVertex.segmentReplicaIndex == startSegmentReplicaIndex);

    html << "<h1>Local assembly graph near segment " << startSegmentId <<
        " replica " << startSegmentReplicaIndex << "</h1>";

    // Create this local assembly graph.
    mode3a::LocalAssemblyGraph localAssemblyGraph(snapshot, startVertexId, maxDistance, minLinkCoverage);
    html << "<p>The local assembly graph has " << num_vertices(localAssemblyGraph) <<
        " vertices and " << num_edges(localAssemblyGraph) << " edges." << endl;

    // Compute its layout.
    localAssemblyGraph.computeLayout(options, timeout);
    localAssemblyGraph.computeSegmentTangents();

    // Display the local assembly graph.
    localAssemblyGraph.writeHtml(html, options);


}



void Assembler::exploreMode3aTangleMatrix(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters for the request.

    uint64_t snapshotIndex = 0;
    getParameterValue(request, "snapshotIndex", snapshotIndex);

    uint64_t segmentId;
    const bool segmentIdIsPresent = getParameterValue(request, "segmentId", segmentId);

    uint64_t segmentReplicaIndex = 0;
    getParameterValue(request, "segmentReplicaIndex", segmentReplicaIndex);


    // Write the form.
    html <<
        "<h2>Tangle matrix</h2>"
        "<form>"
        "<table>"

        "<tr>"
        "<td>Snapshot index"
        "<td class=centered><input type=text required name=snapshotIndex size=8 style='text-align:center'"
        " value='" << snapshotIndex <<
        "'>"

        "<tr>"
        "<td class=left>Segment id<td>"
        "<input type=text required name=segmentId size=8 style='text-align:center'"
        " value='" << (segmentIdIsPresent ? to_string(segmentId) : "") <<
        "'>"
        "<tr><td class=left>Segment replica index<td>"
        "<input type=text required name=segmentReplicaIndex size=8 style='text-align:center'"
        " value='" << segmentReplicaIndex <<
        "'>"
        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";



    if(not segmentIdIsPresent) {
        return;
    }


    // Access the requested snapshot.
    const uint64_t snapshotCount = mode3aAssemblyData.assemblyGraphSnapshots.size();
    if(snapshotIndex >= snapshotCount) {
        html << "<br>Invalid snapshot index. The number of available snapshots is " << snapshotCount <<
            ". Valid snapshot indexes are 0 through " << snapshotCount - 1 << ".";
        return;
    }
    const AssemblyGraphSnapshot& snapshot = *mode3aAssemblyData.assemblyGraphSnapshots[snapshotIndex];



    // Locate the vertexId for this segmentId and segmentCopyIndex.
    if(segmentId >= snapshot.vertexTable.size()) {
        html << "<br>Invalid segment id. Valid segment ids are 0 through " << snapshot.vertexTable.size() - 1 << ".";
        return;
    }
    auto v = snapshot.vertexTable[segmentId];
    if(segmentReplicaIndex >= v.size() or v[segmentReplicaIndex] == invalid<uint64_t> ) {
        html << "<br>Invalid segment replica index.<br>Valid segment replica indexes for this segment are:";
        for(uint64_t segmentReplicaIndex=0; segmentReplicaIndex<v.size(); segmentReplicaIndex++) {
            if(v[segmentReplicaIndex] != invalid<uint64_t>) {
                html << " " << segmentReplicaIndex;
            }
        }
        return;
    }
    const uint64_t vertexId = v[segmentReplicaIndex];
    SHASTA_ASSERT(vertexId < snapshot.vertexVector.size());
    const AssemblyGraphSnapshot::Vertex& vertex = snapshot.vertexVector[vertexId];
    SHASTA_ASSERT(vertex.segmentId == segmentId);
    SHASTA_ASSERT(vertex.segmentReplicaIndex == segmentReplicaIndex);

    html << "<h1>Tangle matrix at segment " << segmentId <<
        " replica " << segmentReplicaIndex << "</h1>";

    // Compute the tangle matrix at this vertex.
    std::map< pair<uint64_t, uint64_t>, uint64_t> tangleMatrix;
    snapshot.computeTangleMatrix(vertexId, tangleMatrix);

    if(tangleMatrix.empty()) {
        html << "The tangle matrix is empty.";
        return;
    }



    // For now write it out in flat form.
    for(const auto& p: tangleMatrix) {
        const auto& vertexIds = p.first;
        const uint64_t coverage = p.second;
        const uint64_t vertexId0 = vertexIds.first;
        const uint64_t vertexId1 = vertexIds.second;
        const auto& vertex0 = snapshot.vertexVector[vertexId0];
        const auto& vertex1 = snapshot.vertexVector[vertexId1];

        html << "<br>From " << vertex0.stringId() << " to " << vertex1.stringId() <<
            " coverage " << coverage;

    }
}

