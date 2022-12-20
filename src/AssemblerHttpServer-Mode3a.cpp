// Shasta.
#include "Assembler.hpp"
#include "AssembledSegment.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "html.hpp"
#include "invalid.hpp"
#include "mode3a-AssemblyGraphSnapshot.hpp"
#include "mode3a-LocalAssemblyGraph.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
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
        "<td>Start segment id"
        "<td class=centered>"
        "<input type=text required name=startSegmentId size=8 style='text-align:center'"
        " value='" << (startSegmentIdIsPresent ? to_string(startSegmentId) : "") <<
        "'>"
        "<tr><td class=left>Start segment replica index<td class=centered>"
        "<input type=text required name=startSegmentReplicaIndex size=8 style='text-align:center'"
        " value='" << startSegmentReplicaIndex <<
        "'>"

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
    string message;
    const uint64_t startVertexId = snapshot.getVertexId(startSegmentId, startSegmentReplicaIndex, message);
    if(startVertexId == invalid<uint64_t>) {
        html << "<br>Invalid combination of start segment id and start segment replicaIndex.<br>" << message;
        return;
    }
    SHASTA_ASSERT(startVertexId < snapshot.vertexVector.size());

    // Locate the corresponding snapshot vertex.
    const AssemblyGraphSnapshot::Vertex& startVertex = snapshot.vertexVector[startVertexId];
    SHASTA_ASSERT(startVertex.segmentId == startSegmentId);
    SHASTA_ASSERT(startVertex.segmentReplicaIndex == startSegmentReplicaIndex);

    html << "<h1>Local assembly graph near segment " << startSegmentId <<
        " replica " << startSegmentReplicaIndex << "</h1>";

    // Create this local assembly graph.
    mode3a::LocalAssemblyGraph localAssemblyGraph =
        (options.layoutMethod == "detailedLimited") ?
        mode3a::LocalAssemblyGraph(snapshot, startVertexId) :
        mode3a::LocalAssemblyGraph(snapshot, startVertexId, maxDistance, minLinkCoverage);


    // Display it.
    if(options.layoutMethod == "detailed" or options.layoutMethod == "detailedLimited") {
        localAssemblyGraph.writeDetailedHtml(html, timeout, options.layoutMethod == "detailedLimited");
    } else {

        // Compute its layout.
        localAssemblyGraph.computeLayout(options, timeout);
        localAssemblyGraph.computeSegmentTangents();

        // Display the local assembly graph.
        localAssemblyGraph.writeHtml(html, options);
    }


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

    uint64_t minLinkCoverage = 4;
    getParameterValue(request, "minLinkCoverage", minLinkCoverage);


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
        "<tr><td class=left>Minimum link coverage<br>(for highlighting)<td>"
        "<input type=text required name=minLinkCoverage size=8 style='text-align:center'"
        " value='" << minLinkCoverage <<
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

    // Analyze the simple tangle at this vertex.
    vector<uint64_t> inVertices;
    vector<uint64_t> outVertices;
    std::map<uint64_t, uint64_t> inCoverage;
    std::map<uint64_t, uint64_t> outCoverage;
    std::map< pair<uint64_t, uint64_t>, uint64_t> tangleMatrix;
    snapshot.analyzeSimpleTangleAtVertex(
        vertexId,
        inVertices,
        outVertices,
        inCoverage,
        outCoverage,
        tangleMatrix);

    // Check if the tangle at this vertex is trivial.
    bool isTrivial = false;
    if(inCoverage.empty()) {
        html << "<br>There are no incoming vertices.";
        isTrivial = true;
    }
    if(outCoverage.empty()) {
        html << "<br>There are no outgoing vertices.";
        isTrivial = true;
    }
    if(isTrivial) {
        SHASTA_ASSERT(tangleMatrix.empty());
        return;
    }
    SHASTA_ASSERT(not tangleMatrix.empty());



    // Display the tangle matrix.

    // Header.
    html <<
        "<table>" <<
        "<tr><td rowspan=2 colspan=2>" <<
        "<th colspan=" << outCoverage.size() << ">Out" <<
        "<td>" <<
        "<tr>";
    for(const auto& p: outCoverage) {
        const uint64_t vertexId = p.first;
        html << "<th style='background-color:PaleGreen'>" << snapshot.vertexStringId(vertexId);
    }
    html << "<th style='background-color:LightPink'>Total";

    // One row for each in-vertex.
    html << "<tr><th rowspan=" << inCoverage.size() << ">In";
    bool isFirstRow = true;
    for(const auto& p: inCoverage) {
        const uint64_t vertexId0 = p.first;
        const uint64_t coverage0 = p.second;
        if(isFirstRow) {
            isFirstRow = false;
        } else {
            html << "<tr>";
        }
        html << "<th style='background-color:PaleGreen'>" << snapshot.vertexStringId(vertexId0);
        for(const auto& p: outCoverage) {
            const uint64_t vertexId2 = p.first;
            const uint64_t coverage2 = p.second;
            const uint64_t m = tangleMatrix[make_pair(vertexId0, vertexId2)];
            html << "<td class=centered";
            if(vertexId0 != invalid<uint64_t> and vertexId2 != invalid<uint64_t>
                and coverage0 >= minLinkCoverage and coverage2 >= minLinkCoverage) {
                html << " style='background-color:LightBlue'";
            }
            html << ">";
            if(m > 0) {
                html << m;
            }
        }
        const string highlightColor = (vertexId0 != invalid<uint64_t> and coverage0 >= minLinkCoverage) ? "LightBlue" : "LightPink";
        html << "<td class=centered style='background-color:" << highlightColor << "'>" << coverage0;
    }

    // Final row.
    html << "<tr><td><th style='background-color:LightPink'>Total";
    uint64_t sum = 0;
    for(const auto& p: outCoverage) {
        const uint64_t vertexId2 = p.first;
        const uint64_t coverage = p.second;
        const string highlightColor = (vertexId2 != invalid<uint64_t> and coverage >= minLinkCoverage) ? "LightBlue" : "LightPink";
        html << "<td class=centered style='background-color:" << highlightColor << "'>" << coverage;
        sum += coverage;
    }
    html << "<td class=centered style='background-color:LightCoral'>" << sum << "</table>";




    // Write out previous vertices and last vertices.
    html << "<h3>Tangle details</h3><p>";
    html << "<table><tr>" <<
        "<th>Oriented<br>read<th>VertexId0<th>VertexId1<th>VertexId2";
    auto journeyEntries = snapshot.vertexJourneyEntries[vertexId];
    SHASTA_ASSERT(inVertices.size()  == journeyEntries.size());
    SHASTA_ASSERT(outVertices.size() == journeyEntries.size());
    for(uint64_t i=0; i<journeyEntries.size(); i++) {
        const JourneyEntry& journeyEntry = journeyEntries[i];
        const uint64_t vertexId0 = inVertices[i];
        const uint64_t vertexId2 = outVertices[i];
        html << "<tr><td class=centered>" << journeyEntry.orientedReadId;
        html << "<td class=centered>" << snapshot.vertexStringId(vertexId0);
        html << "<td class=centered>" << snapshot.vertexStringId(vertexId);
        html << "<td class=centered>" << snapshot.vertexStringId(vertexId2);
    }
    html << "</table>";
}



void Assembler::exploreMode3aAssemblyGraphSegment(
    const vector<string>& request,
    ostream& html)
{
    // Access the PackedMarkerGraph
    auto packedMarkerGraphPointer = mode3aAssemblyData.packedMarkerGraph;
    SHASTA_ASSERT(packedMarkerGraphPointer);
    mode3a::PackedMarkerGraph& packedMarkerGraph = *packedMarkerGraphPointer;

    const auto k = assemblerInfo->k;



    // Get request parameters.
    uint64_t snapshotIndex = 0;
    getParameterValue(request, "snapshotIndex", snapshotIndex);

    uint64_t segmentId;
    const bool segmentIdIsPresent = getParameterValue(request, "segmentId", segmentId);

    uint64_t segmentReplicaIndex = 0;
    getParameterValue(request, "segmentReplicaIndex", segmentReplicaIndex);

    string showOrientedReadsString;
    const bool showOrientedReads = HttpServer::getParameterValue(request,
        "showOrientedReads", showOrientedReadsString);

    string showMarkerGraphPathString;
    const bool showMarkerGraphPath = HttpServer::getParameterValue(request,
        "showMarkerGraphPath", showMarkerGraphPathString);

    string showSequenceString;
    const bool showSequence = HttpServer::getParameterValue(request,
        "showSequence", showSequenceString);

    string showSequenceDetailsString;
    const bool showSequenceDetails = HttpServer::getParameterValue(request,
        "showSequenceDetails", showSequenceDetailsString);



    // Write the form.
    html <<
        "<h2>Display details of an assembly graph segment</h2>"
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

        "<tr>"
        "<td>Show oriented reads"
        "<td class=centered> <input type=checkbox name=showOrientedReads" <<
        (showOrientedReads ? " checked=checked" : "") <<
        ">"

        "<tr>"
        "<td>Show marker graph path"
        "<td class=centered> <input type=checkbox name=showMarkerGraphPath" <<
        (showMarkerGraphPath ? " checked=checked" : "") <<
        ">"

        "<tr>"
        "<td>Show sequence"
        "<td class=centered> <input type=checkbox name=showSequence" <<
        (showSequence ? " checked=checked" : "") <<
        ">"

        "<tr>"
        "<td>Show sequence assembly details"
        "<td class=centered> <input type=checkbox name=showSequenceDetails" <<
        (showSequenceDetails ? " checked=checked" : "") <<
        ">"

        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";

    // If the segmentId was not specified, stop here.
    if(not segmentIdIsPresent) {
        return;
    }

    // Check that we have a valid segmentId.
    if(segmentId >= packedMarkerGraph.segments.size()) {
        html << "Invalid segment id. Maximum valid value is " <<
            packedMarkerGraph.segments.size() - 1 << ".";
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


    // Access some information for this segment.
    const auto path = snapshot.packedMarkerGraph.segments[segmentId];
    const auto sequence = snapshot.packedMarkerGraph.segmentCompleteSequence(segmentId);
    const auto journeyEntries = snapshot.vertexJourneyEntries[vertexId];



    html << "<h1>Snapshot " << snapshotIndex << " segment " << segmentId <<
        " replica " << segmentReplicaIndex << "</h1>";


    // Summary table.
    const auto oldPrecision = html.precision(1);
    const auto oldFlags = html.setf(std::ios_base::fixed, std::ios_base::floatfield);
    html <<
        "<p><table>"
        "<tr><th class=left>Length of marker graph path<td class=centered>" <<
        path.size() <<
        "<tr><th class=left>Length of complete assembled sequence<td class=centered>" <<
        sequence.size() <<
        "<tr><th class=left>Average marker graph edge coverage<td class=centered>" <<
        snapshot.packedMarkerGraph.averageMarkerGraphEdgeCoverage(segmentId) <<
        "<tr><th class=left>Number of journey entries ";
    writeInformationIcon(html, "Number of traversal of this vertex by an oriented read. "
        "Some oriented reads may traverse a vertex more than once.");
    html <<
        "<td class=centered>" <<
        journeyEntries.size() <<
        "</table>";
    html.precision(oldPrecision);
    html.flags(oldFlags);



    // Write the journey entries in a table.
    if(showOrientedReads) {
        html <<
            "<h2>Traversals of this segment by oriented reads</h2>"
            "<table>"
            "<tr>"
            "<th>Oriented<br>read"
            "<th>Journey<br>length"
            "<th>Position<br>in journey";
        for(const JourneyEntry& journeyEntry: journeyEntries) {
            html<<
                "<tr>"
                "<td class=centered>" << journeyEntry.orientedReadId <<
                "<td class=centered>" << snapshot.journeys[journeyEntry.orientedReadId.getValue()].size() <<
                "<td class=centered>" << journeyEntry.position;
        }
        html << "</table>";
    }



    // Write the marker graph path.
    if(showMarkerGraphPath) {
        html <<
            "<h2>Marker graph path</h2>"
            "<table>"
            "<tr>"
            "<th>Position"
            "<th>Edge"
            "<th>Coverage"
            "<th>Source<br>vertex"
            "<th>Target<br>vertex";

        for(uint64_t position=0; position<path.size(); position++) {
            const MarkerGraphEdgeId& edgeId = path[position];
            const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
            const MarkerGraph::VertexId vertexId0 = edge.source;
            const MarkerGraph::VertexId vertexId1 = edge.target;

            html << "<tr>"
                "<td class=centered>" << position <<
                "<td class=centered>" <<
                "<a href='exploreMarkerGraphEdge?edgeId=" << edgeId <<
                "'>" << edgeId << "</a>"
                "<td class=centered>" << markerGraph.edgeMarkerIntervals.size(edgeId) <<
                "<td class=centered>" <<
                "<a href='exploreMarkerGraphVertex?vertexId=" << vertexId0 <<
                "'>" << vertexId0 << "</a>"
                "<td class=centered>" <<
                "<a href='exploreMarkerGraphVertex?vertexId=" << vertexId1 <<
                "'>" << vertexId1 << "</a>"
                "\n";



        }
        html << "</table>";
    }



    // Assembled sequence with details.
    if(showSequence or showSequenceDetails) {


        // Assemble the sequence for this segment.
        AssembledSegment assembledSegment;
        assembleMarkerGraphPath(
            0,
            k,
            markers,
            markerGraph,
            path,
            false,
            assembledSegment);


        // Check that the sequence we have is the same as the stored sequence
        // for this segment.
        SHASTA_ASSERT(std::equal(
            assembledSegment.rawSequence.begin(), assembledSegment.rawSequence.end(),
            sequence.begin(), sequence.end()));

        // Write the sequence.
        if(showSequence) {
            html << "<h1>Assembled sequence</h2>";
            html << "<pre style='font-family:courier'>";
            copy(sequence.begin(), sequence.end(),ostream_iterator<Base>(html));
            html << "</pre>";
        }
        if(showSequenceDetails) {
            html <<
                "<h1>Assembly details</h1>";
            assembledSegment.writeHtml(html, true, true,
                0, uint32_t(assembledSegment.rawSequence.size()));
        }
    }

}


