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

    uint64_t startSegmentId;
    const bool startSegmentIdIsPresent = getParameterValue(request, "startSegmentId", startSegmentId);

    uint64_t startSegmentCopyIndex = 0;
    getParameterValue(request, "startSegmentCopyIndex", startSegmentCopyIndex);

    double timeout = 30.;
    getParameterValue(request, "timeout", timeout);



    // Write the form.
    html <<
        "<h2>Display the local assembly graph near a given vertex</h2>"
        "<form>"
        "<table>"

        "<tr>"
        "<td>Snapshot index"
        "<td class=centered><input type=text required name=snapshotIndex size=8 style='text-align:center'"
        " value='" << snapshotIndex <<
        "'>"

        "<tr>"
        "<td>Start segment id"
        "<td class=centered><input type=text required name=startSegmentId size=8 style='text-align:center'"
        " value='" << (startSegmentIdIsPresent ? to_string(startSegmentId) : "") <<
        "'>"

        "<tr>"
        "<td>Start segment copy index"
        "<td class=centered><input type=text required name=startSegmentCopyIndex size=8 style='text-align:center'"
        " value='" << startSegmentCopyIndex <<
        "'>"

        "<tr>"
        "<td>Maximum distance in the assembly graph (edges)"
        "<td class=centered><input type=text name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance <<
        "'>"

        "<tr>"
        "<td>Timeout for graph layout (seconds)"
        "<td class=centered><input type=text name=timeout size=8 style='text-align:center'"
        " value='" << timeout <<
        "'>";

    options.addFormRows(html);

    html <<
        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";



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
    if(startSegmentCopyIndex >= v.size() or v[startSegmentCopyIndex] == invalid<uint64_t> ) {
        html << "<br>Invalid segment copy index.<br>Valid segment copy indexes for this segment are";
        for(uint64_t segmentCopyIndex=0; segmentCopyIndex<v.size(); segmentCopyIndex++) {
            if(v[segmentCopyIndex] != invalid<uint64_t>) {
                html << " " << segmentCopyIndex;
            }
        }
        return;
    }
    const uint64_t startVertexId = v[startSegmentCopyIndex];
    SHASTA_ASSERT(startVertexId < snapshot.vertexVector.size());
}
