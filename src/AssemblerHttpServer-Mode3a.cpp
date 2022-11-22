// Shasta.
#include "Assembler.hpp"
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

    uint64_t startVertexId;
    const bool startVertexIdIsPresent = getParameterValue(request, "startVertexId", startVertexId);

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
        "<td>Start vertex"
        "<td class=centered><input type=text required name=startVertexId size=8 style='text-align:center'"
        " value='" << (startVertexIdIsPresent ? to_string(startVertexId) : "") <<
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
}
