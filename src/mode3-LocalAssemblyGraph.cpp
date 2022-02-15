// Shasta.
#include "mode3-LocalAssemblyGraph.hpp"
#include "computeLayout.hpp"
#include "writeGraph.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <map>
#include <queue>
#include "tuple.hpp"



// Create the LocalAssemblyGraph using a BFS
// that starts at the specified vertex and moves away
// (in both directions) up to the specified distance
mode3::LocalAssemblyGraph::LocalAssemblyGraph(
    const AssemblyGraph& assemblyGraph,
    uint64_t startSegmentId,
    uint64_t maxDistance)
{
    LocalAssemblyGraph& localAssemblyGraph= *this;

    // The BFS queue.
    std::queue<uint64_t> q;

    // Map segments in the AssemblyGraph to vertices in
    // the LocalAssemblyGraph.
    std::map<uint64_t, vertex_descriptor> segmentMap;

    // Initialize the BFS.
    q.push(startSegmentId);
    const vertex_descriptor vStart = addVertex(startSegmentId, 0, assemblyGraph.paths[startSegmentId]);
    segmentMap.insert(make_pair(startSegmentId, vStart));



    // BFS.
    while(not q.empty()) {

        // Dequeue a segment.
        const uint64_t segmentId0 = q.front();
        q.pop();
        const vertex_descriptor v0 = segmentMap[segmentId0];
        const uint64_t distance0 = localAssemblyGraph[v0].distance;
        const uint64_t distance1 = distance0 + 1;

        // Loop over children.
        for(const uint64_t linkId: assemblyGraph.linksBySource[segmentId0]) {
            const Link& link = assemblyGraph.links[linkId];
            const uint64_t segmentId1 = link.segmentId1;
            if(segmentMap.find(segmentId1) != segmentMap.end()) {
                // We already encountered this segment.
                continue;
            }
            const vertex_descriptor v1 = addVertex(segmentId1, distance1, assemblyGraph.paths[segmentId1]);
            segmentMap.insert(make_pair(segmentId1, v1));
            if(distance1 < maxDistance) {
                q.push(segmentId1);
            }
        }

        // Loop over parents.
        for(const uint64_t linkId: assemblyGraph.linksByTarget[segmentId0]) {
            const Link& link = assemblyGraph.links[linkId];
            const uint64_t segmentId1 = link.segmentId0;
            if(segmentMap.find(segmentId1) != segmentMap.end()) {
                // We already encountered this segment.
                continue;
            }
            const vertex_descriptor v1 = addVertex(segmentId1, distance1, assemblyGraph.paths[segmentId1]);
            segmentMap.insert(make_pair(segmentId1, v1));
            if(distance1 < maxDistance) {
                q.push(segmentId1);
            }
        }
    }



    // Add the edges.
    for(const auto& p: segmentMap) {
        const uint64_t segmentId0 = p.first;
        const vertex_descriptor v0 = p.second;

        for(const uint64_t linkId: assemblyGraph.linksBySource[segmentId0]) {
            const Link& link = assemblyGraph.links[linkId];
            const uint64_t segmentId1 = link.segmentId1;
            const auto it1 = segmentMap.find(segmentId1);
            if(it1 == segmentMap.end()) {
                continue;
            }
            const vertex_descriptor v1 = it1->second;
            boost::add_edge(v0, v1, LocalAssemblyGraphEdge(link.coverage), localAssemblyGraph);
        }
    }
}



mode3::LocalAssemblyGraphVertex::LocalAssemblyGraphVertex(
    uint64_t segmentId,
    uint64_t distance,
    const span<const MarkerGraphEdgeInfo> pathArgument) :
    segmentId(segmentId),
    distance(distance)
{
    copy(pathArgument.begin(), pathArgument.end(), back_inserter(path));
}



mode3::LocalAssemblyGraphVertex::LocalAssemblyGraphVertex() :
    segmentId(0),
    distance(0)
{
}



mode3::LocalAssemblyGraph::vertex_descriptor mode3::LocalAssemblyGraph::addVertex(
    uint64_t segmentId,
    uint64_t distance,
    const span<const MarkerGraphEdgeInfo> path)
{
    return add_vertex(LocalAssemblyGraphVertex(segmentId, distance, path), *this);
}




void mode3::LocalAssemblyGraph::writeGraphviz(const string& fileName) const
{
    ofstream s(fileName);
    writeGraphviz(s);
}



void mode3::LocalAssemblyGraph::writeGraphviz(ostream& s) const
{
    const LocalAssemblyGraph localAssemblyGraph = *this;

    s << "digraph LocalAssemblyGraph {\n";

    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        s << localAssemblyGraph[v].segmentId << ";\n";
    }

    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const vertex_descriptor v0 = source(e, localAssemblyGraph);
        const vertex_descriptor v1 = target(e, localAssemblyGraph);
        s << localAssemblyGraph[v0].segmentId << "->";
        s << localAssemblyGraph[v1].segmentId <<
            "[penwidth=" << 0.2 * double(localAssemblyGraph[e].coverage) << "];\n";
    }

    s << "}\n";
}



// Bandage-style svg output, done using sfdp.
void mode3::LocalAssemblyGraph::writeSvg1(const string& fileName, uint64_t sizePixels)
{
    ofstream svg(fileName);
    writeSvg1(svg, sizePixels);
}
void mode3::LocalAssemblyGraph::writeSvg1(ostream& svg, uint64_t sizePixels)
{
    const double lengthPerMarker = 2.;
    const string svgId = "LocalAssemblyGraph";
    const string internalEdgeColor = "green";
    const double internalEdgeThickness = 10;
    const string externalEdgeColor = "black";
    const double externalEdgeThickness = 3;

    const LocalAssemblyGraph& localAssemblyGraph = *this;

    // We use an auxiliary graph which has 2 or more vertices for
    // each vertex of the LocalAssemblyGraph.
    // The number of auxiliary graph vertices
    // corresponding to each LocalAssemblyGraph vertex
    // is determined by the path length.
    using G = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;
    G g;
    std::map<vertex_descriptor, vector<G::vertex_descriptor> > vertexMap;
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {

        // Decide how many auxiliary vertices we want.
        const uint64_t pathLength = localAssemblyGraph[v].path.size();
        const uint64_t auxiliaryCount = max(uint64_t(2),
            uint64_t(std::round(double(pathLength) * lengthPerMarker)));

        // Add the auxiliary vertices.
        vector<G::vertex_descriptor>& auxiliaryVertices = vertexMap[v];
        for(uint64_t i=0; i<auxiliaryCount; i++) {
            auxiliaryVertices.push_back(boost::add_vertex(g));
        }

        // Add the edges between these auxiliary vertices.
        for(uint64_t i=1; i<auxiliaryCount; i++) {
            boost::add_edge(auxiliaryVertices[i-1], auxiliaryVertices[i], g);
        }
    }

    // Add auxiliary graph edges between vertices corresponding to different
    // LocalAssemblyGraph vertices.
    std::map<LocalAssemblyGraph::edge_descriptor, G::edge_descriptor> edgeMap;
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const vertex_descriptor v0 = source(e, localAssemblyGraph);
        const vertex_descriptor v1 = target(e, localAssemblyGraph);
        G::edge_descriptor ee;
        tie(ee, ignore) = add_edge(
            vertexMap[v0].back(),
            vertexMap[v1].front(),
            g);
        edgeMap.insert(make_pair(e, ee));
    }

    // Compute the layout of the auxiliary graph using graphviz.
    std::map<G::vertex_descriptor, array<double, 2> > positionMap;
    if(computeLayout(g, "sfdp", 30., positionMap, "-Gsmoothing=avg_dist") != ComputeLayoutReturnCode::Success) {
        throw runtime_error("Graph layout failed.");
    }

    // Compute the view box.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::min();
    BGL_FORALL_VERTICES_T(u, g, G) {
        const auto& position = positionMap[u];

        // Update the view box to include this vertex.
        xMin = min(xMin, position[0]);
        xMax = max(xMax, position[0]);
        yMin = min(yMin, position[1]);
        yMax = max(yMax, position[1]);
    }
    // Add a bit of extra space.
    const double dx = 0.05 * (xMax - xMin);
    const double dy = 0.05 * (yMax - yMin);
    xMin -= dx;
    xMax += dx;
    yMin -= dy;
    yMax += dy;


    // Begin the svg.
    svg << "<svg id='" << svgId << "' width='" << sizePixels << "' height='" << sizePixels <<
        "' viewbox='" << xMin << " " << yMin << " " << xMax-xMin << " " << yMax-yMin <<
        "'>\n";



    // Write the edges between segments first, so they don't overwrite the segments.
    // Each of these edges corresponds to a Link.
    svg << "<g id='" << svgId << "-links'>\n";
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {

        // Access the LocalAssemblyGraph vertices corresponding to
        // the two segments of this Link.
        const vertex_descriptor v1 = source(e, localAssemblyGraph);
        const vertex_descriptor v2 = target(e, localAssemblyGraph);

        // Get the vectors of auxiliary graph vertices for each of them.
        // These have size at least 2.
        const vector<G::vertex_descriptor>& V1 = vertexMap[v1];
        const vector<G::vertex_descriptor>& V2 = vertexMap[v2];

        // Access the terminal vertices.
        // These will be the start  and end points of the cubic spline.
        const G::vertex_descriptor u1 = V1.back();
        const G::vertex_descriptor u2 = V2.front();
        const auto& p1 = positionMap[u1];
        const auto& p2 = positionMap[u2];

        // Access the vertices adjacent to the terminal ones.
        // These will be used to construct the control points of
        // the cubic spline, by reflection.
        const G::vertex_descriptor uu1 = V1[V1.size() -2];
        const G::vertex_descriptor uu2 = V2[1];
        const auto& pp1 = positionMap[uu1];
        const auto& pp2 = positionMap[uu2];

        // Create the control points for the spline.
        array<double, 2> q1, q2;
        for(uint64_t i=0; i<2; i++) {
            q1[i] = p1[i] + (p1[i] - pp1[i]);
            q2[i] = p2[i] + (p2[i] - pp2[i]);
        }

        svg <<
            "<path d='M " <<
            p1[0] << " " << p1[1] << ", C " <<
            q1[0] << " " << q1[1] << ", " <<
            q2[0] << " " << q2[1] << ", " <<
            p2[0] << " " << p2[1] << "'" <<
            " stroke='" << externalEdgeColor << "'"
            " stroke-width='" << externalEdgeThickness << "px'"
            " stroke-linecap='round'"
            " fill='transparent'"
            " vector-effect='non-scaling-stroke'"
            "/>\n";

    }
    svg << "</g>\n";



    // Write the edges internal to segments.
    svg << "<g id='" << svgId << "-segments'>\n";
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const vector<G::vertex_descriptor>& auxiliaryVertices = vertexMap[v];

        for(uint64_t i=1; i<auxiliaryVertices.size(); i++) {

            // Get vertex positions.
            const G::vertex_descriptor v1 = auxiliaryVertices[i-1];
            const G::vertex_descriptor v2 = auxiliaryVertices[i];
            const auto& position1 = positionMap[v1];
            const auto& position2 = positionMap[v2];

            svg << "<line x1='" << position1[0] << "' y1='" << position1[1] <<
                "' x2='" << position2[0] << "' y2='" << position2[1] << "'";

            svg <<
                " stroke='" << internalEdgeColor << "'"
                " stroke-width='" << internalEdgeThickness << "px'"
                " stroke-linecap='round'"
                " vector-effect='non-scaling-stroke'"
                ">";

            svg << "</line>\n";
        }
    }
    svg << "</g>\n";



    // End the svg.
    svg << "</svg>\n";
}

