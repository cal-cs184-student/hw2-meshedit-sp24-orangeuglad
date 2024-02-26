#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

    /**
     * Evaluates one step of the de Casteljau's algorithm using the given points and
     * the scalar parameter t (class member).
     *
     * @param points A vector of points in 2D
     * @return A vector containing intermediate points or the final interpolated vector
     */
    std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const& points)
    {
        std::vector<Vector2D> intermediate_points;
        for (int n = 0; n < points.size() - 1; n++) {
            intermediate_points.push_back(Vector2D((1 - this->t) * points[n].x + this->t * points[n + 1].x, (1 - this->t) * points[n].y + this->t * points[n + 1].y));
        }
        return intermediate_points;
    }

    /**
     * Evaluates one step of the de Casteljau's algorithm using the given points and
     * the scalar parameter t (function parameter).
     *
     * @param points    A vector of points in 3D
     * @param t         Scalar interpolation parameter
     * @return A vector containing intermediate points or the final interpolated vector
     */
    std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const& points, double t) const
    {
        std::vector<Vector3D> intermediate_points;
        for (int n = 0; n < points.size() - 1; n++) {
            intermediate_points.push_back(Vector3D((1 - t) * points[n].x + t * points[n + 1].x, (1 - t) * points[n].y + t * points[n + 1].y, (1 - t) * points[n].z + t * points[n + 1].z));
        }
        return intermediate_points;
    }

    /**
     * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
     *
     * @param points    A vector of points in 3D
     * @param t         Scalar interpolation parameter
     * @return Final interpolated vector
     */
    Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const& points, double t) const
    {
        std::vector<Vector3D> interpolated_points = points;
        while (interpolated_points.size() != 1) {
            interpolated_points = evaluateStep(interpolated_points, t);
        }
        return interpolated_points[0];
    }

    /**
     * Evaluates the Bezier patch at parameter (u, v)
     *
     * @param u         Scalar interpolation parameter
     * @param v         Scalar interpolation parameter (along the other axis)
     * @return Final interpolated vector
     */
    Vector3D BezierPatch::evaluate(double u, double v) const
    {
        std::vector<Vector3D> interpolated_row_points;
        for (int n = 0; n < this->controlPoints.size(); n++) {
            interpolated_row_points.push_back(evaluate1D(this->controlPoints[n], u));
        }
        return evaluate1D(interpolated_row_points, v);
    }

    Vector3D Vertex::normal(void) const
    {
        // TODO Part 3.
        // Returns an approximate unit normal at this vertex, computed by
        // taking the area-weighted average of the normals of neighboring
        // triangles, then normalizing.
        // use HalfedgeCIter and not HalfedgeIter
        Vector3D normal(0, 0, 0);
        vector<VertexCIter> indices;
        VertexCIter v = this->halfedge()->vertex();
        indices.push_back(v);
        VertexCIter secondVertex;
        HalfedgeCIter h = this->halfedge();      // get the outgoing half-edge of the vertex
        do {
            HalfedgeCIter h_twin = h->twin(); // get the opposite half-edge
            VertexCIter v = h_twin->vertex(); // vertex is the 'source' of the half-edge, so
            // h->vertex() is v, whereas h_twin->vertex()
            // is the neighboring vertex
            indices.push_back(v);
            if (indices.size() == 3) {
                Vector3D edge1 = -(indices[1]->position - indices[0]->position);
                Vector3D edge2 = indices[2]->position - indices[0]->position;
                Vector3D faceNormal = cross(edge1, edge2);
                double faceArea = faceNormal.norm() / 2;
                normal += faceNormal * faceArea;

                indices.erase(indices.begin() + 1); //remove 2nd vertex
            }
            else if (indices.size() == 2) {
                secondVertex = v;
            }
            h = h_twin->next();               // move to the next outgoing half-edge of the vertex
        } while (h != v->halfedge());          // keep going until we are back where we were

        indices.push_back(secondVertex);
        Vector3D edge1 = indices[0]->position - indices[1]->position;
        Vector3D edge2 = indices[2]->position - indices[0]->position;
        Vector3D faceNormal = cross(edge1, edge2);
        double faceArea = faceNormal.norm() / 2;
        normal += faceNormal * faceArea;

        return normal / normal.norm();
    }

    EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0)
    {
        // TODO Part 4.
        // This method should flip the given edge and return an iterator to the flipped edge.
        if (e0->isBoundary())
            return e0;
        //Phase I
        HalfedgeIter h0 = e0->halfedge();
        HalfedgeIter h1 = h0->next();
        HalfedgeIter h2 = h1->next();
        HalfedgeIter h3 = h0->twin();
        HalfedgeIter h4 = h3->next();
        HalfedgeIter h5 = h4->next();
        HalfedgeIter h6 = h1->twin();
        HalfedgeIter h7 = h2->twin();
        HalfedgeIter h8 = h4->twin();
        HalfedgeIter h9 = h5->twin();

        VertexIter v0 = h0->vertex();
        VertexIter v1 = h1->vertex();
        VertexIter v2 = h2->vertex();
        VertexIter v3 = h5->vertex();

        EdgeIter e1 = h1->edge();
        EdgeIter e2 = h2->edge();
        EdgeIter e3 = h4->edge();
        EdgeIter e4 = h5->edge();

        FaceIter f0 = h0->face();
        FaceIter f1 = h3->face();

        //Phase II
        h0->next() = h1;
        h0->twin() = h3;
        h0->vertex() = v3;
        h0->edge() = e0;
        h0->face() = f0;

        h1->next() = h2;
        h1->twin() = h7;
        h1->vertex() = v2;
        h1->edge() = e2;
        h1->face() = f0;

        h2->next() = h0;
        h2->twin() = h8;
        h2->vertex() = v0;
        h2->edge() = e3;
        h2->face() = f0;

        h3->next() = h4;
        h3->twin() = h0;
        h3->vertex() = v2;
        h3->edge() = e0;
        h3->face() = f1;

        h4->next() = h5;
        h4->twin() = h9;
        h4->vertex() = v3;
        h4->edge() = e4;
        h4->face() = f1;

        h5->next() = h3;
        h5->twin() = h6;
        h5->vertex() = v1;
        h5->edge() = e1;
        h5->face() = f1;

        h6->next() = h6->next(); // didn’t change, but set it anyway!
        h6->twin() = h5;
        h6->vertex() = v2;
        h6->edge() = e1;
        h6->face() = h6->face(); // didn’t change, but set it anyway!

        h7->next() = h7->next(); // didn’t change, but set it anyway!
        h7->twin() = h1;
        h7->vertex() = v0;
        h7->edge() = e2;
        h7->face() = h7->face(); // didn’t change, but set it anyway!

        h8->next() = h8->next(); // didn’t change, but set it anyway!
        h8->twin() = h2;
        h8->vertex() = v3;
        h8->edge() = e3;
        h8->face() = h8->face(); // didn’t change, but set it anyway!

        h9->next() = h9->next(); // didn’t change, but set it anyway!
        h9->twin() = h4;
        h9->vertex() = v1;
        h9->edge() = e4;
        h9->face() = h9->face(); // didn’t change, but set it anyway!

        v0->halfedge() = h2;
        v1->halfedge() = h5;
        v2->halfedge() = h3;
        v3->halfedge() = h0;

        e0->halfedge() = h0;
        e1->halfedge() = h5;
        e2->halfedge() = h1;
        e3->halfedge() = h2;
        e4->halfedge() = h4;

        f0->halfedge() = h0;
        f1->halfedge() = h3;

        return e0;
    }

    VertexIter HalfedgeMesh::splitEdge(EdgeIter e0)
    {
        // TODO Part 5.
        // This method should split the given edge and return an iterator to the newly inserted vertex.
        // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
        if (e0->isBoundary())
            return e0->halfedge()->vertex();

        //Phase I
        HalfedgeIter h0 = e0->halfedge();
        HalfedgeIter h1 = h0->next();
        HalfedgeIter h2 = h1->next();
        HalfedgeIter h3 = h0->twin();
        HalfedgeIter h4 = h3->next();
        HalfedgeIter h5 = h4->next();
        HalfedgeIter h6 = h1->twin();
        HalfedgeIter h7 = h2->twin();
        HalfedgeIter h8 = h4->twin();
        HalfedgeIter h9 = h5->twin();

        VertexIter v0 = h0->vertex();
        VertexIter v1 = h1->vertex();
        VertexIter v2 = h2->vertex();
        VertexIter v3 = h5->vertex();

        EdgeIter e1 = h1->edge();
        EdgeIter e2 = h2->edge();
        EdgeIter e3 = h4->edge();
        EdgeIter e4 = h5->edge();

        FaceIter f0 = h0->face();
        FaceIter f1 = h3->face();

        //Phase Ib
        HalfedgeIter h10 = this->newHalfedge();
        HalfedgeIter h11 = this->newHalfedge();
        HalfedgeIter h12 = this->newHalfedge();
        HalfedgeIter h13 = this->newHalfedge();
        HalfedgeIter h14 = this->newHalfedge();
        HalfedgeIter h15 = this->newHalfedge();

        VertexIter v4 = this->newVertex();
        v4->position = (v0->position + v1->position) / 2;

        EdgeIter e5 = this->newEdge();
        EdgeIter e6 = this->newEdge();
        EdgeIter e7 = this->newEdge();
        e5->isNew = false;
        e6->isNew = true;
        e7->isNew = true;

        FaceIter f2 = this->newFace();
        FaceIter f3 = this->newFace();

        // Phase II

        v0->halfedge() = h12;
        v1->halfedge() = h1;
        v2->halfedge() = h2;
        v3->halfedge() = h5;
        v4->halfedge() = h0;

        h0->next() = h1;
        h0->twin() = h3;
        h0->vertex() = v4;
        h0->edge() = e0;
        h0->face() = f0;

        h1->next() = h2;
        h1->twin() = h6;
        h1->vertex() = v1;
        h1->edge() = e1;
        h1->face() = f0;

        h2->next() = h0;
        h2->twin() = h13;
        h2->vertex() = v2;
        h2->edge() = e7;
        h2->face() = f0;

        h3->next() = h4;
        h3->twin() = h0;
        h3->vertex() = v1;
        h3->edge() = e0;
        h3->face() = f1;

        h4->next() = h5;
        h4->twin() = h10;
        h4->vertex() = v4;
        h4->edge() = e6;
        h4->face() = f1;

        h5->next() = h3;
        h5->twin() = h9;
        h5->vertex() = v3;
        h5->edge() = e4;
        h5->face() = f1;

        h6->next() = h6->next();
        h6->twin() = h1;
        h6->vertex() = v2;
        h6->edge() = e1;
        h6->face() = h6->face();

        h7->next() = h7->next();
        h7->twin() = h14;
        h7->vertex() = v0;
        h7->edge() = e2;
        h7->face() = h7->face();

        h8->next() = h8->next();
        h8->twin() = h12;
        h8->vertex() = v3;
        h8->edge() = e3;
        h8->face() = h8->face();

        h9->next() = h9->next();
        h9->twin() = h5;
        h9->vertex() = v1;
        h9->edge() = e4;
        h9->face() = h9->face();

        h10->next() = h11;
        h10->twin() = h4;
        h10->vertex() = v3;
        h10->edge() = e6;
        h10->face() = f2;

        h11->next() = h12;
        h11->twin() = h15;
        h11->vertex() = v4;
        h11->edge() = e5;
        h11->face() = f2;

        h12->next() = h10;
        h12->twin() = h8;
        h12->vertex() = v0;
        h12->edge() = e3;
        h12->face() = f2;

        h13->next() = h14;
        h13->twin() = h2;
        h13->vertex() = v4;
        h13->edge() = e7;
        h13->face() = f3;

        h14->next() = h15;
        h14->twin() = h7;
        h14->vertex() = v2;
        h14->edge() = e2;
        h14->face() = f3;

        h15->next() = h13;
        h15->twin() = h11;
        h15->vertex() = v0;
        h15->edge() = e5;
        h15->face() = f3;

        e0->halfedge() = h0;
        e1->halfedge() = h1;
        e2->halfedge() = h7;
        e3->halfedge() = h8;
        e4->halfedge() = h5;
        e5->halfedge() = h11;
        e6->halfedge() = h4;
        e7->halfedge() = h2;

        f0->halfedge() = h0;
        f1->halfedge() = h3;
        f2->halfedge() = h10;
        f3->halfedge() = h13;

        return v4;
    }



    void MeshResampler::upsample(HalfedgeMesh& mesh)
    {
        // TODO Part 6.
        // This routine should increase the number of triangles in the mesh using Loop subdivision.
        // One possible solution is to break up the method as listed below.

        // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
        // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
        // a vertex of the original mesh.

        // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.

        // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
        // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
        // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
        // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)

        // 4. Flip any new edge that connects an old and new vertex.

        // 5. Copy the new vertex positions into final Vertex::position.

        for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
            int n = v->degree();
            double u;
            if (n == 3) {
                u = 3.0 / 16;
            }
            else {
                u = 3.0 / (8 * n); //not 8n because centroid is already divided by n
            }
            //Vector3D sum;
            //VertexCIter vc;
            //HalfedgeCIter h = v->halfedge();      // get the outgoing half-edge of the vertex
            //do {
            //    HalfedgeCIter h_twin = h->twin(); // get the opposite half-edge
            //    VertexCIter vc = h_twin->vertex(); // vertex is the 'source' of the half-edge, so
            //    // h->vertex() is v, whereas h_twin->vertex()
            //    // is the neighboring vertex
            //    sum += vc->position;
            //    h = h_twin->next();               // move to the next outgoing half-edge of the vertex
            //} while (h != vc->halfedge());          // keep going until we are back where we were

            Vector3D sum;
            VertexCIter vc;
            HalfedgeCIter h = v->getVertex()->halfedge();      // get the outgoing half-edge of the vertex
            do {
                HalfedgeCIter h_twin = h->twin(); // get the opposite half-edge
                VertexCIter v = h_twin->vertex(); // vertex is the 'source' of the half-edge, so
                sum = sum + v->position;
                h = h_twin->next();               // move to the next outgoing half-edge of the vertex
            } while (h != v->halfedge());          // keep going until we are back where we were


            v->newPosition = (1 - n * u) * v->position + u * sum;
            v->isNew = false;


        }

        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            Vector3D A = e->halfedge()->vertex()->position;
            Vector3D B = e->halfedge()->twin()->vertex()->position;
            Vector3D C = e->halfedge()->next()->next()->vertex()->position;
            Vector3D D = e->halfedge()->twin()->next()->next()->vertex()->position;
            e->newPosition = 3.0 / 8 * (A + B) + 1.0 / 8 * (C + D);
        }

        //HalfedgeMesh meshCopy = mesh;
        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            if (!(e->isNew)) {
                bool v0new = e->halfedge()->vertex()->isNew;
                bool v1new = e->halfedge()->twin()->vertex()->isNew;
                if (!v0new && !v1new) {
                    VertexIter v = mesh.splitEdge(e);
                    v->isNew = true;
                    v->newPosition = e->newPosition;
                }
            }
        }
        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            if (e->isNew) {
                bool v0new = e->halfedge()->vertex()->isNew;
                bool v1new = e->halfedge()->twin()->vertex()->isNew;
                if ((v0new && !v1new) || (!v0new && v1new)) {
                    EdgeIter eFlip = mesh.flipEdge(e);
                }
                std::cout << e->halfedge()->vertex()->position;
            }
        }

        for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
            v->position = v->newPosition;
        }
    }
}
