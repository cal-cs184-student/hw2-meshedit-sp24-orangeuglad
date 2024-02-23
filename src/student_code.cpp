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
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    std::vector<Vector2D> intermediate_points;
    for(int n = 0; n < points.size() - 1; n++) {
      intermediate_points.push_back(Vector2D((1-this->t)*points[n].x + this->t*points[n+1].x, (1-this->t)*points[n].y + this->t*points[n+1].y));
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
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    std::vector<Vector3D> intermediate_points;
    for(int n = 0; n < points.size() - 1; n++) {
      intermediate_points.push_back(Vector3D((1-t)*points[n].x + t*points[n+1].x, (1-t)*points[n].y + t*points[n+1].y, (1-t)*points[n].z + t*points[n+1].z));
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
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    std::vector<Vector3D> interpolated_points = points;
    while(interpolated_points.size() != 1) {
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
    for(int n = 0; n < this->controlPoints.size(); n++) {
      interpolated_row_points.push_back(evaluate1D(this->controlPoints[n], u));
    }
    return evaluate1D(interpolated_row_points, v);
  }

  Vector3D Vertex::normal( void ) const
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
        if(indices.size() == 3) {
          Vector3D edge1 = -(indices[1]->position - indices[0]->position);
          Vector3D edge2 = indices[2]->position - indices[0]->position;
          Vector3D faceNormal = cross(edge1, edge2);
          double faceArea = faceNormal.norm()/2;
          normal += faceNormal * faceArea;

          indices.erase(indices.begin() + 1); //remove 2nd vertex
        } else if (indices.size() == 2) {
          secondVertex = v;
        }
        h = h_twin->next();               // move to the next outgoing half-edge of the vertex
    } while(h != v->halfedge());          // keep going until we are back where we were

    indices.push_back(secondVertex);
    Vector3D edge1 = indices[1]->position - indices[0]->position;
    Vector3D edge2 = indices[2]->position - indices[0]->position;
    Vector3D faceNormal = cross(edge1, edge2);
    double faceArea = faceNormal.norm()/2;
    normal += faceNormal * faceArea;

    return normal / normal.norm();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
    return EdgeIter();
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    return VertexIter();
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
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

  }
}
