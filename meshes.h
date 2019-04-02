#ifndef MESHES_H
#define MESHES_H

#include <QString>
#include <QDateTime>

/// vcg imports
#include "meshtypes.h"

//#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/refine_loop.h>
#include <vcg/complex/algorithms/update/selection.h>
#include <vcg/complex/algorithms/curve_on_manifold.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/complex/algorithms/parametrization/poisson_solver.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include<vcg/complex/algorithms/isotropic_remeshing.h>
/* include the support for half edges */
//#include <vcg/complex/algorithms/update/halfedge_indexed.h>

/// wrapper imports and exports
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <chrono>
#include <thread>
#include <boost/config.hpp>
#include <boost/bimap.hpp>

//#include <wrap/io_trimesh/export_ply.h>
//#include <wrap/io_trimesh/export_obj.h>
/// Others

//#include "selectedfaces.h"
//#include "selectedfacesp.h"
//#include "intersectingcurves.h"
//#include "patches.h"
#include "patch3d.h"
#include <quad_tracer.h>
#include "PatternsTakayama/patchg.h"
#include "orient_faces.h"
#include "utils/smoothintersectioncurve.h"
#include "utils/smooth_intersection.h"
#include <vcg/complex/algorithms/update/quality.h>
using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds

class EdgeSplitFixedList_BooleanMesh
{
  typedef typename face::Pos<CMesh::FaceType> PosType;
  public:
      CMesh::PerFaceAttributeHandle<int> meshnumber;
      CMesh::ScalarType min_edgelen;
      EdgeSplitFixedList_BooleanMesh(){min_edgelen=0;}
      EdgeSplitFixedList_BooleanMesh(CMesh &imesh,CMesh::ScalarType threshold=0) {
          meshnumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(imesh,std::string("MeshNumber"));
          min_edgelen=threshold;
      }
      bool operator()(PosType &ep)
      {
          CMesh::FacePointer facep=ep.F();
          CMesh::FacePointer faceopp=ep.FFlip();
          bool isIntersectionEdge= meshnumber[facep]!=meshnumber[faceopp];
          CMesh::ScalarType distance=vcg::Distance(ep.V()->P(),ep.VFlip()->P());
          if(isIntersectionEdge && distance>min_edgelen) // if this edge belongs to the intersection curve and is not too small
              return true;
          else
              return false;
      }
};


class meshes
{
public:
    meshes();
    vector<vector<face::Pos<CMesh::FaceType>>> getIntersectionCurve();
    vector<face::Pos<CMesh::FaceType>> getIntersectionCurveFromPos(face::Pos<CMesh::FaceType>& startPos);
    void computeSelectedRegions(vector<vector<face::Pos<CMesh::FaceType>>>& curves,set<int> &  selectedPartitions1, set<int> & selectedPartitions2);
    void computePatches(std::set<int> &sc1,std::set<int> &sc2);
    void processBooleanMesh(Eigen::VectorXi &criteria);
    void collapseSmallTrianglesIntersectionCurve1(vector<vector<face::Pos<CMesh::FaceType>>> &edges);
    void collapseSmallTrianglesIntersectionCurve2(vector<vector<face::Pos<CMesh::FaceType>>> &edges);
    void testPatches();
    bool makePatchesEven(vector<tempora::Patch3D<CMesh,PMesh>> & listpatches);
    bool makePatchesEven2(vector<tempora::Patch3D<CMesh,PMesh>> & listpatches);
    void updateBooleanIMarkSelectedPartitions(std::set<int> &sc1,std::set<int> &sc2);        
    void exportBooleanTagged();
    void splitIntersectionCurve();
    void mergeEndpointsPatches();
    void collapsingAfterSplit(vector<CMesh::ScalarType> & aedgelen);
    void conditioningIntersection();
    void selectFacesIntersection(bool selone,bool seltwo,int numberNeighbors);
    void processall();


    /* VISUALIZATION METHODS
     */
    void reset();
    int loadmesh1(QString fileName);
    int loadmesh2(QString fileName);
    void initMesh(PMesh& mesh);
    void initMesh1();
    void initMesh2();
    void changeSeeSelectedQuads();
    vector<pair<size_t,PMesh::CoordType>> getIntersectionCurveQuad();

    /* NECESSARY METHODS FOR VISUALIZATION
     */
    void triangulate(CMesh& tmesh,PMesh& mesh);
    void triangulate1();
    void triangulate2();
    /* Internal function. Select all components from a mesh A, inner or outer to the mesh B according to parameter "inner".
     * This method assumes that a selection over the meshA is already done.
     * @gm is a grid of the meshB
     */
    void selectComponents(bool inner,PMesh &meshA, PMesh &meshB, GridStaticPtr<PMesh::FaceType, PMesh::ScalarType>& gm);
    // Select all components from mesh 2 inner ou outer to the mesh 1 according to the argument
    void selectComponents1(bool inner);
    // Select all components from mesh 1 inner ou outer to the mesh 2 according to the argument
    void selectComponents2(bool inner);
    void visualizePatches( QuadMeshTracer<PMesh>& tracer, PMesh &mesh, CMesh &tmesh);
    void visualizePatches1();
    void visualizePatches2();

    /*UTILS
     */
    int convertIndexT2Q(int tindex);

public:
    vector<Segment3f> edges;

    /// the active mesh instance
    PMesh mesh1;
    CMesh tmesh1;

    /// the second mesh instance
    PMesh mesh2;
    CMesh tmesh2;

    /// the triangle boolean result
    CMesh tboolean;

    /// The edgemesh of the patch partitions (for visualization) of each quad mesh
    CMesh borderpatches1;
    CMesh borderpatches2;

    ///Uniform Grids for the quad meshes
    GridStaticPtr<PMesh::FaceType, PMesh::ScalarType> GRIDQ1,GRIDQ2;

    //Quad layouts tracers
    QuadMeshTracer<PMesh> PTracer1;
    QuadMeshTracer<PMesh> PTracer2;

    //intersectingCurves intcurves;
    smoothIntersectionCurve intcurves;

    // collection of intersection curves in TriBoolean store as a vector of half-edges
    vector<vector<face::Pos<CMesh::FaceType>>> CurvesInBoolean;

    // To store all patches covering the intersection curves
    vector<tempora::Patch3D<CMesh,PMesh>> CollectionPatches;

    // For visualization of the patches
    //Edgemesh
    CMesh seepatches;
    // quad mesh resulting
    PMesh quadboolean;

    bool stopAnimation;
    bool flagEdgeMesh;
    bool flagIntersectedCurve;
    bool flagSmoothIntersectedCurve;
    bool flagTmesh1;
    bool flagTmesh2;
    bool flagTriangleBoolean;
    bool flagQuadBoolean;
    bool flagSeeSelectedQuad;
    bool flagExportIntermediateMeshes;
    bool flagGenerateReport;
    int numbersubdivisions;
};

#endif // MESHES_H
