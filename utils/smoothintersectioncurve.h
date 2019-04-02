#ifndef SMOOTHINTERSECTIONCURVE_H
#define SMOOTHINTERSECTIONCURVE_H

#include "meshtypes.h"
#include "myutils.h"

class smoothIntersectionCurve
{
public:
    smoothIntersectionCurve();

    smoothIntersectionCurve(vector<vector<face::Pos<CMesh::FaceType>>> &icurves);
    smoothIntersectionCurve(vector<vector<face::Pos<CMesh::FaceType>>> &icurves, CMesh::ScalarType ithresh);

    // getting a regular smooth version of the part of the intersection curve starting in begin and ending in end
    void getSmoothVersion(face::Pos<CMesh::FaceType> & begin,face::Pos<CMesh::FaceType> & end,vector<CMesh::CoordType> & result);
    void getSmoothVersion(const face::Pos<CMesh::FaceType> & begin,const face::Pos<CMesh::FaceType> & end,vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double>> & result);
    void getSmoothVersion(const face::Pos<CMesh::FaceType> & begin,const face::Pos<CMesh::FaceType> & end,vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,size_t>> & result,bool invertorder=false);

    // creating renderer for the original curve
    void createRenderOriginal();

    // creating renderer for the smooth curve
    void createRenderSmooth();
    void createRenderSmooth2();

    void createSmoothIntersectionCurve();
    void setThreshold( CMesh::ScalarType ithresh);
    int getNumberEdges();
    int getNumberEdges_over_threshold(CMesh::ScalarType ithresh);
    void initialize( vector<vector<face::Pos<CMesh::FaceType>>> &icurves);
    void initialize( vector<vector<face::Pos<CMesh::FaceType>>> &icurves, CMesh::ScalarType ithresh);
    void getAverageEdgeLen(vector<CMesh::ScalarType> &edgelen);
    void processCurves(vector<CMesh::VertexPointer> indexesPos,vector<face::Pos<CMesh::FaceType>> & edgestosplit,vector<int> &timestosubdivide);

    void reset();

    void checkOriginalCurve();
public:

    CMesh renderOriginal,renderSmooth;

    vector<vector<face::Pos<CMesh::FaceType>>> original_curves;

    // an indexes group by each curve
    vector<set<size_t>> great_indexes;

    // For each curve a collection of smooth blocks. These are to debug. It stores the positions

    vector<vector<vector<CMesh::CoordType>>> blocks;

    // Each block stored as a set of indexes into renderSmooth

    vector<vector<vector<size_t>>> blocks_indexes;

    /* For each curve a collection of smooth blocks. Each blocks store the indexes of the halfedge's endpoints
    * that it belongs, i.e [a,b],  and the interpolation parameter u such at the point positions is
    * ua+(1-u)b
    */
    vector<vector<vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double>>>> blocks_interpolation;

    // Each blocks stored as a set of a position and a face into booleanmesh
    vector<vector<vector<std::pair<CMesh::CoordType,size_t>>>> blocks_interpolation_barycentric;

    CMesh::ScalarType threshold;

    // mapping between great indexes and blocks
    vector<unordered_map<size_t,size_t>> indexes_to_blocks;

};

#endif // SMOOTHINTERSECTIONCURVE_H
