#ifndef INTERSECTINGCURVES_H
#define INTERSECTINGCURVES_H

#include "myutils.h"
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>
#include<vcg/complex/algorithms/isotropic_remeshing.h>
#include<vcg/complex/algorithms/smooth.h>
#include <forward_list>

class intersectingCurves
{
public:
    //constructor from two meshes and corresponding loops

    //intersectingCurves(CMesh& ml,CMesh& mr,vector<pair<vector<CMesh::VertexPointer>,vector<CMesh::VertexPointer>>>& loops);

    intersectingCurves(CMesh& ml,CMesh& mr,vector<pair<vector<CMesh::VertexPointer>,vector<CMesh::VertexPointer>>>& loops,
                       vector<pair<vector<vector<Point3f>>,vector<vector<Point3f>>>>& mapSF);


    //Default constructor

    intersectingCurves();

    void reset();

    void updatethresholds();

    void resetPairPositionLoops(int index=-1);

    void redefine(CMesh& ml,CMesh& mr,vector<pair<vector<CMesh::VertexPointer>,vector<CMesh::VertexPointer>>>& loops,
                  vector<pair<vector<vector<Point3f>>,vector<vector<Point3f>>>>& mapSF);

    
    /*
    //Set of functions to perform the merged
    */

    //Returns the vertex in lopp which is nearest to x under Euclidean distance.		
    vcg::Point3f nearest(const std::forward_list<Point3f>& loop, vcg::Point3f& x);
    vcg::Point3f nearest(const std::forward_list<Point3f*>& loop, vcg::Point3f& x);

    int nearestIndex(const std::forward_list<Point3f>& loop, vcg::Point3f& x);
	
    //Envolve the two input loops both closer
    void moreNear(std::forward_list<Point3f>& left,std::forward_list<Point3f>& right,
                  vector<vector<Point3f>>& toprojectleft,
                  vector<vector<Point3f>>& toprojectright);

    //Refine loops. Return the final number of elements of the loop

    int refineLoop(std::forward_list<Point3f>& loop,bool iscircular=true);

    void mergeLoops();

    void setThrMax(float max);

    int sizefowardlist(const std::forward_list<Point3f>& loop);

    std::unordered_map<int,forward_list<Point3f>::iterator> fowardlistToIndexes(std::forward_list<Point3f>& loop);

    void addloopToMesh(const std::forward_list<Point3f>& loop,bool iscircular=true);
    void addloopToMesh(const vector<CMesh::VertexPointer>& loop,bool iscircular=true);

    bool improveloop(vector<CMesh::VertexPointer>& loop, CMesh& m);

    float averageEdgeLenghtLoop(const std::forward_list<Point3f>& loop,bool iscircular=true);

    void improveInitialLoop(std::forward_list<Point3f>& loop,bool iscircular=true);

    // fill fisrt loop with the intersecting loop(need it be computed before) and store the new vertexpointer added to the mesh
    // in the second argument
    void fillloops(int index,vector<CMesh::VertexPointer>& loop,vector<CMesh::VertexPointer>& newloop,CMesh& m);

    void utilfillloops(vector<CMesh::VertexPointer>& loop,vector<CMesh::VertexPointer>& newloop,CMesh& m,bool inmesh);

    void doIsoRemeshing(CMesh & m,bool selectedFaces=true);

    //methods for update properties from outside

    // delete vertices in intersection curve i
    // indexes are the vertex's indexes from the curve i to be delete
    void deleteVerticesInIntersectionCurve(int indexCurve,set<int>& indexesToDelete);

public:
    //An edgemesh vector that store each intersection curve for each pair of loops
    CMesh IntersectedCurves;

    vector<vector<Point3f>> intersectingLoop;
    vector<pair<vector<CMesh::VertexPointer>,vector<CMesh::VertexPointer>>> intersectingLoopPointers;

    // Pairs of matched loops to be merged

    vector<pair<vector<CMesh::VertexPointer>,vector<CMesh::VertexPointer>>> pairLoops;

    // Pairs of 3D Points corresponding to the pairs of loops

    vector<pair<std::forward_list<Point3f>,std::forward_list<Point3f>>> pairPositionLoops;
		
    // Selected Faces to reproject

    vector<pair<vector<vector<Point3f>>,vector<vector<Point3f>>>> maptoproject;


    //copy of the initial two meshes

    CMesh* cloneml;
    CMesh* clonemr;

    float thrZero;
    float thrMax;
    bool merged;
    int maxIteration;
    int numberIteration;
};

#endif // INTERSECTINGCURVES_H
