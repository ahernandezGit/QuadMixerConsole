#ifndef MYUTILS_H
#define MYUTILS_H

#include "meshtypes.h"
#include <Eigen/Geometry>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/refine_loop.h>
#include <vcg/complex/algorithms/edge_collapse.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/inside.h>
#include <vcg/space/intersection2.h>
#include <vcg/space/triangle3.h>
#include <wrap/io_trimesh/export_ply.h>
#include <typeinfo>
#include <fstream>
namespace utility{

    typedef typename tri::BasicVertexPair<CMesh::VertexType> vertexPairCollapse;
    typedef typename face::Pos<CMesh::FaceType> PosType;

    template <typename ScalarType>
    class ParametricCurve{
    public:
        //Stores original points controls
        vector<vcg::Point3<ScalarType>> _points;
        int numberSegments=0;
        double totalLength=0;
        //Stores length of each segment
        vector<double> segmentLengths;

        //Store the parameter in [0,1] of the original points
        vector<double> parameter;
    public:
        ParametricCurve(vector<Point3<ScalarType>>& points){
            if(points.size()<1){
                cout<<"curve is just a point"<<endl;
                return;
            }
           _points=points;
           numberSegments=points.size()-1;
           segmentLengths.clear();
           parameter.clear();
           // Total length of path
           totalLength = 0.0;
           for (int segIdx=0; segIdx <numberSegments; segIdx++) {
               // calculate length of a segment
               double segLen=vcg::Distance(points[segIdx],points[segIdx+1]);
               totalLength += segLen;
               segmentLengths.push_back(segLen);
           }
           assert(numberSegments==segmentLengths.size());
           double currentLenght=0;
           parameter.push_back(0.0);
           for (int i=0; i<numberSegments-1; i++) {
              currentLenght+=segmentLengths[i];
              parameter.push_back(currentLenght/totalLength);
           }
           parameter.push_back(1.0);
           assert((numberSegments+1)==parameter.size());
        }
        /*
                Function returns segment index of point at distance len
             *	len length since start of path that point end on
             *	returns segment index
        */
        int getSegmentIndexAtArcLength(double& len){
                double currentLength = 0.0;
                for (int i=0; i <numberSegments; i++) {
                    currentLength += segmentLengths[i];
                    if (len <= currentLength) {
                        len=segmentLengths[i]-(currentLength-len);
                        return i;
                    }
                }
                return numberSegments-1;
        }

        // u is the parameter between 0 and 1.
        Point3<ScalarType> interpolatePoints(const Point3<ScalarType>& P0, const Point3<ScalarType>& P1, double u)
        {
             assert(u>=0 && u<=1);
             Point3<ScalarType> cP0,cP1;
             cP0=P0;
             cP1=P1;
             double pu=1-u;
             Point3<ScalarType> result=cP0*pu+cP1*u;
             return result;
        }
        /*
                Function returns point on path that is at length len
                len length since start the point ends on
                returns (x,y,z) point at length len
        */
        Point3<ScalarType> getPointAtArcLength(double len){
              assert(len>=0);
              if(len>totalLength)
                  return _points[numberSegments];
              else{
                int segIdx = getSegmentIndexAtArcLength(len);
                double u=len/segmentLengths[segIdx];
                //cout<<"u after "<<u<<endl;
                return interpolatePoints(_points[segIdx],_points[segIdx+1],u);
              }
        }

        double getParameterAtControlPoint(int i){
              assert(i>=0 && i<=numberSegments);
              return parameter[i];
        }
        Point3<ScalarType> interpolateForT(double t){
           if (t < 0.0) t=0.0;
           else if ( t>1.0) t=1.0;
           double len=t*totalLength;
           return getPointAtArcLength(len);
        }
    };

    struct patchSide{
        bool iscurveSide=false;
        int size=0;
        vector<Point3f> points;
        face::Pos<PMesh::FaceType> begin;
        face::Pos<PMesh::FaceType> end;
        //Default constructor
        patchSide(){
            points.clear();
            begin.SetNull();
            end.SetNull();
        }
        // Constructor from vector of pos
        patchSide(vector<face::Pos<PMesh::FaceType>>& pospoints,bool iscurve=false){
            for(int i=0;i<pospoints.size();i++)
                points.push_back(pospoints[i].V()->P());
            begin=pospoints[0];
            end=pospoints[pospoints.size()-1];
            iscurveSide=iscurve;
            size=pospoints.size();
        }
        // Constructor from initial and final pos
        patchSide(face::Pos<PMesh::FaceType>& beginPos,face::Pos<PMesh::FaceType>& endPos, bool iscurve=false){

        }
    };

    struct UVonQuad{
       int origin=-1;
       float onNext=0.0;
       float onPrev=0.0;
       UVonQuad(){}
       UVonQuad(int a,float b,float c){
        origin=a;
        onNext=b;
        onPrev=c;
       }
       void set(int a,float b,float c){
           origin=a;
           onNext=b;
           onPrev=c;
       }
    };

    void saveAsCurvePly(vector<face::Pos<CMesh::FaceType>> &curve,string  filepath);
    void saveAsCurvePly(vector<CMesh::CoordType> &curve,string  filepath,bool iscircular=true);
    void saveMesh(CMesh &m,string filepath);
    void ComputeRealRootsDegree2(double &root1, double &root2, double a, double b, double c,int & valid);

    //Verify is a pair of two vector are the same at least of a index permutation
    bool ispermuted(const vector<int>& list1,const vector<int>& list2, vector<pair<int,int>>& correspondence);
    //Correct curve of points smoothing  corner
    void correctCornersFromCurve(vector<Point3f>& curve, bool iscircular=true);

    // return blocks of consecutive indexes inside a set of indexes from a array of nd elements
    // blocks is the final set of blocks
    // markindexes is the input set of indexes
    // nd is the size of the array of elements
    void groupIndexes(vector<vector<int>>& blocks,set<int>& markindexes,int nd,bool iscircular=true);

    //Check is the array vec has more consecutive increments than decrements
    bool moreIncrementsArray(vector<int>& vec);

    // Check if a vertex is endpoint of a crease edge
    // Need to precompute the crease edges first
    bool checkVertexIsCrease(const CMesh::VertexType* v);

    // Check if vertex is corner
    // Don't need to precompute the crease edges. Check if two adjacent faces containing the common vertex
    // has a dihedral angle near to the angle threshold;
    bool checkVertexIsCorner(const CMesh::VertexType* v,float anglethr);

    // Get the average edge length of a mesh
    // Edge average approximation iterating on selected faces
    // If not faces selected then return the whole average dge lenght
    float edgeAverage(CMesh& m,bool selectedFaces=false);

    // True if some edge from the selected faces has length less than the threshold
    // If not faces were selected then run over all the edges
    bool someEdgeSmall(CMesh& m,float threshold, bool selectedFaces=false);

    // True if some edge from the selected faces has length bigger than the threshold
    // If not faces were selected then run over all the edges
    bool someEdgeBigger(CMesh& m,float threshold, bool selectedFaces=false);

    vector<CMesh::CoordType> getRegularDivisionsAsBarycentricFromPos(face::Pos<CMesh::FaceType> &pos,int subdivisionsPos);

    //Print a point
    void printPoint3f(const vcg::Point3f& p,string message);

    // Compute the aspect ratio of a triangle
    float aspectRatioTriangleFromFace(CMesh::FacePointer f);
    float aspectRatioTriangle(const Point3f& v0,const Point3f& v1,const Point3f& v2);

    // signed area of a triangle
    float signedAreaTriangle(const Point3f& v0,const Point3f& v1,const Point3f& v2);

    //Split selected faces of the mesh with the midpoint method. Each edge of lenght less than threshold
    void splitMesh(CMesh& m, float threshold,bool selectedFaces=false);

    //Split selected faces of the mesh with the midpoint method. Each edge of lenght less than threshold
    void splitMeshP(PMesh& m, float threshold,bool selectedFaces=false);

    //Split the incident edges of all vertex borders of the mesh with the midpoint method.
    void splitMeshStrip(CMesh& m, float threshold,bool selectedFaces=false);

    //Split the incident edges of all vertex borders except the border edges.
    void splitMeshStripNonBorder(CMesh& m, float threshold,bool selectedFaces);

    //Collapse edges of lenght less than threshold
    void collapseMesh(CMesh& m, float threshold,bool selectedFaces=false);

    //Collapse edges of lenght less than threshold in a selected border strip submesh
    // This method use fortely that m have a coherent border strip mesh selected
    void collapseStripMesh(CMesh& m, float threshold,bool selectedFaces=false);

    // Flip edges of the selected Faces of a mesh following the criterium from tha paper
    bool flipEdgesStripMesh(CMesh& m,bool selectedFaces=false);

    //Fast check  if a point is inside a mesh
    bool isInsideMesh(const CMesh& m, const Point3f& p);

    bool isInsideMesh(CMesh::FacePointer vfp, const Point3f& p);

    // return if the point (testx,testy) is inside the polygon (2D). Is the point is on the border so p is considered as inside
    bool pointIsInsidePoly(vector<CMesh::CoordType> &polygon,CMesh::CoordType point);

    // compute the reflexion of the point P over the line AB, the result is stered in RP
    void computeReflexion(CMesh::CoordType RP,CMesh::CoordType P,CMesh::CoordType A,CMesh::CoordType B);

    void correctPlanarBox(vcg::Box3<float>& auxiliarBox,float offset);

    // Project a point into a colection of faces of a mesh
    void easyProjectToMesh(const std::vector<CMesh::FacePointer>& faces,Point3f& p);

    void easyProjectToMesh(const vector<vector<Point3f>>& faces,Point3f& p);

    // Give a score to a direction change of a point in a boundary loop
    // error is a linear function over the angle between this vectors
    // Based on the discrete version of the GrisMesh paper of Igarashi
    // if angle is 0 so error is 3;
    // if angle is 60º so error is 2;
    // if angle is 120º so error is 1;
    // if angle is 180º so error is 0;
    float errorAB(const Point3f& A,const Point3f& B);

    // check orientation
    bool isCounterClockWise(vector<vcg::face::Pos<CMesh::FaceType>> & polygon);
    bool isCounterClockWise(vector<Point3f> & polygon);
    bool isCounterClockWise(vector<Point2f> & polygon);

    /* Convert a 2D polygon with holes into a single boundary polygon by joing the holes to a outer boundary. Recursively this function
     * joins the closest hole to the outer boundary with a valid diagonal.
     * This function is a brute force implementation onn O(n^3)
     * @param outerB represents the exterior boundary of the polygon with holes
     * @param innerB represents the collection of holes of the polygon with holes
     * @param singleB is the output of the method. It is more convenient for us to return a 3D polygon on the plane XY
    */
    void toSingleBoundary(vector<Point2f>& outerB,vector<vector<Point2f>> &innerB, vector<Point3f>& singleB);

    /* The same as the previous function onlywith a set of feature vertices os the input outer bundary
     * This function is a brute force implementation onn O(n^3)
     * @param outerB represents the exterior boundary of the polygon with holes
     * @param innerB represents the collection of holes of the polygon with holes
     * @param singleB is the output of the method. It is more convenient for us to return a 3D polygon on the plane XY
     * @param featureVertex store the feature vertices of the outer boundary. At the end os the algorithm this argument is update
     * to correctly store the feature vertex of the output
    */

    void toSingleBoundary(vector<Point2f>& outerB,vector<vector<Point2f>> &innerB, vector<Point3f>& singleB,set<int> &featureVertex);

    //This splits regularly a segment starting on A and finising on B, numberDivisions is the final number divisions
    // it changes the splitseg the inner points after the splitting
    void  splitSegment(const Point3f& A,const Point3f& B, vector<Point3f>& splitseg ,int numberDivisions);

    template<class ScalarType>
    void  splitSegment(const Point3<ScalarType>& A,const Point3<ScalarType>& B, vector<Point3<ScalarType>>& splitseg ,int numberDivisions){
        vector<Point3<ScalarType>> result;
        splitseg.clear();
        if(numberDivisions<2)
            return;
        Point3<ScalarType> direction=B-A;
        direction.normalized();
        ScalarType len=vcg::Distance(A,B);
        ScalarType step=len/numberDivisions;
        //result.push_back(A);
        for(int i=1;i<numberDivisions;i++){
            result.push_back(A+direction*(i*step));
        }
        //result.push_back(B);
        splitseg=result;
    }

    //split some pos edges and recompute the feature indexes of a given polygon
    // For each i in poslist this method split the segment of indexes [i,i+1] in polygon
    void splitPosList_onPolygon(vector<Point3f> & polygon,vector<size_t> & poslist,vector<int> &timestosplit, set<int> featureOnPolygon,vector<Point3f> &result);

    // This converts a polygon to a regular a polygon
    template<class ScalarType>
    void polygon_to_regular_polygon(vector<vector<Point3<ScalarType>>> &ipoly,vector<vector<Point3<ScalarType>>> &opoly){
        opoly.clear();
        size_t n=ipoly.size();
        for(size_t i=0;i<n;i++){
            vector<Point3<ScalarType>> nside;
            nside.push_back(ipoly[i][0]);
            vector<Point3<ScalarType>> temp;
            splitSegment<ScalarType>(ipoly[i][0],ipoly[i].back(),temp,ipoly[i].size()-1);
            nside.insert(nside.end(),temp.begin(),temp.end());
            nside.push_back(ipoly[i].back());
            opoly.push_back(nside);
        }
    }
    //Save a polygon represented as a list of points for each side, i.e, P=[ [points side 0],[points side 1],...,[points side n]]
    template <class S>
    void polygon_to_txt(vector<vector<Point3<S>>> & poly,string filename){
        if(filename.empty())
            filename+="Exportedpoly.txt";

        std::ofstream fileCurve(filename);
        if (fileCurve.is_open()){
           for(int i=0;i<poly.size();i++){
               for(int j=0;j<poly[i].size();j++){
                   fileCurve<<i<<" "<<poly[i][j].X()<<" "<<poly[i][j].Y()<<" "<<poly[i][j].Z()<< '\n';
               }
           }
        }
        fileCurve.close();
    }
    //Read a polygon represented as a list of points for each side, i.e, P=[ [points side 0],[points side 1],...,[points side n]]
    template <class S>
    void readtxt_polygon(string filename,std::vector<std::vector<vcg::Point3<S>>> &poly){
        if(filename.empty())
           return;
        std::ifstream in(filename, ios::in);
        if (!in){
            cerr << "Cannot open " << filename << endl;
            return;
        }
        string line;
        std::vector<vcg::Point3<S>> tempside;
        std::vector<int> indexes;
        while (std::getline(in, line))
        {
            istringstream s(line);
            double vx,vy, vz;
            int i;
            s >> i;
            s >> vx;
            s >> vy;
            s >> vz;
            indexes.push_back(i);
            tempside.push_back(Point3<S>(vx,vy,vz));
        }
        poly.clear();
        for(size_t i=0;i<=indexes.back();i++){
            std::vector<vcg::Point3<S>> temp;
            for (size_t j=0;j<tempside.size();j++) {
                if(indexes[j]==i){
                    temp.push_back(tempside[j]);
                }
            }
            poly.push_back(temp);
        }

    }
    template <class MeshType>
    void collapseEdges(MeshType &m, vector<face::Pos<typename MeshType::FaceType>> & edgestocollapse){
        if(edgestocollapse.empty())
            return;
        typedef typename MeshType::CoordType CoordType;
        CoordType np;
        if(edgestocollapse.size()==1){
            np=(edgestocollapse[0].VFlip()->P()+edgestocollapse[0].V()->P())/2;
        }
        else{
            np=(edgestocollapse[0].VFlip()->P()+edgestocollapse.back().V()->P())/2;
        }
        for(size_t i=0;i<edgestocollapse.size();i++){
            vcg::tri::BasicVertexPair<CMesh::VertexType> vertexpair(edgestocollapse[i].VFlip(),edgestocollapse[i].V());
            bool isFeasible=vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>>::LinkConditions(vertexpair);
            if(isFeasible){
                vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>> eco;
                eco.Do(m,vertexpair,np);
            }
        }
    }
    template <class MeshType>
    void collapseEdgesToOne(MeshType &m, vector<face::Pos<typename MeshType::FaceType>> & edgestocollapse){
        if(edgestocollapse.size()<=1)
            return;
        typedef typename MeshType::CoordType CoordType;
        CoordType np=edgestocollapse.back().V()->P();

        for(size_t i=1;i<edgestocollapse.size();i++){
            vcg::tri::BasicVertexPair<CMesh::VertexType> vertexpair(edgestocollapse[i].VFlip(),edgestocollapse[i].V());
            bool isFeasible=vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>>::LinkConditions(vertexpair);
            if(isFeasible){
                vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>> eco;
                eco.Do(m,vertexpair,np);
            }
        }
    }


    template <class MeshType>
    void remeshingCurveSegmentByCollapsing(MeshType &m, vector<face::Pos<typename MeshType::FaceType>> & curve, typename MeshType::ScalarType threshold){
        if(curve.empty() || curve.size()==1) // if curve is empty or has only one segment we can't collasing anything
            return;

        typedef typename MeshType::CoordType CoordType;
        typedef typename MeshType::ScalarType ScalarType;
        typedef typename MeshType::FaceType FaceType;
        typedef typename face::Pos<FaceType> PosType;

        ScalarType length=0;
        size_t ndiv=1;
        //computing length
        curve[0].VFlip()->ClearV();
        for(size_t i=0;i<curve.size();i++){
           curve[i].V()->ClearV();
           length+=vcg::Distance(curve[i].VFlip()->P(),curve[i].V()->P());
        }

        //computing number of divisions
        if(length>1.5*threshold){
            ndiv=int(length/threshold)+1;
        }
        // collapsing
        if(curve.size()>1){
            if(ndiv==1){
                collapseEdgesToOne(m,curve);
            }
            else if(ndiv>1){
                ScalarType step=length/ScalarType(ndiv);
                curve[0].VFlip()->SetV();
                curve.back().V()->SetV();
                ScalarType temp_length=0;
                size_t i=0;

                // marking new extreme vertices
                while(i<curve.size()){
                   PosType temp_pos=curve[i];
                   ScalarType distance=vcg::Distance(temp_pos.VFlip()->P(),temp_pos.V()->P());
                   temp_length+=distance;
                   if(temp_length>step){
                       if(temp_length-step>step-temp_length+distance){
                           temp_pos.VFlip()->SetV();
                           //temp_pos.VFlip()->C()=Color4b::Red;
                       }
                       else{
                           temp_pos.V()->SetV();
                           //temp_pos.V()->C()=Color4b::Red;
                       }
                       temp_length=0;
                   }
                   i++;
                }

                // computing blocks to be collapsed
                vector<vector<PosType>> tocollapse;
                vector<PosType> temp_block;
                i=0;
                while(i<curve.size()){
                    temp_block.push_back(curve[i]);
                    if(curve[i].V()->IsV()){
                        tocollapse.push_back(temp_block);
                        temp_block.clear();
                    }
                    i++;
                }

                // deleting last block if this is too small

                //computing virtual length of the last block
                ScalarType last_length=vcg::Distance(tocollapse.back()[0].VFlip()->P(),tocollapse.back().back().V()->P());
                if(last_length<step/3){
                    vector<PosType> last=tocollapse.back();
                    tocollapse.pop_back();
                    vector<PosType> before_last=tocollapse.back();
                    tocollapse.pop_back();
                    before_last.insert(before_last.end(),last.begin(),last.end());
                    tocollapse.push_back(before_last);
                }
                // collapsing
                for(size_t j=0;j<tocollapse.size();j++){                    
                    collapseEdgesToOne(m,tocollapse[j]);
                }

            }            
        }
    }
    template <class FaceType>
    typename FaceType::ScalarType getAverageEdgeLengthFace(FaceType &f){
        assert(f.VN()>1);
        typename FaceType::ScalarType result;
        for(size_t i=0;i<f.VN();i++){
            result+=vcg::Distance(f.V(i)->P(),f.V(f.Next(i))->P());
        }
        result/=f.VN();
        return result;
    }

    // works for 2D meshes with Z coordinate equal to 0
    template <class TriMeshType,class PolyMeshType>
    void getBarycentrycCoordinatesQ2T(TriMeshType &tmesh,PolyMeshType &pmesh,
                                      vector<Point4<typename TriMeshType::ScalarType>> &coordinates){
        typedef typename TriMeshType::ScalarType ScalarType;
        typedef typename TriMeshType::CoordType CoordType;
        typedef typename TriMeshType::FaceIterator FaceIterator;
        coordinates.clear();
        size_t n=pmesh.vert.size();
        coordinates.resize(n);
        Point3<ScalarType> L(0,0,0);
        for(size_t i=0;i<n;i++){
            coordinates[i]=Point4<ScalarType>(0,0,0,0);
            for(FaceIterator fi=tmesh.face.begin();fi!=tmesh.face.end();fi++){
                CoordType a3=fi->V(0)->P();
                CoordType b3=fi->V(1)->P();
                CoordType c3=fi->V(2)->P();
                Point2<ScalarType> p(pmesh.vert[i].P().X(),pmesh.vert[i].P().Y());
                Point2<ScalarType> a2(a3.X(),a3.Y());
                Point2<ScalarType> b2(b3.X(),b3.Y());
                Point2<ScalarType> c2(c3.X(),c3.Y());
                Triangle2<ScalarType> triangle(a2,b2,c2);
                if(vcg::IsInsideTrianglePoint(triangle,p)){
                      vcg::InterpolationParameters2(triangle.P(0),triangle.P(1),triangle.P(2),p,L);
                      Point4<ScalarType> result(tri::Index(tmesh,&*fi),L.X(),L.Y(),L.Z());
                      coordinates[i]=result;
                      break;
                }
            }
        }
    }
    template <class MeshType>
    vcg::face::Pos<typename MeshType::FaceType> getPosFromIndex(MeshType& mesh,int id1,int id2){
        UpdateTopology<MeshType>::FaceFace(mesh);
        vcg::face::Pos<typename MeshType::FaceType> result;
        result.SetNull();
        for(typename MeshType::FaceIterator fi=mesh.face.begin();fi!=mesh.face.end();fi++){
            if(!fi->IsD()){
                for(int k=0;k<fi->VN();k++){
                    if(vcg::tri::Index(mesh,fi->V(k))==id1 ){
                        result.Set(&*fi,k,fi->V(k));
                        break;
                    }
                }
                if(!result.IsNull())
                    break;
            }
        }
        while(vcg::tri::Index(mesh,result.VFlip())!=id2){
            result.NextE();
        }
        return result;
    }
    // Store many quad meshes into a sole quad mesh. Note that we dont worry with coincident vertices.
    template<typename PolyMeshType>
    double minimum_length_quadmesh(PolyMeshType &m){
        double len=m.bbox.Diag();
        for(typename PolyMeshType::FaceIterator fi=m.face.begin();fi!=m.face.end();fi++){
            if(!fi->IsD()){
                for(int k=0;k<fi->VN();k++){
                    double dist=vcg::Distance(fi->V(k)->P(),fi->V((k+1)%fi->VN())->P());
                    if(dist<len){
                        len=dist;
                    }
                }
                double distd1=vcg::Distance(fi->V(0)->P(),fi->V(2)->P());
                double distd2=vcg::Distance(fi->V(1)->P(),fi->V(3)->P());
                len=std::min(len,distd1);
                len=std::min(len,distd2);
            }
        }
        return len;
    }
    template<typename MeshType>
    void mergeQuadMeshes_LibiglFormat(vector<pair<Eigen::MatrixXd,Eigen::MatrixXi>> &meshes,MeshType &result){
        result.Clear();
        int totalvertices=0;
        int totalfaces=0;
        int n=meshes.size();
        for(int i=0;i<n;i++){
             totalvertices+=meshes[i].first.rows();
             totalfaces+=meshes[i].second.rows();
        }
        Eigen::MatrixXd rV(totalvertices,3);
        Eigen::MatrixXi rF(totalfaces,4);
        int vaccumulated=0;
        int faccumulated=0;
        for(int i=0;i<n;i++){
            int vrows=meshes[i].first.rows();
            int frows=meshes[i].second.rows();
            Eigen::MatrixXi aux=meshes[i].second.array()+vaccumulated;
            rV.block(vaccumulated,0,vrows,3)=meshes[i].first;
            rF.block(faccumulated,0,frows,4)=aux;
            vaccumulated+=vrows;
            faccumulated+=frows;
        }
        Allocator<MeshType>::AddVertices(result,rV.rows());
        auto pfi=Allocator<MeshType>::AddFaces(result,rF.rows());

        for(size_t i=0;i<rV.rows();++i)
        {
            typename MeshType::CoordType vv(rV(i,0),rV(i,1),rV(i,2));
            result.vert[i].P() = vv;
            //cout<<rV.row(i)<<endl;
        }

        for(size_t i=0;i<rF.rows();++i)
        {
            assert( rF(i,0)>=0 );
            assert( rF(i,1)>=0 );
            assert( rF(i,2)>=0 );
            assert( rF(i,3)>=0 );
            assert( rF(i,0)<result.vn );
            assert( rF(i,1)<result.vn );
            assert( rF(i,2)<result.vn );
            assert( rF(i,3)<result.vn );
            result.face[i].Alloc(4);
            result.face[i].V(0) = &result.vert[rF(i,0)];
            result.face[i].V(1) = &result.vert[rF(i,1)];
            result.face[i].V(2) = &result.vert[rF(i,2)];
            result.face[i].V(3) = &result.vert[rF(i,3)];
        }
        tri::UpdateBounding<MeshType>::Box(result);
    }

    // conversion functions from and to libigl format
    template< typename MeshType>
    void convertToLIBIGL(const MeshType& m, Eigen::MatrixXd& V,  Eigen::MatrixXi& F){
        if(m.vert.size()==V.rows() && m.face.size()==F.rows()){
            int i=0;
            for(typename MeshType::ConstVertexIterator vi=m.vert.begin();vi!=m.vert.end();vi++,i++)
                V.row(i)=Eigen::Vector3d(vi->cP().X(),vi->cP().Y(),vi->cP().Z());
            i=0;
            int numberCols=F.cols();
            for(typename  MeshType::ConstFaceIterator fi=m.face.begin();fi!=m.face.end();fi++,i++){
                assert(fi->VN()==numberCols);
                for(int k=0;k<numberCols;k++)
                    F(i,k)=tri::Index(m,fi->cV(k));
                //F.row(i)=Eigen::Vector3i(tri::Index(m,fi->cV(0)),tri::Index(m,fi->cV(1)),tri::Index(m,fi->cV(2)));
            }
        }
        else{
            cout<<"Impossible convert to libIGL format"<<endl;
        }
    }

    template< typename MeshType>
    void libiglToVCG(MeshType& m, const Eigen::MatrixXd& V,const Eigen::MatrixXi& F){
        typedef typename MeshType::CoordType CoordType;
        typedef typename MeshType::ScalarType ScalarType;
        vector<CoordType> coordinates;
        vector<Point3i> faceindexes;
        for(int i=0;i<V.rows();i++)
            coordinates.push_back(CoordType(ScalarType(V(i,0)),ScalarType(V(i,1)),ScalarType(V(i,2))));
        for(int i=0;i<F.rows();i++){
            faceindexes.push_back(Point3i(int(F(i,0)),int(F(i,1)),int(F(i,2))));
        }
        cout<<"vector coord"<<coordinates.size()<<endl;
        cout<<"vector indexes"<<faceindexes.size()<<endl;
        tri::BuildMeshFromCoordVectorIndexVector<MeshType,CoordType,Point3i>(m,coordinates,faceindexes);
        vcg::tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFace(m);
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(m);
        UpdateTopology<MeshType>::FaceFace(m);
    }
    template< typename PolyMeshType>
    void Quad_libiglToVCG(PolyMeshType& m, const Eigen::MatrixXd& V,const Eigen::MatrixXi& F){
        typedef typename PolyMeshType::CoordType CoordType;
        typedef typename PolyMeshType::ScalarType ScalarType;
        m.Clear();
        Allocator<PolyMeshType>::AddVertices(m,V.rows());
        auto pfi=Allocator<PolyMeshType>::AddFaces(m,F.rows());

        for(size_t i=0;i<V.rows();++i)
        {
            CoordType vv(V(i,0),V(i,1),V(i,2));
            m.vert[i].P() = vv;
        }

        for(size_t i=0;i<F.rows();++i)
        {
            assert( F(i,0)>=0 );
            assert( F(i,1)>=0 );
            assert( F(i,2)>=0 );
            assert( F(i,3)>=0 );
            assert( F(i,0)<m.vn );
            assert( F(i,1)<m.vn );
            assert( F(i,2)<m.vn );
            assert( F(i,3)<m.vn );
            m.face[i].Alloc(4);
            m.face[i].V(0) = &m.vert[F(i,0)];
            m.face[i].V(1) = &m.vert[F(i,1)];
            m.face[i].V(2) = &m.vert[F(i,2)];
            m.face[i].V(3) = &m.vert[F(i,3)];
        }
        vcg::tri::UpdateBounding<PolyMeshType>::Box(m);
        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(m);
        vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(m);
    }
    // this returns the value of cos(a) where a is the angle between n0 and n1. (scalar prod is cos(a))
    template <typename ScalarType>
    inline ScalarType fastAngle(Point3<ScalarType> n0, Point3<ScalarType> n1)
    {
        return math::Clamp(n0*n1,(ScalarType)-1.0,(ScalarType)1.0);
    }
    template <typename ScalarType>
    inline ScalarType cosinAngleFromVector(Point3<ScalarType> n0, Point3<ScalarType> n1)
    {
        return math::Clamp(n0*n1/(n0.Norm()*n1.Norm()),(ScalarType)-1.0,(ScalarType)1.0);
    }
    // this returns the value of cos(a) where a is the angle of the vertex Ind
    template <typename ScalarType>
    inline ScalarType cosinAngleFromPoints(Point3<ScalarType> p0, Point3<ScalarType> p1,Point3<ScalarType> p2, int Ind)
    {
        Point3<ScalarType> s1,s2;
        switch(Ind){
            case 0: {
                s1=p1-p0;
                s2=p2-p0;
                break;
            }
            case 1: {
                s1=p0-p1;
                s2=p2-p1;
                break;
            }
            case 2: {
                s1=p1-p2;
                s2=p0-p2;
                break;
            }
            default:break;
        }
        return math::Clamp(s1*s2/(s1.Norm()*s2.Norm()),(ScalarType)-1.0,(ScalarType)1.0);
    }
    template< typename MeshType>
    std::vector<face::Pos<typename MeshType::FaceType>> getBorderCurveFromOnePos(face::Pos<typename MeshType::FaceType>& startPos){
        face::Pos<typename MeshType::FaceType> curPos=startPos;
        std::vector<face::Pos<typename MeshType::FaceType>> result;
        result.push_back(curPos);
        do
        {
          curPos.NextB();
          result.push_back(curPos);
        }
        while(curPos!=startPos);
        result.pop_back();
        return result;
    }
    // Only for quad meshes
    template <typename MeshType>
    inline void walkOnFlow(face::Pos<typename MeshType::FaceType>& pos, std::function<bool (face::Pos<typename MeshType::FaceType> &)> action)
    {
        face::Pos<typename MeshType::FaceType> pi=pos;
        int count=0;
        do{
            bool result=action(pi);
            //cout<<"action done"<<endl;
            if(!result) break;
            pi.FlipE();
            if(pi.IsBorder())
                break;
            pi.FlipF();
            pi.FlipE();
            pi.FlipV();
            //if(pi.V()->Q()!=2 || pi.V()->Q()!=4)
              //  break;
            vector<face::Pos<typename MeshType::FaceType>> posVec;
            face::VFOrderedStarFF<PMesh::FaceType>(pi,posVec);
            if(posVec.size()!=2 && posVec.size()!=4)
               break;
            count++;
        }while(pi!=pos && count<1000);
        cout<<"number of iterations walkOnFlow"<<count<<endl;
    }
    // Only for triangle meshes
    template <typename MeshType>
    inline void forEachFacePos(MeshType &m, std::function<void (face::Pos<typename MeshType::FaceType> &)> action)
    {
        for(auto fi=m.face.begin();fi!=m.face.end();++fi)
            if(!(*fi).IsD())
            {
                for(int i=0;i<3;++i)
                {
                    face::Pos<typename MeshType::FaceType> pi(&*fi,i);
                    action(pi);
                }
            }
    }
    // Only for triangle meshes
    template <typename MeshType>
    void ImproveQuality(MeshType &m,typename MeshType::ScalarType angle,bool selectedOnly)
    {
       typedef typename MeshType::ScalarType ScalarType;
       typedef typename face::Pos<typename MeshType::FaceType> PosType;
       typedef typename MeshType ::VertexPointer VertexPointer;
       ScalarType creaseAngleCosThr=cos(math::ToRad(angle));;
       tri::UpdateTopology<MeshType>::FaceFace(m); //collapser does not update FF
       forEachFacePos<MeshType>(m, [&](PosType &p){
         ScalarType angle=fastAngle<ScalarType>(NormalizedTriangleNormal(*(p.F())), NormalizedTriangleNormal(*(p.FFlip())));
         bool testCrease= (angle <= creaseAngleCosThr && angle >= -creaseAngleCosThr);
         PosType p1=p;
         PosType p2=p;
         p1.FlipE();
         p2.FlipF();
         p2.FlipE();
         /// Connecting opposing vertices and check if this second
         ///  is shorter
         /// than the current edge
         /*  oldD1 ___________ newD1
         *       |\          |
         *       |  \        |
         *       |    \      |
         *       |  f  z\    |
         *       |        \  |
         *       |__________\|
         * newD0               oldD0
         */
         VertexPointer OldDiag0 = p.V();
         VertexPointer OldDiag1 = p.VFlip();
         VertexPointer NewDiag0 = p1.VFlip();
         VertexPointer NewDiag1 = p2.VFlip();
         float currentDiag=vcg::Distance(OldDiag0->P(),OldDiag1->P());
         float newDiag=vcg::Distance(NewDiag0->P(),NewDiag1->P());
         bool newIsShorter=newDiag<currentDiag;
         bool ispossible=(NewDiag1 != NewDiag0) && (NewDiag1 != OldDiag0) && (NewDiag1 != OldDiag1);
         bool testQuality=ispossible && newIsShorter;
         if(p.FFlip() > p.F())
           if(((!selectedOnly) || (p.F()->IsS() && p.FFlip()->IsS())) &&
              !testCrease && !p.IsBorder() && testQuality &&
              face::CheckFlipEdgeNormal(*p.F(), p.E(), math::ToRad(10.f)) &&
              face::CheckFlipEdge(*p.F(), p.E()) )
           {
             face::FlipEdge(*p.F(), p.E());
           }
       });
    }

    //Project Point3f to a Quad mesh
    //Get coodinates
    template <typename MeshType>
    Point2f projectOnMesh(MeshType& m, vcg::GridStaticPtr<typename MeshType::FaceType, typename MeshType::ScalarType>& gm, const Point3f& p,
                          std::unordered_map<typename MeshType::VertexPointer,Point2f>& mapUV){
        //tri::RequireVFAdjacency(m);
        //UpdateTopology<MeshType>::VertexFace(m);
        typedef typename MeshType::CoordType CoordType;
        typedef typename MeshType::FaceType FaceType;
        typedef typename MeshType::VertexType VertexType;
        typedef typename MeshType::ScalarType ScalarType;
        typedef typename MeshType::VertexPointer VertexPointer;
        typedef typename vcg::GridStaticPtr<FaceType,ScalarType> SpatialType;
        ScalarType dist;
        CoordType Norm, ip, nearest;
        FaceType* fac = vcg::tri::GetClosestFaceBase<MeshType,SpatialType>(m,gm, p, m.bbox.Diag(),  dist, nearest, Norm, ip );
        for(int i=0;i<fac->VN();i++)
            if(vcg::SquaredDistance(p,fac->V(i)->P())<0.0000001)
                return mapUV[fac->V(i)];

        /*assert( v != NULL );
        vcg::face::VFIterator<FaceType> vfi(v);
        FaceType *fac=NULL;
        if(vcg::SquaredDistance(p,v->P())<0.0000001)
            return mapUV[v];

        for(;!vfi.End();++vfi){
           FaceType* f = vfi.F();
           CoordType normalBase=f->N();
           bool isthis=true;
           for (size_t i=0;i<(size_t)f->VN();i++){
              CoordType p0=f->P0(i);
              CoordType p1=f->P1(i);
              CoordType  normalTriangle=vcg::Normal<CoordType>(p0,p1,p);
              ScalarType dihedralangle=vcg::AngleN(normalBase,normalTriangle);
              dihedralangle=dihedralangle*180/3.1416;
              isthis=isthis && (dihedralangle<=90);
           }
           if(isthis){
               fac=f;
               break;
           }
        }*/
        assert( fac != NULL );
        // Finding best triangle in the found Face
        std::vector<PMesh::CoordType> pointsToTest;
        for(int k= 0; k <fac->VN(); ++k) {
            PMesh::VertexType * v = fac->V(k);
            pointsToTest.push_back(v->P());
        }
        std::vector<int> facestoTest;
        TessellatePlanarPolygon3(pointsToTest,facestoTest);

        VertexPointer U0=fac->V(facestoTest[0]);
        VertexPointer U1=fac->V(facestoTest[1]);
        VertexPointer U2=fac->V(facestoTest[2]);
        VertexPointer V0=fac->V(facestoTest[1]);
        VertexPointer V1=fac->V(facestoTest[2]);
        VertexPointer V2=fac->V(facestoTest[3]);

        vcg::Ray3<ScalarType> r;
        r.Set( p, fac->N());
        r.Normalize();

        if((p-U0->P()).dot(fac->N())>0 || (p-U1->P()).dot(fac->N())>0 || (p-U2->P()).dot(fac->N())>0){
            r.Flip();
        }
        else if((p-V0->P()).dot(fac->N())>0 || (p-V1->P()).dot(fac->N())>0 || (p-V2->P()).dot(fac->N())>0){
            r.Flip();
        }

        ScalarType t,alfa,beta,gamma;
        int U0U1U2=false;
        if(vcg::IntersectionRayTriangle(r,U0->P(),U1->P(),U2->P(),t,alfa,beta))
            U0U1U2=true;
        else if(vcg::IntersectionRayTriangle(r,V0->P(),V1->P(),V2->P(),t,alfa,beta))
            U0U1U2=false;

        Point2f p1,p2,p3,result;
        if(U0U1U2){
            p1=mapUV.at(U0);
            p2=mapUV.at(U1);
            p3=mapUV.at(U2);
        }
        else{
            p1=mapUV.at(V0);
            p2=mapUV.at(V1);
            p3=mapUV.at(V2);
        }
        gamma=1-alfa-beta;
        result=p1.Scale(gamma,gamma)+p2.Scale(alfa,alfa)+p3.Scale(beta,beta);
        return result;
    }
    //Find the closest point in a edgemesh to a target point
    template <typename MeshType>
    int closestInMesh(const MeshType& m, const Point3f& p){
        float dist=std::numeric_limits<float>::infinity();
        int min=-1;
        for(typename MeshType::ConstVertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi){
            if(!vi->IsD()){
                float d=(vi->cP()-p).SquaredNorm();
                if(d<dist){
                    dist=d;
                    min=tri::Index(m,*vi);
                }
            }
        }
        return min;
    }
    template <typename MeshType>
    float ComputeFaceEdgeLengthAverage(MeshType & m, bool selected=false)
     {
       float sum=0;
       int n=0;
       vector<typename MeshType::CoordType> points;
       for(typename MeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
         if(!(*fi).IsD())
           if(!selected || fi->IsS())
           {
               points.clear();
               for(int i=0;i<fi->VN();++i)
                   points.push_back(fi->V(i)->P());
               for(int j=0;j<points.size();j++){
                  sum+=float(vcg::Distance(points[j],points[(j+1)%points.size()]));
                  n+=fi->VN();
               }
           }
       if(n!=0) sum=sum/n;
       else{
           cout<<"non selected faces"<<endl;
           assert(0);
       }
       return sum;
     }
    template <typename MeshType>
    int closestInMesh(const vector<face::Pos<typename MeshType::FaceType>>& m, const Point3f& p){
        float dist=std::numeric_limits<float>::infinity();
        int min=-1;
        int i=0;
        for(typename vector<face::Pos<typename MeshType::FaceType>>::const_iterator vi=m.begin();vi!=m.end();++vi,i++){
            if(!vi->V()->IsD()){
                float d=(vi->V()->cP()-p).SquaredNorm();
                if(d<dist){
                    dist=d;
                    min=i;
                }
            }
        }
        return min;
    }
    template <typename MeshType>
    int closestIn(const vector<typename MeshType::VertexPointer>& curve, const Point3f& p){
        float dist=std::numeric_limits<float>::infinity();
        int min=-1;
        int i=0;
        for(typename vector<typename MeshType::VertexPointer>::const_iterator vi=curve.begin();vi!=curve.end();++vi,i++){
            if(!(*vi)->IsD()){
                float d=((*vi)->cP()-p).SquaredNorm();
                if(d<dist){
                    dist=d;
                    min=i;
                }
            }
        }
        return min;
    }
    template <typename MeshType>
    bool isFullSelected(const MeshType& m){
        for(typename MeshType::ConstFaceIterator fi=m.face.begin();fi!=m.face.end();++fi){
             if(!fi->IsD() && !fi->IsS()){
                 return false;
             }
        }
        return true;
    }
    //Need a selected region where test the algorithm
    template <typename MeshType>
    bool isInsideMesh(const MeshType& m, const Point3f& p,bool selectedFaces=true){
        Point3f normalPointMesh=m.vert[0].N();
        Point3f PointMesh=m.vert[0].P();
        Plane3f facemesh;
        // Find the nearest face to p in selected faces of the mesh
        float dist=std::numeric_limits<float>::infinity();
        float distplano=std::numeric_limits<float>::infinity();
        for(typename MeshType::ConstFaceIterator fi=m.face.begin();fi!=m.face.end();++fi){
            if(fi->IsS() || !selectedFaces){
                vcg::Plane3f pface;
                pface.Init(fi->cV(0)->P(),fi->cN());
                float dplano=vcg::SignedDistancePlanePoint(pface,p);
                if(dplano*dplano<distplano){
                     facemesh=pface;
                     distplano=dplano*dplano;
                }
                for(int k=0;k<fi->VN();k++){
                    float d=(fi->cV(k)->P()-p).SquaredNorm();
                    if(d<dist){
                        //normalPointMesh=fi->cN();
                        normalPointMesh=fi->cV(k)->cN();
                        PointMesh=fi->cV(k)->cP();
                        dist=d;
                    }
                }
            }
        }
        Point3f vector=p-PointMesh;
        if(vcg::Angle(normalPointMesh,facemesh.Direction())>1.57)
            cout<<"contrarios"<<endl;
        //float distance=vcg::SignedDistancePlanePoint(facemesh,p);
        if(vector.dot(facemesh.Direction())<=0) return true;
        else return false;
        //if(distance<0) return true;
    }
    template <typename PolyMeshType>
    void triangulatePolyMesh(PolyMeshType & pm,bool selectedFaces=false){

        vector<typename PolyMeshType::CoordType> points;
        typedef typename PolyMeshType::FaceIterator PolyFaceIterator;
        typedef typename PolyMeshType::FacePointer PolyFacePointer;

        vector<vector<size_t>> indices;
        //Compute new triangles to add
        for(PolyFaceIterator fi = pm.face.begin(); fi != pm.face.end(); ++fi)
        {
          if(!((*fi).IsD()) && (!selectedFaces || (*fi).IsS())){
            points.clear();
            for(int i  = 0; i < (*fi).VN(); ++i) {
              typename	PolyMeshType::VertexType * v = (*fi).V(i);
              points.push_back(v->P());
            }
            std::vector<int> faces;
            TessellatePlanarPolygon3(points,faces);
            for(size_t i = 0; i<faces.size();i+=3){
              vector<size_t> tfaces;
              tfaces.push_back(tri::Index(pm,(*fi).V( faces[i+0] )));
              tfaces.push_back(tri::Index(pm,(*fi).V( faces[i+1] )));
              tfaces.push_back(tri::Index(pm,(*fi).V( faces[i+2] )));
              indices.push_back(tfaces);
            }
          }
        }
        //delete selected polygon faces
        for(PolyFaceIterator fi = pm.face.begin(); fi != pm.face.end(); ++fi)
            if(!((*fi).IsD()) && (!selectedFaces || (*fi).IsS()))
                tri::Allocator<PMesh>::DeleteFace(pm,*fi);
        //Add new triangle faces
        for(int k=0;k<indices.size();k++){
           PolyFaceIterator f= tri::Allocator<PolyMeshType>::AddFace(pm,indices[k][0],indices[k][1],indices[k][2]);
           f->SetS();
        }
        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(pm);
        vcg::tri::UpdateNormal<PolyMeshType>::PerFaceNormalized(pm);
   }
    //extract triangle mesh from a polymesh. Only selected faces are taking into account if the flag is true.
    // If selectedFaces are true and deleteFaces too so the selected faces from polymesh are deleted after the process
    template <typename TriMeshType,typename PolyMeshType>
    void extractTriangleFromSelectedPolyMesh(TriMeshType & tm,PolyMeshType & pm,bool selectedFaces=false,bool deleteFaces=false){
       if(!selectedFaces) tri::PolygonSupport<TriMeshType,PolyMeshType>::ImportFromPolyMesh(tm,pm);
       else{
           tm.Clear();
           vector<typename PolyMeshType::CoordType> points;
           typedef typename PolyMeshType::FaceIterator PolyFaceIterator;
           typedef typename TriMeshType::FaceIterator TriFaceIterator;
           typedef typename PolyMeshType::FacePointer PolyFacePointer;
           typedef typename TriMeshType::VertexIterator  TriVertexIterator;
           typedef typename TriMeshType::VertexPointer  TriVertexPointer;
           typedef typename TriMeshType::CoordType  TriVertexCoordType;

           vector<vector<TriVertexCoordType>> indices;
           vector<TriVertexCoordType> vertToAdd;
           //Compute new triangles to add
           for(PolyFaceIterator fi = pm.face.begin(); fi != pm.face.end(); ++fi)
           {
             if(!((*fi).IsD()) && (!selectedFaces || (*fi).IsS())){
               points.clear();
               for(int i  = 0; i < (*fi).VN(); ++i) {
                 typename	PolyMeshType::VertexType * v = (*fi).V(i);
                 points.push_back(v->P());
               }
               std::vector<int> faces;
               TessellatePlanarPolygon3(points,faces);
               for(size_t i = 0; i<faces.size();i+=3){
                 vector<TriVertexCoordType> tfaces;
                 tfaces.push_back((*fi).V( faces[i+0])->P());
                 tfaces.push_back((*fi).V( faces[i+1])->P());
                 tfaces.push_back((*fi).V( faces[i+2])->P());
                 indices.push_back(tfaces);

               }
             }
           }
           if(deleteFaces){
               //delete selected polygon faces
               for(PolyFaceIterator fi = pm.face.begin(); fi != pm.face.end(); ++fi)
                   if(!((*fi).IsD()) && (!selectedFaces || (*fi).IsS()))
                       tri::Allocator<PMesh>::DeleteFace(pm,*fi);
           }

           //Add new triangle faces
           for(int k=0;k<indices.size();k++){
              TriFaceIterator f= tri::Allocator<TriMeshType>::AddFace(tm,indices[k][0],indices[k][1],indices[k][2]);
              f->SetS();
           }

           //cout<<"vertices "<<tm.VN()<<endl;
           //cout<<"faces "<<tm.FN()<<endl;
           tri::Clean<TriMeshType>::RemoveDuplicateVertex(tm);
           tri::Allocator<TriMeshType>::CompactEveryVector(tm);
           vcg::tri::UpdateTopology<TriMeshType>::FaceFace(tm);
           vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalized(tm);
           vcg::tri::UpdateNormal<TriMeshType>::PerFaceNormalized(tm);
           //cout<<"vertices "<<tm.VN()<<endl;
           //vcg::tri::UpdateNormal<PolyMeshType>::PerFaceNormalized(pm);
       }

   }
    /*template <typename MeshType>
    bool loopIsIn(MeshType & m,const vector<face::Pos<typename MeshType::FaceType>>& loop,bool selectedFaces=true){
        int score=0;
        for(int i=0;i<loop.size();i++){
            if(utility::isInsideMesh<MeshType>(m,loop[i].V()->P(),selectedFaces))
                score++;
        }
        cout<<"score: "<<score*100/loop.size()<<endl;
        if(score>0.7*loop.size())
            return true;
        else
            return false;
    }*/
    template <class MeshType>
    bool loopIsIn(MeshType & m,GridStaticPtr<typename MeshType::FaceType, typename MeshType::ScalarType> & gm,const vector<face::Pos<typename MeshType::FaceType>>& loop){
        int score=0;
        for(int i=0;i<loop.size();i++){
            //if(utility::isInsideMesh<MeshType>(m,loop[i].V()->P(),selectedFaces))
            if(tri::Inside<GridStaticPtr<typename MeshType::FaceType, typename MeshType::ScalarType>,MeshType>::Is_Inside(m,gm,loop[i].V()->P()))
                score++;
        }
        //cout<<"score: "<<score*100/loop.size()<<endl;
        if(score>0.6*loop.size())
            return true;
        else
            return false;
    }
    template <class MeshType,class OMeshType>
    bool loopIsInM(OMeshType & m,GridStaticPtr<typename OMeshType::FaceType, typename OMeshType::ScalarType> & gm,const vector<face::Pos<typename MeshType::FaceType>>& loop){
        int score=0;
        for(int i=0;i<loop.size();i++){
            //if(utility::isInsideMesh<MeshType>(m,loop[i].V()->P(),selectedFaces))
            if(tri::Inside<GridStaticPtr<typename OMeshType::FaceType, typename OMeshType::ScalarType>,OMeshType>::Is_Inside(m,gm,loop[i].V()->P())){
                score++;
                //loop[i].V()->SetS();
            }
        }
        //cout<<"score: "<<score*100/loop.size()<<endl;
        if(score>0.6*loop.size())
            return true;
        else
            return false;
    }
    // need that tha faces to subdivide doing  n-ring around some vertices
    template<class MeshType>
    void RegularSubdivideFaces(MeshType & m,vector<typename MeshType::FacePointer>& facetoSubdivide ){
        typedef typename MeshType::FaceIterator FaceIterator;
        typedef typename MeshType::VertexIterator VertexIterator;
        typedef typename MeshType::VertexPointer VertexPointer;
        typedef typename MeshType::FacePointer FacePointer;
        if(facetoSubdivide.size()!=0){
        FaceIterator lastf = tri::Allocator<MeshType>::AddFaces(m,facetoSubdivide.size()*3);
        VertexIterator lastv = tri::Allocator<MeshType>::AddVertices(m,facetoSubdivide.size()*3);
        /*
         *               v0
         *              /  \
         *            /  f0  \
         *          /          \
         *        mp01----------mp02
         *       /  \    f3    /   \
         *     / f1   \      /  f2   \
         *   /          \  /           \
         *v1 ---------- mp12------------v2
         *
        */

            for(int i=0;i<facetoSubdivide.size();++i){
                FacePointer f0= facetoSubdivide[i];
                f0->N()=facetoSubdivide[i]->N();
                FacePointer f1= &*lastf; ++lastf;
                FacePointer f2= &*lastf; ++lastf;
                FacePointer f3= &*lastf; ++lastf;
                VertexPointer v0 =facetoSubdivide[i]->V(0);
                VertexPointer v1 =facetoSubdivide[i]->V(1);
                VertexPointer v2 =facetoSubdivide[i]->V(2);
                VertexPointer mp01 = &*lastv; ++lastv;
                VertexPointer mp12 = &*lastv; ++lastv;
                VertexPointer mp02 = &*lastv; ++lastv;
                f0->V(0) = v0;   f0->V(1) = mp01; f0->V(2) = mp02;
                f1->V(0) = v1;   f1->V(1) = mp12; f1->V(2) = mp01;
                f2->V(0) = v2;   f2->V(1) = mp02; f2->V(2) = mp12;
                f3->V(0) = mp12; f3->V(1) = mp02; f3->V(2) = mp01;
                f1->N()=f0->N();
                f2->N()=f0->N();
                f3->N()=f0->N();
                if(facetoSubdivide[i]->IsS()){
                      //f0->SetS();
                      f1->SetS();
                      f2->SetS();
                      f3->SetS();
                }
                mp01->P() = (v0>v1) ? (v0->P()+v1->P())/2.0 : (v1->P()+v0->P())/2.0;
                mp01->N() = (v0>v1) ? (v0->N()+v1->N())/2.0 : (v1->N()+v0->N())/2.0;

                mp12->P() = (v1>v2) ? (v1->P()+v2->P())/2.0 : (v2->P()+v1->P())/2.0;
                mp12->N() = (v1>v2) ? (v1->N()+v2->N())/2.0 : (v2->N()+v1->N())/2.0;

                mp02->P() = (v0>v2) ? (v0->P()+v2->P())/2.0 : (v2->P()+v0->P())/2.0;
                mp02->N() = (v0>v2) ? (v0->N()+v2->N())/2.0 : (v2->N()+v0->N())/2.0;
             }

            int vd = tri::Clean<MeshType>::RemoveDuplicateVertex(m);
            printf("Vertex unification %i\n",vd);
            int vu = tri::Clean<MeshType>::RemoveUnreferencedVertex(m);
            //tri::Clean<MeshType>::RemoveNonManifoldFace(m);
            printf("Vertex unref %i\n",vu);
            UpdateTopology<MeshType>::FaceFace(m);
            //UpdateNormal<MeshType>::PerVertex(m);
            //UpdateNormal<MeshType>::PerFace(m);
            //Allocator<MeshType>::CompactEveryVector(m);
        }
    }
    template<class MeshType>
    void correctBadSelectedVertices(MeshType & m){
        typedef typename MeshType::FaceIterator FaceIterator;
        typedef typename MeshType::VertexIterator VertexIterator;
        typedef typename MeshType::VertexPointer VertexPointer;
        typedef typename MeshType::FacePointer FacePointer;
        typedef typename MeshType::FaceType FaceType;
        UpdateFlags<MeshType>::VertexClearV(m);
        //vector<FacePointer> facetoSubdivide;
        UpdateTopology<MeshType>::FaceFace(m);
        for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++)
        {
            if(!(*fi).IsD() && (*fi).IsS())
            {
                for(int k=0;k<3;k++){
                    if(!fi->V(k)->IsV()){
                        face::Pos<FaceType> startPos(&*fi,k,fi->V(k));
                        face::Pos<FaceType> copy=startPos;
                        do{
                            if(copy.IsManifold()){
                                copy.V()->C()=Color4b::Blue;
                                copy.VFlip()->C()=Color4b::Blue;
                                cout<<"edge maniold in tri?"<<endl;
                                //tri::io::ExporterPLY<MeshType>::Save(m,"Mesh-ManifoldEdge.ply",true,tri::io::Mask::IOM_ALL);
                                return;
                            }
                            copy.NextE();
                        }while(copy!=startPos);
                        vector<face::Pos<FaceType>> posVec;
                        vcg::face::VFOrderedStarFF(startPos,posVec);
                        // we search by 2 non-sucessives unselected triangles
                        int changes=0;
                        int psize=posVec.size();
                        for(int i=0;i<psize;i++){
                            if(posVec[i].F()->IsS()!=posVec[(i+1)%psize].F()->IsS())
                                changes++;
                        }

                        if(changes>2){
                            face::Pos<FaceType> rPos=startPos;
                            do{
                                rPos.F()->SetS();
                                rPos.NextE();
                            }while(rPos!=startPos);

                        }
                        fi->V(k)->SetV();
                      }                      
                 }
             }
         }
    }
}

#endif // MYUTILS_H
