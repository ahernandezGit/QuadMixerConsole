#ifndef VECTORCONE_H
#define VECTORCONE_H
#include "meshtypes.h"
#include <list>
#include <memory>
#include "graph.h"
#include <wrap/io_trimesh/export.h>
enum SideType{S_DIAGONAL,S_FIXED,S_UNDEFINED};

class VectorCone
{        
    const float epsilon=0.00000001;
public:

    VectorCone();
    // Constructor given an Initial vector
    VectorCone(const Point3f &v);

    // Returns a unit vector that bisects this cone
    Point3f bisector();

    // Set the cone to a unique vector
    void set(const Point3f & v);

    // Returns the current aperture of the cone
    float angle();

    // Adds a vector to the cone.
    void add(const Point3f &v);

    void printCone();

Point3f min,max;

};

namespace MonotonicFunctions {
    const float epsilon=0.00000001;

    struct DiagonalDS{
        int i,j;
        float imbalance,angularDeviation,size;
        void setDiag(int ii,int jj,float isize, float imb,float adesviation){
            i=ii;
            j=jj;
            size=isize;
            imbalance=imb;
            angularDeviation=adesviation;
        }
        void copy(DiagonalDS &diag){
            i=diag.i;
            j=diag.j;
            size=diag.size;
            imbalance=diag.imbalance;
            angularDeviation=diag.angularDeviation;
        }
    };
    /*
     * Model a side of a polygon(subpolygon) contained into another bigger polygon(great polygon) resulting of the subdivision by diagonals
     * The side is stored as a list of vertex pointer to the vertices of the point mesh formed by the points of the great polygon.
     * */
    struct PolygonSide{
        vector<Point3f> vertices;
        SideType flag;
        int size;
        int index_3D; // index of the 3d connected component side from what this side comes
        int id; // just used when this side is a diagonal
        PolygonSide(){
            vertices.clear();
            flag=S_UNDEFINED;
            size=0;
            index_3D=-1;
            id=-1;
        }
        PolygonSide(const vector<Point3f> &sd,SideType isv=S_UNDEFINED,int i3d=-1,int iid=-1){
            vertices=sd;
            flag=isv;
            size=sd.size();
            index_3D=i3d;
            id=iid;
        }
        bool getparity(){
            return (size-1)%2;
        }
        float length(){
            float result=0;
            for(int i=0;i<vertices.size()-1;i++)
                result+=vcg::Distance(vertices[i],vertices[i+1]);
            return result;
        }
        //reverse the ordering of this side
        void reverse(){
            std::reverse(vertices.begin(),vertices.end());
        }
    };
    struct Polygon{
        vector<PolygonSide> sides;
        int numberID,parent; //used for processing at higer levels
        vector<Polygon*> adjacency;

        Polygon(){
            sides.clear();
            adjacency.clear();
            numberID=-1;
        }
        Polygon(vector<PolygonSide> &isides,int nid=-1){
            //verify coherence
            int n=isides.size();
            for(int i=0;i<n;i++){
               assert(isides[i].vertices.back()==isides[(i+1)%n].vertices[0]);
            }
            sides=isides;
            numberID=nid;
        }

        vector<Point3f> get_as_points(vector<int> &feature){
            feature.clear();
            vector<Point3f> result;
            result.push_back(sides[0].vertices[0]);
            feature.push_back(0);
            for(int i=0;i<sides.size();i++){
                for(int j=1;j<sides[i].size;j++)
                    result.push_back(sides[i].vertices[j]);
                feature.push_back(result.size()-1);
            }
            result.pop_back();
            feature.pop_back();
            return result;
        }
        vector<Point3f> get_as_points(){
            vector<Point3f> result;
            result.push_back(sides[0].vertices[0]);
            for(int i=0;i<sides.size();i++){
                for(int j=1;j<sides[i].size;j++)
                    result.push_back(sides[i].vertices[j]);
            }
            result.pop_back();
            return result;
        }
        void addConnection(Polygon &poly){
            adjacency.push_back(&poly);
        }
        int size(){
            int nsize=0;
            for(int i=0;i<sides.size();i++)
                nsize+=(sides[i].vertices.size()-1);
            return nsize;
        }
        int numbersides(){
            return sides.size();
        }
        bool getConstant(){
            int result=0;
            for(int i=0;i<sides.size();i++)
               if(sides[i].flag!=S_DIAGONAL)
                    result+=(sides[i].vertices.size()-1);

            return result%2;
        }
        float averageEdgeLength_fixed(){
            float result=0;
            int count=0;
            for(int i=0;i<sides.size();i++)
               if(sides[i].flag ==S_FIXED){
                    count+=(sides[i].vertices.size()-1);
                    result+=sides[i].length();
               }
            if(count==0)
                return 0;
            return result/count;
        }
        float averageEdgeLength(){
            float result=0;
            int count=0;
            for(int i=0;i<sides.size();i++){
                count+=(sides[i].vertices.size()-1);
                result+=sides[i].length();
            }
            if(count==0)
                return 0;
            return result/count;
        }

        //reverse side j
        void reverse(int j){
            assert(j>=0 && j<sides.size());
            sides[j].reverse();
        }
        void saveAsPLY(string name){
            vector<Point3f> data=get_as_points();
            CMesh mesh;
            int n=data.size();
            Allocator<CMesh>::AddVertices(mesh,n);
            CMesh::VertexIterator vi=mesh.vert.begin();
            for(int i=0;i<n;i++){
               (*vi).P()=data[i];
               if(i<n-1)
                  ++vi;
            }
            for(int i=0;i<n;i++){
               Allocator<CMesh>::AddEdge(mesh,&mesh.vert[i],&mesh.vert[(i+1)%n]);
            }
            string filename="../exported/polygonVectorCone_"+QTime::currentTime().toString().toStdString()+"-"+name+".ply";
            tri::io::ExporterPLY<CMesh>::Save(mesh,filename.c_str(),tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
        }
        Polygon subdivideDiagonals(set<int> &added){
            vector<PolygonSide> newsides;
            float edgelen=averageEdgeLength_fixed();
            edgelen/=4;   // minimum edge lenght to sample the diagonals
            int count=0;   // store point index in the new polygon
            for(int i=0;i<sides.size();i++){
               if(sides[i].flag ==S_DIAGONAL){
                   float lenght_diagonal=vcg::Distance(sides[i].vertices[0],sides[i].vertices[1]);
                   int number_divisions=int(lenght_diagonal/edgelen+1);
                   float step=lenght_diagonal/number_divisions;
                   count++; // vertices[0]
                   if(number_divisions>1){
                     int number_newpoints=number_divisions-1;
                     vector<Point3f> newpts;
                     Point3f direction=sides[i].vertices[1]-sides[i].vertices[0];
                     direction.normalized();
                     for(int j=1;j<=number_newpoints;j++){
                         Point3f aux=sides[i].vertices[0]+direction*(j*step);
                         aux.Z()=0;
                         newpts.push_back(aux);
                         count++;
                         added.insert(count);
                     }
                     count++;// vertices[1]
                    newpts.push_back(sides[i].vertices[1]);
                    newpts.insert(newpts.begin(),sides[i].vertices[0]);
                    PolygonSide tempside(newpts,sides[i].flag,sides[i].index_3D,sides[i].id);
                    newsides.push_back(tempside);

                   }
                   else{
                       newsides.push_back(sides[i]);
                       count+=(sides[i].vertices.size()-1);
                   }
               }
               else{
                  newsides.push_back(sides[i]);
                  count+=(sides[i].vertices.size()-1);
               }
            }
            Polygon result(newsides,numberID);
            return result;
        }
    };
    // Compute a monotone partition of a curve so that
    // each curve segment does not span a cone with angle greater than ang
    vector<int> monotonePartition(vector<Point3f> & pts, float ang=60, int wherestart=0);

    // Partitions a curve in monotonic chains and returns them as an array of curves
    vector<vector<int>> monotonicPartitionCurves(vector<Point3f> curvePoints, float ang, int wherestart=0);

    /**
     * Finds a set of line segments that splits the polygon given by curve
     * using the vertices of partition as starting points.
     * @param {Curve} curve
     * @param {Number[]} partition
     */
    DiagonalDS diagonals(vector<Point3f> &curve,vector<int> &partition,float sw=4,float aw=6,float iw=3, int index=-1);

    /**
     * Check if a diagonal of a polygon(closed) is completely in the polygon.
     */
    bool isDiagonalInterior(int id1,int id2,vector<Point3f> &polygon);

    vector<vector<Point3f>> splitCurve(vector<Point3f> &polygon, float limitAngle=60, int maxPartition=5,float sw=4,float aw=6,float iw=3);

    vector<vector<Point3f>> recursiveSplitCurve(vector<Point3f> &polygon, float limitAngle=60, int maxPartition=5, float sw=4,float aw=6,float iw=3);

    // This class represents a 2D polygon that we want to subdivide into subpatches feasible to be quadrangulated with the Takayama's method
    // We trace certain indexes of the input polygon to robustly mark side's endpoints of the final quadrangulating patches.
    class PrePatch2D
    {
    public:
        //Vertex Indexes corresponding to the endpoints of diagonals
        std::set<int> precorners;

        // Diagonals
        vector<std::pair<int,int>> diagonalset;
        // Initial polygon to subdivide. Mesh only with vertices
        //This is used in order to get a fixed id for the precorners and diagonal
        //CMesh initial_polygon;

        // Graph storing the adjacent between the final subpatches
        myutils::Graph graph;

        //list of polygons contains into this prepatch
        vector<Polygon> allpolygons;

        //Add a vertex to the list of precorners
        void add_precorner(vector<int>& ipc);

        PrePatch2D();        
        PrePatch2D(Polygon& polygon);
        void set(Polygon& polygon);
        void setnull();
        vector<vector<CMesh::VertexPointer>> splitCurve(vector<CMesh::VertexPointer> &polygon, float limitAngle=60, int maxPartition=4,float sw=4,float aw=6,float iw=3);
        vector<vector<CMesh::VertexPointer>> recursiveSplitCurve(vector<CMesh::VertexPointer>& polygon,float limitAngle=60, int maxPartition=5, float sw=4,float aw=6,float iw=3);

        vector<Polygon> splitCurve(Polygon &poly, float limitAngle=60, int maxPartition=4,float sw=4,float aw=6,float iw=3);
        vector<Polygon> recursiveSplitCurve(Polygon& polygon,float limitAngle=45, int maxPartition=4, float sw=4,float aw=6,float iw=6);


        // return the final polygon subdivision as a collection of ordered list of sides, i.e , each prepatch is the form
        // poly=[points in side 1, points in side 2,...,points in side n].
        vector<vector<vector<Point3f>>> getPolygonsAsPrePatch();
        void processPrePatch(Polygon& poly);

        void constructGraph(vector<vector<CMesh::VertexPointer>> &polygons);

        int isOnFeatures(const vector<int> &feature,int ind);
    };
}

#endif // VECTORCONE_H
