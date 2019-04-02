#ifndef PARAMETERIZATIONLSCM_H
#define PARAMETERIZATIONLSCM_H

#include "meshtypes.h"
#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <wrap/igl/lscm_parametrization.h>
#include <vcg/space/intersection2.h>

class LSCMparameterization
{
    typedef typename CMesh::FaceType FaceTypeT;
    typedef typename CMesh::ScalarType ScalarType;
    typedef typename CMesh::FacePointer FacePointer;
    typedef typename CMesh::FaceIterator FaceIterator;
    typedef typename CMesh::VertexPointer VertexPointer;
    typedef typename CMesh::FaceType FaceType;
    typedef typename CMesh::VertexIterator VertexIterator;
    typedef typename CMesh::CoordType CoordType;
    //typedef typename CMesh:: template PerFaceAttributeHandle<int> AttributeFaceType;

    const size_t invalidindex=std::numeric_limits<size_t>::max();

private:

    //CMesh mesh; //it does not work. Problem related to ill-formed constructor, this class is implicitily deleted if we add this member

    // reference mesh where we got the portion to be parameterized
    CMesh *reference_mesh;

    // mapping between the parameterized vertices list and the vertices indexes of the reference mesh
    vector<size_t> map_mesh_to_reference;

    // mapping between the vertices indexes of the reference mesh and the parameterized vertices list
    unordered_map<size_t,size_t> map_reference_to_mesh;

    // uv coordinates for the list of parameterized vertices
    vector<vcg::Point2<ScalarType>> uv;

    // faces into the portion of the mesh to be parameterized
    vector<size_t> faces;


public:

    LSCMparameterization(){
        reset();
    }
    void reset(){
        reference_mesh=NULL;
        map_mesh_to_reference.clear();
        map_reference_to_mesh.clear();
        uv.clear();
        faces.clear();
    }
    // minuv and maxuv are the indexes in the reference mesh that we will set as (0,0) and (1,1)
    LSCMparameterization(CMesh & ref,vector<size_t> &ifaces,size_t minuv, size_t maxuv){
        initialize(ref,ifaces,minuv,maxuv);
    }
    void initialize(CMesh &ref,vector<size_t> &ifaces,size_t minuv, size_t maxuv){
        reference_mesh=&ref;
        map_mesh_to_reference.clear();
        map_reference_to_mesh.clear();
        uv.clear();
        faces=ifaces;
        //
        // Creating the portion of the reference mesh corresponding to the portion
        //

        // storing initial selected vertices and faces
        vector<size_t> initial_selected_vertices,initial_selected_faces;
        selectionVertices(*reference_mesh,initial_selected_vertices,false);
        selectionFaces(*reference_mesh,initial_selected_faces,false);
        vcg::tri::UpdateSelection<CMesh>::VertexClear(*reference_mesh);
        vcg::tri::UpdateSelection<CMesh>::FaceClear(*reference_mesh);

        // selecting the faces and vertices corresponding to the portion of the mesh to be parameterize
        for(size_t i=0;i<faces.size();i++){
            reference_mesh->face[faces[i]].SetS();
        }
        vcg::tri::UpdateSelection<CMesh>::VertexFromFaceStrict(*reference_mesh);
        CMesh mesh;
        vcg::tri::Append<CMesh,CMesh>::Mesh(mesh,*reference_mesh,true);

        //creating the mapping between the new mesh and the reference
        for(VertexIterator vi=reference_mesh->vert.begin(); vi!=reference_mesh->vert.end(); ++vi){
          if(!(*vi).IsD() && (*vi).IsS()){
            size_t ind=vcg::tri::Index(*reference_mesh,*vi);
            map_mesh_to_reference.push_back(ind);
            size_t i=map_mesh_to_reference.size()-1;
            map_reference_to_mesh.emplace(ind,i);
            if(vcg::Distance(vi->P(),mesh.vert[i].P())>0.0001){
                cout<<"LSCMParameterization : points in list are different from 3D"<<endl;
            }
          }
        }
        //applying the least square conformal map parameterization
        vcg::tri::UpdateFlags<CMesh>::VertexClearS(mesh);
        mesh.vert[map_reference_to_mesh[minuv]].T().P()=vcg::Point2<ScalarType>(0,0);
        mesh.vert[map_reference_to_mesh[maxuv]].T().P()=vcg::Point2<ScalarType>(1.0,1.0);
        mesh.vert[map_reference_to_mesh[minuv]].SetS();
        mesh.vert[map_reference_to_mesh[maxuv]].SetS();
        vcg::tri::OptimizeUV_LSCM<CMesh>(mesh,CMesh::VertexType::SELECTED);

        //storing the uv coordinates
        size_t numberVertices=size_t(mesh.VN());
        uv.resize(numberVertices);
        for(size_t i=0;i<numberVertices;i++){
            uv[i]=mesh.vert[i].T().P();
        }
        selectionVertices(*reference_mesh,initial_selected_vertices,true);
        selectionFaces(*reference_mesh,initial_selected_faces,true);

    }
    CoordType reproject(ScalarType u,ScalarType v){
        // search in what triangle this (u,v) coordinate is
        Point2<ScalarType> p(u,v);
        Point3<ScalarType> L(0,0,0);
        for(size_t i=0;i<faces.size();i++){
            FaceType fi=reference_mesh->face[faces[i]];
            size_t i0=map_reference_to_mesh[tri::Index(*reference_mesh,fi.V(0))];
            size_t i1=map_reference_to_mesh[tri::Index(*reference_mesh,fi.V(1))];
            size_t i2=map_reference_to_mesh[tri::Index(*reference_mesh,fi.V(2))];
            Triangle2<ScalarType> triangle(uv[i0],uv[i1],uv[i2]);
            if(vcg::IsInsideTrianglePoint(triangle,p)){
                  vcg::InterpolationParameters2(triangle.P(0),triangle.P(1),triangle.P(2),p,L);
                  CoordType v0_3D=fi.V(0)->P();
                  CoordType v1_3D=fi.V(1)->P();
                  CoordType v2_3D=fi.V(2)->P();
                  CoordType result=v0_3D*L.X()+v1_3D*L.Y()+v2_3D*L.Z();
                  return result;
            }
        }
        cout<<"reprojecting fails"<<endl;
        return CoordType(0,0,0);
    }
    void reprojecting_and_update(Eigen::MatrixXd &points){
        assert(points.cols()==3);
        for(int i=0;i<points.rows();i++){
            /*if(i==8){
                cout<<"point 8 "<<aux.X()<<" "<<aux.Y()<<endl;
            }*/
            CMesh::CoordType projected=reproject(ScalarType(points(i,0)),ScalarType(points(i,1)));
            /*if(i==99 || i==306){
                cout<<"i :"<<i<<", value :"<<aux.X()<<" , "<<aux.Y()<<endl;
                cout<<"projected :"<<projected.X()<<" , "<<projected.Y()<<" , "<<projected.Y()<<endl;
                projected=unproject2D(aux);
            }*/
            if(projected.X()==NAN || projected.Y()==NAN || projected.Z()==NAN){
                cout<<"NAN appears "<<endl;
            }
            points.row(i)=Eigen::Vector3d(projected.X(),projected.Y(),projected.Z());
        }
    }
    void reprojecting_and_update(PMesh &mesh){
        for(size_t i=0;i<mesh.vert.size();i++){
            CMesh::CoordType projected=reproject(mesh.vert[i].P().X(),mesh.vert[i].P().Y());

            if(projected.X()==NAN || projected.Y()==NAN || projected.Z()==NAN){
                cout<<"NAN appears "<<endl;
            }
            mesh.vert[i].P()=projected;
        }
    }
    Point2<ScalarType> project(face::Pos<FaceType> &hedge){
        size_t index=tri::Index(*reference_mesh,hedge.V());
        assert(index<size_t(reference_mesh->vert.size()) && index>=0);
        size_t id=map_reference_to_mesh[index];
        return uv[id];
    }
    void exportMeshParametrization(CMesh &result){
       assert(!uv.empty());
       vector<Point3i> faceindexes;
       vector<CoordType> coordinates2D;
       for(size_t i=0;i<faces.size();i++){
           FaceType fi=reference_mesh->face[faces[i]];
           size_t v0=map_reference_to_mesh[tri::Index(*reference_mesh,fi.V(0))];
           size_t v1=map_reference_to_mesh[tri::Index(*reference_mesh,fi.V(1))];
           size_t v2=map_reference_to_mesh[tri::Index(*reference_mesh,fi.V(2))];
           faceindexes.push_back(Point3i(v0,v1,v2));
       }
       for(size_t i=0;i<uv.size();i++){
          vcg::Point2<ScalarType> point=uv[i];
          coordinates2D.push_back(CoordType(point.X(),point.Y(),0.0));
       }
       tri::BuildMeshFromCoordVectorIndexVector<CMesh,CoordType,Point3i>(result,coordinates2D,faceindexes);
       vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(result);
       vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(result);
       UpdateTopology<CMesh>::FaceFace(result);
    }
    void selectionVertices(CMesh &m, vector<size_t> &selected, bool select=false){
        if(select){
            assert(!selected.empty());
            vcg::tri::UpdateSelection<CMesh>::VertexClear(m);
            for(size_t i=0;i<selected.size();i++){
                m.vert[selected[i]].SetS();
            }

        }
        else{
            selected.clear();
            for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi){
                if(!vi->IsD() && vi->IsS()){
                    selected.push_back(vcg::tri::Index(m,&*vi));
                }
            }
        }
    }
    void selectionFaces(CMesh &m, vector<size_t> &selected, bool select=false){
        if(select){
            assert(!selected.empty());
            vcg::tri::UpdateSelection<CMesh>::FaceClear(m);
            for(size_t i=0;i<selected.size();i++){
                m.face[selected[i]].SetS();
            }

        }
        else{
            selected.clear();
            for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi){
                if(!fi->IsD() && fi->IsS()){
                    selected.push_back(vcg::tri::Index(m,&*fi));
                }
            }
        }
    }


};

#endif // PARAMETERIZATIONLSCM_H
