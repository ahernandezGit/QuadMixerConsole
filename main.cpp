#include <iostream>
#include <ctime>
#include <chrono>
#include <CGAL/config.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include "meshes.h"

using namespace std;
enum ActiveMesh{MESHONE=0,MESHTWO};

void initPolyMesh(string message,ActiveMesh activemesh,meshes & twomeshes){
    if(activemesh==MESHONE){
        twomeshes.initMesh1();
        //Check if the mesh is valid
        bool isPureQuad=true;
        for(PMesh::FaceIterator fi=twomeshes.mesh1.face.begin();fi!=twomeshes.mesh1.face.end();++fi) if(!fi->IsD())
        {
         if(fi->VN()!=4)
           isPureQuad=false;
           break;
        }
        if(!isPureQuad){
            twomeshes.mesh1.Clear();
            cout<<"Mesh1 is not pure quad"<<endl;
            return;
        }
        int nonmanifoldedges=tri::Clean<PMesh>::CountNonManifoldEdgeFF(twomeshes.mesh1);
        if(nonmanifoldedges>0){
            twomeshes.mesh1.Clear();
            cout<<"Mesh1 has non manifold edges"<<endl;
            return;
        }
        twomeshes.visualizePatches1();
        CMesh temp;
        vcg::tri::PolygonSupport<CMesh,PMesh>::ImportFromPolyMesh(temp,twomeshes.mesh1);
        vcg::tri::UpdateBounding<CMesh>::Box(temp);
        vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(temp);
        vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(temp);
        vector<CMesh::FaceType*> retVec;
        tri::Clean<CMesh>::SelfIntersections(temp, retVec);
        if(retVec.size()>0){
            cout<<message+" Mesh1 has selfintersecting faces"<<endl;
            return;
        }

    }
    else{
        twomeshes.initMesh2();
        bool isPureQuad=true;
        for(PMesh::FaceIterator fi=twomeshes.mesh2.face.begin();fi!=twomeshes.mesh2.face.end();++fi) if(!fi->IsD())
        {
         if(fi->VN()!=4)
           isPureQuad=false;
           break;
        }
        if(!isPureQuad){
            twomeshes.mesh2.Clear();
            cout<<"Mesh2 is not pure quad"<<endl;
            return;
        }
        int nonmanifoldedges=tri::Clean<PMesh>::CountNonManifoldEdgeFF(twomeshes.mesh2);
        if(nonmanifoldedges>0){
            twomeshes.mesh2.Clear();
            cout<<"Mesh2 has non manifold edges"<<endl;
            return;
        }
        twomeshes.visualizePatches2();
        CMesh temp;
        vcg::tri::PolygonSupport<CMesh,PMesh>::ImportFromPolyMesh(temp,twomeshes.mesh2);
        vcg::tri::UpdateBounding<CMesh>::Box(temp);
        vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(temp);
        vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(temp);
        vector<CMesh::FaceType*> retVec;
        tri::Clean<CMesh>::SelfIntersections(temp, retVec);
        if(retVec.size()>0){
            cout<<message+" Mesh2 has selfintersecting faces"<<endl;
            return;
        }
    }
}
void loadPolyMesh(QString fileName,ActiveMesh activemesh, meshes & twomeshes){
    int err=0;
    // try to load the mesh from the file
    if(activemesh==MESHONE)
       err=twomeshes.loadmesh1(fileName);
    else
       err=twomeshes.loadmesh2(fileName);

    //cout<<"err= "<<err<<endl;
    // check if successfully loaded
    if ((err == 0) ||(err%2!=0))
        initPolyMesh("Loaded \""+fileName.toStdString()+"\n",activemesh,twomeshes);
}
int main(int argc, char *argv[])
{
    QString i1="../../../models/sphere2.obj";
    QString i2="../../../models/duck.obj";

    //assert(argc==3);

    meshes twomeshes;
    twomeshes.reset();
    loadPolyMesh(i1,MESHTWO,twomeshes);
    loadPolyMesh(i2,MESHONE,twomeshes);
    twomeshes.triangulate1();
    twomeshes.triangulate2();

    if(twomeshes.mesh1.FN()==0 || twomeshes.mesh2.FN()==0 ){
        cout<<"Missing one of the input quad meshes"<<endl;
        return 0;
    }
    clock_t c2, c1;
    twomeshes.triangulate1();
    twomeshes.triangulate2();
    if(tri::Clean<CMesh>::CountHoles(twomeshes.tmesh1)>0 ||
            tri::Clean<CMesh>::CountHoles(twomeshes.tmesh2)>0){
        cout<<"Input is not watertight"<<endl;
        return 0;
    }

    tri::Clean<CMesh>::RemoveDuplicateVertex(twomeshes.tmesh1);
    tri::Clean<CMesh>::RemoveDuplicateVertex(twomeshes.tmesh2);
    tri::Clean<CMesh>::RemoveUnreferencedVertex(twomeshes.tmesh1);
    tri::Clean<CMesh>::RemoveUnreferencedVertex(twomeshes.tmesh2);
    tri::Allocator<CMesh>::CompactVertexVector(twomeshes.tmesh1);
    tri::Allocator<CMesh>::CompactVertexVector(twomeshes.tmesh2);
    tri::Allocator<CMesh>::CompactFaceVector(twomeshes.tmesh1);
    tri::Allocator<CMesh>::CompactFaceVector(twomeshes.tmesh2);
    UpdateTopology<CMesh>::FaceFace(twomeshes.tmesh1);
    UpdateTopology<CMesh>::FaceFace(twomeshes.tmesh2);
    c1 = clock();
    //gettimeofday(&c1, NULL);
    Eigen::MatrixXd V1(twomeshes.tmesh1.vert.size(),3);
    Eigen::MatrixXd V2(twomeshes.tmesh2.vert.size(),3);
    Eigen::MatrixXi F1(twomeshes.tmesh1.face.size(),3);
    Eigen::MatrixXi F2(twomeshes.tmesh2.face.size(),3);
    Eigen::MatrixXd V3;
    Eigen::MatrixXi F3;
    Eigen::VectorXi J;
    utility::convertToLIBIGL<CMesh>(twomeshes.tmesh1,V1,F1);
    utility::convertToLIBIGL<CMesh>(twomeshes.tmesh2,V2,F2);
    //igl::readOFF("cheburashka.off",V1,F1);
    //igl::readOFF("decimated-knight.off",V2,F2);
    igl::copyleft::cgal::mesh_boolean(V1,F1,V2,F2,igl::MESH_BOOLEAN_TYPE_UNION,V3,F3,J);
    utility::libiglToVCG<CMesh>(twomeshes.tboolean,V3,F3);
    //tri::io::ExporterPLY<CMesh>::Save(twomeshes.tboolean,"parameterization/boolean.ply",tri::io::Mask::IOM_WEDGTEXCOORD,false);
    twomeshes.processBooleanMesh(J);
    twomeshes.visualizePatches1();
    twomeshes.visualizePatches2();
    twomeshes.processall();

    c2 = clock();
    float tempo = (c2-c1)*1000/CLOCKS_PER_SEC;
    cout<<"Processed in "+QString::number(tempo).toStdString()+" ms"<<endl;

    std::cout<<"finishing..."<<std::endl;
    return 1;
}
