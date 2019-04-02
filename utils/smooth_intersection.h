#ifndef SMOOTH_INTERSECTION_H
#define SMOOTH_INTERSECTION_H

#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/complex/algorithms/smooth.h>

template <class MeshType>
class SmoothIntersection
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;

    static void ReprojectStep(MeshType &mesh,
                              MeshType &target,
                              vcg::GridStaticPtr<FaceType,ScalarType> &Grid)
    {
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (!mesh.vert[i].IsS())continue;
            CoordType ClosePt;
            ScalarType maxD,minD;
            maxD=mesh.bbox.Diag();
            FaceType *f=NULL;
            f=vcg::tri::GetClosestFaceBase(target,Grid,mesh.vert[i].P(),maxD,minD,ClosePt);
            assert(f!=NULL);
            mesh.vert[i].P()=ClosePt;
        }
    }

    static void LaplacianStep(MeshType &mesh,
                              bool OnlyIntersection,
                              ScalarType Damp)
    {
        assert(Damp>=0);
        assert(Damp<=1);
        //cumulate
        std::vector<CoordType> AvgPos(mesh.vert.size(),CoordType(0,0,0));
        std::vector<size_t> NumPos(mesh.vert.size(),0);

        for (size_t i=0;i<mesh.face.size();i++)
        {
            std::vector<size_t> IndexVF;
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                size_t IndexV=vcg::tri::Index(mesh,mesh.face[i].V(j));
                IndexVF.push_back(IndexV);
            }

            for (size_t j=0;j<IndexVF.size();j++)
            {
                size_t IndexV=vcg::tri::Index(mesh,mesh.face[i].V(j));

                int Q0=mesh.face[i].V(j)->Q();
                CoordType Pos0=mesh.vert[IndexV].P();

                for (size_t j=0;j<IndexVF.size();j++)
                {

                    if (IndexVF[j]==IndexV)continue;//same one

                    int Q1=mesh.face[i].V(j)->Q();
                    if ((OnlyIntersection)&&((Q0<0)||(Q1<0)))continue;

                    //if smooth the line check the qualities should be one after the other and !=-1;
                    if ((OnlyIntersection)&&(fabs(Q0-Q1)!=1))continue;
                    AvgPos[IndexVF[j]]+=Pos0;
                    NumPos[IndexVF[j]]++;
                }
            }
        }
        //average
        for (size_t i=0;i<AvgPos.size();i++)
        {
            int Q=mesh.vert[i].Q();
            if ((OnlyIntersection)&&(Q<0))continue;//internal
            if ((!OnlyIntersection)&&(Q>=0))continue;//border
            if (NumPos[i]==0)continue;
            if (!mesh.vert[i].IsS())continue;
            mesh.vert[i].P()=AvgPos[i]/NumPos[i]*(1-Damp)+mesh.vert[i].P()*Damp;
        }
    }

 public:

    static void Smooth(MeshType &mesh,ScalarType Kernel=0.05,int Steps=10)
    {
        //save the original quality
        std::vector<int> QInit;
        std::vector<VertexType*> seeds;
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            QInit.push_back(mesh.vert[i].Q());
            if (mesh.vert[i].Q()==-1)continue;
            seeds.push_back(&mesh.vert[i]);
            //mesh.vert[i].SetS();
        }
        vcg::tri::Geodesic<MeshType>::Compute(mesh,seeds);

        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);
        vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);

        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (mesh.vert[i].Q()>Kernel*mesh.bbox.Diag())continue;
            if (mesh.vert[i].IsB())continue;
            mesh.vert[i].SetS();
        }
        for (size_t i=0;i<mesh.vert.size();i++)
            mesh.vert[i].Q()=QInit[i];

        vcg::tri::Smooth<MeshType>::VertexCoordLaplacian(mesh,Steps,true);

        //make a copy of the mesh
        MeshType TempMesh;
        vcg::tri::Append<MeshType,MeshType>::Mesh(TempMesh,mesh);
        vcg::GridStaticPtr<FaceType,ScalarType> Grid;
        vcg::tri::UpdateBounding<MeshType>::Box(TempMesh);
        vcg::tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFaceNormalized(TempMesh);
        Grid.Set(TempMesh.face.begin(),TempMesh.face.end());

        //smooth the intersection
        for (size_t i=0;i<Steps;i++)
        {
            LaplacianStep(mesh,true,0.1);
            ReprojectStep(mesh,TempMesh,Grid);
            LaplacianStep(mesh,false,0.1);
            ReprojectStep(mesh,TempMesh,Grid);
        }

    }
    static void Smooth(MeshType &mesh, vector<vector<vcg::face::Pos<typename MeshType::FaceType>>> &curves,ScalarType Kernel=0.05,int Steps=10)
    {
        // Storing vertices in the intersection curves
        std::vector<VertexType*> seeds;
        size_t numberCurves=curves.size();
        for (size_t i=0;i<numberCurves;i++)
        {
            size_t curveSize=curves[i].size();
            for (size_t j=0;j<curveSize;j++){
                seeds.push_back(curves[i][j].V());
            }
            //mesh.vert[i].SetS();
        }
        vcg::tri::Geodesic<MeshType>::Compute(mesh,seeds);

        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);
        vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);

        // Selecting vertices in close to the intersection curve
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (mesh.vert[i].Q()>Kernel*mesh.bbox.Diag())continue;
            if (mesh.vert[i].IsB())continue;
            mesh.vert[i].SetS();
        }

        vcg::tri::Smooth<MeshType>::VertexCoordLaplacian(mesh,Steps,true);

        //make a copy of the mesh
        MeshType TempMesh;
        vcg::tri::Append<MeshType,MeshType>::Mesh(TempMesh,mesh);
        vcg::GridStaticPtr<FaceType,ScalarType> Grid;
        vcg::tri::UpdateBounding<MeshType>::Box(TempMesh);
        vcg::tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFaceNormalized(TempMesh);
        Grid.Set(TempMesh.face.begin(),TempMesh.face.end());

        //smooth the intersection
        for (size_t i=0;i<Steps;i++)
        {
            LaplacianStep(mesh,true,0.1);
            ReprojectStep(mesh,TempMesh,Grid);
            LaplacianStep(mesh,false,0.1);
            ReprojectStep(mesh,TempMesh,Grid);
        }

    }
};

#endif
