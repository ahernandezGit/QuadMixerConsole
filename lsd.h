#ifndef LSD_H
#define LSD_H

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <limits>
#include <memory>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/update/selection.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <wrap/io_trimesh/export.h>
using namespace std;

template< class MeshType>
class LaplacianSurfaceDeformation{

    typedef typename MeshType::VertexIterator VertexIterator;
    typedef typename MeshType::VertexPointer VertexPointer;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef Eigen::SparseMatrix<ScalarType> SpMat;
    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> DMat;
    typedef Eigen::Triplet<ScalarType> Triplet;
    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> Vec;


    const size_t invalidindex=std::numeric_limits<size_t>::max();

public:
    enum TypeDeformation{CorrectLaplacian=0,OnlyLaplacian,ZeroLaplacian};

    //Default constructor
    LaplacianSurfaceDeformation(){
        mesh.Clear();
        IndexL.clear();
        ValuesL.clear();
        Ts.clear();
        basemesh=NULL;
        constrain_weight=100;
    }

    LaplacianSurfaceDeformation(MeshType & m, const std::vector<size_t> &roi, const std::vector<size_t> &handle){
        set(m,roi,handle);
    }

    void set(MeshType &m, const std::vector<size_t> &roi, const std::vector<size_t> &handle){
        mesh.Clear();
        map_vertices_base_to_mesh.clear();
        Handle.clear();
        IndexL.clear();
        ValuesL.clear();
        Ts.clear();
        basemesh=&m;
        ROI_base=roi;
        Handle_base=handle;
        constrain_weight=10000;
    }
    void InitSparse(const std::vector<std::pair<int,int> > &Index,
                               const std::vector<ScalarType> &Values,
                               const int m,
                               const int n,
                               Eigen::SparseMatrix<ScalarType>& X)
        {
            assert(Index.size()==Values.size());

            std::vector<Eigen::Triplet<ScalarType> > IJV;
            IJV.reserve(Index.size());

            for(size_t i= 0;i<Index.size();i++)
            {
                int row=Index[i].first;
                int col=Index[i].second;
                ScalarType val=Values[i];

                assert(row<m);
                assert(col<n);

                IJV.push_back(Eigen::Triplet<ScalarType>(row,col,val));
            }
            X.resize(m,n);
            X.setFromTriplets(IJV.begin(),IJV.end());
    }
    void computeTi(){
        // computing T_i for each vertex in ROI
        //see slide 51 of https://www.cse.wustl.edu/~taoju/cse554/lectures/lect08_Deformation.pd

        Ts.clear();
        Neighbours.clear();
        Ts.resize(mesh.vert.size());
        Neighbours.resize(mesh.vert.size());
        for (int i = 0; i <mesh.vert.size(); ++i) {
            // set of {i} and the neigbbours of i.

            std::vector<size_t> iAndNeighbours;
            vector<VertexPointer> startvec;
            vcg::face::VVStarVF<FaceType>(&mesh.vert[i],startvec);
            iAndNeighbours.push_back(i);
            for (int j = 0; j < startvec.size(); ++j) {
                iAndNeighbours.push_back(vcg::tri::Index(mesh,startvec[j]));
            }
            Neighbours[i]=iAndNeighbours;

            DMat At;
            At.setZero(7, iAndNeighbours.size() * 3);

            for (int j = 0; j < iAndNeighbours.size(); ++j) {
                int k = iAndNeighbours[j];
                CoordType point=mesh.vert[k].P();

                double vk[3];
                vk[0] = ScalarType(point.X());
                vk[1] = ScalarType(point.Y());
                vk[2] = ScalarType(point.Z());

                const int x = 0;
                const int y = 1;
                const int z = 2;

                At(0, j * 3 + 0) = +vk[x];
                At(1, j * 3 + 0) = 0;
                At(2, j * 3 + 0) = +vk[z];
                At(3, j * 3 + 0) = -vk[y];
                At(4, j * 3 + 0) = +1;
                At(5, j * 3 + 0) = 0;
                At(6, j * 3 + 0) = 0;

                At(0, j * 3 + 1) = +vk[y];
                At(1, j * 3 + 1) = -vk[z];
                At(2, j * 3 + 1) = 0;
                At(3, j * 3 + 1) = +vk[x];
                At(4, j * 3 + 1) = 0;
                At(5, j * 3 + 1) = +1;
                At(6, j * 3 + 1) = 0;

                At(0, j * 3 + 2) = +vk[z];
                At(1, j * 3 + 2) = +vk[y];
                At(2, j * 3 + 2) = -vk[x];
                At(3, j * 3 + 2) = 0;
                At(4, j * 3 + 2) = 0;
                At(5, j * 3 + 2) = 0;
                At(6, j * 3 + 2) = 1;
            }

            DMat invprod = (At * At.transpose()).inverse();
            DMat pseudoinv = invprod * At;
            Ts[i] = pseudoinv;
            // Ts[i] now contains (A^T A ) A^T (see equation 12 from paper.)
        }

    }
    void computeRoiDelta(){
        DMat b;
        b.setZero(mesh.vert.size(),3);
        vector<size_t> index_zero;
        for(size_t i=0;i<mesh.vert.size();i++){
            CoordType point=mesh.vert[i].P();
            b.row(i)=Eigen::Matrix<ScalarType,3,1>(point.X(),point.Y(),point.Z());
            if(point.X()==0 && point.Y()==0)
                index_zero.push_back(i); // to correct bad sparse multiplication: it gives nan values
        }
        // computing &_i for each vertex in ROI
        //cout<<"computing &_i for each vertex in ROI"<<endl;
        cout<<"computing roidelta"<<endl;
        //cout<<b<<endl;
        roiDelta=L*b;

    }
    void prepareDeform(){
        //we extract submesh from original
        vcg::tri::UpdateSelection<MeshType>::Clear(*basemesh);
        for(size_t i=0;i<ROI_base.size();i++){
            if(!basemesh->vert[ROI_base[i]].IsD()){
                basemesh->vert[ROI_base[i]].SetS();
            }
        }
        map_vertices_base_to_mesh.resize(basemesh->vert.size(),invalidindex);
        size_t i=0;
        for(VertexIterator vi=basemesh->vert.begin(); vi!=basemesh->vert.end(); ++vi){
           if(!(*vi).IsD() && (*vi).IsS()){
             size_t ind=vcg::tri::Index(*basemesh,*vi);
             map_vertices_base_to_mesh[ind]=i;
             i++;
           }
        }
        vcg::tri::UpdateSelection<MeshType>::FaceFromVertexStrict(*basemesh);
        vcg::tri::Append<MeshType,MeshType>::Mesh(mesh,*basemesh,true);
        vcg::tri::Allocator<MeshType>::CompactEveryVector(mesh);
        vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
        vcg::tri::UpdateTopology<MeshType>::VertexFace(mesh);
        vcg::tri::UpdateBounding<MeshType>::Box(mesh);
        vcg::tri::UpdateFlags<MeshType>::FaceBorderFromFF(mesh);

        /*for(size_t i=0;i<map_vertices_base_to_mesh.size();i++){
            if(map_vertices_base_to_mesh[i]!=invalidindex){
                size_t aux=map_vertices_base_to_mesh[i];
                cout<<" basemesh: "<<basemesh->vert[i].P().X()<<" "<<basemesh->vert[i].P().Y()<<" "<<basemesh->vert[i].P().Z()<<endl;
                cout<<" mesh: "<<mesh.vert[aux].P().X()<<" "<<mesh.vert[aux].P().Y()<<" "<<mesh.vert[aux].P().Z()<<endl;
            }
        }*/

        size_t matr_size=mesh.vert.size(); // number of elements in the ROI

        //get the entries for cotangen laplacian matrix

        bool usecotanWeight=true;
        MeshToMatrix<MeshType>::GetLaplacianMatrix(mesh,IndexL,ValuesL,usecotanWeight,1,false);

        //initialize sparse laplacian matrix
        //cout<<"initialize sparse laplacian matrix"<<endl;
        InitSparse(IndexL,ValuesL,matr_size,matr_size,L);
        /*computeRoiDelta();
        if(correctLaplacian)
            computeTi();
        */
        //cout<<"L matrix "<<endl;
        //cout<<L<<endl;
        // Computing boundary handles and adding this to the total handles
        //cout<<"Computing boundary handles and adding this to the total handles"<<endl;
        size_t constrain_size=Handle_base.size();
        Handle.clear();
        for(size_t i=0;i<constrain_size;i++)
            Handle.push_back(map_vertices_base_to_mesh[Handle_base[i]]);

        //cout<<"computing boundary start"<<endl;
            //boundary handles

            vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);
            vcg::tri::UpdateFlags<MeshType>::FaceBorderFromFF(mesh);
            vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
            vcg::tri::UpdateSelection<MeshType>::VertexFromBorderFlag(mesh);
            vcg::tri::UpdateSelection<MeshType>::FaceFromVertexLoose(mesh);
            vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);
            vcg::tri::UpdateSelection<MeshType>::FaceFromVertexStrict(mesh);
            /*for(int k=0;k<2;k++)
                vcg::tri::UpdateSelection<MeshType>::FaceDilate(mesh);
            vcg::tri::UpdateSelection<MeshType>::VertexFromFaceStrict(mesh);
            for(VertexIterator vi=mesh.vert.begin();vi!=mesh.vert.end();vi++){
                if(!vi->IsD() && vi->IsS()){
                    Handle.push_back(vcg::tri::Index(mesh,&*vi));
                    vi->C()=vcg::Color4b::Blue;
                }
            }*/
            //cout<<"computing boundary finish"<<endl;
        //cout<<constrain_size<<" "<<Handle.size()<<endl;
        //vcg::tri::io::ExporterPLY<MeshType>::Save(mesh,"roi.ply",vcg::tri::io::Mask::IOM_VERTCOLOR+vcg::tri::io::Mask::IOM_VERTQUALITY,false);

    }
    void doDeform(vector<CoordType> &handle_positions,TypeDeformation correctLaplacian){
        assert(handle_positions.size()==Handle_base.size());        

        DMat positions,newpositions;
        newpositions.setZero(mesh.vert.size(),3);
        for(size_t i=0;i<mesh.vert.size();i++){
            CoordType point=mesh.vert[i].P();
            newpositions.row(i)<<point.X(),point.Y(),point.Z();
        }
        positions=newpositions;
        for(size_t i=0;i<handle_positions.size();i++){
            newpositions(Handle[i],0)=handle_positions[i].X();
            newpositions(Handle[i],1)=handle_positions[i].Y();
            newpositions(Handle[i],2)=handle_positions[i].Z();
        }
        //<<"new positions"<<endl;
        //cout<<newpositions<<endl;
        // computing Ti(V')δi
        computeRoiDelta();
        //cout<<"roidelta"<<endl;
        //cout<<roiDelta<<endl;

        if(correctLaplacian==CorrectLaplacian)
            computeTi();
        //cout<<"computing Ti(V')δi"<<endl;
        size_t ROI_size=mesh.vert.size();
        size_t handle_size=Handle.size();
        DMat b(ROI_size+handle_size,3);
        for (int i = 0; i <ROI_size; ++i) {
            ScalarType dx=roiDelta(i,0),dy=roiDelta(i,1),dz=roiDelta(i,2);
            if(correctLaplacian==CorrectLaplacian){
                // set of {i} and the neigbbours of i.
                std::vector<size_t> iAndNeighbours=Neighbours[i];

                Vec pi(iAndNeighbours.size()*3);
                for (int j = 0; j < iAndNeighbours.size(); ++j) {
                    int k = iAndNeighbours[j];
                    pi[ j * 3 + 0] =newpositions(k,0);
                    pi[ j * 3 + 1] =newpositions(k,1);
                    pi[ j * 3 + 2] =newpositions(k,2);
                }
                Vec sht(7);
                sht=Ts[i]*pi;
                //sht[0]=1.0; // in order not to consider scaling

                DMat D(3,7);

                D<<dx,0,dz,-dy,1.0,0.0,0.0,
                   dy,-dz,0.0,dx,0.0,1.0,0.0,
                   dz,dy,-dx,0.0,0.0,0.0,1.0;
                Vec mul=D*sht;
                Vec aux(3);
                aux<<dx,dy,dz;
                ScalarType scale=aux.norm()/mul.norm();
                /* cout<<"vertex "<<i<<endl;
                 * cout<<" scale "<<1/scale<<" "<<sht[0]<<endl;
                */
                mul*=scale;
                b.row(i)=mul;
            }
            if(correctLaplacian==OnlyLaplacian){
                b.row(i)<<dx,dy,dz;
            }
            if(correctLaplacian==ZeroLaplacian){
                b.row(i)<<0,0,0;
            }
        }
        // computing augmented  matrix and vector
        std::vector<std::pair<int,int>> IndexM=IndexL;
        std::vector<ScalarType> ValuesM=ValuesL;
        for(size_t i=0;i<Handle.size();i++){
            size_t idv=Handle[i];
            IndexM.push_back(make_pair(ROI_size+i,idv));
            ValuesM.push_back(constrain_weight);
            b(ROI_size+i,0)=constrain_weight*newpositions(idv,0);
            b(ROI_size+i,1)=constrain_weight*newpositions(idv,1);
            b(ROI_size+i,2)=constrain_weight*newpositions(idv,2);
        }
        SpMat M;
        InitSparse(IndexM,ValuesM,ROI_size+handle_size,ROI_size,M);

        // Now we solve
        // Ax = b
        // where A is the augmented matrix, and the value of b depends on whether we are optimizing with the correctLaplacian flag
        // by solving, we obtain the deformed surface coordinates.
        SpMat Mt=M.transpose();
        DMat y=Mt*b;
        DMat X;
        //cout<<y.block(0,0,10,3)<<endl<<endl;
        if(ROI_size<=20000){
            std::shared_ptr<Eigen::SimplicialCholesky<SpMat>> solver;
            solver.reset(new Eigen::SimplicialCholesky<SpMat>(Mt*M));
            //Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> solver;
            //solver.compute(M);
             //X=solver.solve(b);
            X = solver->solve(y);
        }
        else{
            cout<<"Large sparse solver"<<endl;
            //Eigen::SparseQR<SpMat,Eigen::COLAMDOrdering<int>> solver;
            //Eigen::ConjugateGradient<SpMat,Eigen::Lower> solver;
            //Eigen::SimplicialLDLT<SpMat> solver;
            Eigen::LeastSquaresConjugateGradient<SpMat> solver;
            solver.compute(M);
            if(solver.info() != Eigen::Success) {
                // decomposition failed
                cout<<"decomposition failed" <<endl;
                return;
                //X=solver->solve(y);
            }
            else{
                X=solver.solve(b);
            }
        }


        // if minimizing with correctLaplacian, a local scaling is introduced by the solver.
        // so we need to normalize the delta coordinates of the deformed vertices back to their
        // original lengths.
        // otherwise, the mesh will increase in size when manipulating the mesh, which is not desirable.

        /*DMat solutionDelta = L * X;

        for (int i = 0; i < ROI_size; ++i) {
            Eigen::Vector3d deltai,roiDeltai;
            roiDeltai<<roiDelta(i,0),roiDelta(i,1),roiDelta(i,2);
            deltai<<solutionDelta(i,0),solutionDelta(i,1),solutionDelta(i,2);
            double len =deltai.norm();
            double originalLength = roiDeltai.norm();
            double scale = originalLength / len;
            b(i,0)= scale * solutionDelta(i,0);
            b(i,1)= scale * solutionDelta(i,1);
            b(i,2)= scale * solutionDelta(i,2);
        }

        y = Mt * b;
        //solver.reset(new Eigen::SimplicialCholesky<SpMat>(Mt*M));
        DMat normalizedSolution = solver->solve(y);
        cout<<X.block(0,0,5,3)<<endl;
        cout<<endl;
        cout<<normalizedSolution.block(0,0,5,3)<<endl;
        X=normalizedSolution;*/
        for(size_t i=0;i<mesh.vert.size();i++){
            mesh.vert[i].P()=CoordType(X(i,0),X(i,1),X(i,2));
        }
        for(size_t i=0;i<basemesh->vert.size();i++){
            size_t index=map_vertices_base_to_mesh[i];
            if(index!=invalidindex)
                basemesh->vert[i].P()=CoordType(X(index,0),X(index,1),X(index,2));
        }
        UpdateBounding<MeshType>::Box(mesh);
        UpdateNormal<MeshType>::PerVertexNormalizedPerFaceNormalized(mesh);

        /*Eigen::SparseQR<SpMat,Eigen::COLAMDOrdering<int>> solversp;
        solver.reset(new Eigen::SimplicialCholesky<SpMat>(Mt*M));
        solversp.compute(M);
        DMat X;
        if(solversp.info() != Eigen::Success) {
            // decomposition failed
            X=solver->solve(y);
        }
        else{
            X=solversp.solve(b);
        }
        //DMat X = A.colPivHouseholderQr().solve(b);*/
    }
public:

    MeshType mesh;
    SpMat L; // pure ROI laplacian matrix
    MeshType* basemesh;
    vector<size_t> ROI_base;
    vector<size_t> Handle_base;
    vector<size_t> Handle;

    // array of entries for the augmented matrix
    std::vector<std::pair<int,int>> IndexL;
    std::vector<ScalarType> ValuesL;
    std::vector<DMat> Ts;
    DMat roiDelta; // this is the laplacian &_i
    vector<vector<size_t>> Neighbours;

    /* Mapping between the basemesh indexes and the ROI indexes, i.e ,
    *  basemesh.vert[i]=mesh.vert[map[i]]
    */
    vector<size_t> map_vertices_base_to_mesh;

    ScalarType constrain_weight;
};

#endif // LSD_H
