#ifndef PATCH3D_H
#define PATCH3D_H
#include <QDir>
#include "cdt2d.h"
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/mesh_assert.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include <wrap/io_trimesh/export_obj.h>
//#include "naturalparametrization.h"
#include "parameterizationlscm.h"
#include "vectorcone.h"
#include "myutils.h"
#include "PatternsTakayama/patchg.h"
#include "lsd.h"
#include "orient_faces.h"


//#include <wrap/igl/lscm_parametrization.h>

namespace  tempora {

template<typename MeshType>
struct POS_Side{
    vcg::face::Pos<typename MeshType::FaceType> pos;
    int sideIndex;
};

// Model a side of the patch
template<class MeshType>
class Side{
public:
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::VertexPointer VertexPointer;
    Side(){
        directedside.clear();
        quality.clear();
        size=0;
        marked_index.clear();
        iscurve=false;
        isclosed=false;
        ismarked=false;
        origin.SetNull();
        end.SetNull();
    }
    void setPOS(std::vector<vcg::face::Pos<FaceType>> &idside){
         quality.clear();
         directedside=idside;         
         origin=idside[0];
         end=idside[idside.size()-1];
         size=idside.size();
         marked_index.clear();
         iscurve=false;
         ismarked=false;
         vector<int> default_vec(5,-1);
         quality.resize(size,default_vec);
         if(origin.VFlip()->P()==end.V()->P())
               isclosed=true;
    }
    void setCurveFlag(bool cf){
        iscurve=cf;
    }
    void clearQuality(){
        for(int i=0;i<quality.size();i++){
            vector<int> default_vec(5,-1);
            quality[i]=default_vec;
        }
    }
    size_t getK(VertexPointer V0,VertexPointer V1, bool tryopp=true){
        size_t k=0;
        for(typename std::vector<vcg::face::Pos<FaceType>>::iterator it=directedside.begin();it!=directedside.end();it++,k++){

            if((it->V()==V1 && it->VFlip()==V0) || (tryopp && (it->V()==V0 && it->VFlip()==V1))){
                return k;
            }
        }
        return k;
    }
    size_t getK(MeshType &ref,size_t i0,size_t i1, bool tryopp=true){
        size_t k=0;
        for(typename std::vector<vcg::face::Pos<FaceType>>::iterator it=directedside.begin();it!=directedside.end();it++,k++){

            if((tri::Index(ref,it->V())==i0 && tri::Index(ref,it->VFlip())==i1) || (tryopp && (tri::Index(ref,it->V())==i1 && tri::Index(ref,it->VFlip())==i0))){
                return k;
            }
        }
        return k;
    }
    ScalarType length(){
       ScalarType result=0;
       for(int i=0;i<directedside.size();i++)
           result+=vcg::Distance(directedside[i].VFlip()->P(),directedside[i].V()->P());

       return result;
    }
    std::vector<vcg::face::Pos<FaceType>> directedside;    
    int size;
    set<size_t> marked_index;
    bool iscurve;
    bool isclosed;
    bool ismarked;

    // Properties used for procesing at higer levels of the code
    vector<vector<int>> quality;

    vcg::face::Pos<FaceType> origin,end;

};
/* Class that models a patch for our purposes and gives us useful methods in order to get a quadrangulation of it.
 * The type of patches that we want to manage is the result of a boolean operation between a 4-sided surface patch
 * and a collection of free-form surface curves (curves over the surface). The cut process is done in a triangular
 * version of the quad meshes. Then, the patch to manage is a surface piece of a triangle mesh, associated to a quad
 * surface patch. It may have several connected components.
 * We set up our representation of the patch with three meshes and an integer: A quad mesh and its default triangular version and a triangle mesh result
 * of the boolean operation. This last contains perface and pervertex attributes that allow us to go from one mesh to the other. The two triangle meshes
 * are necessary because the boolean operation is performed into triangles (in order to deal with non-planar quads). The integer is an ID that identifies all
 * faces belonging to the patch. This integer is present in all three meshes.
*/
template<class TriMeshType,class PMeshType>
class Patch3D
{
    typedef vcg::face::Pos<typename TriMeshType::FaceType> FacePos;
    typedef typename PMeshType::FaceType FaceTypeQ;
    typedef typename TriMeshType::FaceType FaceTypeT;
    typedef typename TriMeshType::ScalarType ScalarType;
    typedef typename TriMeshType::FacePointer FacePointer;
    typedef typename TriMeshType::VertexPointer VertexPointer;
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexIterator VertexIterator;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType:: template PerFaceAttributeHandle<int> AttributeFaceType;

    // Model a connected component
    struct DSCC{

        // Property for processing at higher levels of the code.
        int quality;
        vector<Side<TriMeshType>> sidesT;
        vector<size_t> faces;
        MonotonicFunctions::PrePatch2D boundaries2d;
        LSCMparameterization parameterization;
        /* Store the adjacency with other patches
         * Each pair are in the format : (component side index, adjacent patch number ) Set in a high level of the code
         */
        vector<std::pair<int,int>> adjacency;
        bool hasholes;
        DSCC(){
          sidesT.clear();
          parameterization.reset();
          faces.clear();
          hasholes=false;
          boundaries2d.setnull();
          quality=-1;
        }
        size_t getFaceIndex(TriMeshType &ref,typename TriMeshType::VertexPointer &V0,typename TriMeshType::VertexPointer &V1, bool tryopp=true){
            size_t k=0;
            for(size_t j=0;j<sidesT.size();j++){
                for(size_t i=0;i<sidesT[j].directedside.size();i++){
                    if((sidesT[j].directedside[i].V()==V1 && sidesT[j].directedside[i].VFlip()==V0) || (tryopp && (sidesT[j].directedside[i].V()==V0 && sidesT[j].directedside[i].VFlip()==V1))){
                        return tri::Index(ref,sidesT[j].directedside[i].F());
                    }
                }
            }
            return k;
        }
        size_t getFaceIndex(TriMeshType &ref,size_t i0,size_t i1, bool tryopp=true){
            size_t k=0;
            for(size_t j=0;j<sidesT.size();j++){
                for(size_t i=0;i<sidesT[j].directedside.size();i++){
                    if((tri::Index(ref,sidesT[j].directedside[i].V())==i0 && tri::Index(ref,sidesT[j].directedside[i].VFlip())==i1) ||
                       (tryopp && (tri::Index(ref,sidesT[j].directedside[i].V())==i1 && tri::Index(ref,sidesT[j].directedside[i].VFlip())==i0))){
                        return tri::Index(ref,sidesT[j].directedside[i].F());
                    }
                }
            }
            return k;
        }
    };


public:

    PMeshType *quadmesh;        //quad mesh
    TriMeshType *booleanmesh;
    TriMeshType *tmesh;         //triangulation of the quad mesh

    // Manage the parameterization of the patch over the 4-sided patch over the quad mesh

    //NaturalParametrization parameterization;

    //Face indexes belonging to the patch (in booleanmesh)
    std::vector<size_t> PatchFaces;

    // Defaul color of the patch
    vcg::Color4b color;

    // ID identifying the patch's faces in the reference's meshes
    int numberpartitionQuadMesh;

    // List of connected components
    vector<struct DSCC> components;

    int numberHoles;
    int numberCurveSides;
    int numberConnectedComponents;

    bool exportIntermediateSteps;

    string patchName;
    string sessionName;

public:
    /*
     * Defaul constructor
     */
    Patch3D(){
        quadmesh=NULL;
        tmesh=NULL;
        booleanmesh=NULL;
        PatchFaces.clear();
        components.clear();
        //color=vcg::Color4b::DarkGreen;
        color=vcg::Color4b::White;
        numberpartitionQuadMesh=-1;
        numberHoles=0;
        numberCurveSides=0;
        numberConnectedComponents=0;
        exportIntermediateSteps=false;
        patchName="DefaultPatch3D";
        sessionName="DefaultSession";
    }
    /*
     * Constructor with three meshes and an ID. The attributes are passed because convenience
     */
    Patch3D(TriMeshType &tbool,TriMeshType &trmesh,PMeshType &mesh,std::vector<size_t> &faces, int inumberPart=-1, string iname="DefaultPatch3D", string isession="DefaultSession"){
        quadmesh=&mesh;
        tmesh=&trmesh;
        booleanmesh=&tbool;
        PatchFaces=faces;
        components.clear();       
        patchName=iname;
        sessionName=isession;
        color=vcg::Color4b::White;
        numberpartitionQuadMesh=inumberPart;
        numberHoles=0;
        numberCurveSides=0;
        numberConnectedComponents=0;
        exportIntermediateSteps=false;
        computeBorders(false);
    }
    /*
     * To set the quality property
    */
    void setQuality(size_t i,int q){
        assert(i<components.size());
        components[i].quality=q;
    }
    /*
     * Reset set of faces belonging to the patch
    */
    void setfaces(std::vector<size_t> &faces){
        PatchFaces=faces;
    }
    /*
     * Reset color of the patch
    */
    void setcolor(vcg::Color4b &icolor){
        color=icolor;
    }
    /*
     * Mark all vertices of the boolean mesh patch as non-visited
    */
    void cleanVisitedVertices(){
        for(int i=0;i<PatchFaces.size();i++){
             booleanmesh->face[PatchFaces[i]].V(0)->ClearV();
             booleanmesh->face[PatchFaces[i]].V(1)->ClearV();
             booleanmesh->face[PatchFaces[i]].V(2)->ClearV();
        }
    }
    /*
     * Mark all faces of the boolean mesh patch as non-visited
    */
    void cleanVisitedFaces(){
        for(int i=0;i<PatchFaces.size();i++)
             booleanmesh->face[PatchFaces[i]].ClearV();
    }
    /*
     * Paint patch, i.e, all faces of the boolean mesh patch
    */
    void paint(){
        //vcg::tri::UpdateColor<CMesh>::PerFaceConstant(*booleanmesh,vcg::Color4b::White);
        for(int i=0;i<PatchFaces.size();i++)
            booleanmesh->face[PatchFaces[i]].C()=color;
    }
    /*
     * Setting the flag ExportIntermediateSteps. This flag is used for saving the intermediate patches and other steps of the algorithm
     */
    void setExportIntermediateSteps(bool flag_eis){
       exportIntermediateSteps=flag_eis;
    }
    /*
     * Compute the connected components of a set of faces in the patch
     * @param mesh is the triangular mesh after the boolean operation. This argument is unnecessary but because of convenience...
    */
    void connectedComponents(TriMeshType &mesh, std::vector<std::vector<size_t>> &CCV){
            AttributeFaceType hMeshNumber= vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(*booleanmesh,std::string("MeshNumber"));
            vcg::tri::RequireFFAdjacency(mesh);
            CCV.clear();
            cleanVisitedFaces();
            std::stack<FacePointer> sf;
            FacePointer fpt=&mesh.face[PatchFaces[0]];
            for(int i=0;i<PatchFaces.size();i++)
            {
              FacePointer fi=&mesh.face[PatchFaces[i]];
              if(!fi->IsV())
              {
                fi->SetV();
                std::vector<size_t> tempCCV;
                sf.push(fi);
                while (!sf.empty())
                {
                  fpt=sf.top();
                  tempCCV.push_back(vcg::tri::Index(mesh,fpt));
                  sf.pop();
                  for(int j=0; j<fpt->VN(); ++j)
                  {
                    vcg::face::Pos<FaceType> pf(fpt,j,fpt->V(j));
                    bool isborder=(pf.FFlip()->IMark()!=numberpartitionQuadMesh);
                    isborder=isborder || (hMeshNumber[pf.F()]!=hMeshNumber[pf.FFlip()]);
                    if(!isborder)
                    {
                      FacePointer l = fpt->FFp(j);
                      if( !(*l).IsV() )
                      {
                        (*l).SetV();
                        sf.push(l);
                      }
                    }
                  }
                }
                CCV.push_back(tempCCV);
              }
            }
    }
    /* Compute the total number of edges that forms the boundary of the i-th connected component of this patch, i.e,
     *  the sum of the edges of each side of this connected component
     * */
    int numberDivisions(size_t i){
        if(numberConnectedComponents==0 || components.size()<=i){
            return 0;
        }
        else{
            int result=0;
            for(size_t j=0;j<components[i].sidesT.size();j++){
                result+=components[i].sidesT[j].size;
                //cout<<result<<endl;
            }
            return result;
        }
    }

    /*
     * Compute the connected components and the data structure as such.
     *
    */
    void computeBorders(bool computeparameterization=true){
        AttributeFaceType hMeshNumber= vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(*booleanmesh,std::string("MeshNumber"));
        AttributeFaceType hCorrespondence= vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(*booleanmesh,std::string("Correspondence"));

        //Coherence
        for(int i=0;i<PatchFaces.size();i++){
             assert(booleanmesh->face[PatchFaces[i]].IMark()==numberpartitionQuadMesh);
        }

        // Obtaining connected components
        components.clear();
        cleanVisitedVertices();
        vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromFF(*booleanmesh);
        std::vector<std::vector<size_t>> CCV;
        connectedComponents(*booleanmesh,CCV);
        numberConnectedComponents=CCV.size();
        if(numberConnectedComponents>1){
            for(int i=0;i<CCV.size();i++){
                //cout<<"connected component "<<CCV[i].size()<<endl;
            }
            for(int i=0;i<CCV[1].size();i++)
                booleanmesh->face[CCV[1][i]].C()=color;
        }
        for(int i=0;i<numberConnectedComponents;i++){
            // Making a DSCC
            DSCC sidesCB;
            sidesCB.faces=CCV[i];
            // Making the sides
            // First we compute the border holes of the component
            for(int k=0;k<CCV[i].size();k++){
                  for(int j=0;j<booleanmesh->face[CCV[i][k]].VN();++j)
                  {
                      vcg::face::Pos<FaceTypeT> pf(&(booleanmesh->face[CCV[i][k]]),j,booleanmesh->face[CCV[i][k]].V(j));                    
                      bool isIntersectionCurveEdge=(hMeshNumber[pf.F()]!=hMeshNumber[pf.FFlip()]);

                      if(isIntersectionCurveEdge && !pf.V()->IsV()){
                          std::vector<vcg::face::Pos<FaceTypeT>> protoside;
                          protoside.push_back(pf);
                          vcg::face::Pos<FaceTypeT> pfbegin=pf;
                          bool isfrontier;
                          do{
                              do{
                                 pf.NextE();
                                 isfrontier=(pf.FFlip()->IMark()!=numberpartitionQuadMesh) || (hMeshNumber[pf.F()]!=hMeshNumber[pf.FFlip()]);
                              }while(!isfrontier);
                              pf.FlipF();
                              pf.FlipV();
                              pf.V()->SetV();
                              protoside.push_back(pf);
                              isIntersectionCurveEdge= (hMeshNumber[pf.F()]!=hMeshNumber[pf.FFlip()]);
                          }while(pf!=pfbegin && isIntersectionCurveEdge);

                          if(pf==pfbegin){
                              /*pf.F()->C()=vcg::Color4b::DarkBlue;
                              int triIndex=hCorrespondence[pf.F()];
                              int qIndex=(triIndex%2==0)? int(triIndex/2):int((triIndex-1)/2);
                              quadmesh->face[qIndex].C()=vcg::Color4b::DarkBlue;*/
                              pf.V()->SetV();
                              protoside.pop_back();
                              Side<TriMeshType> tempSide;
                              tempSide.setPOS(protoside);
                              tempSide.setCurveFlag(true);
                              sidesCB.sidesT.push_back(tempSide);
                              sidesCB.hasholes=true;
                              numberHoles++;
                              numberCurveSides++;
                              cout<<"A curved loop found"<<endl;
                          }
                      }
                  }
            }
            // Second, we compute the border curve compose of quad patch sides and intersection curve sides
               cleanVisitedVertices();
               for(int k=0;k<sidesCB.sidesT.size();k++){
                   for(int r=0;r<sidesCB.sidesT[k].size;r++)
                       sidesCB.sidesT[k].directedside[r].V()->SetV();
               }
               // we start  finding a border half-edge that is not in a hole.
               vcg::face::Pos<FaceTypeT> pbegin;
               pbegin.SetNull();
               for(int k=0;k<CCV[i].size();k++){
                   if(!booleanmesh->face[CCV[i][k]].IsD()){
                       for(int j=0;j<booleanmesh->face[CCV[i][k]].VN();++j)
                       {
                           vcg::face::Pos<FaceTypeT> pf(&(booleanmesh->face[CCV[i][k]]),j,booleanmesh->face[CCV[i][k]].V(j));
                           bool isborder=(pf.FFlip()->IMark()!=numberpartitionQuadMesh);
                           //bool isNonIntersectionCurve=!pf.F()->IsS();
                           if(isborder && !pf.V()->IsV()){
                               pbegin=pf;
                               break;
                           }
                       }
                       if(!pbegin.IsNull())
                           break;
                   }
               }

            if(pbegin.IsNull()){ //Degenerate case. Just one border curve that is a intersection curve.
               sidesCB.hasholes=false;               
               numberHoles--;

               //choosing a good pair of extreme points in order to parameterize
               size_t zerozero,oneone;
               assert(sidesCB.sidesT[0].size>2);
               int middle= int(sidesCB.sidesT[0].size/2);

               zerozero=tri::Index(*booleanmesh,sidesCB.sidesT[0].origin.VFlip());
               oneone=tri::Index(*booleanmesh,sidesCB.sidesT[0].directedside[middle].V());
               sidesCB.parameterization.initialize(*booleanmesh,CCV[i],zerozero,oneone);
               components.push_back(sidesCB);
               CMesh projected;
               sidesCB.parameterization.exportMeshParametrization(projected);
               if(exportIntermediateSteps){
                   string stringcomponent=checkPath("");
                   stringcomponent+="/component"+std::to_string(i)+".ply";
                   utility::saveMesh(projected,stringcomponent);
               }
               return;
            }
            // We walk over all neighboring triangles after we find a border half-edge in the direction of a 4-side patch into the quad mesh,
            // in other words, if the half-edge is over a non intersection curve edge.

            bool isNonIntersectionCurveEdge= (hMeshNumber[pbegin.F()]==hMeshNumber[pbegin.FFlip()]);
            int count=0;
            while(!isNonIntersectionCurveEdge && count<1000){
                bool isfrontier=false;
                do{
                   pbegin.NextE();
                   isfrontier=(pbegin.F()->IMark()!=numberpartitionQuadMesh);
                }while(!isfrontier);
                pbegin.FlipF();
                pbegin.FlipV();
                //pbegin.F()->C()=vcg::Color4b::DarkRed;
                isNonIntersectionCurveEdge=  (hMeshNumber[pbegin.F()]==hMeshNumber[pbegin.FFlip()]);
                count++;
            }

            //Third, we go back to the begining of the quad patch side or to the first cut with a intersection curve.

             vcg::face::Pos<FaceTypeT> pbefore=pbegin;
             //POS_Side<PMeshType> quadbegin=getQuadBorderFromTriBorder(pbegin);
             isNonIntersectionCurveEdge=true;
             bool hasSamePatchNumber=true;
             do{
                  pbefore=pbegin;
                  //pbegin.FlipV();
                  bool isfrontier=true;
                  do{
                     pbegin.NextE();
                     isfrontier=(pbegin.F()->IMark()!=numberpartitionQuadMesh) || (hMeshNumber[pbegin.F()]!=hMeshNumber[pbegin.FFlip()]);
                  }while(!isfrontier);
                  pbegin.FlipF();
                  pbegin.FlipV();
                  isNonIntersectionCurveEdge=(hMeshNumber[pbegin.F()]==hMeshNumber[pbegin.FFlip()]);
                  if(!isNonIntersectionCurveEdge)
                      hasSamePatchNumber=true;
                  else{
                      hasSamePatchNumber=(pbefore.FFlip()->IMark()==pbegin.FFlip()->IMark());
                  }
             }while(isNonIntersectionCurveEdge && hasSamePatchNumber);

             // Now, we store in sequence the sides of this border curve

             //pbegin=pbefore;
             //pbegin.F()->C()=color;
             vcg::face::Pos<FaceTypeT> pstart=pbegin;
             isNonIntersectionCurveEdge= (hMeshNumber[pbegin.F()]==hMeshNumber[pbegin.FFlip()]);
             std::vector<vcg::face::Pos<FaceTypeT>> posside;
             do{
                 posside.push_back(pbegin);
                 pbefore=pbegin;
                 bool isfrontier=false;
                 do{
                     pbegin.NextE();
                     isfrontier=(pbegin.F()->IMark()!=numberpartitionQuadMesh) || (hMeshNumber[pbegin.F()]!=hMeshNumber[pbegin.FFlip()]);
                 }while(!isfrontier);
                 pbegin.FlipF();
                 pbegin.FlipV();
                 isNonIntersectionCurveEdge= (hMeshNumber[pbegin.F()]==hMeshNumber[pbegin.FFlip()]);
                 bool isNonIntersectionCurveEdgeBefore= (hMeshNumber[pbefore.F()]==hMeshNumber[pbefore.FFlip()]);
                 bool hasSamePatchNumber;
                 if(isNonIntersectionCurveEdge && isNonIntersectionCurveEdgeBefore){
                     hasSamePatchNumber=(pbefore.FFlip()->IMark()==pbegin.FFlip()->IMark());
                 }
                 else
                     hasSamePatchNumber=true;

                 //finishing one side originated from a quad side
                 if(isNonIntersectionCurveEdge && !hasSamePatchNumber && isNonIntersectionCurveEdgeBefore){
                     Side<TriMeshType> tempSide;
                     tempSide.setPOS(posside);
                     /*for(auto ip=posside.begin();ip!=posside.end();ip++)
                         ip->F()->C()=Color4b::DarkBlue;*/
                     posside.clear();
                     sidesCB.sidesT.push_back(tempSide);
                 }
                 //finishing one side as a cutting quad side
                 if(!isNonIntersectionCurveEdge && isNonIntersectionCurveEdgeBefore){
                     Side<TriMeshType> tempSide;
                     tempSide.setPOS(posside);
                     /*for(auto ip=posside.begin();ip!=posside.end();ip++)
                         ip->F()->C()=Color4b::DarkGray;*/
                     posside.clear();
                     sidesCB.sidesT.push_back(tempSide);
                 }
                 //finishing one side as a part of some intersection curve
                 if(isNonIntersectionCurveEdge && !isNonIntersectionCurveEdgeBefore)
                 {
                     Side<TriMeshType> tempSide;
                     tempSide.setPOS(posside);
                     tempSide.setCurveFlag(true);
                     numberCurveSides++;
                     /*for(auto ip=posside.begin();ip!=posside.end();ip++)
                         ip->F()->C()=Color4b::DarkGreen;*/
                     posside.clear();
                     sidesCB.sidesT.push_back(tempSide);
                 }
             }while(pbegin!=pstart);
             //choosing a good pair of extreme points in order to parameterize
             size_t numbersides=sidesCB.sidesT.size();             
             int middle=int(numbersides/2.0);
             size_t zerozero,oneone;
             if(numbersides>2){
                 zerozero=tri::Index(*booleanmesh,sidesCB.sidesT[0].origin.VFlip());
                 oneone=tri::Index(*booleanmesh,sidesCB.sidesT[middle].end.V());
             }
             else{ // sidesCB has 2 sides ,then the code above gives only one index
                 zerozero=tri::Index(*booleanmesh,sidesCB.sidesT[0].origin.VFlip());
                 oneone=tri::Index(*booleanmesh,sidesCB.sidesT[middle].origin.VFlip());
             }
             sidesCB.parameterization.initialize(*booleanmesh,CCV[i],zerozero,oneone);

             CMesh projected;
             sidesCB.parameterization.exportMeshParametrization(projected);
             if(exportIntermediateSteps){
                 string stringcomponent=checkPath("");
                 stringcomponent+="/component"+std::to_string(i)+".ply";
                 utility::saveMesh(projected,stringcomponent);
             }
             components.push_back(sidesCB);
        }



        /*if(computeparameterization)
            setBoundaries2D();*/
    }    
    /*
     * Get a list of 2D boundaries formed by all boundaries of all connected components projected to the unit square.
     * The algorithm takes into account that certain edges from certain connected components will be splitted and some sides will be complete resampling
     * in order to satisfy all requirements of parity and good quadrangulation
     *
     * @param lisfFrontiers is the collection of single boundary 2D polygons representing the projection of all boundaries of all
     * connected components of the patch
     * @feature: list of features indexes corresponding to each single polygon in listFronties
     * great_indexes: for each curve in the collection of intersection curves this argument lists fixed indexes
     * blocks: for each curve in the collection of intersection curves this argument  lists a complete resampling of the curve in blocks
     * blocks_interpolation: for each curve in the collection of intersection curves this argument  lists a complete resampling of the curve in blocks as barycentric coordinates
     * indexes_to_blocks: for each curve in the collection of intersection curves this maps great_indexes to blocks_interpolation values
     */
    /*template<class CurveClass>
    void get2DFrontiersLowLevel(vector<vector<Point3f>> & listFrontiers,vector<set<int>> &feature,
                                CurveClass &intersection_curve){

        typedef std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,size_t> interpolation_data;
        typedef pair<FacePos,interpolation_data> frontier_point;
        if(numberpartitionQuadMesh==-1) return;
        listFrontiers.clear();
        feature.clear();
        assert(numberConnectedComponents==components.size());
        for(auto cc=components.begin();cc!=components.end();cc++){
            set<int> featureVertices;
            vector<frontier_point> boundary;
            vector<Point3f> temp3DFrontier;
            // If this component has holes then we store these in a array and join properly with the exterior boundary in a later step
            if(cc->hasholes){
                vector<vector<Point2f>> holes2D;
                vector<Point2f> eboundary2D;
                // Here we store the projected holes
                for(int j=0;j<cc->sidesT.size();j++){
                    if(cc->sidesT[j].isclosed && cc->sidesT[j].iscurve){
                        vector<Point2f> tempholes;
                        vector<interpolation_data> tempholes_barycentric;
                        //utility::saveAsCurvePly(cc->sidesT[j].directedside,"parameterization/tempHoleOriginal.ply");
                        intersection_curve.getSmoothVersion(cc->sidesT[j].origin,cc->sidesT[j].end,tempholes_barycentric);
                        tempholes_barycentric.pop_back(); // because this is a hole
                        vector<CMesh::CoordType> tosave;
                        vector<CMesh::CoordType> tosave2;
                        for(size_t i=0;i<tempholes_barycentric.size();i++){
                            CMesh::VertexPointer a=std::get<0>(tempholes_barycentric[i]);
                            CMesh::VertexPointer b=std::get<1>(tempholes_barycentric[i]);
                            size_t u=std::get<2>(tempholes_barycentric[i]);
                            size_t ia=tri::Index(*booleanmesh,a);
                            size_t ib=tri::Index(*booleanmesh,b);
                            //tosave.push_back(a->P()*u+b->P()*(1.0-u));
                            //tosave.push_back(intersection_curve.renderSmooth.vert[u]);
                            size_t indexPos=cc->sidesT[j].getK(*booleanmesh,ia,ib);
                            CMesh::CoordType position=intersection_curve.renderSmooth.vert[u].P();
                            size_t indexFace=tri::Index(*booleanmesh,cc->sidesT[j].directedside[indexPos].F());
                            Point2d projection=parameterization.getUVFromFace_Position(*booleanmesh,indexFace,position);
                            tempholes.push_back(Point2f(projection.X(),projection.Y()));
                        }
                        //utility::saveAsCurvePly(tosave,"parameterization/tempHoleInterpolation.ply");
                        //utility::saveAsCurvePly(tosave2,"parameterization/tempHoleVertex.ply");
                        if(utility::isCounterClockWise(tempholes))
                            std::reverse(tempholes.begin(),tempholes.end());

                        holes2D.push_back(tempholes);
                    }
                    else{
                       featureVertices.insert(vcg::tri::Index(*booleanmesh,cc->sidesT[j].end.V()));
                    }
                }
                // Here we store the exterior boundary as a collection of enhaced edges (pos + interpolation data)
                // This because in intersection curve edges we complete resampling and in some other edges
                // we split them in order to solve the parity of the network of patches.
                // After considering the correct split and interpolation information we project and create the 2d
                // exterior boundary.
                for(int j=0;j<cc->sidesT.size();j++){ // this works because we first store all holes and after the other sides
                    if(!(cc->sidesT[j].isclosed && cc->sidesT[j].iscurve)){
                         for(size_t i=0;i<cc->sidesT[j].size;i++){
                            FacePos tpos=cc->sidesT[j].directedside[i];
                            std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,size_t> tinter(0,0,0);
                            boundary.push_back(make_pair(tpos,tinter));
                         }
                    }
                    if(cc->sidesT[j].iscurve){
                        vector<interpolation_data> tempholes_barycentric;
                        intersection_curve.getSmoothVersion(cc->sidesT[j].origin,cc->sidesT[j].end,tempholes_barycentric);
                        tempholes_barycentric.pop_back(); // because in the circulation we consider the VFlip pos position
                        for(size_t i=0;i<tempholes_barycentric.size();i++){
                            FacePos tpos;
                            tpos.SetNull();
                            boundary.push_back(make_pair(tpos,tempholes_barycentric[i]));
                        }
                    }
                }

                set<int> featureThis;
                for(int j=0;j<boundary.size();j++){
                    frontier_point temp=boundary[j];
                    if(!temp.first.IsNull()){
                        FacePos tempos=temp.first;
                        tempos.FlipV();  // in order to begin on a corner vertex
                        if(featureVertices.find(tri::Index(*booleanmesh,tempos.V()))!=featureVertices.end())
                            featureThis.insert(j);
                        Point2d tempoint=parameterization.addSegmentPOS(*booleanmesh,tempos);
                        eboundary2D.push_back(Point2f(tempoint.X(),tempoint.Y()));
                    }
                    else{
                        CMesh::VertexPointer a=std::get<0>(temp.second);
                        CMesh::VertexPointer b=std::get<1>(temp.second);
                        size_t u=std::get<2>(temp.second);
                        size_t ia=tri::Index(*booleanmesh,a);
                        size_t ib=tri::Index(*booleanmesh,b);
                        size_t indexFace=cc->getFaceIndex(*booleanmesh,ia,ib);
                        CMesh::CoordType position=intersection_curve.renderSmooth.vert[u].P();
                        Point2d tempoint=parameterization.getUVFromFace_Position(*booleanmesh,indexFace,position);
                        eboundary2D.push_back(Point2f(tempoint.X(),tempoint.Y()));
                    }
                }

                if(!utility::isCounterClockWise(eboundary2D)){
                    std::reverse(eboundary2D.begin(),eboundary2D.end());
                    int sizeThis=eboundary2D.size();
                    set<int> newfeature;
                    for(set<int>::iterator its=featureThis.begin();its!=featureThis.end();its++)
                        newfeature.insert(sizeThis-(*its)-1);
                    utility::toSingleBoundary(eboundary2D,holes2D,temp3DFrontier,newfeature);
                    featureVertices=newfeature;
                }
                else{
                    utility::toSingleBoundary(eboundary2D,holes2D,temp3DFrontier,featureThis);
                    featureVertices=featureThis;
                }

            }
            else{ // If it has no holes then we store all sides in sequence
                vector<FacePos> edgestosplit;
                vector<int> timestosub;
                for(int j=0;j<cc->sidesT.size();j++){
                    //storing edges of this side
                    // if this side is not part of the intersection curve, so we store the pos edges
                    if(!cc->sidesT[j].iscurve){
                        for(size_t i=0;i<cc->sidesT[j].size;i++){
                            FacePos tpos=cc->sidesT[j].directedside[i];
                            std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,size_t> tinter(0,0,0);
                            boundary.push_back(make_pair(tpos,tinter));
                        }
                        // storing sides to be splitted
                        set<size_t> & marked=cc->sidesT[j].marked_index;
                        for(const size_t &k:marked){
                            edgestosplit.push_back(cc->sidesT[j].directedside[k]);
                            timestosub.push_back(cc->sidesT[j].quality[k][4]);
                        }
                    }
                    else{ // if not so we store the interpolation data
                        vector<interpolation_data> temp_barycentric;
                        intersection_curve.getSmoothVersion(cc->sidesT[j].origin,cc->sidesT[j].end,temp_barycentric);
                        temp_barycentric.pop_back(); // because in the circulation we consider the VFlip pos position
                        for(size_t i=0;i<temp_barycentric.size();i++){
                            FacePos tpos;
                            tpos.SetNull();
                            boundary.push_back(make_pair(tpos,temp_barycentric[i]));
                        }
                    }
                    featureVertices.insert(vcg::tri::Index(*booleanmesh,cc->sidesT[j].end.V()));

                }
                if(!boundary.empty()){
                    set<int> featurethis=featureVertices;
                    vector<size_t> poslist;
                    featureVertices.clear();
                    for(int j=0;j<boundary.size();j++){
                        FacePos tempos=boundary[j].first;
                        std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double> tinter=boundary[j].second;
                        if(!tempos.IsNull()){
                            for(size_t k=0;k<edgestosplit.size();k++) {
                                if(edgestosplit[k]==tempos){
                                    poslist.push_back(j);  // to after add the points into the marked edges
                                    break;
                                }
                            }
                            //tempos.FlipV(); // in order to begin on a corner vertex

                            if(featurethis.find(tri::Index(*booleanmesh,tempos.V()))!=featurethis.end())
                                featureVertices.insert(j);

                            Point2d tempoint=parameterization.addSegmentPOS(*booleanmesh,tempos); //2D projection
                            Point2f temp2DFrontier=Point2f(tempoint.X(),tempoint.Y());
                            temp3DFrontier.push_back(Point3f(temp2DFrontier.X(),temp2DFrontier.Y(),0.0));
                        }
                        else{
                            CMesh::VertexPointer a=std::get<0>(tinter);
                            CMesh::VertexPointer b=std::get<1>(tinter);
                            size_t u=std::get<2>(tinter);
                            size_t ia=tri::Index(*booleanmesh,a);
                            size_t ib=tri::Index(*booleanmesh,b);
                            size_t indexFace=cc->getFaceIndex(*booleanmesh,ia,ib);
                            CMesh::CoordType position=intersection_curve.renderSmooth.vert[u].P();

                            if(indexFace==0)
                                return;
                            Point2d tempoint=parameterization.getUVFromFace_Position(*booleanmesh,indexFace,position);
                            temp3DFrontier.push_back(Point3f(tempoint.X(),tempoint.Y(),0.0));
                        }
                    }

                    vector<Point3f> tempF;
                    utility::splitPosList_onPolygon(temp3DFrontier,poslist,timestosub,featureVertices,tempF);
                    if(!tempF.empty())
                        temp3DFrontier=tempF;
                }
            }
            if(!temp3DFrontier.empty()){
                if(!utility::isCounterClockWise(temp3DFrontier)){
                    std::reverse(temp3DFrontier.begin(),temp3DFrontier.end());
                    set<int> newfeaturevertices;
                    int sizetemp=temp3DFrontier.size();
                    for(const int &fid : featureVertices)
                        newfeaturevertices.insert(sizetemp-1-fid);
                    featureVertices=newfeaturevertices;
                }
                listFrontiers.push_back(temp3DFrontier);
                feature.push_back(featureVertices);
            }
        }
        // to create the directory where we store the polygons
        string directory=checkPath("PolygonsToSplit");
        for(size_t i=0;i<listFrontiers.size();i++){
            string path=directory+"/poly0"+std::to_string(i)+".ply";
            utility::saveAsCurvePly(listFrontiers[i],path.c_str());
        }
    }*/
    /*
     * Check if there exist the corresponding string passed as an argument into the parent directory ../exported/patchName
     * If not exists the directory then this function creates it and return the total path
     * If the passed string is empty then it returns the parent directory ../exported/patchName
     */
    string checkPath(string idirectory){
        string directory;
        if(idirectory.empty())
            directory="../exported/"+sessionName+"/"+patchName;
        else
            directory="../exported/"+sessionName+"/"+patchName+"/"+idirectory;
        QDir dir(directory.c_str());
        if (!dir.exists()) {
            dir.mkpath(".");
        }
        return directory;
    }
    void get2DFrontiers3(vector<vector<Point3f>> & listFrontiers,vector<set<int>> &feature){
        if(numberpartitionQuadMesh==-1) return;
        listFrontiers.clear();
        assert(numberConnectedComponents==components.size());
        for(int i=0;i<components.size();i++){
            vector<Point3f> temp3DFrontier;
            set<int> featureVertices;
            // If this component has holes then we store these in a array and join properly with the exterior boundary
            if(components[i].hasholes){
                vector<vector<Point2<ScalarType>>> holes2D;
                vector<FacePos> eboundary;
                vector<Point2<ScalarType>> eboundary2D;
                for(int j=0;j<components[i].sidesT.size();j++){
                    if(components[i].sidesT[j].isclosed && components[i].sidesT[j].iscurve){
                        vector<Point2<ScalarType>> tempholes;
                        set<size_t> & marked=components[i].sidesT[j].marked_index;

                        for(int k=0;k<components[i].sidesT[j].size;k++){
                            Point2<ScalarType> tempoint=components[i].parameterization.project(components[i].sidesT[j].directedside[k]);
                            if(!marked.empty() && marked.find(k)!=marked.end()){
                                int timestosubdivide=components[i].sidesT[j].quality[k][4]+1;
                                FacePos temp=components[i].sidesT[j].directedside[k];
                                temp.FlipV();
                                Point2<ScalarType> last=components[i].parameterization.project(temp);
                                Point3<ScalarType> origin(last.X(),last.Y(),0);
                                Point3<ScalarType> endp(tempoint.X(),tempoint.Y(),0);
                                vector<Point3<ScalarType>> newpoints;
                                utility::splitSegment<ScalarType>(origin,endp,newpoints,timestosubdivide);
                                for (size_t r=0;r<newpoints.size();r++) {
                                    tempholes.push_back(Point2<ScalarType>(newpoints[i].X(),newpoints[i].Y()));
                                }
                            }
                            tempholes.push_back(tempoint);
                        }

                        if(utility::isCounterClockWise(tempholes))
                            std::reverse(tempholes.begin(),tempholes.end());

                        holes2D.push_back(tempholes);
                    }
                    else{
                       featureVertices.insert(vcg::tri::Index(*booleanmesh,components[i].sidesT[j].end.V()));
                    }
                }
                set<int> featureThis;
                for(int j=0;j<components[i].sidesT.size();j++){ // this works because we first store all holes and after the other sides
                    if(!(components[i].sidesT[j].isclosed && components[i].sidesT[j].iscurve)){

                        set<size_t> & marked=components[i].sidesT[j].marked_index;
                        for(size_t r=0;r<components[i].sidesT[j].size;r++){
                           FacePos temp=components[i].sidesT[j].directedside[r];
                           temp.FlipV();
                           Point2<ScalarType> tempoint=components[i].parameterization.project(temp);
                           eboundary2D.push_back(tempoint);
                           if(featureVertices.find(tri::Index(*booleanmesh,temp.V()))!=featureVertices.end())
                               featureThis.insert(eboundary2D.size()-1);
                           if(marked.find(r)!=marked.end()){
                               int timestosubdivide=components[i].sidesT[j].quality[r][4]+1;
                               Point3<ScalarType> origin(tempoint.X(),tempoint.Y(),0);
                               temp.FlipV();
                               Point2<ScalarType> endp2=components[i].parameterization.project(temp);
                               Point3<ScalarType> endp(endp2.X(),endp2.Y(),0);
                               vector<Point3<ScalarType>> newpoints;
                               utility::splitSegment<ScalarType>(origin,endp,newpoints,timestosubdivide);
                               for (size_t s=0;s<newpoints.size();s++) {
                                   eboundary2D.push_back(Point2<ScalarType>(newpoints[i].X(),newpoints[i].Y()));
                               }
                           }
                        }
                    }
                }


                if(!utility::isCounterClockWise(eboundary2D)){
                    std::reverse(eboundary2D.begin(),eboundary2D.end());
                    int sizeThis=eboundary2D.size();
                    set<int> newfeature;
                    for(set<int>::iterator its=featureThis.begin();its!=featureThis.end();its++)
                        newfeature.insert(sizeThis-(*its)-1);
                    utility::toSingleBoundary(eboundary2D,holes2D,temp3DFrontier,newfeature);
                    featureVertices=newfeature;
                }
                else{
                    utility::toSingleBoundary(eboundary2D,holes2D,temp3DFrontier,featureThis);
                    featureVertices=featureThis;
                }

            }
            else{ // If it has no holes then we store all sides in sequence

                vector<FacePos> tempFrontier;
                vector<FacePos> edgestosplit;
                vector<int> timestosub;
                for(int j=0;j<components[i].sidesT.size();j++){
                    tempFrontier.insert(tempFrontier.end(),components[i].sidesT[j].directedside.begin(),components[i].sidesT[j].directedside.end());
                    featureVertices.insert(vcg::tri::Index(*booleanmesh,components[i].sidesT[j].end.V()));
                    set<size_t> & marked=components[i].sidesT[j].marked_index;
                    for(const size_t &k:marked){
                        edgestosplit.push_back(components[i].sidesT[j].directedside[k]);
                        timestosub.push_back(components[i].sidesT[j].quality[k][4]);
                    }
                    if(components[i].sidesT[j].iscurve){

                    }
                }
                if(!tempFrontier.empty()){
                    set<int> featurethis=featureVertices;
                    vector<size_t> poslist;
                    featureVertices.clear();
                    for(int j=0;j<tempFrontier.size();j++){
                        FacePos temp=tempFrontier[j];

                        for (size_t k=0;k<edgestosplit.size();k++) {
                            if(edgestosplit[k]==temp){
                                poslist.push_back(j);  // to after add the points into the marked edges
                                break;
                            }
                        }
                        temp.FlipV(); // in order to begin on a corner vertex

                        if(featurethis.find(tri::Index(*booleanmesh,temp.V()))!=featurethis.end())
                            featureVertices.insert(j);

                        Point2<ScalarType> tempoint=components[i].parameterization.project(temp); //2D projection
                        Point2f temp2DFrontier=Point2f(tempoint.X(),tempoint.Y());
                        temp3DFrontier.push_back(Point3f(temp2DFrontier.X(),temp2DFrontier.Y(),0.0));
                    }
                    vector<Point3f> tempF;
                    utility::splitPosList_onPolygon(temp3DFrontier,poslist,timestosub,featureVertices,tempF);
                    if(!tempF.empty())
                        temp3DFrontier=tempF;
                }
            }

            if(!temp3DFrontier.empty()){
                if(!utility::isCounterClockWise(temp3DFrontier)){
                    std::reverse(temp3DFrontier.begin(),temp3DFrontier.end());
                    set<int> newfeaturevertices;
                    int sizetemp=temp3DFrontier.size();
                    for(const int &fid : featureVertices)
                        newfeaturevertices.insert(sizetemp-1-fid);
                    featureVertices=newfeaturevertices;
                }
                assert(temp3DFrontier.size()%2==0);
                listFrontiers.push_back(temp3DFrontier);
                feature.push_back(featureVertices);
            }

        }
        // to create the directory where we store the polygons
        if(exportIntermediateSteps){
            string directory=checkPath("PolygonsToSplit");
            for(size_t i=0;i<listFrontiers.size();i++){
                string path=directory+"/poly0"+std::to_string(i)+".ply";
                utility::saveAsCurvePly(listFrontiers[i],path.c_str());
            }
        }
    }
    /*
     * Get a list of 2D boundaries formed by all boundaries of all connected components projected to the unit square
     * @param lisfFrontiers is the collection of single boundary 2D polygons representing the projection of all boundaries of all
     * connected components of the patch
     */
    /*void get2DFrontiers2(vector<vector<Point3f>> & listFrontiers,vector<set<int>> &feature){
        if(numberpartitionQuadMesh==-1) return;
        listFrontiers.clear();
        assert(numberConnectedComponents==components.size());
        for(int i=0;i<components.size();i++){
            vector<Point3f> temp3DFrontier;
            set<int> featureVertices;
            // If this component has holes then we store these in a array and join properly with the exterior boundary
            if(components[i].hasholes){
                vector<vector<Point2f>> holes2D;
                vector<FacePos> eboundary;
                vector<Point2f> eboundary2D;
                for(int j=0;j<components[i].sidesT.size();j++){
                    if(components[i].sidesT[j].isclosed && components[i].sidesT[j].iscurve){
                        vector<Point2f> tempholes;
                        for(int k=0;k<components[i].sidesT[j].size;k++){
                            Point2d tempoint=parameterization.addSegmentPOS(*booleanmesh,components[i].sidesT[j].directedside[k]);
                            tempholes.push_back(Point2f(tempoint.X(),tempoint.Y()));
                        }

                        if(utility::isCounterClockWise(tempholes))
                            std::reverse(tempholes.begin(),tempholes.end());

                        holes2D.push_back(tempholes);
                    }
                    else{
                       featureVertices.insert(vcg::tri::Index(*booleanmesh,components[i].sidesT[j].end.V()));
                    }
                }
                for(int j=0;j<components[i].sidesT.size();j++){ // this works because we first store all holes and after the other sides
                    if(!(components[i].sidesT[j].isclosed && components[i].sidesT[j].iscurve)){
                            eboundary.insert(eboundary.end(),components[i].sidesT[j].directedside.begin(),components[i].sidesT[j].directedside.end());
                    }
                }
                assert(eboundary[0].VFlip()==eboundary.back().V());
                set<int> featureThis;
                for(int j=0;j<eboundary.size();j++){
                    FacePos temp=eboundary[j];
                    temp.FlipV(); // in order to begin on a corner vertex
                    if(featureVertices.find(tri::Index(*booleanmesh,temp.V()))!=featureVertices.end())
                        featureThis.insert(j);
                    Point2d tempoint=parameterization.addSegmentPOS(*booleanmesh,temp);
                    eboundary2D.push_back(Point2f(tempoint.X(),tempoint.Y()));
                }

                if(!utility::isCounterClockWise(eboundary2D)){
                    std::reverse(eboundary2D.begin(),eboundary2D.end());
                    int sizeThis=eboundary2D.size();
                    set<int> newfeature;
                    for(set<int>::iterator its=featureThis.begin();its!=featureThis.end();its++)
                        newfeature.insert(sizeThis-(*its)-1);
                    utility::toSingleBoundary(eboundary2D,holes2D,temp3DFrontier,newfeature);
                    featureVertices=newfeature;
                }
                else{
                    utility::toSingleBoundary(eboundary2D,holes2D,temp3DFrontier,featureThis);
                    featureVertices=featureThis;
                }

            }
            else{ // If it has no holes then we store all sides in sequence

                vector<FacePos> tempFrontier;
                vector<FacePos> edgestosplit;
                vector<int> timestosub;
                for(int j=0;j<components[i].sidesT.size();j++){
                    tempFrontier.insert(tempFrontier.end(),components[i].sidesT[j].directedside.begin(),components[i].sidesT[j].directedside.end());
                    featureVertices.insert(vcg::tri::Index(*booleanmesh,components[i].sidesT[j].end.V()));
                    set<size_t> & marked=components[i].sidesT[j].marked_index;
                    for(const size_t &k:marked){
                        edgestosplit.push_back(components[i].sidesT[j].directedside[k]);
                        timestosub.push_back(components[i].sidesT[j].quality[k][4]);
                    }
                    if(components[i].sidesT[j].iscurve){

                    }
                }
                if(!tempFrontier.empty()){                   
                    set<int> featurethis=featureVertices;                    
                    vector<size_t> poslist;
                    featureVertices.clear();
                    for(int j=0;j<tempFrontier.size();j++){
                        FacePos temp=tempFrontier[j];

                        for (size_t k=0;k<edgestosplit.size();k++) {
                            if(edgestosplit[k]==temp){
                                poslist.push_back(j);  // to after add the points into the marked edges
                                break;
                            }
                        }
                        temp.FlipV(); // in order to begin on a corner vertex

                        if(featurethis.find(tri::Index(*booleanmesh,temp.V()))!=featurethis.end())
                            featureVertices.insert(j);

                        Point2d tempoint=parameterization.addSegmentPOS(*booleanmesh,temp); //2D projection
                        Point2f temp2DFrontier=Point2f(tempoint.X(),tempoint.Y());
                        temp3DFrontier.push_back(Point3f(temp2DFrontier.X(),temp2DFrontier.Y(),0.0));
                    }
                    vector<Point3f> tempF;
                    utility::splitPosList_onPolygon(temp3DFrontier,poslist,timestosub,featureVertices,tempF);
                    if(!tempF.empty())
                        temp3DFrontier=tempF;
                }
            }

            if(!temp3DFrontier.empty()){
                if(!utility::isCounterClockWise(temp3DFrontier)){
                    std::reverse(temp3DFrontier.begin(),temp3DFrontier.end());
                    set<int> newfeaturevertices;
                    int sizetemp=temp3DFrontier.size();
                    for(const int &fid : featureVertices)
                        newfeaturevertices.insert(sizetemp-1-fid);
                    featureVertices=newfeaturevertices;
                }
                listFrontiers.push_back(temp3DFrontier);
                feature.push_back(featureVertices);
            }

        }
        // to create the directory where we store the polygons
        string directory=checkPath("PolygonsToSplit");
        for(size_t i=0;i<listFrontiers.size();i++){
            string path=directory+"/poly0"+std::to_string(i)+".ply";
            utility::saveAsCurvePly(listFrontiers[i],path.c_str());
        }
    }
    void get2DFrontiers(vector<vector<Point3f>> & listFrontiers,vector<set<int>> &feature){
        if(numberpartitionQuadMesh==-1) return;
        listFrontiers.clear();
        assert(numberConnectedComponents==components.size());
        for(int i=0;i<components.size();i++){
            vector<Point3f> temp3DFrontier;
            set<int> featureVertices;
            // If this component has holes then we store these in a array and join properly with the exterior boundary
            if(components[i].hasholes){
                vector<vector<Point2f>> holes2D;
                vector<FacePos> eboundary;
                vector<Point2f> eboundary2D;
                for(int j=0;j<components[i].sidesT.size();j++){
                    if(components[i].sidesT[j].isclosed && components[i].sidesT[j].iscurve){
                        vector<Point2f> tempholes;
                        for(int k=0;k<components[i].sidesT[j].size;k++){
                            Point2d tempoint=parameterization.addSegmentPOS(*booleanmesh,components[i].sidesT[j].directedside[k]);
                            tempholes.push_back(Point2f(tempoint.X(),tempoint.Y()));
                        }

                        if(utility::isCounterClockWise(tempholes))
                            std::reverse(tempholes.begin(),tempholes.end());

                        holes2D.push_back(tempholes);
                    }
                    else{
                       featureVertices.insert(vcg::tri::Index(*booleanmesh,components[i].sidesT[j].end.V()));
                    }
                }
                for(int j=0;j<components[i].sidesT.size();j++){ // this works because we first store all holes and after the other sides
                    if(!(components[i].sidesT[j].isclosed && components[i].sidesT[j].iscurve)){
                            eboundary.insert(eboundary.end(),components[i].sidesT[j].directedside.begin(),components[i].sidesT[j].directedside.end());
                    }
                }
                assert(eboundary[0].VFlip()==eboundary.back().V());
                set<int> featureThis;
                for(int j=0;j<eboundary.size();j++){
                    FacePos temp=eboundary[j];
                    temp.FlipV(); // in order to begin on a corner vertex
                    if(featureVertices.find(tri::Index(*booleanmesh,temp.V()))!=featureVertices.end())
                        featureThis.insert(j);
                    Point2d tempoint=parameterization.addSegmentPOS(*booleanmesh,temp);
                    eboundary2D.push_back(Point2f(tempoint.X(),tempoint.Y()));
                }

                if(!utility::isCounterClockWise(eboundary2D)){
                    std::reverse(eboundary2D.begin(),eboundary2D.end());
                    int sizeThis=eboundary2D.size();
                    set<int> newfeature;
                    for(set<int>::iterator its=featureThis.begin();its!=featureThis.end();its++)
                        newfeature.insert(sizeThis-(*its)-1);
                    utility::toSingleBoundary(eboundary2D,holes2D,temp3DFrontier,newfeature);
                    featureVertices=newfeature;
                }
                else{
                    utility::toSingleBoundary(eboundary2D,holes2D,temp3DFrontier,featureThis);
                    featureVertices=featureThis;
                }

            }
            else{ // If it has no holes then we store all sides in sequence

                vector<FacePos> tempFrontier;
                for(int j=0;j<components[i].sidesT.size();j++){
                    tempFrontier.insert(tempFrontier.end(),components[i].sidesT[j].directedside.begin(),components[i].sidesT[j].directedside.end());
                    featureVertices.insert(vcg::tri::Index(*booleanmesh,components[i].sidesT[j].end.V()));
                }
                if(!tempFrontier.empty()){
                    vector<Point2f> temp2DFrontier;
                    set<int> featurethis=featureVertices;
                    featureVertices.clear();
                    for(int j=0;j<tempFrontier.size();j++){
                        FacePos temp=tempFrontier[j];
                        temp.FlipV(); // in order to begin on a corner vertex

                        if(featurethis.find(tri::Index(*booleanmesh,temp.V()))!=featurethis.end())
                            featureVertices.insert(j);

                        Point2d tempoint=parameterization.addSegmentPOS(*booleanmesh,temp); //2D projection
                        temp2DFrontier.push_back(Point2f(tempoint.X(),tempoint.Y()));
                        temp3DFrontier.push_back(Point3f(temp2DFrontier[j].X(),temp2DFrontier[j].Y(),0.0));
                    }
                }
            }
            if(!temp3DFrontier.empty()){
                if(!utility::isCounterClockWise(temp3DFrontier)){
                    std::reverse(temp3DFrontier.begin(),temp3DFrontier.end());
                    set<int> newfeaturevertices;
                    int sizetemp=temp3DFrontier.size();
                    for(const int &fid : featureVertices)
                        newfeaturevertices.insert(sizetemp-1-fid);
                    featureVertices=newfeaturevertices;
                }
                listFrontiers.push_back(temp3DFrontier);
                feature.push_back(featureVertices);
            }

        }
    }*/
    /*
     * Version 2, for polygons
     * */
    void setBoundaries2D(){
        //exportProjectedTriangulation();
        vector<vector<Point3f>> polygons;  //list of polygons to subdivide
        vector<set<int>> features_on_polygons;  //list of corresponding feature vertices of the polygons to subdivide
        get2DFrontiers3(polygons,features_on_polygons);
        /*for(size_t i=0;i<polygons.size();i++){
            utility::saveAsCurvePly(polygons[i],checkPath("PolygonsToSplit")+"/polytosee"+std::to_string(i)+".ply",true);
        }*/
        //get2DFrontiers2(polygons,features_on_polygons);
        //get2DFrontiersLowLevel<CurveClass>(polygons,features_on_polygons,intersection_curve);
        //assert(components.size()==polygons.size());
        QString auc=QTime::currentTime().toString();
        for(size_t i=0;i<polygons.size();i++){

            //projecting the boundaries of the 3d patch. This in order to do a match after
            vector<vector<Point2<CMesh::ScalarType>>> projected;
            for(int j=0;j<components[i].sidesT.size();j++){                
                vector<Point2<CMesh::ScalarType>> aux;
                FacePos auxpos=components[i].sidesT[j].directedside[0];
                auxpos.FlipV();
                aux.push_back(components[i].parameterization.project(auxpos));
                for(int k=0;k<components[i].sidesT[j].size;k++){
                    aux.push_back(components[i].parameterization.project(components[i].sidesT[j].directedside[k]));
                }                
                projected.push_back(aux);
            }
            //creating the sides of this polygon
            set<int> featuresthis=features_on_polygons[i];
            vector<MonotonicFunctions::PolygonSide> sidesthis;
            int start=0;            
            int psize=polygons[i].size();
            for(int j=0;j<psize;j++){
                if(featuresthis.find(j)!=featuresthis.end()){
                    start=j;
                    break;
                }
            }
            int end=start;
            vector<Point3f> aux;
            do{
               aux.push_back(polygons[i][start]);
               if(featuresthis.find(start)!=featuresthis.end() && aux.size()>1){ // to avoid sides with one point
                   sidesthis.push_back(MonotonicFunctions::PolygonSide(aux));
                   /*cout<<"side "<<endl;
                   for(int r=0;r<aux.size();r++)
                       std::cout << "[ "<<aux[r].X() <<" "<<aux[r].Y()<<" ] ,";
                   cout<<endl;
                   */
                   aux.clear();
                   aux.push_back(polygons[i][start]);
               }
               start=(start+1)%psize;
            }while(start!=end);
            aux.push_back(polygons[i][end]);            
            sidesthis.push_back(MonotonicFunctions::PolygonSide(aux));
            assert(sidesthis.back().vertices.back()==sidesthis[0].vertices[0]);
            for(int j=0;j<sidesthis.size();j++){
                int sidessize=sidesthis[j].size;
                Point2<CMesh::ScalarType> totest2d;
                int totest;
                if(sidessize>2)
                    totest=int(sidessize/2);
                else
                    totest=1;

                Point3f totest3d=sidesthis[j].vertices[totest];
                totest2d=Point2<CMesh::ScalarType>(totest3d.X(),totest3d.Y());
                //searching where this point belongs
                bool found=false;
                for(int k=0;k<projected.size();k++){
                    for(int r=0;r<projected[k].size();r++)
                    {
                        if(vcg::Distance(projected[k][r],totest2d)<1e-6){
                           if(sidessize==2 && projected[k].size()!=1)
                                break;
                           sidesthis[j].index_3D=k;
                           sidesthis[j].flag=S_FIXED;
                        }
                    }
                    if(found)
                        break;
                }
            }
            //Note that the polygon is in counterclockwise order (polygons in counterclockwise). This is important for the quadrangulator
            MonotonicFunctions::Polygon polygonthis(sidesthis);
            if(polygonthis.size()==92)
                cout<<"size before 92"<< endl;
            if(polygonthis.size()==93)
                cout<<"size before 93"<< endl;
            polygonthis.compacting();

            components[i].boundaries2d.set(polygonthis);
        }
    }
    //getting the average edge length of the component c considering only the boundary-non intetsecting edges
    ScalarType getLengthOnlyMeshEdges(size_t c){
        ScalarType result=0;
        size_t count=0;
        for(size_t i=0;i<components[c].sidesT.size();i++){
            if(!components[c].sidesT[i].iscurve){
                result+=components[c].sidesT[i].length();
                count+=components[c].sidesT[i].size;
            }
        }
        if(count!=0)
            result/=count;
        else{
            // if the unique border is an intersection curve we take the average edgelength
            for(size_t i=0;i<components[c].faces.size();i++){
                result+=vcg::Perimeter(booleanmesh->face[components[c].faces[i]]);
            }
            result/=(3*components[c].faces.size());
        }
        return result;
    }
    /*
     * To debug. Get a list of 3D boundaries formed by all boundaries of all connected components projected to the unit square
     * @param lisfFrontiers is the collection of all boundaries of all  connected components of the patch
     */
    void get3DFrontiers(vector<vector<Point3f>> & listFrontiers){
        if(numberpartitionQuadMesh==-1) return;
        listFrontiers.clear();
        assert(numberConnectedComponents==components.size());
        for(int i=0;i<components.size();i++){
            vector<Point3f> temp3DFrontier;

            // If this component has holes then we store these in a array and join properly with the exterior boundary
            if(components[i].hasholes){
                vector<vector<Point3f>> holes3D;
                vector<FacePos> eboundary;
                vector<Point3f> eboundary3D;
                for(int j=0;j<components[i].sidesT.size();j++){
                    if(components[i].sidesT[j].isclosed && components[i].sidesT[j].iscurve){
                        vector<Point3f> tempholes;
                        for(int k=0;k<components[i].sidesT[j].size;k++)
                            tempholes.push_back(components[i].sidesT[j].directedside[k].V()->P());
                        holes3D.push_back(tempholes);
                    }
                }
                for(int j=0;j<components[i].sidesT.size();j++){ // this works because we first store all holes and after the other sides
                    if(!(components[i].sidesT[j].isclosed && components[i].sidesT[j].iscurve)){
                            eboundary.insert(eboundary.end(),components[i].sidesT[j].directedside.begin(),components[i].sidesT[j].directedside.end());
                    }
                }
                assert(eboundary[0].VFlip()==eboundary.back().V());
                for(int j=0;j<eboundary.size();j++)
                    eboundary3D.push_back(eboundary[j].V()->P());

                listFrontiers.insert(listFrontiers.end(),holes3D.begin(),holes3D.end());
                temp3DFrontier=eboundary3D;
            }
            else{ // If it has no holes then we store all sides in sequence

                vector<FacePos> tempFrontier;
                for(int j=0;j<components[i].sidesT.size();j++)
                    tempFrontier.insert(tempFrontier.end(),components[i].sidesT[j].directedside.begin(),components[i].sidesT[j].directedside.end());
                if(!tempFrontier.empty()){
                    for(int j=0;j<tempFrontier.size();j++){
                        temp3DFrontier.push_back(tempFrontier[j].VFlip()->P());
                    }
                }
            }
            if(!temp3DFrontier.empty())
                listFrontiers.push_back(temp3DFrontier);

        }
    }
    vector<vector<vector<Point3f>>> processAsPolygons(){
        vector<vector<vector<Point3f>>> patches;
        for(int i=0;i<components.size();i++){
            vector<vector<vector<Point3f>>> temppatches=components[i].boundaries2d.getPolygonsAsPrePatch();
            patches.insert(patches.end(),temppatches.begin(),temppatches.end());
        }
        return patches;
    }

    /*void processAsPatches(vector<std::pair<Eigen::MatrixXd,Eigen::MatrixXi>> &result){
        result.clear();
        vector<vector<vector<Point3f>>> patches;
        for(int i=0;i<components.size();i++){
            vector<vector<vector<Point3f>>> temppatches=components[i].boundaries2d.getPolygonsAsPrePatch();
            patches.insert(patches.end(),temppatches.begin(),temppatches.end());
        }
        vector<std::pair<Eigen::MatrixXd,Eigen::MatrixXi>> patches_as_mesh;
        //int totalvertices=0;
        //int totalfaces=0;
        UpdateTopology<CMesh>::FaceFace(*tmesh);
        for(int i=0;i<patches.size();i++){
            PatchG<PMeshType,patchgen::PatchParam> patchDS;
            if(patches[i].size()<2 || patches[i].size()>6 ){
                cout<<"problem with patches number sides"<<endl;
                continue;
            }
            else{
                patchDS.create_and_process(patches[i]);
                int Mask=0;
                int msc;
                QTime::fromMSecsSinceStartOfDay(msc);
                string aux="../exported/patchmesh2D"+std::to_string(patchDS.numberSides)+" - "+std::to_string(i) + "- "+QString::number(msc).toStdString()+".obj";
                vcg::tri::io::ExporterOBJ<PMeshType>::Save(patchDS.mesh,aux.c_str(), Mask);
                if(patchDS.finish){
                    vcg::tri::UpdateColor<PMeshType>::PerFaceConstant(patchDS.mesh,vcg::Color4b::Scatter(patches.size(),i));
                    //vcg::tri::Append<PMeshType,PMeshType>::Mesh(result,patchDS.mesh); //don't work
                    Eigen::MatrixXd V(patchDS.mesh.vert.size(),3);
                    Eigen::MatrixXi F(patchDS.mesh.face.size(),4);
                    utility::convertToLIBIGL<PMeshType>(patchDS.mesh,V,F);
                    parameterization.updateThreshold(V);
                    parameterization.unproject2D_and_update(V);
                    //totalvertices+=V.rows();
                    //totalfaces+=F.rows();
                    patches_as_mesh.push_back(make_pair(V,F));
                }
            }
        }
        result=patches_as_mesh;
        // creating a unique mesh
        //utility::mergeQuadMeshes_LibiglFormat(patches_as_mesh,result);
    }
    void processAsPatches(PMeshType &result){
        result.Clear();
        UpdateTopology<CMesh>::FaceFace(*tmesh);
        for(int i=0;i<components.size();i++){
            vector<vector<vector<Point3f>>> patches=components[i].boundaries2d.getPolygonsAsPrePatch();
            for(int i=0;i<patches.size();i++){
                PatchG<PMeshType,patchgen::PatchParam> patchDS;
                if(patches[i].size()<2 || patches[i].size()>6 ){
                    cout<<"problem with patches number sides"<<endl;
                    continue;
                }
                else{
                    patchDS.create_and_process(patches[i]);
                    //string aux="../exported/patchmesh2D"+std::to_string(patchDS.numberSides)+" - "+std::to_string(i) + "- "+QString::number(msc).toStdString()+".obj";
                    //vcg::tri::io::ExporterOBJ<PMeshType>::Save(patchDS.mesh,aux.c_str(), Mask);
                    if(patchDS.finish){
                        //vcg::tri::UpdateColor<PMeshType>::PerFaceConstant(patchDS.mesh,vcg::Color4b::Scatter(patches.size(),i));
                        vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormalByFitting(patchDS.mesh);
                        vcg::tri::UpdateBounding<PMeshType>::Box(patchDS.mesh);
                        vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(patchDS.mesh);
                        vcg::tri::UpdateTopology<PMeshType>::FaceFace(patchDS.mesh);
                        vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(patchDS.mesh);

                        vcg::tri::Append<PMeshType,PMeshType>::Mesh(result,patchDS.mesh);

                        vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormalByFitting(result);
                        vcg::tri::UpdateBounding<PMeshType>::Box(result);
                        vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(result);                        
                        vcg::tri::UpdateTopology<PMeshType>::FaceFace(result);
                        vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(result);

                    }
                }
            }
        }
        Clean<PMeshType>::RemoveDuplicateVertex(result);

        double len=utility::minimum_length_quadmesh<PMeshType>(result);
        tri::MeshAssert<PMesh>::OnlyQuadFace(result);
            cout<<"is quad before merging "<<endl;
        Clean<PMeshType>::MergeCloseVertex(result,len/10);
        tri::MeshAssert<PMesh>::OnlyQuadFace(result);
            cout<<"is quad after"<<endl;
        //vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormalByFitting(result);
        vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormals(result);
        vcg::tri::UpdateBounding<PMeshType>::Box(result);
        vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(result);        
        vcg::tri::UpdateTopology<PMeshType>::FaceFace(result);
        vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(result);
        int Mask=0;
        string aux2D=checkPath("");
        aux2D+="/quadrangulated2D.obj";
        vcg::tri::io::ExporterOBJ<PMeshType>::Save(result,aux2D.c_str(), Mask);

        Eigen::MatrixXd V(result.vert.size(),3);
        Eigen::MatrixXi F(result.face.size(),4);
        utility::convertToLIBIGL<PMeshType>(result,V,F);
        parameterization.updateThreshold(V);
        parameterization.unproject2D_and_update(V);
        result.Clear();
        utility::Quad_libiglToVCG<PMeshType>(result,V,F);
        vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormalByFitting(result);
        vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(result);

        Mask=0;
        string aux3D=checkPath("");
        aux3D+="/quadrangulated3D.obj";
        vcg::tri::io::ExporterOBJ<PMeshType>::Save(result,aux3D.c_str(), Mask);
    }*/
    void processAsPatchesLSCM(PMeshType &result){
        result.Clear();
        UpdateTopology<CMesh>::FaceFace(*tmesh);
        for(int i=0;i<components.size();i++){
            PMeshType tempresult;
            vector<vector<vector<Point3f>>> patches=components[i].boundaries2d.getPolygonsAsPrePatch();
            for(size_t j=0;j<patches.size();j++){
                PatchG<PMeshType,patchgen::PatchParam> patchDS;
                if(patches[j].size()<2 || patches[j].size()>6 ){
                    cout<<"problem with patches number sides"<<endl;
                    continue;
                }
                else{                    
                    //Making a quadrangulation of this patch
                    patchDS.create_and_process(patches[j]);
                    vcg::PolygonalAlgorithm<PMesh>::UpdateFaceNormals(patchDS.mesh);
                    vcg::tri::UpdateBounding<PMesh>::Box(patchDS.mesh);
                    vcg::tri::UpdateTopology<PMesh>::FaceFace(patchDS.mesh);
                    vcg::tri::UpdateFlags<PMesh>::VertexBorderFromFaceAdj(patchDS.mesh);
                    vcg::PolygonalAlgorithm<PMesh>::LaplacianReproject(patchDS.mesh);
                    OrientFaces<PMesh>::AutoOrientFaces(patchDS.mesh);
                    if(patchDS.finish){
                        //vcg::tri::UpdateColor<PMeshType>::PerFaceConstant(patchDS.mesh,vcg::Color4b::Scatter(patches.size(),i));
                        //vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormalByFitting(patchDS.mesh);
                        vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormals(patchDS.mesh);
                        vcg::tri::UpdateBounding<PMeshType>::Box(patchDS.mesh);
                        vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(patchDS.mesh);
                        vcg::tri::UpdateTopology<PMeshType>::FaceFace(patchDS.mesh);
                        vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(patchDS.mesh);
                        vcg::tri::MeshAssert<PMesh>::OnlyQuadFace(patchDS.mesh);
                        if(exportIntermediateSteps){
                            int Mask=0;
                            string aux2D=checkPath("PolygonsToSplit");
                            aux2D+="/quadrangulated2D-subpatch-"+std::to_string(i)+"-"+std::to_string(j)+".obj";
                            vcg::tri::io::ExporterOBJ<PMeshType>::Save(patchDS.mesh,aux2D.c_str(), Mask);
                        }
                        vcg::tri::Append<PMeshType,PMeshType>::Mesh(tempresult,patchDS.mesh);

                        vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormals(tempresult);
                        vcg::tri::UpdateBounding<PMeshType>::Box(tempresult);
                        vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(tempresult);
                        vcg::tri::UpdateTopology<PMeshType>::FaceFace(tempresult);
                        vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(tempresult);
                    }
                }
            }


            // cleaning steps
            PMesh::ScalarType min_edge_len=utility::minimum_length_quadmesh(tempresult);
            vcg::tri::Clean<PMesh>::RemoveDuplicateVertex(tempresult);
            vcg::tri::Clean<PMesh>::MergeCloseVertex(tempresult,min_edge_len/100);
            // saving 2D models
            if(exportIntermediateSteps){
                int maskaux=0;
                string aux2D=checkPath("PolygonsToSplit");
                aux2D+="/quadrangulated2D-total.obj";
                vcg::tri::io::ExporterOBJ<PMeshType>::Save(tempresult,aux2D.c_str(),maskaux);
            }

            // getting the triangle version
            CMesh tri_tempresult;
            components[i].parameterization.exportMeshParametrization(tri_tempresult);
            if(exportIntermediateSteps){
                 string aux3Dtriangle=checkPath("PolygonsToSplit");
                 aux3Dtriangle+="/triangulatedPatch.ply";
                 components[i].parameterization.exportMeshAsTrianglePatch(aux3Dtriangle);
            }
            components[i].parameterization.reprojecting_and_update(tempresult);
            vcg::tri::Append<PMeshType,PMeshType>::Mesh(result,tempresult);

            vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormals(result);
            vcg::tri::UpdateBounding<PMeshType>::Box(result);
            vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(result);
            vcg::tri::UpdateTopology<PMeshType>::FaceFace(result);
            vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(result);

        }
        Clean<PMeshType>::RemoveDuplicateVertex(result);

        double len=utility::minimum_length_quadmesh<PMeshType>(result);
        tri::MeshAssert<PMesh>::OnlyQuadFace(result);
            cout<<"is quad before merging "<<endl;
        Clean<PMeshType>::MergeCloseVertex(result,len/10);
        tri::MeshAssert<PMesh>::OnlyQuadFace(result);
            cout<<"is quad after"<<endl;
        //vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormalByFitting(result);
        vcg::tri::UpdateBounding<PMeshType>::Box(result);
        vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(result);
        //vcg::tri::UpdateTopology<PMeshType>::FaceFace(result);
        //vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(result);
        if(exportIntermediateSteps){
            int Mask=0;
            string aux3D=checkPath("");
            aux3D+="/quadrangulated3D.obj";
            vcg::tri::io::ExporterOBJ<PMeshType>::Save(result,aux3D.c_str(), Mask);
        }
        //setting mesh number property
        AttributeFaceType hMeshNumberTri= vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(*booleanmesh,std::string("MeshNumber"));
        int numberOfMesh=hMeshNumberTri[components[0].faces[0]];
        for(PMesh::FaceIterator qfi=result.face.begin();qfi!=result.face.end();qfi++){
            qfi->Q()=numberOfMesh;
        }

    }
    void processAsPatchesLSCM2(PMeshType &result){
        result.Clear();
        UpdateTopology<CMesh>::FaceFace(*tmesh);
        for(int i=0;i<components.size();i++){
            PMeshType tempresult;
            vector<vector<vector<Point3f>>> patches=components[i].boundaries2d.getPolygonsAsPrePatch();
            for(size_t j=0;j<patches.size();j++){
                PatchG<PMeshType,patchgen::PatchParam> patchDS;
                if(patches[j].size()<2 || patches[j].size()>6 ){
                    cout<<"problem with patches number sides"<<endl;
                    continue;
                }
                else{

                    string auxpoly=checkPath("PolygonsToSplit");
                    string auxpoly_flattened=auxpoly;
                    auxpoly+="/polyborder-subpatch-"+std::to_string(i)+"-"+std::to_string(j)+".txt";
                    auxpoly_flattened+="/flattened_polyborder-subpatch-"+std::to_string(i)+"-"+std::to_string(j)+".txt";
                    vector<vector<Point3f>> flattened_poly;
                    //utility::polygon_to_regular_polygon(patches[j],flattened_poly);
                    //utility::polygon_to_txt<CMesh::ScalarType>(patches[j],auxpoly);
                    //utility::polygon_to_txt<CMesh::ScalarType>(flattened_poly,auxpoly_flattened);

                    //Creating a triangulation of this patch
                    CMesh triangle,boundary;
                    vector<CMesh::CoordType> points,inner;
                    vector<pair<size_t,size_t>> indexes;
                    vector<CMesh::CoordType> positions;
                    for(size_t k=0;k<patches[j].size();k++){
                        for(size_t r=0;r<patches[j][k].size()-1;r++){
                            tri::Allocator<CMesh>::AddVertex(boundary,patches[j][k][r]);
                            points.push_back(patches[j][k][r]);
                            positions.push_back(flattened_poly[k][r]);
                            indexes.push_back(make_pair(k,r));
                        }
                    }
                    cdt2d<CMesh::CoordType> triangle_cdt(points,inner);
                    vector<CMesh::CoordType> triangle_vertices;
                    vector<vector<unsigned int>> triangle_faces_temp;
                    vector<Point3i> triangle_faces_final;
                    triangle_cdt.get_triangulation(triangle_vertices,triangle_faces_temp);
                    for(size_t k=0;k<triangle_faces_temp.size();k++){
                        Point3i thisp(triangle_faces_temp[k][0],triangle_faces_temp[k][1],triangle_faces_temp[k][2]);
                        triangle_faces_final.push_back(thisp);
                    }
                    tri::BuildMeshFromCoordVectorIndexVector(triangle,triangle_vertices,triangle_faces_final);
                    tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFaceNormalized(triangle);
                    tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFaceNormalized(triangle);
                    UpdateTopology<CMesh>::FaceFace(triangle);
                    UpdateFlags<CMesh>::FaceBorderFromFF(triangle);
                    CMesh::ScalarType edgeaverage=0;
                    size_t s=0;
                    for(CMesh::FaceIterator fi=triangle.face.begin();fi!=triangle.face.end();fi++){
                        for(size_t k=0;k<3;k++){
                            if(vcg::face::IsBorder(*fi,k)){
                                edgeaverage+=vcg::Distance(fi->V(k)->P(),fi->V(fi->Next(k))->P());
                                s++;
                            }
                        }
                    }
                    edgeaverage/=s;

                    tri::UpdateBounding<CMesh>::Box(triangle);
                    CMesh::ScalarType samplingRadius =2*edgeaverage;
                    vector<CMesh::CoordType> pointVec;
                    int dx=triangle.bbox.DimX()/samplingRadius;
                    int dy=triangle.bbox.DimY()/samplingRadius;
                    tri::PoissonSampling<CMesh>(triangle,pointVec,0,samplingRadius);
                    CMesh BasicPoissonMesh;
                    tri::BuildMeshFromCoordVector(BasicPoissonMesh,pointVec);

                    cdt2d<CMesh::CoordType> triangulation(points,pointVec);
                    vector<CMesh::CoordType> vertices;
                    vector<vector<unsigned int>> faces_temp;
                    vector<Point3i> faces_final;
                    triangulation.get_triangulation(vertices,faces_temp);
                    for(size_t k=0;k<faces_temp.size();k++){
                        Point3i thisp(faces_temp[k][0],faces_temp[k][1],faces_temp[k][2]);
                        faces_final.push_back(thisp);
                    }
                    tri::BuildMeshFromCoordVectorIndexVector(triangle,vertices,faces_final);
                    utility::saveMesh(triangle,checkPath("PolygonsToSplit")+"/triangulation-subpatch"+std::to_string(i)+"-"+std::to_string(j)+".ply");

                    UpdateTopology<CMesh>::FaceFace(triangle);
                    UpdateFlags<CMesh>::VertexBorderFromFaceAdj(triangle);
                    CMesh projected;
                    tri::Append<CMesh,CMesh>::MeshCopy(projected,triangle);

                    GridStaticPtr<CMesh::VertexType, CMesh::ScalarType> GRID;
                    GRID.Set(boundary.vert.begin(),boundary.vert.end());
                    vector<size_t> handle;
                    vector<CMesh::CoordType> positions2;
                    for(CMesh::VertexIterator vi=projected.vert.begin();vi!=projected.vert.end();vi++){
                        if(vi->IsB()){
                          CMesh::ScalarType maxdist=3000,mindist=0;
                          CMesh::VertexType* vert = tri::GetClosestVertex<CMesh,GridStaticPtr<CMesh::VertexType,CMesh::ScalarType>>(boundary,GRID,vi->P(),maxdist,mindist);
                          assert(vert!=NULL);
                          pair<size_t,size_t> values=indexes[tri::Index(boundary,vert)];
                          handle.push_back(tri::Index(projected,&*vi));
                          positions2.push_back(flattened_poly[values.first][values.second]);
                        }
                    }
                    std::vector<size_t> roi(projected.vert.size()) ;
                    std::iota(std::begin(roi), std::end(roi), 0);
                    LaplacianSurfaceDeformation<CMesh> lsd(projected,roi,handle);
                    lsd.prepareDeform();
                    lsd.doDeform(positions2,LaplacianSurfaceDeformation<CMesh>::TypeDeformation(0));

                    //Making a quadrangulation of the flattened version
                    patchDS.create_and_process(flattened_poly);

                    // Going back to the original shape

                    vcg::tri::UpdateNormal<PMesh>::PerVertexNormalized(patchDS.mesh);
                    vcg::tri::UpdateBounding<PMesh>::Box(patchDS.mesh);
                    vcg::tri::UpdateTopology<PMesh>::FaceFace(patchDS.mesh);
                    vcg::tri::UpdateFlags<PMesh>::VertexBorderFromFaceAdj(patchDS.mesh);

                    vcg::PolygonalAlgorithm<PMesh>::LaplacianReproject(patchDS.mesh);
                    if(exportIntermediateSteps){
                        string auxfla=checkPath("PolygonsToSplit")+"/quadrangulation-flattened-subpatch-"+std::to_string(i)+"-"+std::to_string(j)+".obj";
                        vcg::tri::io::ExporterOBJ<PMesh>::Save(patchDS.mesh,auxfla.c_str(),vcg::tri::io::Mask::IOM_BITPOLYGONAL);
                    }
                        //vcg::PolygonalAlgorithm<PMesh>::LaplacianReproject<CMesh>(patches.mesh,triangle,2);
                        //vcg::PolygonalAlgorithm<PMesh>::SmoothReprojectPCA(patches.mesh);
                        //OrientFaces<PMesh>::AutoOrientFaces(patchDS.mesh);
                        //vcg::PolygonalAlgorithm<PMesh>::UpdateFaceNormalByFitting(patches.mesh);

                    OrientFaces<PMesh>::AutoOrientFaces(patchDS.mesh);
                    OrientFaces<PMesh>::AutoOrientFaces(patchDS.mesh);

                    vector<Point4<CMesh::ScalarType>> coordinates;
                    utility::getBarycentrycCoordinatesQ2T<CMesh,PMesh>(projected,patchDS.mesh,coordinates);
                     for(size_t k=0;k<patchDS.mesh.vert.size();k++){
                         if(!patchDS.mesh.vert[k].IsB()){
                             size_t index_face=size_t(coordinates[k][0]);
                             CMesh::CoordType bc(coordinates[k][1],coordinates[k][2],coordinates[k][3]);
                             CMesh::FaceType fi=triangle.face[index_face];
                             PMesh::CoordType new_P=fi.V(0)->P()*bc.X()+fi.V(1)->P()*bc.Y()+fi.V(2)->P()*bc.Z();
                             patchDS.mesh.vert[k].P()=new_P;
                         }
                         else{
                             int index_poly=-1;
                             PMesh::CoordType pp=patchDS.mesh.vert[k].P();
                             CMesh::CoordType ppc(pp.X(),pp.Y(),pp.Z());
                             CMesh::ScalarType mindist=triangle.bbox.Diag();
                             for(size_t ii=0;ii<positions.size();ii++){
                                CMesh::ScalarType dist=vcg::Distance(ppc,positions[ii]);
                                if(dist<mindist){
                                    mindist=dist;
                                    index_poly=ii;
                                }
                             }
                             //assert(vert->IsB());
                             pair<size_t,size_t> values=indexes[index_poly];
                             CMesh::CoordType npc=patches[j][values.first][values.second];
                             PMesh::CoordType np(npc.X(),npc.Y(),npc.Z());
                             patchDS.mesh.vert[k].P()=np;
                         }
                     }
                    vcg::PolygonalAlgorithm<PMesh>::UpdateFaceNormals(patchDS.mesh);
                    vcg::tri::UpdateNormal<PMesh>::PerVertexNormalizedPerFaceNormalized(patchDS.mesh);
                    if(exportIntermediateSteps){
                        string auxcurrent=checkPath("PolygonsToSplit")+"/quadrangulation-current-subpatch-"+std::to_string(i)+"-"+std::to_string(j)+".obj";
                        vcg::tri::io::ExporterOBJ<PMesh>::Save(patchDS.mesh,auxcurrent.c_str(),vcg::tri::io::Mask::IOM_BITPOLYGONAL);
                    }
                    //vcg::PolygonalAlgorithm<PMesh>::LaplacianReproject(patchDS.mesh);

                    if(patchDS.finish){
                        //vcg::tri::UpdateColor<PMeshType>::PerFaceConstant(patchDS.mesh,vcg::Color4b::Scatter(patches.size(),i));
                        //vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormalByFitting(patchDS.mesh);
                        vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormals(patchDS.mesh);
                        vcg::tri::UpdateBounding<PMeshType>::Box(patchDS.mesh);
                        vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(patchDS.mesh);
                        vcg::tri::UpdateTopology<PMeshType>::FaceFace(patchDS.mesh);
                        vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(patchDS.mesh);
                        vcg::tri::MeshAssert<PMesh>::OnlyQuadFace(patchDS.mesh);
                        if(exportIntermediateSteps){
                            int Mask=0;
                            string aux2D=checkPath("PolygonsToSplit");
                            aux2D+="/quadrangulated2D-subpatch"+std::to_string(i)+"-"+std::to_string(j)+".obj";
                            vcg::tri::io::ExporterOBJ<PMeshType>::Save(patchDS.mesh,aux2D.c_str(), Mask);
                        }

                        vcg::tri::Append<PMeshType,PMeshType>::Mesh(tempresult,patchDS.mesh);

                        vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormals(tempresult);
                        vcg::tri::UpdateBounding<PMeshType>::Box(tempresult);
                        vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(tempresult);
                        vcg::tri::UpdateTopology<PMeshType>::FaceFace(tempresult);
                        vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(tempresult);
                    }
                }
            }
            string aux2D=checkPath("PolygonsToSplit");
            aux2D+="/quadrangulated2D-total.obj";
            // cleaning steps
            PMesh::ScalarType min_edge_len=utility::minimum_length_quadmesh(tempresult);
            vcg::tri::Clean<PMesh>::RemoveDuplicateVertex(tempresult);
            vcg::tri::Clean<PMesh>::MergeCloseVertex(tempresult,min_edge_len/100);

            // saving 2D models
            if(exportIntermediateSteps){
                int maskaux=0;
                vcg::tri::io::ExporterOBJ<PMeshType>::Save(tempresult,aux2D.c_str(),maskaux);
            }
            // getting the triangle version
            CMesh tri_tempresult;
            components[i].parameterization.exportMeshParametrization(tri_tempresult);


            components[i].parameterization.reprojecting_and_update(tempresult);

            vcg::tri::Append<PMeshType,PMeshType>::Mesh(result,tempresult);

            vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormals(result);
            vcg::tri::UpdateBounding<PMeshType>::Box(result);
            vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(result);
            //vcg::tri::io::ExporterOBJ<PMeshType>::Save(result,"last.obj", vcg::tri::io::Mask::IOM_BITPOLYGONAL);
            vcg::tri::UpdateTopology<PMeshType>::FaceFace(result);
            vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(result);

        }
        Clean<PMeshType>::RemoveDuplicateVertex(result);

        double len=utility::minimum_length_quadmesh<PMeshType>(result);
        tri::MeshAssert<PMesh>::OnlyQuadFace(result);
            cout<<"is quad before merging "<<endl;
        Clean<PMeshType>::MergeCloseVertex(result,len/10);
        tri::MeshAssert<PMesh>::OnlyQuadFace(result);
            cout<<"is quad after"<<endl;
        //vcg::PolygonalAlgorithm<PMeshType>::UpdateFaceNormalByFitting(result);
        vcg::tri::UpdateBounding<PMeshType>::Box(result);
        vcg::tri::UpdateNormal<PMeshType>::PerVertexNormalized(result);
        //vcg::tri::UpdateTopology<PMeshType>::FaceFace(result);
        //vcg::tri::UpdateFlags<PMeshType>::VertexBorderFromFaceAdj(result);
        if(exportIntermediateSteps){
            int Mask=0;
            string aux3D=checkPath("");
            aux3D+="/quadrangulated3D.obj";
            vcg::tri::io::ExporterOBJ<PMeshType>::Save(result,aux3D.c_str(), Mask);
        }
    }
    /*
     * test if a Pos belongs to some side of some connected component of this patch
     * Returns the necessary info in the affirmative case
     */
    bool isPosOnThisPatch(FacePos &pos,Eigen::Vector4i &info,bool testopposite=false){
        info.setConstant(-1);
        for(int i=0;i<components.size();i++){
            for (int j=0;j<components[i].sidesT.size();j++) {
                for(int k=0;k<components[i].sidesT[j].size;k++){
                    if(pos==components[i].sidesT[j].directedside[k]){
                        info[0]=i;info[1]=j;info[2]=k;info[3]=1;
                        return true;
                    }
                    if(testopposite){
                        FacePos copy=pos;
                        copy.FlipV();
                        if(copy==components[i].sidesT[j].directedside[k]){
                            info[0]=i;info[1]=j;info[2]=k;info[3]=0;
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }
    /*
     * LowLevel. Return the first commom set of pos edges belonging to a side of the connected component with index "idc" and having the quality[][0] property equal to quality.
     * @param idc: index of the connected component to be tested
     * @param quality: quality value to be tested with all sides of this component
     * @param result: resulting vector of pos edges
     */
    int  getCommonPosByQuality(size_t idc,int quality, vector<size_t> &result){
        size_t ns=components[idc].sidesT.size();
        int foundside=-1;
        result.clear();
        //for each side
        for(size_t s=0;s<ns;s++){
            size_t nts=components[idc].sidesT[s].size;
            // for each edge in this side
            for(size_t k=0;k<nts;k++){
                if(components[idc].sidesT[s].quality[k][0]==quality){
                    result.push_back(k);
                }
            }
            if(!result.empty()){
                foundside=s;
                break;
            }
        }
        return foundside;
    }
    /*void exportProjectedTriangulation(){
        //UpdateFlags<CMesh>::VertexClearV(*booleanmesh);

        // Getting the min and max indexes of the vertices of this patch
        int min=std::max(int(booleanmesh->VN()),int(booleanmesh->vert.size()));
        int max=0;
        for(vector<size_t>::iterator it=PatchFaces.begin();it!=PatchFaces.end();it++){
            for(int i=0;i<3;i++){
                booleanmesh->face[*it].V(i)->ClearV();
                int index=tri::Index(*booleanmesh,booleanmesh->face[*it].V(i));
                if(index>max)
                    max=index;
                if(index<min)
                    min=index;
            }
        }

        //creating a 2d mesh of the triangulation
        TriMeshType projected;
        TriMeshType original;
        vector<Point3i> faceindexes;
        vector<CoordType> coordinates2D;
        vector<CoordType> coordinates3D;
        int numberVertices=max-min+1;
        for(int i=0;i<numberVertices;i++){
            coordinates2D.push_back(CoordType(0,0,0));
            coordinates3D.push_back(CoordType(0,0,0));
        }

        for(vector<size_t>::iterator it=PatchFaces.begin();it!=PatchFaces.end();it++){
            for(int i=0;i<3;i++){
                if(!booleanmesh->face[*it].V(i)->IsV()){
                    int index=tri::Index(*booleanmesh,booleanmesh->face[*it].V(i));
                    index-=min;
                    FacePos thispos(&(booleanmesh->face[*it]),i,booleanmesh->face[*it].V(i));
                    Point2d aux=parameterization.addSegmentPOS(*booleanmesh,thispos);
                    coordinates2D[index]=CoordType(aux.X(),aux.Y(),0);
                    coordinates3D[index]=CoordType(thispos.V()->P());
                    booleanmesh->face[*it].V(i)->SetV();
                }
            }
            int index0=tri::Index(*booleanmesh,booleanmesh->face[*it].V(0));
            int index1=tri::Index(*booleanmesh,booleanmesh->face[*it].V(1));
            int index2=tri::Index(*booleanmesh,booleanmesh->face[*it].V(2));
            index0-=min;
            index1-=min;
            index2-=min;
            faceindexes.push_back(Point3i(index0,index1,index2));
        }
        tri::BuildMeshFromCoordVectorIndexVector<TriMeshType,CoordType,Point3i>(projected,coordinates2D,faceindexes);
        tri::BuildMeshFromCoordVectorIndexVector<TriMeshType,CoordType,Point3i>(original,coordinates3D,faceindexes);
        vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalizedPerFace(projected);
        vcg::tri::UpdateNormal<TriMeshType>::PerFaceNormalized(projected);
        UpdateTopology<TriMeshType>::FaceFace(projected);
        vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalizedPerFace(original);
        vcg::tri::UpdateNormal<TriMeshType>::PerFaceNormalized(original);
        UpdateTopology<TriMeshType>::FaceFace(original);


        std::string meshstring=checkPath("");
        std::string patchstring="/projectedTriangulation.ply";
        std::string filename=meshstring+patchstring;
        std::string filename3D=meshstring+"/originalPatch3D.ply";
        tri::io::ExporterPLY<TriMeshType>::Save(projected,filename.c_str(),tri::io::Mask::IOM_WEDGTEXCOORD,false);
        tri::io::ExporterPLY<TriMeshType>::Save(original,filename3D.c_str(),tri::io::Mask::IOM_WEDGTEXCOORD,false);
    }*/
    static int convert34(int ti){
        if(ti%2==0){
            return ti/2;
        }
        else{
            return (ti-1)/2;
        }
    }
    /*
     * Print information about the patch
     */
    void printInformationPatch(){
        cout<<"-----------About the Patch ---------------"<<endl;
        cout<<"Number holes "<< numberHoles<<endl;
        cout<<"Number curved sides "<< numberCurveSides<<endl;
        cout<<"About components "<<endl;
        cout<<numberConnectedComponents<<" connected components: "<<endl;
        for(int i=0;i<components.size();i++){
            cout<<" number sides "<<i<<"-cc : "<<components[i].sidesT.size()<<endl;
            for(int j=0;j<components[i].sidesT.size();j++)
                if(components[i].sidesT[j].iscurve)
                    cout<<"  curved side with "<< components[i].sidesT[j].size<<" edges" <<endl;
            cout<<" number faces "<<i<<"-cc : " <<components[i].faces.size()<<endl<<endl;
        }
        cout<<"----------------------------------------- "<<endl;
    }
};
}
#endif // PATCH3D_H
