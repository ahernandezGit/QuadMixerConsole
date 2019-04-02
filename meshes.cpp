#include "meshes.h"
#include "aabb.h"
//#include "naturalparametrization.h"
#include "vectorcone.h"
#include "utils/unwgraph.h"
#include <algorithm>

meshes::meshes()
{
    reset();
}
void meshes::processall(){
    vector<vector<face::Pos<CMesh::FaceType>>> intersection=getIntersectionCurve();
    set<int> selectedPartitions1,selectedPartitions2;
    computeSelectedRegions(intersection,selectedPartitions1,selectedPartitions2);
    computePatches(selectedPartitions1,selectedPartitions2);
}

void meshes::conditioningIntersection(){
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    CMesh::PerFaceAttributeHandle<int> hCorrespondence = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("Correspondence"));

    //computing average edgelen of the quads around the intersection curves
    vector<vector<face::Pos<CMesh::FaceType>>> &curves=intcurves.original_curves;
    vector<CMesh::ScalarType> ael_ic(curves.size());
    for(size_t i=0;i<curves.size();i++){
        PMesh::ScalarType average=0;
        for(size_t j=0;j<curves[i].size();j++){
            CMesh::FacePointer fp=curves[i][j].F();
            CMesh::FacePointer fop=curves[i][j].FFlip();
            PMesh::FacePointer qf,qfop;;
            if(hMeshNumber[fp]==1){
                qf=&mesh1.face[convertIndexT2Q(hCorrespondence[fp])];
                qfop=&mesh2.face[convertIndexT2Q(hCorrespondence[fop])];
            }
            else{
                qf=&mesh2.face[convertIndexT2Q(hCorrespondence[fp])];
                qfop=&mesh1.face[convertIndexT2Q(hCorrespondence[fop])];
            }
            average+=utility::getAverageEdgeLengthFace<PMesh::FaceType>(*qf)/2;
            average+=utility::getAverageEdgeLengthFace<PMesh::FaceType>(*qf)/2;
        }
        average/=curves[i].size();
        ael_ic[i]=average;
    }

    mergeEndpointsPatches(); // some collapse operations
    getIntersectionCurve();
    splitIntersectionCurve();
    splitIntersectionCurve();
    //Allocator<CMesh>::CompactEveryVector(tboolean);
    UpdateColor<CMesh>::PerFaceConstant(tboolean);
    collapsingAfterSplit(ael_ic);
    selectFacesIntersection(true,false,2);
    tri::Smooth<CMesh>::VertexColorLaplacian(tboolean,1,true);
    selectFacesIntersection(false,true,2);
    tri::Smooth<CMesh>::VertexColorLaplacian(tboolean,1,true);
}
// It oonly uses the mesh number attribute
vector<vector<face::Pos<CMesh::FaceType>>> meshes::getIntersectionCurve(){
    // First pass: Finding all edges of the intersection curve
    tri::UpdateFlags<CMesh>::FaceClearV(tboolean);
    tri::UpdateFlags<CMesh>::VertexClearV(tboolean);
    vector<vector<face::Pos<CMesh::FaceType>>> curves;
    vector<face::Pos<CMesh::FaceType>> PosSoup;
    face::Pos<CMesh::FaceType> p;
    p.SetNull();
    //bool found=false;
    // Remember: if J(i) < FA.rows() then face i of the output mesh is a sub-triangle of the J(i)th face of A.
    // If J(i)>FA.rows() then face i of the output mesh is a sub-triangle of the J(i)-FA.rows()th face of B.
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    //CMesh::PerFaceAttributeHandle<int> hCorrespondence = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("Correspondence"));
    for(auto fi=tboolean.face.begin();fi!=tboolean.face.end();++fi)
        if(!(*fi).IsD()){
            for(int i=0;i<3;++i){
                face::Pos<CMesh::FaceType> pi(&*fi,i,fi->V(i));
                CMesh::FacePointer facep=pi.F();
                CMesh::FacePointer faceopp=pi.FFlip();
                //if(faceopp->IsV()) continue;
                //bool isIntersectionEdge=(criteria(tri::Index(tboolean,facep))<numberfacesA) && (criteria(tri::Index(tboolean,faceopp))>=numberfacesA);
                //isIntersectionEdge=isIntersectionEdge || ((criteria(tri::Index(tboolean,facep))>=numberfacesA) && (criteria(tri::Index(tboolean,faceopp))<numberfacesA));
                bool isIntersectionEdge= hMeshNumber[facep]!=hMeshNumber[faceopp];
                //faceopp->SetV();
                bool wasvisited= pi.V()->IsV() && pi.VFlip()->IsV();
                if(isIntersectionEdge && !wasvisited){
                    pi.V()->SetV();
                    pi.VFlip()->SetV();
                    PosSoup.push_back(pi);
                    //p=pi;
                    //found=true;
                    //break;
                }
            }
            //fi->SetV();
            //if(found) break;
        }
    tri::UpdateFlags<CMesh>::VertexClearV(tboolean);
    tri::UpdateFlags<CMesh>::VertexClearS(tboolean);
    tri::UpdateFlags<CMesh>::FaceClearS(tboolean);
    CMesh::VertexPointer vpt;
    //Second pass: Obtaining a circulation of the intersection curves.
    int count=0;
    while(PosSoup.size()!=0 && count<1000){
        face::Pos<CMesh::FaceType> currPos=PosSoup.back();
        vpt=currPos.V();
        if(!vpt->IsV() && !currPos.VFlip()->IsV()){
            curves.push_back(getIntersectionCurveFromPos(currPos));
        }
        PosSoup.pop_back();
        count++;
    }
    CurvesInBoolean=curves;
    intcurves.reset();
    intcurves.initialize(curves);
    return curves;

}
void meshes::computeSelectedRegions(vector<vector<face::Pos<CMesh::FaceType>>>& curves, set<int> &  selectedPartitions1, set<int> & selectedPartitions2){
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    CMesh::PerFaceAttributeHandle<int> hCorrespondence = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("Correspondence"));
    // Store to vizualization
    for(int i=0;i<curves.size();i++){
        for(int j=0;j<curves[i].size();j++){
            curves[i][j].V()->SetS();
            curves[i][j].FFlip()->SetS();
            curves[i][j].F()->SetS();
        }
    }

    //Recovering selected quads from selected triangles
    for(CMesh::FaceIterator fi = tboolean.face.begin(); fi != tboolean.face.end(); ++fi)
    {
        if(!(*fi).IsD() && (*fi).IsS()){
            int indexq=convertIndexT2Q(hCorrespondence[fi]);
            if(hMeshNumber[fi]==1){
                mesh1.face[indexq].SetS();
                selectedPartitions1.emplace(PTracer1.FacePatch[indexq]);
            }
            else{
                mesh2.face[indexq].SetS();
                selectedPartitions2.emplace(PTracer2.FacePatch[indexq]);
            }
        }
    }
    // Selecting quads from selected partitions
    // Initialize the  Mark component with the partition number
    for (size_t i=0;i<PTracer1.FacePatch.size();i++){
        mesh1.face[i].IMark()=PTracer1.FacePatch[i];
        tmesh1.face[2*i].IMark()=PTracer1.FacePatch[i];
        tmesh1.face[2*i+1].IMark()=PTracer1.FacePatch[i];
        if(selectedPartitions1.find(PTracer1.FacePatch[i])!=selectedPartitions1.end()){
           mesh1.face[i].SetS();
        }
    }

    for (size_t i=0;i<PTracer2.FacePatch.size();i++){
        mesh2.face[i].IMark()=PTracer2.FacePatch[i];
        tmesh2.face[2*i].IMark()=PTracer2.FacePatch[i];
        tmesh2.face[2*i+1].IMark()=PTracer2.FacePatch[i];
        if(selectedPartitions2.find(PTracer2.FacePatch[i])!=selectedPartitions2.end()){
           mesh2.face[i].SetS();
        }
    }

    // For "deleting" interior quads to the other mesh
    selectComponents1(true);
    selectComponents2(true);

}
void meshes::computePatches(std::set<int> &sc1,std::set<int> &sc2)
{
    QString as=QDateTime::currentDateTime().toString("dd-MM-yyyy__hh:mm:ss");
    string session = "session__"+as.toStdString();
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    CMesh::PerFaceAttributeHandle<int> hCorrespondence = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("Correspondence"));
    CollectionPatches.clear();

    //tri::UpdateSelection<CMesh>::FaceDilate(tboolean);
    //tri::UpdateColor<CMesh>::PerFaceConstant(tboolean,Color4b::LightGreen,true);
    /*UpdateQuality<CMesh>::VertexConstant(tboolean,-1.0);
    size_t index=0;
    for(size_t i=0;i<CurvesInBoolean.size();i++){
        for(size_t j=0;j<CurvesInBoolean[i].size();j++){
            CurvesInBoolean[i][j].VFlip()->Q()=float(index);
            index++;
        }
    }
    SmoothIntersection<CMesh>::Smooth(tboolean);
    UpdateQuality<CMesh>::VertexConstant(tboolean,-1.0);*/


    //collapseSmallTrianglesIntersectionCurve2(CurvesInBoolean);

    //splitIntersectionCurve();
    //return;

    //Transferring patches number from the trimeshes to the boolean mesh. This requires the Correspondence attribute correctly setted
    // It is mandatory had not made any split operation    
    updateBooleanIMarkSelectedPartitions(sc1,sc2);
    conditioningIntersection();

    // paint patches
    size_t npatches=sc1.size()+sc2.size()+1;
    size_t r=0;
    for (auto it=sc1.begin(); it != sc1.end(); ++it,r++){
        vector<size_t> partitionfaces;
        for(CMesh::FaceIterator fi = tboolean.face.begin(); fi != tboolean.face.end(); ++fi)
        {
            size_t findex=tri::Index(tboolean,*fi);
            if(!fi->IsD() && (hMeshNumber[fi]==1) && (fi->IMark()==*it)){
                //cout<<"hcor "<<hCorrespondence[fi]<<endl;
                if(hCorrespondence[fi]>-1)
                    tmesh1.face[hCorrespondence[fi]].C()=vcg::Color4b::Scatter(npatches,r);
                fi->C()=vcg::Color4b::Scatter(npatches,r);
            }
        }
    }
    for (auto it=sc2.begin(); it != sc2.end(); ++it,r++){
        vector<size_t> partitionfaces;
        for(CMesh::FaceIterator fi = tboolean.face.begin(); fi != tboolean.face.end(); ++fi)
        {
            size_t findex=tri::Index(tboolean,*fi);
            if(!fi->IsD() && (hMeshNumber[fi]==2) && (fi->IMark()==*it)){
                if(hCorrespondence[fi]>-1)
                    tmesh2.face[hCorrespondence[fi]].C()=vcg::Color4b::Scatter(npatches,r);

                fi->C()=vcg::Color4b::Scatter(npatches,r);
            }
        }
    }

    // Storing the faces belonging to the patches. Brute search for while
    std::unordered_map<int,vector<size_t>> partition2FacesMesh1; //store the pair (numberPartition,faces)
    for (auto it=sc1.begin(); it != sc1.end(); ++it){
        vector<size_t> partitionfaces;
        for(CMesh::FaceIterator fi = tboolean.face.begin(); fi != tboolean.face.end(); ++fi)
        {
            size_t findex=tri::Index(tboolean,*fi);
            if(!fi->IsD() && (hMeshNumber[fi]==1) && (fi->IMark()==*it)){
                partitionfaces.push_back(findex);
            }
        }
        partition2FacesMesh1.emplace(*it,partitionfaces);
    }
    std::unordered_map<int,vector<size_t>> partition2FacesMesh2;
    for (auto it=sc2.begin(); it != sc2.end(); ++it){
        vector<size_t> partitionfaces;

        for(CMesh::FaceIterator fi = tboolean.face.begin(); fi != tboolean.face.end(); ++fi)
        {
            size_t findex=tri::Index(tboolean,*fi);
            if(!fi->IsD() && (hMeshNumber[fi]==2) && (fi->IMark()==*it)){
                partitionfaces.push_back(findex);
            }
        }
        partition2FacesMesh2.emplace(*it,partitionfaces);
    }


    /*UpdateQuality<CMesh>::VertexConstant(tboolean,-1.0);
    size_t index=0;
    for(size_t i=0;i<CurvesInBoolean.size();i++){
        for(size_t j=0;j<CurvesInBoolean[i].size();j++){
            CurvesInBoolean[i][j].VFlip()->Q()=float(index);
            CurvesInBoolean[i][j].F()->C()=Color4b::Blue;
            index++;
        }
    }
    SmoothIntersection<CMesh>::Smooth(tboolean);
    UpdateQuality<CMesh>::VertexConstant(tboolean,-1.0);*/


    for (auto pit=sc1.begin(); pit != sc1.end(); ++pit){
        string name="mesh1_patch_id="+std::to_string(*pit);
        tempora::Patch3D<CMesh,PMesh> patch1(tboolean,tmesh1,mesh1,partition2FacesMesh1[*pit],*pit,name,session);
        //fpatch.paint();
        //fpatch.exportProjectedTriangulation();
        CollectionPatches.push_back(patch1);
        cout<<"patch with "<<patch1.numberConnectedComponents<<" conected components "<<endl;
    }
    for (auto pit=sc2.begin(); pit != sc2.end(); ++pit){
        string name="mesh2_patch_id="+std::to_string(*pit);
        tempora::Patch3D<CMesh,PMesh> patch2(tboolean,tmesh2,mesh2,partition2FacesMesh2[*pit],*pit,name,session);
        CollectionPatches.push_back(patch2);
        cout<<"patch with "<<patch2.numberConnectedComponents<<" conected components "<<endl;
    }

    // Painting the patches
    for(size_t i=0;i<CollectionPatches.size();i++){
        Color4b temp=vcg::Color4b::Scatter(CollectionPatches.size()+1,i);
        CollectionPatches[i].setcolor(temp);
        CollectionPatches[i].paint();
    }

    // Doing all patches as even patches, i.e, we add a vertices if necessary
    if(!makePatchesEven(CollectionPatches))
        return;

    edges.clear();
    for(size_t i=0;i<CollectionPatches.size();i++){
        for(size_t j=0;j<CollectionPatches[i].components.size();j++){
           cout<<"components with "<<CollectionPatches[i].numberDivisions(j)<<" divisions "<<endl;
           for(size_t k=0;k<CollectionPatches[i].components[j].sidesT.size();k++){
               edges.push_back(Segment3f(CollectionPatches[i].components[j].sidesT[k].origin.VFlip()->P(),
                                         CollectionPatches[i].components[j].sidesT[k].origin.V()->P()));
               if(CollectionPatches[i].components[j].sidesT[k].ismarked){
                  set<size_t> index=CollectionPatches[i].components[j].sidesT[k].marked_index;
                  //cout<<"side "<<k<<" is marked"<<endl;
                  for(const size_t & r:index) {
                   cout<<"   faces "<<tri::Index(tboolean,CollectionPatches[i].components[j].sidesT[k].directedside[r].F())<<" "<<tri::Index(tboolean,CollectionPatches[i].components[j].sidesT[k].directedside[r].FFlip())<<" times "<<CollectionPatches[i].components[j].sidesT[k].quality[r][4]<<endl;
                  }
               }
           }
        }
    }
    //computing threshold to sample the intersection curves. Computing edges in the intersection curve to be splitted
    /*vector<CMesh::VertexPointer> fixed_indexes;
    vector<face::Pos<CMesh::FaceType>> edgestosplit;
    vector<int> timestosplit;
    CMesh::ScalarType thres=0;
    size_t count=0;
    for(size_t i=0;i<CollectionPatches.size();i++){
        for(size_t j=0;j<CollectionPatches[i].components.size();j++){
           thres+= CollectionPatches[i].getLengthOnlyMeshEdges(j);
           count++;
           for(size_t k=0;k<CollectionPatches[i].components[j].sidesT.size();k++){
               if(CollectionPatches[i].components[j].sidesT[k].iscurve){// if it is part of the intersection curve
                    fixed_indexes.push_back(CollectionPatches[i].components[j].sidesT[k].origin.VFlip());
                    fixed_indexes.push_back(CollectionPatches[i].components[j].sidesT[k].end.V());
                    set<size_t> &tosplit=CollectionPatches[i].components[j].sidesT[k].marked_index;
                    for (const size_t & ind:tosplit) {
                        edgestosplit.push_back(CollectionPatches[i].components[j].sidesT[k].directedside[ind]);
                        timestosplit.push_back(CollectionPatches[i].components[j].sidesT[k].quality[ind][4]);
                    }

               }
           }
        }
    }
    thres/=count;
    cout<<"trhesloh "<<thres<<endl;*/
    //intcurves.setThreshold(thres);
    //intcurves.processCurves(fixed_indexes,edgestosplit,timestosplit);
    //glWrapSIC.m=&intcurves.renderSmooth;

    for(size_t i=0;i<CollectionPatches.size();i++){
        CollectionPatches[i].setBoundaries2D();   // after have computed the relationships and corrected parity we compute the 2D polygons
    }

    //exportBooleanTagged();
    quadboolean.Clear();

    for(size_t i=0;i<CollectionPatches.size();i++){
        PMesh tempmesh0;
        CollectionPatches[i].processAsPatchesLSCM(tempmesh0);
        /*if(i==1){
            UpdateSelection<PMesh>::FaceAll(tempmesh0);
            cout<<"selected faces "<<UpdateSelection<PMesh>::FaceCount(tempmesh0)<<endl;
        }*/
        vcg::tri::Append<PMesh,PMesh>::Mesh(quadboolean,tempmesh0);
        vcg::PolygonalAlgorithm<PMesh>::UpdateFaceNormalByFitting(quadboolean);
        vcg::tri::UpdateNormal<PMesh>::PerVertexNormalized(quadboolean);
        vcg::tri::UpdateBounding<PMesh>::Box(quadboolean);
        //vcg::tri::UpdateTopology<PMesh>::FaceFace(quadboolean);
        //vcg::tri::UpdateFlags<PMesh>::VertexBorderFromFaceAdj(quadboolean);

    }

    int Mask=0;
    string auxorigin="../exported/"+session+"/result3D-original-before.obj";
    if(flagExportIntermediateMeshes)
        vcg::tri::io::ExporterOBJ<PMesh>::Save(quadboolean,auxorigin.c_str(), Mask);
    //cout<<"selected faces "<<UpdateSelection<PMesh>::FaceCount(quadboolean)<<endl;
    Clean<PMesh>::RemoveDuplicateVertex(quadboolean);
    tri::MeshAssert<PMesh>::OnlyQuadFace(quadboolean);
        cout<<"is quad before"<<endl;
    //double len=utility::minimum_length_quadmesh<PMesh>(quadboolean);
    //tri::Clean<PMesh>::MergeCloseVertex(quadboolean,len/5);
    //cout<<"selected faces "<<UpdateSelection<PMesh>::FaceCount(quadboolean)<<endl;
    //UpdateTopology<PMesh>::FaceFace(quadboolean);
    tri::MeshAssert<PMesh>::OnlyQuadFace(quadboolean);
        cout<<"is quad after"<<endl;
    //Clean<PMesh>::MergeCloseVertex(quadboolean,0.0001);
    //tri::UpdateTopology<PMesh>::FaceFace(quadboolean);
    //tri::UpdateTopology<PMesh>::EdgeEdge(quadboolean);
    vcg::PolygonalAlgorithm<PMesh>::UpdateFaceNormals(quadboolean);
    if(flagExportIntermediateMeshes){
        string auxorigin1="../exported/"+session+"/result3D-original.obj";
        vcg::tri::io::ExporterOBJ<PMesh>::Save(quadboolean,auxorigin1.c_str(), vcg::tri::io::Mask::IOM_BITPOLYGONAL);
    }
    //vcg::PolygonalAlgorithm<PMesh>::UpdateFaceNormalByFitting(quadboolean);
    vcg::tri::UpdateNormal<PMesh>::PerVertexNormalized(quadboolean);
    vcg::tri::UpdateBounding<PMesh>::Box(quadboolean);
    vcg::tri::UpdateTopology<PMesh>::FaceFace(quadboolean);
    vcg::tri::UpdateFlags<PMesh>::VertexBorderFromFaceAdj(quadboolean);
    OrientFaces<PMesh>::AutoOrientFaces(quadboolean);
    if(flagExportIntermediateMeshes)
        vcg::tri::io::ExporterOBJ<PMesh>::Save(quadboolean,auxorigin.c_str(), Mask);
    CMesh temptri;
    vcg::PolygonalAlgorithm<PMesh>::TriangulateToTriMesh<CMesh>(quadboolean,temptri);
    vcg::tri::UpdateBounding<CMesh>::Box(temptri);
    vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(temptri);
    vcg::tri::UpdateTopology<CMesh>::FaceFace(temptri);
    vcg::tri::UpdateFlags<CMesh>::FaceBorderFromFF(temptri);

    vector<pair<size_t,PMesh::CoordType>> intcurveQuad=getIntersectionCurveQuad();
    for(PMesh::FaceIterator qfi=quadboolean.face.begin();qfi!=quadboolean.face.end();qfi++){
        if(int(qfi->Q())==1){
            qfi->SetS();
        }
        else{
            qfi->ClearS();
        }
    }
    vcg::PolygonalAlgorithm<PMesh>::LaplacianReproject(quadboolean,temptri,50,0.5,0.5,true);
    /*for(size_t k=0;k<intcurveQuad.size();k++){
        quadboolean.vert[intcurveQuad[k].first].P()=intcurveQuad[k].second;
    }*/
    for(PMesh::FaceIterator qfi=quadboolean.face.begin();qfi!=quadboolean.face.end();qfi++){
        if(int(qfi->Q())==2){
            qfi->SetS();
        }
        else{
            qfi->ClearS();
        }
    }
    vcg::PolygonalAlgorithm<PMesh>::LaplacianReproject(quadboolean,temptri,50,0.5,0.5,true);
    /*for(size_t k=0;k<intcurveQuad.size();k++){
        quadboolean.vert[intcurveQuad[k].first].P()=intcurveQuad[k].second;
    }*/
    //vcg::PolygonalAlgorithm<PMesh>::SmoothReprojectPCA(quadboolean);
    OrientFaces<PMesh>::AutoOrientFaces(quadboolean);
    UpdateFlags<PMesh>::FaceClearS(quadboolean);
    //vcg::PolygonalAlgorithm<PMesh>::LaplacianReproject(quadboolean,tmesh1);
    //vcg::PolygonalAlgorithm<PMesh>::Laplacian(quadboolean);
    //vcg::PolygonalAlgorithm<PMesh>::SmoothPCA(quadboolean);
    //bool isoriented,orientable;
    tri::UpdateTopology<PMesh>::FaceFace(quadboolean);
    //tri::UpdateTopology<PMesh>::TestFaceFace(quadboolean);
    //vcg::PolygonalAlgorithm<PMesh>::UpdateFaceNormals(quadboolean);
    //Clean<PMesh>::OrientCoherentlyMesh(quadboolean,isoriented,orientable);

    //vcg::PolygonalAlgorithm<PMesh>::LaplacianReprojectBorder(quadboolean);
    Mask=0;
    string aux="../exported/"+session+"/result3D.obj";
    vcg::tri::io::ExporterOBJ<PMesh>::Save(quadboolean,aux.c_str(), vcg::tri::io::Mask::IOM_BITPOLYGONAL);
    flagQuadBoolean=true;

}
void meshes::testPatches(){


}
//simulate the boolean operation between two meshes when we draw a curve over the tmesh1

bool meshes::makePatchesEven(vector<tempora::Patch3D<CMesh,PMesh>> & listpatches){
    // Populating  quality of all connected components in sequence
    int count=0;
    vector<pair<size_t,size_t>> count_to_patch; // for each integer it give us the corresponding patch and component index
    for(size_t i=0;i<listpatches.size();i++){
        //for each connected component
        size_t ncc=listpatches[i].components.size();
        for(size_t c=0;c<ncc;c++){
            listpatches[i].setQuality(c,count);
            count_to_patch.push_back(make_pair(i,c));
            count++;
        }
    }

    // Populating relations between sides

    // clearing flags
    for(size_t i=0;i<listpatches.size();i++){
        //for each connected component
        size_t ncc=listpatches[i].components.size();
        for(size_t c=0;c<ncc;c++){
            size_t ns=listpatches[i].components[c].sidesT.size();
            //for each side
            for(size_t s=0;s<ns;s++){
                auto side=listpatches[i].components[c].sidesT[s];
                side.clearQuality();
            }
        }
    }
    for(size_t i=0;i<listpatches.size();i++){
        //for each connected component
        size_t ncc=listpatches[i].components.size();
        for(size_t c=0;c<ncc;c++){
            size_t ns=listpatches[i].components[c].sidesT.size();
            //for each side
            for(size_t s=0;s<ns;s++){
                int nts=listpatches[i].components[c].sidesT[s].size;
                // for each edge in this side
                for(size_t k=0;k<nts;k++){
                    vcg::face::Pos<CMesh::FaceType> edge=listpatches[i].components[c].sidesT[s].directedside[k];
                    edge.FlipF();
                    //for all other patches
                    for(size_t j=0;j<listpatches.size();j++){
                        if(j==i){ continue;}
                        else{
                            Eigen::Vector4i info;                 
                            if(listpatches[j].isPosOnThisPatch(edge,info,true)){
                                //cout<<"info sender "<<i<<" "<<c<<endl;
                                listpatches[i].components[c].sidesT[s].quality[k][0]=listpatches[j].components[info[0]].quality;
                                listpatches[i].components[c].sidesT[s].quality[k][1]=info[0];
                                listpatches[i].components[c].sidesT[s].quality[k][2]=info[1];
                                listpatches[i].components[c].sidesT[s].quality[k][3]=info[2];
                                listpatches[i].components[c].sidesT[s].quality[k][4]=0; // number of divisions of this edge
                                //cout<<"info enter "<<j<<" "<<c<<endl;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    // Making the set of odd patches
    set<int> odd_set;
    for(size_t i=0;i<listpatches.size();i++){
        //for each connected component
        size_t ncc=listpatches[i].components.size();
        for(size_t c=0;c<ncc;c++){
            if((listpatches[i].numberDivisions(c)%2)!=0){
                odd_set.insert(listpatches[i].components[c].quality);
            }
        }
    }
    cout<<"number odd "<<odd_set.size()<<endl;
    if(odd_set.size()%2!=0){
        cout<<"odd number of odd polygons"<<endl;
        cout<<odd_set.size()<<endl;
        return false;
    }

    // Making a graph
    std::vector<int> adjmatrix[count_to_patch.size()];
    for(size_t i=0;i<count_to_patch.size();i++){
       auto cc=listpatches[count_to_patch[i].first].components[count_to_patch[i].second]; // component
       for(size_t j=0;j<cc.sidesT.size();j++){
           if(cc.sidesT[j].quality[0][0]!=-1){
              int correspond=cc.sidesT[j].quality[0][0];
              unwgraph::add_edge(adjmatrix,i,correspond);
           }
       }
    }
    cout<<"number nodes"<<count_to_patch.size()<<endl;
    for(int i=0;i<count_to_patch.size();i++){
        set<int> neighbors(adjmatrix[i].begin(),adjmatrix[i].end());
        std::vector<int> simplified(neighbors.begin(),neighbors.end());
        adjmatrix[i]=simplified;
        //Eigen::VectorXi aux=Eigen::VectorXi::Map(simplified.data(),simplified.size());
        //cout<<aux.transpose()<<endl;
    }
    cout<<"before"<<endl;
    for(const int &ind: odd_set)
        std::cout<<ind<<" ";
    cout<<endl;

    //Discarding easy neighbors odd patches

        set<int>::iterator i=odd_set.begin();
        while(i!=odd_set.end() && !odd_set.empty()){

                vector<int> &neighbors=adjmatrix[*i];
                for(size_t j=0;j<neighbors.size();j++){
                    if(neighbors[j]==(*i)) continue;
                    set<int>::iterator itj=odd_set.find(neighbors[j]);
                    if(itj!=odd_set.end()){

                        // marking sides as solved
                        size_t inlist_i=count_to_patch[size_t(*i)].first;
                        size_t incomponent_i=count_to_patch[size_t(*i)].second;
                        size_t inlist_j=count_to_patch[size_t(neighbors[j])].first;
                        size_t incomponent_j=count_to_patch[size_t(neighbors[j])].second;
                        vector<size_t> commonedges;
                        int sidei=listpatches[inlist_i].getCommonPosByQuality(incomponent_i,neighbors[j],commonedges);
                        assert(sidei!=-1);
                        //assert(listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].quality[0][0]==neighbors[j]);

                        size_t insidesT_j=listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].quality[commonedges[0]][2];

                        listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].ismarked=true;
                        listpatches[inlist_j].components[incomponent_j].sidesT[insidesT_j].ismarked=true;

                        // updating pos edges  to be subdivided with the correct number of times

                        CMesh::ScalarType biggeredgelen=0;
                        size_t tobemarked=0;
                        for(size_t s=0;s<commonedges.size();s++){
                            CMesh::CoordType p1=listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].directedside[commonedges[s]].V()->P();
                            CMesh::CoordType p2=listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].directedside[commonedges[s]].VFlip()->P();
                            CMesh::ScalarType thislen=vcg::Distance(p1,p2);
                            // if this side has no subdivision we prefer more than others
                            if(listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].quality[commonedges[s]][4]==0){
                                thislen*=2;
                            }
                            if(thislen>biggeredgelen){
                                biggeredgelen=thislen;
                                tobemarked=commonedges[s];
                            }
                        }

                        size_t indirectedside_j=listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].quality[tobemarked][3]; // corresponding tobemarked into the other patch

                        assert(listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].quality[tobemarked][0]==neighbors[j]);
                        assert(listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].quality[tobemarked][1]==int(incomponent_j));
                        assert(listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].quality[tobemarked][2]==int(insidesT_j));

                        listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].marked_index.insert(tobemarked);
                        listpatches[inlist_j].components[incomponent_j].sidesT[insidesT_j].marked_index.insert(indirectedside_j);
                        listpatches[inlist_i].components[incomponent_i].sidesT[size_t(sidei)].quality[tobemarked][4]+=1;   // subdivide one more time
                        listpatches[inlist_j].components[incomponent_j].sidesT[insidesT_j].quality[indirectedside_j][4]+=1;

                        //erasing
                        odd_set.erase(i);
                        odd_set.erase(itj);
                        i=odd_set.begin();
                        break;
                    }
                }

            cout<<"all nodes"<<endl;
            for(const int &ind: odd_set)
                std::cout<<ind<<" ";
            cout<<endl;
            if(!odd_set.empty())
                i++;
        }


    cout<<"after"<<endl;
    for(const int &ind: odd_set)
        std::cout<<ind<<" ";
    cout<<endl;
    assert(odd_set.size()%2==0);

    // Removing hard odd patches
    i=odd_set.begin();
    int v=count_to_patch.size();
    while(!odd_set.empty()){
        size_t len=count_to_patch.size();
        std::vector<int> path;
        int node=-1;
        // checking the shortes path in all possible paths
        for(const int &j : odd_set){
            if(j==(*i)) continue;
            int pred[v], dist[v];
            if (unwgraph::BFS(adjmatrix, (*i), j, v, pred, dist) == true){ // if there is a path between *i and j
                if(dist[j]<int(len)){
                    len=size_t(dist[j]);
                    node=j;
                    path.clear();
                    //storing path
                    int crawl = j;
                    path.push_back(crawl);
                    while (pred[crawl] != -1) {
                        path.push_back(pred[crawl]);
                        crawl = pred[crawl];
                    }
                }
            }
        }
        assert(node!=-1);

        // printing path from source to destination
        std::cout << "\nPath is::\n"<<endl;
        for (int s = path.size() - 1; s >= 0; s--)
            std::cout << path[s] << " ";

        // Resolving parity along the path of patches. Here we find a path over the sides of the patches
       for(size_t r=0;r<path.size()-1;r++){
            size_t next=r; next++;
            size_t inlist_r=count_to_patch[size_t(path[r])].first;
            size_t incomponent_r=count_to_patch[size_t(path[r])].second;
            size_t inlist_n=count_to_patch[size_t(path[next])].first;
            size_t incomponent_n=count_to_patch[size_t(path[next])].second;

            // we find the common side between this patches
            vector<size_t> commonedges;
            size_t side_s=0;
            size_t ns=listpatches[inlist_r].components[incomponent_r].sidesT.size();
            //for each side
            for(size_t s=0;s<ns;s++){
                size_t nts=size_t(listpatches[inlist_r].components[incomponent_r].sidesT[s].size);
                // for each edge in this side
                for(size_t k=0;k<nts;k++){
                    if(listpatches[inlist_r].components[incomponent_r].sidesT[s].quality[k][0]==path[next])
                          commonedges.push_back(k);
                }
                if(!commonedges.empty()){
                    side_s=s;
                    break;
                }
            }
            // ranking edges
            CMesh::ScalarType biggeredgelen=0;
            size_t tobemarked=0;
            for(size_t s=0;s<commonedges.size();s++){
                CMesh::CoordType p1=listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].directedside[commonedges[s]].V()->P();
                CMesh::CoordType p2=listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].directedside[commonedges[s]].VFlip()->P();
                CMesh::ScalarType thislen=vcg::Distance(p1,p2);
                // if this side has no subdivision we prefer more than others
                if(listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].quality[commonedges[s]][4]==0){
                    thislen*=2;
                }
                if(thislen>biggeredgelen){
                    biggeredgelen=thislen;
                    tobemarked=commonedges[s];
                }
            }

            size_t insidesT_n=listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].quality[tobemarked][2];

            listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].ismarked=true;
            listpatches[inlist_n].components[incomponent_n].sidesT[insidesT_n].ismarked=true;

            size_t indirectedside_n=listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].quality[tobemarked][3]; // corresponding tobemarked into the other patch

            assert(listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].quality[tobemarked][0]==path[next]);
            assert(listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].quality[tobemarked][1]==int(incomponent_n));
            assert(listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].quality[tobemarked][2]==int(insidesT_n));

            listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].marked_index.insert(tobemarked);
            listpatches[inlist_n].components[incomponent_n].sidesT[insidesT_n].marked_index.insert(indirectedside_n);
            listpatches[inlist_r].components[incomponent_r].sidesT[size_t(side_s)].quality[tobemarked][4]+=1;   // subdivide one more time
            listpatches[inlist_n].components[incomponent_n].sidesT[insidesT_n].quality[indirectedside_n][4]+=1;

       }
        odd_set.erase(i);
        odd_set.erase(node);
        if(!odd_set.empty())
            i=odd_set.begin();
    }

    return true;

}
bool meshes::makePatchesEven2(vector<tempora::Patch3D<CMesh,PMesh>> & listpatches){
    // Populating  quality of all connected components in sequence
    int count=0;
    vector<pair<size_t,size_t>> count_to_patch;
    for(size_t i=0;i<listpatches.size();i++){
        //for each connected component
        size_t ncc=listpatches[i].components.size();
        for(size_t c=0;c<ncc;c++){
            listpatches[i].setQuality(c,count);
            count_to_patch.push_back(make_pair(i,c));
            count++;
        }
    }

    // Populating relations between sides
    for(size_t i=0;i<listpatches.size();i++){
        //for each connected component
        size_t ncc=listpatches[i].components.size();
        for(size_t c=0;c<ncc;c++){
            size_t ns=listpatches[i].components[c].sidesT.size();
            //for each side
            for(size_t s=0;s<ns;s++){
                auto side=listpatches[i].components[c].sidesT[s];
                side.clearQuality();
            }
        }
    }
    for(size_t i=0;i<listpatches.size();i++){
        //for each connected component
        size_t ncc=listpatches[i].components.size();
        for(size_t c=0;c<ncc;c++){
            size_t ns=listpatches[i].components[c].sidesT.size();
            //for each side
            for(size_t s=0;s<ns;s++){
                int nts=listpatches[i].components[c].sidesT[s].size;
                // for each edge in this side
                for(size_t k=0;k<nts;k++){
                    vcg::face::Pos<CMesh::FaceType> edge=listpatches[i].components[c].sidesT[s].directedside[k];
                    edge.FlipF();
                    //for all other patches
                    for(size_t j=0;j<listpatches.size();j++){
                        if(j==i){ continue;}
                        else{
                            Eigen::Vector4i info;
                            if(listpatches[j].isPosOnThisPatch(edge,info,true)){
                                //cout<<"info sender "<<i<<" "<<c<<endl;
                                listpatches[i].components[c].sidesT[s].quality[k][0]=listpatches[j].components[info[0]].quality;
                                listpatches[i].components[c].sidesT[s].quality[k][1]=info[0];
                                listpatches[i].components[c].sidesT[s].quality[k][2]=info[1];
                                listpatches[i].components[c].sidesT[s].quality[k][3]=info[2];
                                listpatches[i].components[c].sidesT[s].quality[k][4]=0; // number of divisions of this edge
                                //cout<<"info enter "<<j<<" "<<c<<endl;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    // Making the set of odd patches
    set<int> odd_set;
    for(size_t i=0;i<listpatches.size();i++){
        //for each connected component
        size_t ncc=listpatches[i].components.size();
        for(size_t c=0;c<ncc;c++){
            if((listpatches[i].numberDivisions(c)%2)!=0){
                odd_set.insert(listpatches[i].components[c].quality);  
            }
        }
    }
    cout<<"number odd "<<odd_set.size()<<endl;
    if(odd_set.size()%2!=0){
        cout<<"odd number of odd polygons"<<endl;
        cout<<odd_set.size()<<endl;
        return false;
    }
    // Deleting odd parity
    /*size_t idpatch;
    size_t idcount;
    int iterations=0;
    while(!odd_set.empty() && iterations<100){
      idcount=*odd_set.begin();
      idpatch=count_to_patch[idcount].first;
      size_t i=count_to_patch[idcount].second;
      tempora::Patch3D<CMesh,PMesh>& patch=listpatches[idpatch];
      bool finish=false;

      //for each side
      for(size_t j=0;j<patch.components[i].sidesT.size();j++) {
          tempora::Side<CMesh>& side=patch.components[i].sidesT[j];
          //first we try non market sides
          if(!side.ismarked && side.quality[0][0]!=-1){
              for(int k=0;k<side.directedside.size();k++){
                  vector<int> correspondence=side.quality[size_t(k)];
                  size_t correspond_count=size_t(correspondence[0]);
                  size_t correspond_patch=count_to_patch[correspond_count].first;
                  size_t correspond_i=count_to_patch[correspond_count].second;
                  bool correspondence_isodd=(odd_set.find(correspond_count)!=odd_set.end());
                  if(correspondence_isodd){
                      //Mark these two sides
                      side.ismarked=true;
                      //side.marked_index=size_t(k);
                      side.marked_index.push_back(size_t(k));
                      side.quality[k][4]=side.quality[k][4]+1;
                      listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].ismarked=true;
                      //listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].marked_index=correspondence[3];
                      listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].marked_index.push_back(correspondence[3]);
                      listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].quality[correspondence[3]][4]=side.quality[k][4];
                      // deleting these patches from the odd components list
                      //cout<<"set with "<<odd_set.size()<<" elements"<<endl;
                      set<int>::iterator it=odd_set.find(idcount);
                      odd_set.erase(it);
                      it=odd_set.find(correspond_count);
                      odd_set.erase(it);

                      finish=true;
                  }
                  if(finish)
                      break;
              }
          }
          if(finish)
              break;
      }
      if(!finish){
          //we try marked sides
          for(size_t j=0;j<patch.components[i].sidesT.size();j++) {
              tempora::Side<CMesh>& side=patch.components[i].sidesT[j];              
              if(side.ismarked && side.quality[0][0]!=-1){
                  for(int k=0;k<side.directedside.size();k++){
                      vector<int> correspondence=side.quality[size_t(k)];
                      size_t correspond_count=size_t(correspondence[0]);
                      size_t correspond_patch=count_to_patch[correspond_count].first;
                      size_t correspond_i=count_to_patch[correspond_count].second;
                      bool correspondence_isodd=(odd_set.find(correspond_count)!=odd_set.end());
                      if(correspondence_isodd){
                          //Mark these two sides
                          side.ismarked=true;
                          //side.marked_index=size_t(k);
                          side.marked_index.push_back(size_t(k));
                          side.quality[k][4]=side.quality[k][4]+1;
                          listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].ismarked=true;
                          //listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].marked_index=correspondence[3];
                          listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].marked_index.push_back(correspondence[3]);
                          listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].quality[correspondence[3]][4]=side.quality[k][4];
                          // deleting these patches from the odd components list
                          //cout<<"set with "<<odd_set.size()<<" elements"<<endl;
                          set<int>::iterator it=odd_set.find(idcount);
                          odd_set.erase(it);
                          it=odd_set.find(correspond_count);
                          odd_set.erase(it);

                          finish=true;
                      }
                      if(finish)
                          break;
                  }
              }
              if(finish)
                  break;
          }
          //If not exists an adjacent odd patch we attempt to transfer the parity to the first neighbor
          for(size_t j=0;j<patch.components[i].sidesT.size();j++) {
              tempora::Side<CMesh>& side=patch.components[i].sidesT[j];
              if(side.quality[0][0]!=-1){
                  for(size_t k=0;k<side.directedside.size();k++){
                          vector<int> correspondence=side.quality[k];
                          size_t correspond_count=size_t(correspondence[0]);
                          size_t correspond_patch=count_to_patch[correspond_count].first;
                          size_t correspond_i=count_to_patch[correspond_count].second;

                          //Mark these two sides
                          side.ismarked=true;
                          side.marked_index.push_back(k);
                          side.quality[k][4]=side.quality[k][4]+1;
                          listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].ismarked=true;
                          listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].marked_index.push_back(correspondence[3]);
                          listpatches[correspond_patch].components[correspond_i].sidesT[correspondence[2]].quality[correspondence[3]][4]=side.quality[k][4];
                          // deleting this patch from the odd components list
                          set<int>::iterator it=odd_set.find(idcount);
                          odd_set.erase(it);

                          //Adding the adjacent component to the list of odd components
                          odd_set.insert(correspond_count);
                          finish=true;
                          break;
                  }
              }
              if(finish)
                  break;
          }
      }
      iterations++;
    }*/
    return true;
}

vector<face::Pos<CMesh::FaceType>> meshes::getIntersectionCurveFromPos(face::Pos<CMesh::FaceType> &startPos){
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    face::Pos<CMesh::FaceType> curPos=startPos;
    vector<face::Pos<CMesh::FaceType>> result;
    result.push_back(curPos);
    curPos.V()->SetV();
    //cout<<"gerorder"<<endl;
    //cout<<" "<<tri::Index(*my,curPos.V());
    int count=0;
    bool isIntersectionEdge=false;
    do
    { //next edge of the intersection curve
      do{
         curPos.NextE();
         isIntersectionEdge=((hMeshNumber[curPos.F()]==1) && (hMeshNumber[curPos.FFlip()]==2))
                                      || ((hMeshNumber[curPos.F()]==2) && (hMeshNumber[curPos.FFlip()]==1));
      }while(!isIntersectionEdge);
      curPos.FlipF();
      curPos.FlipV();
      //cout<<" "<<tri::Index(*my,curPos.V());
      //curPos.V()->C()=Color4b::Blue;
      curPos.V()->SetV();
      //curPos.V()->SetS();
      result.push_back(curPos);
      //utility::saveAsCurvePly(result, "parameterization/curvetemporal.ply");
      count++;
    }while(curPos.V()!=startPos.V() && count<100000);
    //cout<<endl;
    result.pop_back();
    return result;
}

void meshes::processBooleanMesh(Eigen::VectorXi &criteria){    
    //query if these attribute are present or not
    bool hasMeshNumber = vcg::tri::HasPerFaceAttribute(tboolean,"MeshNumber");
    bool hasTriMeshCorrespondence = vcg::tri::HasPerFaceAttribute(tboolean,"Correspondence");
    bool hasTriMeshVertexMap = vcg::tri::HasPerVertexAttribute(tboolean,"VertexMap");
    if(hasMeshNumber)
        vcg::tri::Allocator<CMesh>::DeletePerFaceAttribute(tboolean,"MeshNumber");
    if(hasTriMeshCorrespondence)
        vcg::tri::Allocator<CMesh>::DeletePerFaceAttribute(tboolean,"Correspondence");
    if(hasTriMeshVertexMap)
        vcg::tri::Allocator<CMesh>::DeletePerVertexAttribute(tboolean,"VertexMap");

    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    CMesh::PerFaceAttributeHandle<int> hCorrespondence = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("Correspondence"));
    CMesh::PerVertexAttributeHandle<int> hVertexMap = vcg::tri::Allocator<CMesh>:: GetPerVertexAttribute<int>(tboolean,std::string("VertexMap"));
    size_t numberfacesA=tmesh1.face.size();
    /* Populate the Meshnumber and Correspondence attributes.
     * The second give the correspondence between faces from the boolean mesh to the triangle mesh. It is a useful attribute for the first stages of the algorithm
     * Remember: if J(i) < FA.rows() then face i of the output mesh is a sub-triangle of the J(i)th face of A.
     * If J(i)>FA.rows() then face i of the output mesh is a sub-triangle of the J(i)-FA.rows()th face of B.
     */
    for(auto fi=tboolean.face.begin();fi!=tboolean.face.end();++fi)
        if(!(*fi).IsD()){
            int booleanindex=vcg::tri::Index(tboolean,&*fi);
            if(criteria(booleanindex)<int(numberfacesA)){
                hMeshNumber[fi]=1;
                hCorrespondence[fi]=criteria(booleanindex);
            }
            else{
                hMeshNumber[fi]=2;
                hCorrespondence[fi]=criteria(booleanindex)-int(numberfacesA);
            }
        }

    //Cleaning steps. The boolean mesh could have zero area faces and other bad features
    Clean<CMesh>::RemoveZeroAreaFace(tboolean);
    Clean<CMesh>::RemoveDuplicateVertex(tboolean);
    Clean<CMesh>::RemoveUnreferencedVertex(tboolean);
    Clean<CMesh>::MergeCloseVertex(tboolean,0.000001);
    Clean<CMesh>::RemoveUnreferencedVertex(tboolean);
    Allocator<CMesh>::CompactEveryVector(tboolean);
    UpdateTopology<CMesh>::FaceFace(tboolean);
    UpdateTopology<CMesh>::VertexFace(tboolean);
    UpdateNormal<CMesh>::PerFaceNormalized(tboolean);
    UpdateNormal<CMesh>::PerVertexFromCurrentFaceNormal(tboolean);

    //Populate the VertexMap attribute that give us a correspondence between vertices from the boolean to the triangle mesh
    UpdateFlags<CMesh>::VertexClearV(tboolean);
    for(auto fi=tboolean.face.begin();fi!=tboolean.face.end();++fi){
        if(!fi->IsD()){
            int faceI=hCorrespondence[fi];
            int meshI=hMeshNumber[fi];
            CMesh::FacePointer fp;
            if(meshI==1)
                fp=&(tmesh1.face[size_t(faceI)]);
            else
                fp=&(tmesh2.face[size_t(faceI)]);
            for(int j=0;j<3;j++){
                if(!fi->V(j)->IsV()){
                    for(int k=0;k<3;k++){
                        double diff=vcg::Distance(fi->V(j)->P(),fp->V(k)->P());
                        if(diff<1e-6){
                            if(meshI==1)
                                hVertexMap[fi->V(j)]=tri::Index(tmesh1,fp->V(k));
                            else
                                hVertexMap[fi->V(j)]=tri::Index(tmesh2,fp->V(k));
                            fi->V(j)->SetV();
                            break;
                        }
                    }
                    if(!fi->V(j)->IsV()){
                        hVertexMap[fi->V(j)]=-1;
                        fi->V(j)->SetV();
                     }
                }
            }
        }
    }
}
void meshes::collapseSmallTrianglesIntersectionCurve1(vector<vector<face::Pos<CMesh::FaceType>>>& edges){
    vcg::tri::UpdateSelection<CMesh>::FaceClear(tboolean);
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    CMesh::PerFaceAttributeHandle<int> hCorrespondence = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("Correspondence"));
    CMesh::ScalarType thr=0;
    for(size_t i=0;i<edges.size();i++) {
       for(size_t j=0;j<edges[i].size();j++) {
           CMesh::FacePointer fi=edges[i][j].F();
           CMesh::FacePointer fio=edges[i][j].FFlip();
           size_t tfi=hCorrespondence[fi];
           size_t tfio=hCorrespondence[fio];
           PMesh::FacePointer qfi,qfio;
           fi->SetS();
           fio->SetS();
           if(hMeshNumber[fi]==1){
               qfi=&(mesh1.face[convertIndexT2Q(tfi)]);
           }
           else{
               qfi=&(mesh2.face[convertIndexT2Q(tfi)]);
           }
           if(hMeshNumber[fio]==1){
               qfio=&(mesh1.face[convertIndexT2Q(tfio)]);
           }
           else{
               qfio=&(mesh2.face[convertIndexT2Q(tfio)]);
           }

           int fisize=qfi->VN();
           int fiosize=qfio->VN();
           for(int k=0;k<qfi->VN();k++)
              thr+=vcg::Distance(qfi->V(k)->P(),qfi->V((k+1)%fisize)->P());
           for(int k=0;k<qfio->VN();k++)
              thr+=vcg::Distance(qfio->V(k)->P(),qfio->V((k+1)%fiosize)->P());
           thr/=(fisize+fiosize);
       }
    }
    vcg::tri::UpdateSelection<CMesh>::FaceDilate(tboolean);
    vcg::tri::UpdateSelection<CMesh>::FaceDilate(tboolean);

    IsotropicRemeshing<CMesh>::Params params0;
    params0.SetTargetLen(thr/2);
    params0.SetFeatureAngleDeg(90);
    params0.iter=3;
    params0.selectedOnly=true;
    params0.smoothFlag=true;
    params0.splitFlag=false;
    params0.projectFlag=true;
    IsotropicRemeshing<CMesh>::Do(tboolean,params0);
}
void meshes::collapseSmallTrianglesIntersectionCurve2(vector<vector<face::Pos<CMesh::FaceType>>>& edges){
    typedef typename tri::BasicVertexPair<CMesh::VertexType> vertexPairCollapse;
    // First pass: marking the edges that will be collapsed
    vector<vector<size_t>> tocollapse(edges.size());
    for(size_t i=0;i<edges.size();i++) {
       vector<size_t> index;
       for(size_t j=0;j<edges[i].size();j++) {
          face::Pos<CMesh::FaceType> copy=edges[i][j];          
          // opposite vertex in each adjacent triangle
          copy.FlipE();
          CMesh::CoordType r=copy.VFlip()->P();
          copy.FlipE(); // go back
          copy.FlipF(); //adjacent face
          copy.FlipE();
          CMesh::CoordType s=copy.VFlip()->P();
          copy.FlipE(); // go back
          copy.FlipF();
          // vector defining opposite angles
          CMesh::CoordType a=copy.V()->P()-r;
          CMesh::CoordType b=copy.VFlip()->P()-r;
          CMesh::CoordType c=copy.V()->P()-s;
          CMesh::CoordType d=copy.VFlip()->P()-s;
          // angles
          CMesh::ScalarType angle1=vcg::Angle(a,b);
          CMesh::ScalarType angle2=vcg::Angle(c,d);
          CMesh::ScalarType angle=angle1+angle2;
          if(angle<CMesh::ScalarType(M_PI/6)){
             index.push_back(j);
             copy.F()->C()=Color4b::Red;
             copy.FFlip()->C()=Color4b::Red;
          }
       }
       tocollapse[i]=index;
    }
    for(size_t i=0;i<tocollapse.size();i++){
        vector<size_t> indexes=tocollapse[i];
        size_t n=edges[i].size();
        size_t nn=indexes.size();
        // finding a index to start
        size_t start=tocollapse[i][0];
        for(size_t j=1;j<indexes.size();j++){
            size_t next=(start+1)%n;
            if(indexes[j]!=next){
               start=j;
               break;
            }
        }
        //collapsing iteratively
        for(size_t j=0;j<nn;j++){
            size_t id=(j+start)%nn;
            vertexPairCollapse vertexpair(edges[i][indexes[id]].VFlip(),edges[i][indexes[id]].V());
            CMesh::CoordType np=(edges[i][indexes[id]].VFlip()->P()+edges[i][indexes[id]].V()->P())/2;
            bool isFeasible=vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>>::LinkConditions(vertexpair);
            if(isFeasible){
                vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>> eco;
                eco.Do(tboolean,vertexpair,np);
            }
            else{
                cout<<"edge to collapse but not collapsed"<<endl;
            }
        }
    }
    // updating the edges list
    vcg::tri::UpdateTopology<CMesh>::FaceFace(tboolean);
    // First pass: Finding all edges of the intersection curve
    tri::UpdateFlags<CMesh>::FaceClearV(tboolean);
    tri::UpdateFlags<CMesh>::VertexClearV(tboolean);
    vector<vector<face::Pos<CMesh::FaceType>>> curves;
    vector<face::Pos<CMesh::FaceType>> PosSoup;
    face::Pos<CMesh::FaceType> p;
    p.SetNull();
    //bool found=false;
    // Remember: if J(i) < FA.rows() then face i of the output mesh is a sub-triangle of the J(i)th face of A.
    // If J(i)>FA.rows() then face i of the output mesh is a sub-triangle of the J(i)-FA.rows()th face of B.
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    for(auto fi=tboolean.face.begin();fi!=tboolean.face.end();++fi)
        if(!(*fi).IsD()){
            for(int i=0;i<3;++i){
                face::Pos<CMesh::FaceType> pi(&*fi,i,fi->V(i));
                CMesh::FacePointer facep=pi.F();
                CMesh::FacePointer faceopp=pi.FFlip();
                //if(faceopp->IsV()) continue;
                //bool isIntersectionEdge=(criteria(tri::Index(tboolean,facep))<numberfacesA) && (criteria(tri::Index(tboolean,faceopp))>=numberfacesA);
                //isIntersectionEdge=isIntersectionEdge || ((criteria(tri::Index(tboolean,facep))>=numberfacesA) && (criteria(tri::Index(tboolean,faceopp))<numberfacesA));
                bool isIntersectionEdge= hMeshNumber[facep]!=hMeshNumber[faceopp];
                //faceopp->SetV();
                bool wasvisited= pi.V()->IsV() && pi.VFlip()->IsV();
                if(isIntersectionEdge && !wasvisited){
                    pi.V()->SetV();
                    pi.VFlip()->SetV();
                    PosSoup.push_back(pi);
                    //p=pi;
                    //found=true;
                    //break;
                }
            }
            //fi->SetV();
            //if(found) break;
        }
    tri::UpdateFlags<CMesh>::VertexClearV(tboolean);
    tri::UpdateFlags<CMesh>::VertexClearS(tboolean);
    tri::UpdateFlags<CMesh>::FaceClearS(tboolean);
    CMesh::VertexPointer vpt;
    //Second pass: Obtaining a circulation of the intersection curves.
    int count=0;
    while(PosSoup.size()!=0 && count<1000){
        face::Pos<CMesh::FaceType> currPos=PosSoup.back();
        vpt=currPos.V();
        if(!vpt->IsV() && !currPos.VFlip()->IsV()){
            curves.push_back(getIntersectionCurveFromPos(currPos));
        }
        PosSoup.pop_back();
        count++;
    }
    for(int i=0;i<curves.size();i++){
        for(int j=0;j<curves[i].size();j++){
            curves[i][j].V()->SetS();
            curves[i][j].FFlip()->SetS();
            curves[i][j].F()->SetS();
        }
    }
    edges=curves;
}
void meshes::splitIntersectionCurve(){
    // storing meshnumber and hcorrespondence properties
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    CMesh::PerFaceAttributeHandle<int> hCorrespondence = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("Correspondence"));
    CMesh::PerVertexAttributeHandle<int> hVertexMap = vcg::tri::Allocator<CMesh>:: GetPerVertexAttribute<int>(tboolean,std::string("VertexMap"));

    // Transferring the mesh number and correspondence attributes. Storing original IMark
    unordered_map<CMesh::FacePointer,int> imark;
    for(CMesh::FaceIterator fi=tboolean.face.begin();fi!=tboolean.face.end();fi++){
        if(!fi->IsD()){
           fi->Q()=float(hMeshNumber[fi]);
           imark.emplace(&*fi,fi->IMark());
           //fi->IMark()=hCorrespondence[fi];
        }
    }

    // Getting number of new vertices
    vector<CMesh::ScalarType> edgelen;
    intcurves.getAverageEdgeLen(edgelen);
    CMesh::ScalarType min_edge_length=*std::min_element(edgelen.begin(),edgelen.end());
    CMesh::ScalarType max_edge_length=*std::max_element(edgelen.begin(),edgelen.end());
    if(double(min_edge_length)<0.000001)
        min_edge_length=max_edge_length;

    min_edge_length/=2;
    //int vertices_number=intcurves.getNumberEdges_over_threshold(min_edge_length);
    size_t nv_before=0;
    for(CMesh::VertexIterator vi=tboolean.vert.begin();vi!=tboolean.vert.end();vi++){
        if(!vi->IsD())
            nv_before++;
    }
    EdgeSplitFixedList_BooleanMesh ep(tboolean,min_edge_length);
    tri::MidPoint<CMesh> midFunctor(&tboolean);
    tri::RefineE(tboolean,midFunctor,ep,false);
    size_t nv_after=0;
    for(CMesh::VertexIterator vi=tboolean.vert.begin();vi!=tboolean.vert.end();vi++){
        if(!vi->IsD())
            nv_after++;
    }
    //assert(size_t(vertices_number)==(nv_after-nv_before));
    // recovering meshnumber and correspondence properties
    for(CMesh::FaceIterator fi=tboolean.face.begin();fi!=tboolean.face.end();fi++){
            hMeshNumber[fi]=int(fi->Q());
            //hCorrespondence[fi]=int(fi->IMark());
    }
    // setting vertexMap property for the new vertices
    size_t vertices_number=nv_after-nv_before;
    for(size_t z=0;z<vertices_number;z++){
            hVertexMap[nv_before+z]=-1;
    }
    //Allocator<CMesh>::CompactEveryVector(tboolean);
    /*for(CMesh::VertexIterator vii=tboolean.vert.begin();vii!=tboolean.vert.end();vii++){
        cout<<" % "<<hVertexMap[vii];
    }*/
    //painting for debugging
    /*tri::UpdateColor<CMesh>::PerFaceConstant(tboolean,Color4b::White);
    for(CMesh::FaceIterator fi=tboolean.face.begin();fi!=tboolean.face.end();fi++){
       if(!fi->IsD() && hMeshNumber[fi]==1){
            fi->C()=Color4b::Red;
       }
       else if(!fi->IsD() && hMeshNumber[fi]==2){
            fi->C()=Color4b::Blue;
       }
       else{
            fi->C()=Color4b::Green;
       }
       if(hMeshNumber[fi]!=1 && hMeshNumber[fi]!=2){
          cout<<"mesh number "<<hMeshNumber[fi]<<endl;
          cout<<"correspondence "<<hCorrespondence[fi]<<endl;
       }
    }*/

    // recovering original IMark
    /*for(auto & pairs: imark){
        if(!pairs.first->IsD()){
            pairs.first->C()=Color4b::Black;
            pairs.first->IMark()=pairs.second;
            //cout<<pairs.second;
        }
    }*/
    // recovering intersection curves as POS
    getIntersectionCurve();
}
// It requires the intersection curve already computed
void meshes::mergeEndpointsPatches(){
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    vector<CMesh::ScalarType> edgelen;
    intcurves.getAverageEdgeLen(edgelen);
    //computing extremes of patches
    vector<vector<face::Pos<CMesh::FaceType>>> & original_curves=intcurves.original_curves;
    vector<vector<size_t>> extremes;
    // for each intersection curve
    for(size_t i=0;i<original_curves.size();i++){
        set<size_t> endpoints_set;
        size_t n=original_curves[i].size();
        // testing patch number flag
        for(size_t j=0;j<n;j++){
            original_curves[i][j].VFlip()->ClearV();
            CMesh::FacePointer f1=(hMeshNumber[original_curves[i][j].F()]==1)?original_curves[i][j].F():original_curves[i][j].FFlip();
            CMesh::FacePointer f2=(hMeshNumber[original_curves[i][j].F()]==2)?original_curves[i][j].F():original_curves[i][j].FFlip();
            CMesh::FacePointer nf1=(hMeshNumber[original_curves[i][(j+1)%n].F()]==1)?original_curves[i][(j+1)%n].F():original_curves[i][(j+1)%n].FFlip();
            CMesh::FacePointer nf2=(hMeshNumber[original_curves[i][(j+1)%n].F()]==2)?original_curves[i][(j+1)%n].F():original_curves[i][(j+1)%n].FFlip();
            if((f1->IMark()!=nf1->IMark()) || (f2->IMark()!=nf2->IMark())){
                endpoints_set.insert(j);
            }
        }
        vector<size_t> endpoints(endpoints_set.begin(),endpoints_set.end());
        extremes.push_back(endpoints);
    }
    if(extremes.empty()) return;

    //cleaning wrong marked extremes, i.e , those thar belong to tha same ṕatch in the same quad mesh
    for(size_t i=0;i<extremes.size();i++){        
        size_t nn=extremes[i].size();
        if(nn==0) continue;
        size_t n=original_curves[i].size();
        set<size_t> endpoints_set;
        for(size_t j=0;j<nn;j++){
            CMesh::FacePointer bf1=(hMeshNumber[original_curves[i][extremes[i][j]].F()]==1)?original_curves[i][extremes[i][j]].F():original_curves[i][extremes[i][j]].FFlip();
            CMesh::FacePointer bf2=(hMeshNumber[original_curves[i][extremes[i][j]].F()]==2)?original_curves[i][extremes[i][j]].F():original_curves[i][extremes[i][j]].FFlip();
            CMesh::FacePointer f1=(hMeshNumber[original_curves[i][(extremes[i][j]+1)%n].F()]==1)?original_curves[i][(extremes[i][j]+1)%n].F():original_curves[i][(extremes[i][j]+1)%n].FFlip();
            CMesh::FacePointer f2=(hMeshNumber[original_curves[i][(extremes[i][j]+1)%n].F()]==2)?original_curves[i][(extremes[i][j]+1)%n].F():original_curves[i][(extremes[i][j]+1)%n].FFlip();
            int patchnumber1=f1->IMark();
            int patchnumber2=f2->IMark();
            CMesh::FacePointer nf1=(hMeshNumber[original_curves[i][(extremes[i][(j+1)%nn]+1)%n].F()]==1)?original_curves[i][(extremes[i][(j+1)%nn]+1)%n].F():original_curves[i][(extremes[i][(j+1)%nn]+1)%n].FFlip();
            CMesh::FacePointer nf2=(hMeshNumber[original_curves[i][(extremes[i][(j+1)%nn]+1)%n].F()]==2)?original_curves[i][(extremes[i][(j+1)%nn]+1)%n].F():original_curves[i][(extremes[i][(j+1)%nn]+1)%n].FFlip();
            bool isendpoint1=(bf1->IMark()!=patchnumber1);
            bool isendpoint2=(bf2->IMark()!=patchnumber2);
            bool isnextendpoint1=(nf1->IMark()!=patchnumber1) ;
            bool isnextendpoint2=(nf2->IMark()!=patchnumber2);
            if((isendpoint1!=isnextendpoint1) && (isendpoint2!=isnextendpoint2)){
                if((isendpoint1 && isnextendpoint2) || (isendpoint2 && isnextendpoint1)){
                    endpoints_set.insert(extremes[i][j]);
                    endpoints_set.insert(extremes[i][(j+1)%nn]);
                    original_curves[i][extremes[i][j]].V()->SetV();
                    original_curves[i][extremes[i][(j+1)%nn]].V()->SetV();
                    //original_curves[i][extremes[i][j]].V()->C()=Color4b::Red;
                    //original_curves[i][extremes[i][(j+1)%nn]].V()->C()=Color4b::Red;
                }
            }
            //original_curves[i][extremes[i][j]].V()->C()=Color4b::Red;
        }
        vector<size_t> endpoints(endpoints_set.begin(),endpoints_set.end());
        extremes[i]=endpoints;
    }
    // for each intersection curve
    for(size_t i=0;i<original_curves.size();i++){
        // computing edges sets
        size_t n=original_curves[i].size();
        vector<vector<face::Pos<CMesh::FaceType>>> blocks;
        if(extremes[i].empty()){
            blocks.push_back(original_curves[i]);
        }
        else{
            vector<face::Pos<CMesh::FaceType>> temp;
            face::Pos<CMesh::FaceType> start=original_curves[i][extremes[i][0]];
            int count=1;
            face::Pos<CMesh::FaceType> movel=original_curves[i][(extremes[i][0]+count)%n];
            while(movel!=start) {
                temp.push_back(movel);
                count++;
                movel=original_curves[i][(extremes[i][0]+count)%n];
                if(movel.V()->IsV()){
                    temp.push_back(movel);
                    blocks.push_back(temp);
                    temp.clear();
                    if((movel!=start)){
                        count++;
                        movel=original_curves[i][(extremes[i][0]+count)%n];
                    }
                }
                if(count>1000){
                    cout<<"while"<<endl;
                }
            }
        }
        // computing blocks to be collapsed
        vector<size_t> tobecollapse;
        size_t size_blocks=blocks.size();
        for(size_t j=0;j<size_blocks;j++){
            CMesh::ScalarType thislen=0;
            for(size_t k=0;k<blocks[j].size();k++){
                thislen+=vcg::Distance(blocks[j][k].VFlip()->P(),blocks[j][k].V()->P());
            }            
            if(thislen<3*edgelen[i]){
                tobecollapse.push_back(j);
            }
        }
        // collapsing
        for(size_t j=0;j<tobecollapse.size();j++){
            for(size_t k=0;k<blocks[tobecollapse[j]].size();k++){
                blocks[tobecollapse[j]][k].F()->C()=Color4b::Red;
                blocks[tobecollapse[j]][k].FFlip()->C()=Color4b::Red;
            }
            utility::collapseEdges<CMesh>(tboolean,blocks[tobecollapse[j]]);
            cout<<"collapsed close patch endpoints in the intersection curve"<<endl;
        }
    }
    UpdateTopology<CMesh>::FaceFace(tboolean);
}
void meshes::collapsingAfterSplit(vector<CMesh::ScalarType> & aedgelen){
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    vector<vector<face::Pos<CMesh::FaceType>>> & curves=intcurves.original_curves;
    UpdateTopology<CMesh>::VertexFace(tboolean);
    //computing average edge length of each intersection curve

    size_t nc=curves.size();

    for(size_t i=0;i<nc;i++){
        //For each intersection curve
        vector<size_t> extremes;
        set<size_t> endpoints_set;
        vector<vector<face::Pos<CMesh::FaceType>>> blocks;
        size_t n=curves[i].size();
        //computing blocks to remesh
        for(size_t j=0;j<n;j++)
            curves[i][j].V()->ClearV();

        for(size_t j=0;j<n;j++){
            CMesh::FacePointer f1=(hMeshNumber[curves[i][j].F()]==1)?curves[i][j].F():curves[i][j].FFlip();
            CMesh::FacePointer f2=(hMeshNumber[curves[i][j].F()]==2)?curves[i][j].F():curves[i][j].FFlip();
            CMesh::FacePointer nf1=(hMeshNumber[curves[i][(j+1)%n].F()]==1)?curves[i][(j+1)%n].F():curves[i][(j+1)%n].FFlip();
            CMesh::FacePointer nf2=(hMeshNumber[curves[i][(j+1)%n].F()]==2)?curves[i][(j+1)%n].F():curves[i][(j+1)%n].FFlip();
            //cout<<" f1 "<<f1->IMark()<<" "<<nf1->IMark()<<endl;
            //cout<<" f2 "<<f2->IMark()<<" "<<nf2->IMark()<<endl;
            if((f1->IMark()!=nf1->IMark()) || (f2->IMark()!=nf2->IMark())){
                //curves[i][j].V()->C()=Color4b::Red;
                curves[i][j].V()->SetV();
                endpoints_set.insert(j);
            }
        }
        vector<size_t> endpoints(endpoints_set.begin(),endpoints_set.end());
        if(!endpoints.empty())
            extremes=endpoints;
        else
            extremes.push_back(0);
        vector<face::Pos<CMesh::FaceType>> temp;
        for(size_t j=0;j<n;j++){
            face::Pos<CMesh::FaceType> pos=curves[i][(j+extremes[0])%n];
            temp.push_back(pos);
            if(pos.V()->IsV() || j==(n-1)){
                blocks.push_back(temp);
                temp.clear();
            }
        }
        assert(temp.empty());
        for(size_t j=0;j<blocks.size();j++){
            // for each block

            // computing edge length of the quad around these edges
            CMesh::ScalarType threshold=aedgelen[i];

            utility::remeshingCurveSegmentByCollapsing(tboolean,blocks[j],threshold);
        }
    }
    tri::Smooth<CMesh>::VertexCoordLaplacian(tboolean,1,false,true);
    UpdateTopology<CMesh>::FaceFace(tboolean);
    getIntersectionCurve();
}
// It requires the Correspondence attribute correctly setted
void meshes::updateBooleanIMarkSelectedPartitions(std::set<int> &sc1,std::set<int> &sc2){
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    CMesh::PerFaceAttributeHandle<int> hCorrespondence = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("Correspondence"));
    for(CMesh::FaceIterator fi = tboolean.face.begin(); fi != tboolean.face.end(); ++fi){
        if(!(*fi).IsD()){
            fi->IMark()=-1;
            int indexq=convertIndexT2Q(hCorrespondence[fi]);
            if(hMeshNumber[fi]==1){
                if(sc1.find(PTracer1.FacePatch[indexq])!=sc1.end()){
                   //fi->C()=mesh1.face[indexq].C();
                   fi->IMark()=PTracer1.FacePatch[indexq];
                }
            }
            else{
                if(sc2.find(PTracer2.FacePatch[indexq])!=sc2.end()){
                   //fi->C()=mesh2.face[indexq].C();
                   fi->IMark()=PTracer2.FacePatch[indexq];
                }
            }
        }
    }
}
vector<pair<size_t,PMesh::CoordType>> meshes::getIntersectionCurveQuad(){
    vector<pair<size_t,PMesh::CoordType>> result;
    set<size_t> indexes;
    for(auto fi=quadboolean.face.begin();fi!=quadboolean.face.end();++fi)
        if(!(*fi).IsD()){
            for(int i=0;i<fi->VN();++i){
                face::Pos<PMesh::FaceType> pi(&*fi,i,fi->V(i));
                PMesh::FacePointer facep=pi.F();
                PMesh::FacePointer faceopp=pi.FFlip();
                bool isIntersectionEdge= int(facep->Q())!= int(faceopp->Q());
                bool wasvisited= pi.V()->IsV() && pi.VFlip()->IsV();
                if(isIntersectionEdge && !wasvisited){
                    pi.V()->SetV();
                    pi.VFlip()->SetV();
                    indexes.insert(tri::Index(quadboolean,pi.V()));
                    indexes.insert(tri::Index(quadboolean,pi.VFlip()));
                }
            }
        }
    for(auto &ind: indexes){
        result.push_back(make_pair(ind,quadboolean.vert[ind].P()));
    }
    return result;
}
void meshes::selectFacesIntersection(bool selone,bool seltwo,int numberNeighbors){
    assert(numberNeighbors>0);
    vector<vector<face::Pos<CMesh::FaceType>>> &curves=intcurves.original_curves;
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    UpdateSelection<CMesh>::FaceClear(tboolean);
    for(size_t i=0;i<curves.size();i++){
      for(size_t j=0;j<curves[i].size();j++){
          if(selone){
              if(hMeshNumber[curves[i][j].F()]==1)
                  curves[i][j].F()->SetS();
              else
                  curves[i][j].FFlip()->SetS();

          }
          if(seltwo){
              if(hMeshNumber[curves[i][j].F()]==2)
                  curves[i][j].F()->SetS();
              else
                  curves[i][j].FFlip()->SetS();

          }
      }
    }
    for(size_t i=0;i<size_t(numberNeighbors-1);i++)
        UpdateSelection<CMesh>::FaceDilate(tboolean);
    if(!selone || !seltwo){
        for(CMesh::FaceIterator fi=tboolean.face.begin();fi!=tboolean.face.end();fi++){
            if(!fi->IsD() && fi->IsS()){
                if(!selone && hMeshNumber[fi]==1)
                    fi->ClearS();
                if(!seltwo && hMeshNumber[fi]==2)
                    fi->ClearS();
            }
        }
    }
}
void meshes::exportBooleanTagged(){
    UpdateQuality<CMesh>::VertexConstant(tboolean,-1.0);
    CMesh intersection;
    size_t index=0;
    for(size_t i=0;i<CurvesInBoolean.size();i++){
        for(size_t j=0;j<CurvesInBoolean[i].size();j++){
            CurvesInBoolean[i][j].VFlip()->Q()=float(index);
            index++;
        }
    }
    UpdateSelection<CMesh>::Clear(tboolean);
    CMesh::PerFaceAttributeHandle<int> hMeshNumber = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("MeshNumber"));
    CMesh::PerFaceAttributeHandle<int> hCorrespondence = vcg::tri::Allocator<CMesh>:: GetPerFaceAttribute<int>(tboolean,std::string("Correspondence"));
    for(CMesh::FaceIterator fi=tboolean.face.begin();fi!=tboolean.face.end();fi++){
        if(!fi->IsD()){
            size_t ftri=hCorrespondence[fi];
            size_t fquad=convertIndexT2Q(ftri);
            if(hMeshNumber[fi]==1){
                if(mesh1.face[fquad].IsS()){
                    fi->SetS();
                }
            }
            else{
                if(mesh2.face[fquad].IsS()){
                    fi->SetS();
                }
            }
        }
    }
    UpdateSelection<CMesh>::VertexFromFaceStrict(tboolean);
    UpdateColor<CMesh>::PerFaceConstant(tboolean,Color4b::Red,true);
    vcg::tri::Append<CMesh,CMesh>::Mesh(intersection,tboolean,true);
    tri::io::ExporterPLY<CMesh>::Save(intersection,"forNicoIntersection.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTQUALITY);
    tri::io::ExporterPLY<CMesh>::Save(tboolean,"forNico.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTQUALITY);
}





/* VISUALIZATION METHODS */




void meshes::reset(){
    flagTmesh1=false;
    flagTmesh2=false;
    flagEdgeMesh=false;
    flagIntersectedCurve=false;
    flagSmoothIntersectedCurve=false;
    flagTriangleBoolean=false;
    flagSeeSelectedQuad=false;
    flagQuadBoolean=false;
    numbersubdivisions=0;
    stopAnimation=false;
    tmesh1.Clear();
    tmesh2.Clear();
    mesh1.Clear();
    mesh2.Clear();
    tboolean.Clear();
    intcurves.reset();
    seepatches.Clear();
    quadboolean.Clear();
    borderpatches1.Clear();
    borderpatches2.Clear();
    CollectionPatches.clear();
    CurvesInBoolean.clear();
    edges.clear();
}
int meshes::loadmesh1(QString fileName){
    mesh1.Clear();
    tmesh1.Clear();
    int err = vcg::tri::io::Importer<PMesh>::Open(mesh1, fileName.toStdString().c_str());
    return err;
}
int meshes::loadmesh2(QString fileName){
    mesh2.Clear();
    tmesh2.Clear();
    int  err = vcg::tri::io::Importer<PMesh>::Open(mesh2, fileName.toStdString().c_str());    
    return err;
}
void meshes::initMesh(PMesh &mesh){
    tboolean.Clear();
    intcurves.reset();
    CollectionPatches.clear();
    vcg::tri::Allocator<CMesh>::DeletePerFaceAttribute(tboolean, "numberPartition");
    // update bounding box
    // update topology
    vcg::tri::UpdateTopology<PMesh>::FaceFace(mesh);
    vcg::tri::Clean<PMesh>::RemoveNonManifoldFace(mesh);
    vcg::tri::UpdateBounding<PMesh>::Box(mesh);
    tri::UpdatePosition<PMesh>::Translate(mesh,mesh.bbox.Center().Scale(-1,-1,-1));
    tri::UpdatePosition<PMesh>::Scale(mesh,1/mesh.bbox.Diag());
    vcg::tri::UpdateBounding<PMesh>::Box(mesh);
    vcg::tri::UpdateTopology<PMesh>::FaceFace(mesh);
    //assert(vcg::tri::Clean<PMesh>::CountNonManifoldEdgeFF(mesh1));

    // compute the half edges
    //vcg::tri::UpdateHalfEdges<PMesh>::FromIndexed(mesh);
    // check for consistency
    //assert(vcg::tri::UpdateHalfEdges<PMesh>::CheckConsistency(mesh));

    // update normals TODO: compute average normal in the polygon
    vcg::tri::UpdateNormal<PMesh>::PerFaceNormalized(mesh);
    // update flags
    vcg::tri::UpdateFlags<PMesh>::Clear(mesh);
    UpdateFlags<PMesh>::VertexBorderFromFaceAdj(mesh); // necessary for tracing patch partition
    // create a triangle mesh from a polygon mesh
    /*vcg::tri::PolygonSupport<CMesh,PMesh>::ImportFromPolyMesh(tmesh1,mesh1);
    // update bounding box
    vcg::tri::UpdateBounding<CMesh>::Box(tmesh1);
    // update Normals
    vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(tmesh1);
    vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(tmesh1);
    */
    /*
    int size =  mesh1.face.size();

        // add a face to each face with more than 3 vertices ( just one pass)

        for(int i = 0; i < size; ++i)
            if(!(mesh1.face[i].IsD()))
            if(mesh1.face[i].VN()>3){
                PMesh::HEdgePointer ef =  mesh1.face[i].FHp();
                PMesh::HEdgePointer ef1 = ef -> HNp();
                ef1 = ef1->HNp();
                vcg::tri::UpdateHalfEdges<PMesh>::AddHEdge(mesh1, ef, ef1 );
            }
        assert(vcg::tri::UpdateHalfEdges<PMesh>::CheckConsistency(mesh1));
        size =  mesh1.face.size();

        // remove an edge for each face

        for(int i = 0; i < size; ++i)
            if(!(mesh1.face[i].IsD() ))
            {
                PMesh::HEdgePointer ef =  mesh1.face[i].FHp();
                if( ef->HOp()->HFp() !=NULL){
                    vcg::tri::UpdateHalfEdges<PMesh>::RemoveHEdge(mesh1,ef);
                }
            }

        // check for consistency
        assert(vcg::tri::UpdateHalfEdges<PMesh>::CheckConsistency(mesh1));
    */
}
void meshes::initMesh1(){
    initMesh(mesh1);
    PTracer1.setModel(mesh1);
    PTracer1.MotorCycle=true;
    PTracer1.TracePartitions();
    PTracer1.SaveColoredMesh();
    tri::UpdateFlags<PMesh>::FaceClearS(mesh1);
    GRIDQ1.Set(mesh1.face.begin(),mesh1.face.end());
    tmesh1.Clear();
    //vcg::tri::UpdateColor<PMesh>::PerFaceConstant(mesh1, vcg::Color4b::LightGray);
}
void meshes::initMesh2(){
    initMesh(mesh2);
    PTracer2.setModel(mesh2);
    PTracer2.MotorCycle=true;
    PTracer2.TracePartitions();
    PTracer2.SaveColoredMesh();
    tri::UpdateFlags<PMesh>::FaceClearS(mesh2);
    GRIDQ2.Set(mesh2.face.begin(),mesh2.face.end());
    tmesh2.Clear();
    //vcg::tri::UpdateColor<PMesh>::PerFaceConstant(mesh2, vcg::Color4b::LightGray);
    /*// create a triangle mesh from a polygon mesh
    vcg::tri::PolygonSupport<CMesh,PMesh>::ImportFromPolyMesh(tmesh2,mesh2);
    // update bounding box
    vcg::tri::UpdateBounding<CMesh>::Box(tmesh2);
    // update Normals
    vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(tmesh2);
    vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(tmesh2);
    */


}
void meshes::changeSeeSelectedQuads(){
    flagSeeSelectedQuad=!flagSeeSelectedQuad;
}








/* NECESSARY METHODS FOR VISUALIZATION */


void meshes::triangulate(CMesh & tmesh,PMesh & mesh){
    tmesh.Clear();
    // create a triangle mesh from a polygon mesh
    vcg::tri::PolygonSupport<CMesh,PMesh>::ImportFromPolyMesh(tmesh,mesh);
    // update bounding box
    vcg::tri::UpdateBounding<CMesh>::Box(tmesh);
    // update Normals
    vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(tmesh);
    vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(tmesh);
}
void meshes::triangulate1(){
    triangulate(tmesh1,mesh1);
    vcg::tri::UpdateColor<CMesh>::PerFaceConstant(tmesh1, vcg::Color4b::LightGray);
}
void meshes::triangulate2(){
    triangulate(tmesh2,mesh2);
    vcg::tri::UpdateColor<CMesh>::PerFaceConstant(tmesh2, vcg::Color4b::LightGray);
}

void meshes::selectComponents(bool inner, PMesh &meshA,PMesh &meshB, GridStaticPtr<PMesh::FaceType, PMesh::ScalarType> &gm){
    UpdateSelection<PMesh>::FaceInvert(meshA);
    vector<PMesh::FacePointer> faces;
    // computing an edge flow from each component representative (missing to finish , for whle a interior face is tested)
    auto go_border = [](face::Pos<PMesh::FaceType> &p) // funtor that takes an edge and goes to the border edge in the edge flow of the first one
    {
        face::Pos<PMesh::FaceType> copy=p;
        bool isborder=false;
        int count=0;
        while(!isborder && count<10000){
            p.F()->C()=Color4b::DarkRed;
            p.FlipE();
            p.FlipV();
            p.FlipE();
            p.FlipF();
            isborder=(p.F()->IsS() && !p.FFlip()->IsS()) || (p.FFlip()->IsS() && !p.F()->IsS());
            count++;
            if(p==copy)
                break;
        }
    };


    for(PMesh::FaceIterator fi=meshA.face.begin();fi!=meshA.face.end();fi++){
        if(!fi->IsD() && fi->IsS()){
            faces.push_back(&*fi);
        }
    }
    if(faces.empty()){
        UpdateSelection<PMesh>::FaceInvert(meshA);
        return;
    }
    else{
        //computing connected components (of the selected faces )
        std::vector< std::pair<int,PMesh::FacePointer> > connectedComp;
        vector<vector<PMesh::FacePointer>> cc_faces;
        tri::UpdateFlags<PMesh>::FaceClearV(meshA);
        std::stack<PMesh::FacePointer> sf;
        PMesh::FacePointer fpt=faces[0];
        for(std::vector<PMesh::FacePointer>::iterator fi=faces.begin();fi!=faces.end();fi++)
        {
          if(!((*fi)->IsD()) && !(*fi)->IsV())
          {
            vector<PMesh::FacePointer> ccm;
            (*fi)->SetV();
            connectedComp.push_back(std::make_pair(0,*fi));
            sf.push(*fi);
            while (!sf.empty())
            {
              fpt=sf.top();
              ccm.push_back(fpt);
              ++connectedComp.back().first;
              sf.pop();
              for(int j=0; j<fpt->VN(); ++j)
              {
                bool isnotborderselected=fpt->IsS() && fpt->FFp(j)->IsS();
                if( !face::IsBorder(*fpt,j) && isnotborderselected )
                {
                  PMesh::FacePointer l = fpt->FFp(j);
                  if( !(*l).IsV() )
                  {
                    (*l).SetV();
                    sf.push(l);
                  }
                }
              }
            }
            cc_faces.push_back(ccm);
          }
        }

        //vector<vector<PMesh::FacePointer>> edge_flows;
        // computing connected components to be selected
        PMesh::CoordType test;
        vector<size_t> selected_components;
        for(size_t i=0;i<connectedComp.size();i++){
            PMesh::FacePointer fi=connectedComp[i].second;
            if(connectedComp[i].first==1){
                test= vcg::PolygonalAlgorithm<PMesh>::GetFaceGetBary(*fi);
            }
            else{
                for(size_t j=0;j<size_t(fi->VN());j++){
                    bool isborder=(fi->IsS() && !fi->FFp(int(j))->IsS()) || (!fi->IsS() && fi->FFp(int(j))->IsS());
                    if(!isborder){
                       test=(fi->V(j)->P()+fi->V(fi->Next(j))->P())/2;
                       break;
                    }
                }
            }
            if(inner){
                if(tri::Inside<GridStaticPtr<PMesh::FaceType, PMesh::ScalarType>,PMesh>::Is_Inside(meshB,gm,test)){
                        selected_components.push_back(i);
                        //fi->C()=Color4b::DarkRed;
                }
            }
            else{
                if(!tri::Inside<GridStaticPtr<PMesh::FaceType, PMesh::ScalarType>,PMesh>::Is_Inside(meshB,gm,test)){
                        selected_components.push_back(i);
                }
            }
            //vector<PMesh::FacePointer> ef;

            //going to and border edge
            //face::Pos<PMesh::FaceType> e(&*fi,0,fi->V(0));
            //go_border(e);

        }
        //UpdateColor<PMesh>::PerFaceConstant(meshA,Color4b::DarkRed);
        // come back to the original selection
        UpdateSelection<PMesh>::FaceInvert(meshA);

        // selecting the correct components
        for(size_t i=0;i<selected_components.size();i++){
            for(size_t k=0;k<cc_faces[selected_components[i]].size();k++)
                cc_faces[selected_components[i]][k]->SetS();
        }
    }
}
void meshes::selectComponents1(bool inner){
    GRIDQ2.Set(mesh2.face.begin(),mesh2.face.end());
    selectComponents(inner,mesh1,mesh2,GRIDQ2);
}
void meshes::selectComponents2(bool inner){
    GRIDQ1.Set(mesh1.face.begin(),mesh1.face.end());
    selectComponents(inner,mesh2,mesh1,GRIDQ1);
}
void meshes::visualizePatches(QuadMeshTracer<PMesh> &tracer, PMesh &mesh, CMesh &tmesh){
    std::set<std::pair<size_t,size_t>> tedges=tracer.getTraceEdges();
    std::set<std::pair<size_t,size_t>>::iterator iteE;
    for (iteE=tedges.begin();iteE!=tedges.end();iteE++){
        CMesh::CoordType pos0=mesh.vert[(*iteE).first].P();
        CMesh::CoordType pos1=mesh.vert[(*iteE).second].P();
        vcg::tri::Allocator<CMesh>::AddEdge(tmesh,pos0,pos1);
    }
}
void meshes::visualizePatches1(){
    borderpatches1.Clear();
    visualizePatches(PTracer1,mesh1,borderpatches1);

}
void meshes::visualizePatches2(){
    borderpatches2.Clear();
    visualizePatches(PTracer2,mesh2,borderpatches2);
}





/*UTILS */


int meshes::convertIndexT2Q(int tindex){
    if(tindex%2==0)
        return int(tindex/2);
    else
        return int((tindex-1)/2);
}
