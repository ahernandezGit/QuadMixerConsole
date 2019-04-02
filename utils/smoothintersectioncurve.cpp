#include "smoothintersectioncurve.h"


smoothIntersectionCurve::smoothIntersectionCurve()
{
  reset();
}
smoothIntersectionCurve::smoothIntersectionCurve(vector<vector<face::Pos<CMesh::FaceType>>> &icurves){
    initialize(icurves);
}
smoothIntersectionCurve::smoothIntersectionCurve(vector<vector<face::Pos<CMesh::FaceType>>> &icurves, CMesh::ScalarType ithresh){
    initialize(icurves,ithresh);
}
void smoothIntersectionCurve::initialize(vector<vector<face::Pos<CMesh::FaceType>>> &icurves, CMesh::ScalarType ithresh){
    renderOriginal.Clear();
    renderSmooth.Clear();
    original_curves=icurves;
    great_indexes.clear();
    blocks_indexes.clear();
    indexes_to_blocks.clear();
    blocks_interpolation.clear();
    blocks.clear();
    threshold=ithresh;
    createRenderOriginal();
    createRenderSmooth();
}
void smoothIntersectionCurve::initialize( vector<vector<face::Pos<CMesh::FaceType>>> &icurves){
    renderOriginal.Clear();
    renderSmooth.Clear();
    original_curves=icurves;
    checkOriginalCurve();
    indexes_to_blocks.clear();
    blocks_indexes.clear();
    blocks_interpolation.clear();
    great_indexes.clear();
    blocks.clear();
    threshold=0;
    createRenderOriginal();
}
void smoothIntersectionCurve::setThreshold( CMesh::ScalarType ithresh){
    threshold=ithresh;
    checkOriginalCurve();
}
void smoothIntersectionCurve::checkOriginalCurve(){
    for(size_t j=0;j<original_curves.size();j++){
        for(size_t i=0;i<original_curves[j].size()-1;i++)
            assert(original_curves[j][i+1].VFlip()==original_curves[j][i].V());
       assert(original_curves[j][0].VFlip()==original_curves[j].back().V());
    }
}
int smoothIntersectionCurve::getNumberEdges(){
    int result=0;
    for(size_t j=0;j<original_curves.size();j++){
        result+=original_curves[j].size()-1;
    }
    return result;
}
int smoothIntersectionCurve::getNumberEdges_over_threshold(CMesh::ScalarType ithresh){
    int result=0;
    for(size_t j=0;j<original_curves.size();j++){
        for(size_t k=0;k<original_curves[j].size();k++){
            CMesh::ScalarType distance=vcg::Distance(original_curves[j][k].V()->P(),original_curves[j][k].VFlip()->P());
            if(distance>ithresh)
                result++;
        }
    }
    return result;
}

void smoothIntersectionCurve::getAverageEdgeLen(vector<CMesh::ScalarType> &edgelen){
    edgelen.clear();
    for(size_t g=0;g<original_curves.size();g++){
        size_t n=original_curves[g].size();
        CMesh::ScalarType len=0;
        for(size_t i=0;i<n;i++){
            CMesh::FacePointer f=original_curves[g][i].F();
            CMesh::VertexPointer v=original_curves[g][i].V();
            int z=original_curves[g][i].E();
            if(!(f->cV(f->Prev(z))!=v && (f->cV(f->Next(z))==v || f->cV(z)==v))){
                f->C()=Color4b::Blue;
                return;
            }
            len+=vcg::Distance(original_curves[g][i].VFlip()->P(),original_curves[g][i].V()->P());
        }
        len/=n;
        edgelen.push_back(len);
    }
}
void smoothIntersectionCurve::createRenderOriginal(){
    renderOriginal.Clear();
    CMesh tempOriginal;
    for(size_t g=0;g<original_curves.size();g++){
        size_t n=original_curves[g].size();
        for(size_t i=0;i<n;i++){
           tri::Allocator<CMesh>::AddEdge(tempOriginal,original_curves[g][i].VFlip()->P(),original_curves[g][i].V()->P());
        }
        tri::Clean<CMesh>::RemoveDuplicateVertex(tempOriginal);
        tri::Clean<CMesh>::RemoveUnreferencedVertex(tempOriginal);
        tri::Allocator<CMesh>::CompactEveryVector(tempOriginal);
        tri::UpdateTopology<CMesh>::VertexEdge(tempOriginal);
        tri::UpdateTopology<CMesh>::EdgeEdge(tempOriginal);
        tri::UpdateColor<CMesh>::PerVertexConstant(tempOriginal,Color4b::Black);
        tri::io::ExporterPLY<CMesh>::Save(tempOriginal,("parameterization/renderOriginal"+std::to_string(original_curves[g].size())+std::to_string(g)+".ply").c_str(),tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
        Append<CMesh,CMesh>::Mesh(renderOriginal,tempOriginal);
    }
}
void smoothIntersectionCurve::reset(){
    renderOriginal.Clear();
    renderSmooth.Clear();
    original_curves.clear();
    great_indexes.clear();
    blocks_indexes.clear();
    blocks_interpolation.clear();
    blocks.clear();
    indexes_to_blocks.clear();
    threshold=0;
}
void smoothIntersectionCurve::createRenderSmooth(){
    for(size_t g=0;g<blocks_indexes.size();g++){
        CMesh tempSmooth;
        for(size_t i=0;i<blocks_indexes[g].size();i++){
            for(size_t j=0;j<blocks_indexes[g][i].size()-1;j++){
               tri::Allocator<CMesh>::AddEdge(renderSmooth,blocks_indexes[g][i][j],blocks_indexes[g][i][j+1]);
               tri::Allocator<CMesh>::AddEdge(tempSmooth,blocks[g][i][j],blocks[g][i][j+1]);
            }
        }
        tri::Clean<CMesh>::RemoveDuplicateVertex(tempSmooth);
        tri::Clean<CMesh>::RemoveUnreferencedVertex(tempSmooth);
        tri::Allocator<CMesh>::CompactEveryVector(tempSmooth);
        tri::UpdateTopology<CMesh>::VertexEdge(tempSmooth);
        tri::UpdateTopology<CMesh>::EdgeEdge(tempSmooth);
        tri::io::ExporterPLY<CMesh>::Save(tempSmooth,("parameterization/renderSmooth"+std::to_string(original_curves[g].size())+std::to_string(g)+".ply").c_str(),tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
    }
    tri::UpdateTopology<CMesh>::EdgeEdge(renderSmooth);
    tri::io::ExporterPLY<CMesh>::Save(renderSmooth,"parameterization/renderSmoothFinal.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
}
void smoothIntersectionCurve::createRenderSmooth2(){
    renderSmooth.Clear();

    for(size_t g=0;g<blocks.size();g++){
        CMesh tempSmooth;
        for(size_t i=0;i<blocks[g].size();i++){
            for(size_t j=0;j<blocks[g][i].size()-1;j++){
               tri::Allocator<CMesh>::AddEdge(tempSmooth,blocks[g][i][j],blocks[g][i][j+1]);
            }
        }
        tri::Clean<CMesh>::RemoveDuplicateVertex(tempSmooth);
        tri::Clean<CMesh>::RemoveUnreferencedVertex(tempSmooth);
        tri::Allocator<CMesh>::CompactEveryVector(tempSmooth);
        tri::UpdateTopology<CMesh>::VertexEdge(tempSmooth);
        tri::UpdateTopology<CMesh>::EdgeEdge(tempSmooth);
        //tri::UpdateColor<CMesh>::PerVertexConstant(tempSmooth,Color4b::Black);
        tri::io::ExporterPLY<CMesh>::Save(tempSmooth,("parameterization/renderSmooth"+std::to_string(original_curves[g].size())+std::to_string(g)+".ply").c_str(),tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);

        vector<CMesh::CoordType> aux;
        for (const size_t &ind: great_indexes[g]){
            aux.push_back(original_curves[g][ind].VFlip()->P());
        }
        utility::saveAsCurvePly(aux,("parameterization/greatindexes"+std::to_string(g)+".ply").c_str());
        Append<CMesh,CMesh>::Mesh(renderSmooth,tempSmooth);
    }

    /*CMesh temp;
    for(size_t g=0;g<blocks.size();g++){
        for(size_t i=0;i<blocks[g].size();i++){
            temp.Clear();
            for(size_t j=0;j<blocks[g][i].size()-1;j++){
               tri::Allocator<CMesh>::AddEdge(temp,blocks[g][i][j],blocks[g][i][j+1]);
            }
            string path=("parameterization/block"+std::to_string(i)+".ply");
            tri::io::ExporterPLY<CMesh>::Save(temp,path.c_str(),tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
        }
    }*/
}
void smoothIntersectionCurve::createSmoothIntersectionCurve(){
    renderSmooth.Clear();

    // add total number of vertices
    blocks_indexes.resize(blocks.size());
    size_t ind=0;
    for(size_t g=0;g<blocks.size();g++){
        blocks_indexes[g].resize(blocks[g].size());
        for(size_t i=0;i<blocks[g].size();i++){
            blocks_indexes[g][i].resize(blocks[g][i].size());
            for(size_t j=0;j<blocks[g][i].size()-1;j++){
               tri::Allocator<CMesh>::AddVertex(renderSmooth,blocks[g][i][j]);
               blocks_indexes[g][i][j]=ind;
               ind++;
            }
            blocks_indexes[g][i].back()=ind;
        }
        blocks_indexes[g].back().back()=blocks_indexes[g][0][0];
    }
    tri::io::ExporterPLY<CMesh>::Save(renderSmooth,"parameterization/renderSmooth000.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
}
void smoothIntersectionCurve::processCurves(vector<CMesh::VertexPointer> indexesPos, vector<face::Pos<CMesh::FaceType>> & edgestosplit,vector<int> &timestosubdivide){
    assert(threshold!=0.0);
    great_indexes.clear();
    blocks.clear();
    blocks_interpolation.clear();
    indexes_to_blocks.clear();
    great_indexes.resize(original_curves.size());
    blocks.resize(original_curves.size());
    blocks_interpolation.resize(original_curves.size());
    indexes_to_blocks.resize(original_curves.size());
        for(size_t j=0;j<original_curves.size();j++){
            great_indexes[j].clear();
        }
    for(size_t i=0;i<indexesPos.size();i++){
        for(size_t j=0;j<original_curves.size();j++){
            for(size_t k=0;k<original_curves[j].size();k++){
                if(original_curves[j][k].VFlip()==indexesPos[i]){
                    great_indexes[j].insert(k);
                }
            }
        }
    }
    for(size_t j=0;j<original_curves.size();j++){
        cout<<endl;
        for(const size_t & ind:great_indexes[j])
            cout<<" "<<ind;
        cout<<endl;
    }
    // generating blocks

    for(size_t j=0;j<original_curves.size();j++){
        vector<size_t> greatindexes(great_indexes[j].begin(),great_indexes[j].end());
        size_t n=original_curves[j].size();
        //for each block
        size_t numberSplits=0;
        vector<vector<CMesh::CoordType>> jblocks;
        vector< vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double>>> jblocks_interpolation;
        unordered_map<size_t,size_t> map;
        for(size_t r=0;r<greatindexes.size();r++){
            size_t start=greatindexes[r];
            size_t end=greatindexes[(r+1)%greatindexes.size()];
            CMesh::ScalarType len=0;
            //checking if some edge must be splited and computing length
            for(size_t l=start;l!=end;l=(l+1)%n){
                for(size_t s=0;s<edgestosplit.size();s++){
                    if(edgestosplit[s]==original_curves[j][l]){
                        cout<<" edge "<<l<<" in curve "<<j<<" to be splitted "<<timestosubdivide[s]<<" times"<<endl;
                        numberSplits+=size_t(timestosubdivide[s]);
                    }
                }
                //utility::printPoint3f(original_curves[j][l].VFlip()->P(),"source ");
                //utility::printPoint3f(original_curves[j][l].V()->P(),"target");
                len+=vcg::Distance(original_curves[j][l].VFlip()->P(),original_curves[j][l].V()->P());
            }
           // computing number of segments to be sampled mainteining the parity
           size_t numberSegments=0;
           bool iseven;
           if(end>start)
               iseven=((end-start+numberSplits)%2)==0;
           else
               iseven=((end+n-start+numberSplits)%2)==0;

           numberSegments=int(std::max<float>(std::round(len/threshold),1.0));
           if(iseven && numberSegments%2!=0){
               numberSegments++;
           }
           else if(!iseven && numberSegments%2==0){
               numberSegments++;
           }

           vector<Point3<CMesh::ScalarType>> ipoints;
           vector<face::Pos<CMesh::FaceType>> ipoints_interpolation;
           for(size_t l=start;l!=end;l=(l+1)%n){
               ipoints.push_back(original_curves[j][l].VFlip()->P());
               ipoints_interpolation.push_back(original_curves[j][l]); // these are half-edges
           }           
           ipoints.push_back(original_curves[j][end].VFlip()->P());
           //utility::saveAsCurvePly(ipoints,"parameterization/blockbefore.ply",false);
           utility::ParametricCurve<CMesh::ScalarType> parameterization(ipoints);
           vector<CMesh::CoordType> newpoints;
           vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double>> newpoints_interpolation;
           newpoints_interpolation.push_back(make_tuple(ipoints_interpolation[0].VFlip(),ipoints_interpolation[0].V(),1.0));
           newpoints.push_back(ipoints[0]);
           double step=len/numberSegments;
           for(int id=1;id<int(numberSegments);id++){
               newpoints.push_back(parameterization.getPointAtArcLength(id*step));
               // getting interpolation parameters
               double idlen=id*step;
               int segIdx = parameterization.getSegmentIndexAtArcLength(idlen);
               double u=idlen/parameterization.segmentLengths[size_t(segIdx)];
               newpoints_interpolation.push_back(make_tuple(ipoints_interpolation[size_t(segIdx)].VFlip(),ipoints_interpolation[size_t(segIdx)].V(),1.0-u));
           }
           newpoints_interpolation.push_back(make_tuple(ipoints_interpolation.back().VFlip(),ipoints_interpolation.back().V(),0.0));
           newpoints.push_back(ipoints.back());
           //utility::saveAsCurvePly(newpoints,"parameterization/blockafter.ply",false);
           assert(((ipoints.size()+numberSplits)%2)==(newpoints.size()%2));
           jblocks.push_back(newpoints);
           jblocks_interpolation.push_back(newpoints_interpolation);
           map.emplace(start,jblocks.size()-1);
        }
        blocks[j]=jblocks;
        blocks_interpolation[j]=jblocks_interpolation;
        indexes_to_blocks[j]=map;
    }
    createSmoothIntersectionCurve();
    createRenderSmooth();
}
void smoothIntersectionCurve::getSmoothVersion(face::Pos<CMesh::FaceType> & begin,face::Pos<CMesh::FaceType> & end,vector<CMesh::CoordType> & result){
    // finding the corresponding curve and endpoints indexes
    vector<size_t> thisendpoints;
    for(size_t j=0;j<original_curves.size();j++){
        for(size_t i=0;i<original_curves[j].size();i++){
            if(begin==original_curves[j][i]){
                assert(great_indexes[j].find(i)!=great_indexes[j].end());
                thisendpoints.push_back(i);
                for(size_t k=0;k<original_curves[j].size();k++){
                    if(end==original_curves[j][k]){
                        assert(great_indexes[j].find(k)!=great_indexes[j].end());
                        thisendpoints.push_back(k);
                    }
                }
            }
        }
    }
}
void smoothIntersectionCurve::getSmoothVersion(const face::Pos<CMesh::FaceType> & begin,const  face::Pos<CMesh::FaceType> & end,vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double>> & result){
    // Finding the corresponding curve and endpoints indexes
    result.clear();
    face::Pos<CMesh::FaceType> ibegin=begin,iend=end;
    size_t curve_index;
    vector<size_t> block_indexes;
    bool found=false;
    for(size_t j=0;j<original_curves.size();j++){
        for(size_t i=0;i<original_curves[j].size();i++){
            face::Pos<CMesh::FaceType> beginopp=ibegin;
            beginopp.FlipF();
            // finding fisrt block in this path
            if(ibegin==original_curves[j][i] || beginopp==original_curves[j][i]){
                assert(great_indexes[j].find(i)!=great_indexes[j].end()); // this index must be in the great_indexes array by construction
                found=true;
                curve_index=j;
                vector<size_t> greatindexes(great_indexes[j].begin(),great_indexes[j].end());
                unordered_map<size_t,size_t> &map=indexes_to_blocks[j];
                size_t thisindex_blocks=map[i];
                vector<vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double>>> &jblocks=blocks_interpolation[j];                ;
                // finding the set of blocks until the end vertexpointer
                size_t sizecurve=jblocks.size();
                bool allwell=false;
                for(size_t k=0;k<sizecurve;k++){
                    size_t id=(thisindex_blocks+k)%sizecurve;
                    CMesh::VertexPointer b=std::get<1>(jblocks[id].back());
                    block_indexes.push_back(id);
                    if(iend.V()==b){
                        //assert(great_indexes[j].find(id)!=great_indexes[j].end());
                        allwell=true;
                        break;
                    }
                }
                assert(allwell);
            }
        }
    }
    // if not found so we need to invert the orientation
    if(!found){
      iend.FlipV();
      ibegin.FlipV();
      face::Pos<CMesh::FaceType> aux=iend;
      iend=ibegin;
      ibegin=aux;
      getSmoothVersion(ibegin,iend,result);
    }
    cout<<"block indexes "<<endl;
    for(size_t i=0;i<block_indexes.size();i++){
        cout<<" "<<block_indexes[i];
    }
    //Creating la resulting path
    vector<CMesh::CoordType> path;
    for(size_t i=0;i<block_indexes.size();i++){
        vector<CMesh::CoordType> temp=blocks[curve_index][block_indexes[i]];
        path.insert(path.end(),temp.begin(),temp.end());
    }
    utility::saveAsCurvePly(path,"parameterization/blocktest.ply",false);
    // populating resulting vector
    vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double>> temp=blocks_interpolation[curve_index][block_indexes[0]];
    result.insert(result.end(),temp.begin(),temp.end());
    for(size_t i=1;i<block_indexes.size();i++){
        vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double>> block_in=blocks_interpolation[curve_index][block_indexes[i]];
        block_in.erase(block_in.begin());
        result.insert(result.end(),block_in.begin(),block_in.end());
    }
}
void smoothIntersectionCurve::getSmoothVersion(const face::Pos<CMesh::FaceType> & begin,const face::Pos<CMesh::FaceType> & end,vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,size_t>> & result,bool invertorder){
    // Finding the corresponding curve and endpoints indexes
    result.clear();
    face::Pos<CMesh::FaceType> ibegin=begin,iend=end;
    size_t curve_index;
    vector<size_t> block_indexes;
    bool found=false;
    for(size_t j=0;j<original_curves.size();j++){
        for(size_t i=0;i<original_curves[j].size();i++){
            face::Pos<CMesh::FaceType> beginopp=ibegin;
            beginopp.FlipF();
            // finding fisrt block in this path
            if(ibegin==original_curves[j][i] || beginopp==original_curves[j][i]){
                assert(great_indexes[j].find(i)!=great_indexes[j].end()); // this index must be in the great_indexes array by construction
                found=true;
                curve_index=j;
                vector<size_t> greatindexes(great_indexes[j].begin(),great_indexes[j].end());
                unordered_map<size_t,size_t> &map=indexes_to_blocks[j];
                size_t thisindex_blocks=map[i];
                vector<vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double>>> &jblocks=blocks_interpolation[j];                ;
                // finding the set of blocks until the end vertexpointer
                size_t sizecurve=jblocks.size();
                bool allwell=false;
                for(size_t k=0;k<sizecurve;k++){
                    size_t id=(thisindex_blocks+k)%sizecurve;
                    CMesh::VertexPointer b=std::get<1>(jblocks[id].back());
                    block_indexes.push_back(id);
                    if(iend.V()==b){
                        //assert(great_indexes[j].find(id)!=great_indexes[j].end());
                        allwell=true;
                        break;
                    }
                }
                assert(allwell);
            }
        }
    }
    // if not found so we need to invert the orientation
    if(!found){
      iend.FlipV();
      ibegin.FlipV();
      face::Pos<CMesh::FaceType> aux=iend;
      iend=ibegin;
      ibegin=aux;
      getSmoothVersion(ibegin,iend,result,true);
    }
    /*cout<<"block indexes "<<endl;
    for(size_t i=0;i<block_indexes.size();i++){
        cout<<" "<<block_indexes[i];
    }
    //Creating la resulting path
    vector<CMesh::CoordType> path;
    for(size_t i=0;i<block_indexes.size();i++){
        vector<CMesh::CoordType> temp=blocks[curve_index][block_indexes[i]];
        path.insert(path.end(),temp.begin(),temp.end());
    }*/
    //utility::saveAsCurvePly(path,"parameterization/blocktest.ply",false);
    // populating resulting vector
    for(size_t i=0;i<block_indexes.size();i++){
        vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,double>> temp_vp=blocks_interpolation[curve_index][block_indexes[i]];
        vector<size_t> temp_ind=blocks_indexes[curve_index][block_indexes[i]];
        vector<std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,size_t>> temp;
        assert(temp_ind.size()==temp_vp.size());
        for(size_t j=0;j<temp_vp.size();j++){
            CMesh::VertexPointer a=std::get<0>(temp_vp[j]);
            CMesh::VertexPointer b=std::get<1>(temp_vp[j]);
            std::tuple<CMesh::VertexPointer,CMesh::VertexPointer,size_t> data(a,b,temp_ind[j]);
            temp.push_back(data);
        }
        if(i!=0)
            temp.erase(temp.begin());
        result.insert(result.end(),temp.begin(),temp.end());
    }
    if(invertorder){
        std::reverse(result.begin(),result.end());
    }
}
