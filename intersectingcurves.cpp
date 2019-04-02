#include "intersectingcurves.h"
#include<wrap/io_trimesh/export.h>
intersectingCurves::intersectingCurves()
{
    cloneml=NULL;
    clonemr=NULL;
    thrZero=0;
    thrMax=0;
    merged=false;
    maxIteration=30;
    numberIteration=0;
    IntersectedCurves.Clear();
    intersectingLoopPointers.clear();
    pairLoops.clear();
    pairPositionLoops.clear();
    maptoproject.clear();
    intersectingLoop.clear();
}

intersectingCurves::intersectingCurves(CMesh &ml, CMesh &mr, vector<pair<vector<CMesh::VertexPointer>, vector<CMesh::VertexPointer> > > &loops,
                                       vector<pair<vector<vector<Point3f>>,vector<vector<Point3f>>>>& mapSF){

    redefine(ml,mr,loops,mapSF);
}
void intersectingCurves::reset(){
    IntersectedCurves.Clear();
    intersectingLoop.clear();
    intersectingLoopPointers.clear();
}
void intersectingCurves::updatethresholds(){
    thrMax=(utility::edgeAverage(*cloneml)+utility::edgeAverage(*clonemr))/2;
    thrZero=thrMax/10;
    thrZero=(thrZero>0.01)?0.01:thrZero;
}
void intersectingCurves::setThrMax(float max){
    thrMax=max;
    thrZero=thrMax/10;
    thrZero=(thrZero>0.01)?0.01:thrZero;
}

void intersectingCurves::redefine(CMesh &ml, CMesh &mr, vector<pair<vector<CMesh::VertexPointer>, vector<CMesh::VertexPointer> > >&loops,
                                  vector<pair<vector<vector<Point3f>>,vector<vector<Point3f>>>>& mapSF){

    cloneml=&ml;
    clonemr=&mr;
    maxIteration=30;
    numberIteration=0;
    merged=false;
    IntersectedCurves.Clear();
    intersectingLoopPointers.clear();
    pairPositionLoops.clear();
    intersectingLoop.clear();
    pairLoops.clear();
    maptoproject.clear();
    pairLoops.assign(loops.begin(),loops.end());
    updatethresholds();

    for(int i=0;i<loops.size();i++){
        std::forward_list<Point3f> first;
        std::forward_list<Point3f> second;
        std::forward_list<Point3f>::iterator itf=first.before_begin();
        std::forward_list<Point3f>::iterator its=second.before_begin();
        for(int j=0;j<loops[i].first.size();j++)
            itf=first.insert_after(itf,loops[i].first[j]->P());
        for(int j=0;j<loops[i].second.size();j++)
            its=second.insert_after(its,loops[i].second[j]->P());
        pairPositionLoops.push_back(make_pair(first,second));
    }
    maptoproject=mapSF;

}

void intersectingCurves::resetPairPositionLoops(int index){

    if( pairLoops.size()!=0 && index==-1){
        pairPositionLoops.clear();
        for(int i=0;i<pairLoops.size();i++){
            std::forward_list<Point3f> first;
            std::forward_list<Point3f> second;
            std::forward_list<Point3f>::iterator itf=first.before_begin();
            std::forward_list<Point3f>::iterator its=second.before_begin();
            for(int j=0;j<pairLoops[i].first.size();j++)
                itf=first.insert_after(itf,pairLoops[i].first[j]->P());
            for(int j=0;j<pairLoops[i].second.size();j++)
                its=second.insert_after(its, pairLoops[i].second[j]->P());
            pairPositionLoops.push_back(make_pair(first,second));
        }
    }
    else if(pairLoops.size()!=0 && index>=0 && index<pairLoops.size()){
        std::forward_list<Point3f> first;
        std::forward_list<Point3f> second;
        std::forward_list<Point3f>::iterator itf=first.before_begin();
        std::forward_list<Point3f>::iterator its=second.before_begin();
        for(int j=0;j<pairLoops[index].first.size();j++)
            itf=first.insert_after(itf,pairLoops[index].first[j]->P());
        for(int j=0;j<pairLoops[index].second.size();j++)
            its=second.insert_after(its, pairLoops[index].second[j]->P());
        pairPositionLoops[index]=make_pair(first,second);
    }
}

vcg::Point3f intersectingCurves::nearest(const std::forward_list<Point3f>& loop, vcg::Point3f& x){
       Point3f min;
       float dist=std::numeric_limits<float>::infinity();
       forward_list<Point3f>::const_iterator it=loop.begin();
       while(it!=loop.end())
       {
            float d=((*it)-x).SquaredNorm();

            if(d<dist){
                min=(*it);
                dist=d;
            }
            it++;
       }

       return min;
}
vcg::Point3f intersectingCurves::nearest(const std::forward_list<Point3f*>& loop, vcg::Point3f& x){
       Point3f min;
       float dist=std::numeric_limits<float>::infinity();
       forward_list<Point3f*>::const_iterator it=loop.begin();
       while(it!=loop.end())
       {
            float d=((**it)-x).SquaredNorm();

            if(d<dist){
                min=(**it);
                dist=d;
            }
            it++;
       }

       return min;
}

int intersectingCurves::nearestIndex(const std::forward_list<Point3f>& loop, vcg::Point3f& x){
       int n=0;
       int min=-1;
       float dist=std::numeric_limits<float>::infinity();
       forward_list<Point3f>::const_iterator it=loop.begin();
       while(it!=loop.end())
       {
            float d=((*it)-x).SquaredNorm();
            if(d<dist){
                min=n;
                dist=d;
            }
            n++;
            it++;
       }
       return min;
}

void intersectingCurves::moreNear(std::forward_list<Point3f>& left,std::forward_list<Point3f>& right,
                                  vector<vector<Point3f>>& toprojectleft,
                                  vector<vector<Point3f>>& toprojectright){
        std::forward_list<Point3f>::iterator it;
        std::forward_list<Point3f*> copyleft;
        std::forward_list<Point3f*> copyright;
        std::vector<vcg::Point3f> newLeftLoop;
        std::vector<vcg::Point3f> newRightLoop;

        float t=0.3;
        //cout<<" size loopleft "<<leftLoop.size()<<endl;
        //cout<<" size loopRight "<<rightLoop.size()<<endl;

        //Preprocessing
        //doing copy of the initial lists
        std::forward_list<Point3f*>::iterator itp=copyleft.before_begin();
        for(auto& e:left)
            itp=copyleft.insert_after(itp,&e);
        itp=copyright.before_begin();
        for(auto& e:right)
            itp=copyright.insert_after(itp,&e);

        //computing fixed points
        itp=copyleft.before_begin();
        while(itp!=copyleft.end())
        {
            std::forward_list<Point3f*>::iterator past=itp;
            itp++;
            if(itp==copyleft.end()) break;
            vcg::Point3f vn=nearest(right,**itp);
            // clean lists copies of fixed points
            if((vn-(**itp)).SquaredNorm()<thrZero*thrZero)
                itp=copyleft.erase_after(past);
        }
        itp=copyright.before_begin();
        while(itp!=copyright.end())
        {
            std::forward_list<Point3f*>::iterator past=itp;
            itp++;
            if(itp==copyright.end()) break;
            vcg::Point3f vn=nearest(left,**itp);
            // clean lists copies of fixed points
            if((vn-(**itp)).SquaredNorm()<thrZero*thrZero)
                itp=copyright.erase_after(past);
        }


        //walk more near only in non fixed points
        for(auto& e:copyleft)
        {
            vcg::Point3f vn;
            vn=nearest(copyright,*e);
            //tri::Allocator<CMesh>::AddEdge(IntersectedCurves,*e,vn);
            //vn.Scale(t,t,t);
            vcg::Point3f vm((1-t)*(e->X())+t*vn[0],(1-t)*(e->Y())+t*vn[1],(1-t)*(e->Z())+t*vn[2]);
            utility::easyProjectToMesh(toprojectleft,vm);
            newLeftLoop.push_back(vm);
        }

        for(auto& e:copyright)
        {
            vcg::Point3f vn;
            vn=nearest(copyleft,*e);
            //tri::Allocator<CMesh>::AddEdge(IntersectedCurves,*e,vn);
            //vn.Scale(t,t,t);
            vcg::Point3f vm((1-t)*(e->X())+t*vn[0],(1-t)*(e->Y())+t*vn[1],(1-t)*(e->Z())+t*vn[2]);
            utility::easyProjectToMesh(toprojectright,vm);
            newRightLoop.push_back(vm);
        }
        itp=copyleft.begin();
        for(int i=0;i<newLeftLoop.size();i++){
            **itp=newLeftLoop[i];
            itp++;
        }
        itp=copyright.begin();
        for(int i=0;i<newRightLoop.size();i++){
            **itp=newRightLoop[i];
            itp++;
        }
}

int intersectingCurves::refineLoop(std::forward_list<Point3f>& loop,bool iscircular){

       std::forward_list<Point3f>::iterator li;
       float average=thrMax;
       float edgemin=0.3*average;
       float minimo=(edgemin>thrZero)?thrZero:edgemin;
       //float minimo=(edgemin<thrZero)?thrZero:edgemin;
       int n=0;// number of points;

       //Here n count number of segments. In the final we add if necessary +1;
       // Do the split step
       li=loop.begin();

       while(li!=loop.end()){
         std::forward_list<Point3f>::iterator past=li;
         li++;
         if(li==loop.end()) break;
         n++;
         float d=vcg::Distance(*li,*past);
         if(d>average){
             Point3f newp=((*li)+(*past))/2;
             loop.insert_after(past,newp);
             n++;
         }
       }
       if(iscircular){
         n++;
         li=loop.before_begin();
         std::advance(li,n);
         float d=vcg::Distance(*li,*loop.begin());
         if(d>average){
           Point3f newp=((*li)+(*loop.begin()))/2;
           loop.insert_after(li,newp);
           n++;
         }
       }
       //Do the collapse step
       li=loop.begin();
       while(li!=loop.end()){
         std::forward_list<Point3f>::iterator past=li;
         li++;
         if(li==loop.end()) break;
         float d=vcg::Distance(*li,*past);
         if(d<edgemin){
             Point3f newp=((*li)+(*past))/2;
             *past=newp;
             li=loop.erase_after(past);
             n--;
         }
       }
       if(iscircular){
         li=loop.before_begin();
         std::advance(li,n);
         float d=vcg::Distance(*loop.begin(),*li);
         if(d<edgemin){
           Point3f newp=((*loop.begin())+(*li))/2;
           *li=newp;
           loop.pop_front();
           n--;
         }
       }
       else{
           n++; //now n is the number of points of this polyline
       }
       return n;
}
void intersectingCurves::mergeLoops(){
    if(!merged){
       //For each loops to merge, we iterate until converge
       bool premerged=true;
       float oldThrMax=thrMax;
       for(int i=0;i<pairPositionLoops.size();i++){
           cout<<"iteration pair "<<i<<" "<<endl;
           numberIteration=0;
           bool thisismerged=false;
           float averageleft=averageEdgeLenghtLoop(pairPositionLoops[i].first);
           float averageright=averageEdgeLenghtLoop(pairPositionLoops[i].second);
           //thrMax=std::min(averageleft,averageright);
           //setThrMax(std::min(averageleft,averageright));
           setThrMax((averageleft+averageright)/2);
           //improveInitialLoop(pairPositionLoops[i].first);
           //improveInitialLoop(pairPositionLoops[i].second);
           addloopToMesh(pairPositionLoops[i].first);
           addloopToMesh(pairPositionLoops[i].second);
           int nleft=sizefowardlist(pairPositionLoops[i].first);
           int nright=sizefowardlist(pairPositionLoops[i].second);
           //cout<<"lefn  "<<nleft<<endl;
           //cout<<"rightn "<<nright<<endl;
           int  epsilon=-1; // if the size difference of the polylines is less than it so this polylines must be merged
           while(numberIteration<maxIteration && !thisismerged){
              if(std::abs(nleft-nright)>epsilon && nleft>3 && nright>3){
                  moreNear(pairPositionLoops[i].first,pairPositionLoops[i].second,maptoproject[i].first,maptoproject[i].second);
                  nleft=refineLoop(pairPositionLoops[i].first);
                  nright=refineLoop(pairPositionLoops[i].second);
                  //cout<<"left size  "<<nleft<<" "<<sizefowardlist(pairPositionLoops[i].first)<<endl;
                  //cout<<"right size  "<<nright<<" "<<sizefowardlist(pairPositionLoops[i].second)<<endl;
              }
              else{
                  cout<<"number interation pair "<<i<<" is "<<numberIteration<<endl;
                  cout<<"final size left polyline"<<nleft<<endl;
                  cout<<"final size right polyline"<<nright<<endl;
                  //addloopToMesh(pairPositionLoops[i].first);
                  //addloopToMesh(pairPositionLoops[i].second);
                  vector<Point3f> tempIntersectingLoop;
                  std::unordered_map<int,forward_list<Point3f>::iterator>  indexLeft=fowardlistToIndexes(pairPositionLoops[i].first);
                  std::unordered_map<int,forward_list<Point3f>::iterator>  indexRight=fowardlistToIndexes(pairPositionLoops[i].second);
                  vector<std::pair<int,int>> indexPairs;
                  forward_list<Point3f> final;
                  std::forward_list<Point3f>::iterator itf=final.before_begin();
                  if(nleft<= nright){
                       //computing best match points
                      for(int j=0;j<nleft;j++){
                            int nIndex=nearestIndex(pairPositionLoops[i].second,*indexLeft[j]);
                            indexPairs.push_back(make_pair(j,nIndex));
                      }
                      //Know for each element in the second loop how many match are in the first loop
                      std::unordered_map<int,std::vector<int>> mymap;
                      for(int k=0;k<nright;k++){
                          bool isin=false;
                          vector<int> leftAssociated;
                          for(int j=0;j<indexPairs.size();j++){
                              if(k==indexPairs[j].second){
                                  isin=true;
                                  leftAssociated.push_back(j);
                              }
                          }
                          if(isin)
                               mymap.emplace(k, leftAssociated);
                      }
                      for (auto& x: mymap)
                      {
                           if(x.second.size()==1){
                               Point3f pl=*indexLeft[x.second.front()];
                               utility::easyProjectToMesh(maptoproject[i].first,pl);
                               Point3f pr=*indexRight[x.first];
                               utility::easyProjectToMesh(maptoproject[i].second,pr);
                               Point3f med=(pl+pr)/2;
                               //a last projection to get better results
                               pl=med;
                               pr=med;
                               utility::easyProjectToMesh(maptoproject[i].first,pl);
                               utility::easyProjectToMesh(maptoproject[i].second,pr);
                               med=(pl+pr)/2;
                               *indexLeft[x.second.front()]=med;
                               *indexRight[x.first]=med;
                           }
                           else{  //Compute indexes of leftLoop that have the same right pair
                              Point3f averagep=*indexRight[x.first];
                              utility::easyProjectToMesh(maptoproject[i].second,averagep);
                              int min=0;
                              float d=vcg::SquaredDistance(*indexRight[x.first],*indexLeft[x.second[0]]);
                              for(int j=1;j<x.second.size();j++){
                                  float dd=vcg::SquaredDistance(averagep,*indexLeft[x.second[j]]);
                                  if(dd<d){
                                      d=dd;
                                      min=j;
                                  }
                                  //averagep+=(*indexLeft[x.second[j]]);
                              }
                              Point3f pl=*indexLeft[x.second[min]];
                              utility::easyProjectToMesh(maptoproject[i].first,pl);
                              averagep+=pl;
                              averagep/=2;
                              //a last projection to get better results
                              pl=averagep;
                              utility::easyProjectToMesh(maptoproject[i].second,averagep);
                              utility::easyProjectToMesh(maptoproject[i].first,pl);
                              averagep+=pl;
                              averagep/=2;

                              *indexRight[x.first]=averagep;
                              for(int j=0;j<x.second.size();j++){
                                  *indexLeft[x.second[j]]=averagep;
                              }
                           }
                      }
                      std::set<int> toskip;
                      for(int k=0;k<nleft;k++){
                          if(mymap[indexPairs[k].second].size()!=1)
                              for(int rr=1;rr<mymap[indexPairs[k].second].size();rr++)
                                toskip.insert(mymap[indexPairs[k].second][rr]);
                          if(!toskip.empty()){
                            bool isrepeat=false;
                            for(std::set<int>::iterator it=toskip.begin(); it!=toskip.end(); ++it)
                              if(k==*it)
                                  isrepeat=true;
                            if(!isrepeat){
                                itf=final.insert_after(itf,*indexLeft[indexPairs[k].first]);
                                tempIntersectingLoop.push_back(*indexLeft[indexPairs[k].first]);
                                //utility::printPoint3f((*indexLeft[indexPairs[k].first]),"gola");
                            }

                          }
                          else{
                              itf=final.insert_after(itf,*indexLeft[indexPairs[k].first]);
                              tempIntersectingLoop.push_back(*indexLeft[indexPairs[k].first]);
                              //utility::printPoint3f((*indexLeft[indexPairs[k].first]),"gola");
                          }
                      }

                  }
                  else{
                      //computing best match points
                     for(int j=0;j<nright;j++){
                           int nIndex=nearestIndex(pairPositionLoops[i].first,*indexRight[j]);
                           indexPairs.push_back(make_pair(j,nIndex));
                     }

                     std::unordered_map<int,std::vector<int>> mymap;
                     for(int k=0;k<nleft;k++){
                         bool isin=false;
                         vector<int> rightAssociated;
                         for(int j=0;j<indexPairs.size();j++){
                             if(k==indexPairs[j].second){
                                 isin=true;
                                 rightAssociated.push_back(j);
                             }
                         }
                         if(isin)
                              mymap.emplace(k, rightAssociated);
                     }
                     for (auto& x: mymap)
                     {
                          if(x.second.size()==1){
                              Point3f pl=*indexLeft[x.first];
                              utility::easyProjectToMesh(maptoproject[i].first,pl);
                              Point3f pr=*indexRight[x.second.front()];
                              utility::easyProjectToMesh(maptoproject[i].second,pr);
                              Point3f med=(pl+pr)/2;
                              //a last projection to get better results
                              pl=med;
                              pr=med;
                              utility::easyProjectToMesh(maptoproject[i].first,pl);
                              utility::easyProjectToMesh(maptoproject[i].second,pr);
                              med=(pl+pr)/2;
                              //utility::easyProjectToMesh(maptoproject[i].first,med);
                              *indexRight[x.second.front()]=med;
                              *indexLeft[x.first]=med;
                          }
                          else{  //Compute indexes of rightLoop that have the same left pair
                             Point3f averagep=*indexLeft[x.first];
                             utility::easyProjectToMesh(maptoproject[i].first,averagep);
                             int min=0;
                             float d=vcg::SquaredDistance(*indexLeft[x.first],*indexRight[x.second[0]]);
                             for(int j=1;j<x.second.size();j++){
                                 float dd=vcg::SquaredDistance(averagep,*indexRight[x.second[j]]);
                                 if(dd<d){
                                     d=dd;
                                     min=j;
                                 }
                                 //averagep+=(*indexRight[x.second[j]]);
                             }
                             Point3f pr=*indexRight[x.second[min]];
                             utility::easyProjectToMesh(maptoproject[i].second,pr);
                             averagep+=pr;
                             averagep/=2;
                             //a last projection to get better results
                             pr=averagep;
                             utility::easyProjectToMesh(maptoproject[i].first,averagep);
                             utility::easyProjectToMesh(maptoproject[i].second,pr);
                             averagep+=pr;
                             averagep/=2;
                             *indexLeft[x.first]=averagep;
                             for(int j=0;j<x.second.size();j++){
                                 *indexRight[x.second[j]]=averagep;
                             }
                          }
                     }
                     std::set<int> toskip;
                     for(int k=0;k<nright;k++){
                         if(mymap[indexPairs[k].second].size()!=1)
                             for(int rr=1;rr<mymap[indexPairs[k].second].size();rr++)
                               toskip.insert(mymap[indexPairs[k].second][rr]);
                         if(!toskip.empty()){
                           bool isrepeat=false;
                           for(std::set<int>::iterator it=toskip.begin(); it!=toskip.end(); ++it)
                             if(k==*it)
                                 isrepeat=true;
                           if(!isrepeat){
                               itf=final.insert_after(itf,*indexRight[indexPairs[k].first]);
                               tempIntersectingLoop.push_back(*indexRight[indexPairs[k].first]);
                           }
                         }
                         else{
                             itf=final.insert_after(itf,*indexRight[indexPairs[k].first]);
                             tempIntersectingLoop.push_back(*indexRight[indexPairs[k].first]);
                         }
                     }
                  }
                  //addloopToMesh(final);
                  //tri::Smooth<CMesh>::VertexCoordLaplacian(IntersectedCurves,50);
                  utility::correctCornersFromCurve(tempIntersectingLoop);
                  intersectingLoop.push_back(tempIntersectingLoop);
                  thisismerged=true;
                  cout<<"FINAL SIZE INTERSECTION CURVE"<<sizefowardlist(final)<<endl;
                  /*cout<<"n left "<<sizefowardlist(final)<<endl;
                  cout<<"n right "<<sizefowardlist(final)<<endl;
                  */
              }
              numberIteration++;
              if(numberIteration==maxIteration){
                  if(epsilon<3){
                      numberIteration=0;
                      epsilon++;
                      resetPairPositionLoops(i);
                      averageleft=averageEdgeLenghtLoop(pairPositionLoops[i].first);
                      averageright=averageEdgeLenghtLoop(pairPositionLoops[i].second);
                      //setThrMax(std::min(averageleft,averageright));
                      setThrMax((averageleft+averageright)/2);
                      //improveInitialLoop(pairPositionLoops[i].first);
                      //improveInitialLoop(pairPositionLoops[i].second);
                      nleft=sizefowardlist(pairPositionLoops[i].first);
                      nright=sizefowardlist(pairPositionLoops[i].second);
                  }
              }
           }
           premerged=premerged && thisismerged;
       }
       merged=premerged;
       if(merged){
           setThrMax(oldThrMax);
           /*for(int i=0;i<pairLoops.size();i++){
                vector<CMesh::VertexPointer> first(pairLoops[i].first);
                vector<CMesh::VertexPointer> second(pairLoops[i].second);
                bool ip1=improveloop(first,(*cloneml));
                bool ip2=improveloop(second,(*clonemr));
                if(ip1 && ip2)
                   pairLoops[i]=make_pair(first,second);
           }*/

           for(int i=0;i<intersectingLoop.size();i++){
               vector<CMesh::VertexPointer> newll,newlll;
               vector<CMesh::VertexPointer> first(pairLoops[i].first);
               vector<CMesh::VertexPointer> second(pairLoops[i].second);
               fillloops(i,second,newll,(*clonemr));
               fillloops(i,first,newlll,(*cloneml));
               intersectingLoopPointers.push_back(make_pair(newlll,newll));
               pairLoops[i]=make_pair(first,second);
           }
           vector<CMesh::FaceType*> retVec;
           tri::Clean<CMesh>::SelfIntersections(*cloneml, retVec);
           cout<<"Mesh had "<<retVec.size()<< " selfintersecting faces\n"<<endl;
           tri::Clean<CMesh>::SelfIntersections(*clonemr, retVec);
           cout<<"Mesh had "<<retVec.size()<< " selfintersecting faces\n"<<endl;
           //vcg::tri::Append< CMesh,CMesh>::Mesh(*cloneml,*clonemr);
           //vcg::tri::Clean<CMesh>::MergeCloseVertex(*cloneml,thrZero/10);
           //doIsoRemeshing(*cloneml);
           //doIsoRemeshing(*clonemr);
           tri::UpdateTopology<CMesh>::FaceFace(*cloneml);
           tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFaceNormalized(*cloneml);
           tri::UpdateTopology<CMesh>::FaceFace(*clonemr);
           tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFaceNormalized(*clonemr);
           //vcg::tri::io::ExporterPLY<CMesh>::Save(*cloneml, "exported.ply");
           //vcg::tri::UpdateColor<CMesh>::PerFaceConstant(*cloneml, vcg::Color4b::Red,true);
           //vcg::tri::UpdateColor<CMesh>::PerFaceConstant(*clonemr, vcg::Color4b::Red,true);
           //doIsoRemeshing(*cloneml);
           //tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFaceNormalized(*cloneml);
           //doIsoRemeshing(*clonemr);
       }
       else{
           cout<<"no entrei"<<endl;
       }
    }
    else cout<<"just merged"<<endl;
}
int intersectingCurves::sizefowardlist(const std::forward_list<Point3f>& loop){
    int n=0;
    for(auto& e:loop)
      n++;
    return n;
}
std::unordered_map<int,forward_list<Point3f>::iterator> intersectingCurves::fowardlistToIndexes(std::forward_list<Point3f>& loop){
    std::unordered_map<int,forward_list<Point3f>::iterator> result;
    int n=0;
    forward_list<Point3f>::iterator it=loop.begin();
    while(it!=loop.end()){
       result.emplace(n,it);
       n++;
       it++;
    }
    return result;
}
float intersectingCurves::averageEdgeLenghtLoop(const std::forward_list<Point3f> &loop,bool iscircular){
    int n=0;
    float length=0;
    for(forward_list<Point3f>::const_iterator li=loop.begin();li!=loop.end();li++){
        forward_list<Point3f>::const_iterator past=li;
        li++;
        if(li==loop.end()) break;
        length+=vcg::Distance(*li,*past);
        n++;
        li=past;
    }
    if(iscircular){
        forward_list<Point3f>::const_iterator li=loop.before_begin();
        n++;
        std::advance(li,n);
        length+=vcg::Distance(*li,*loop.begin());
    }
    length=length/n;
    //cout<<"number of point"<<sizefowardlist(loop)<<endl;
    //cout<<"number of segments"<<n<<endl;
    return length;
}
void intersectingCurves::addloopToMesh(const std::forward_list<Point3f>& loop,bool iscircular){
    tri::Allocator<CMesh>::CompactEveryVector(IntersectedCurves);
    std::forward_list<Point3f>::const_iterator it=loop.begin();
    std::forward_list<Point3f>::const_iterator past=loop.before_begin();
    int n=0;
    while(it!=loop.end()){
           past=it;
           it++;
           if(it!=loop.end()) {
               tri::Allocator<CMesh>::AddEdge(IntersectedCurves,*past,*it);
               n++;
           }
    }
    it=loop.before_begin();
    std::advance(it,n+1);
    if(iscircular)
           tri::Allocator<CMesh>::AddEdge(IntersectedCurves,*it,*loop.begin());

    tri::Clean<CMesh>::RemoveDuplicateVertex(IntersectedCurves);
    tri::Clean<CMesh>::RemoveUnreferencedVertex(IntersectedCurves);
    tri::Allocator<CMesh>::CompactEveryVector(IntersectedCurves);
    tri::UpdateTopology<CMesh>::VertexEdge(IntersectedCurves);
    tri::UpdateTopology<CMesh>::EdgeEdge(IntersectedCurves);
    tri::UpdateColor<CMesh>::PerVertexConstant(IntersectedCurves,Color4b::Black);
}
void intersectingCurves::addloopToMesh(const vector<CMesh::VertexPointer>& loop,bool iscircular){
    tri::Allocator<CMesh>::CompactEveryVector(IntersectedCurves);
    int n=loop.size();
    for(int i=0;i<loop.size()-1;i++){
       tri::Allocator<CMesh>::AddEdge(IntersectedCurves,loop[i]->P(),loop[(i+1)%n]->P());
    }
    if(iscircular)
       tri::Allocator<CMesh>::AddEdge(IntersectedCurves,loop[n-1]->P(),loop[0]->P());

    tri::Clean<CMesh>::RemoveDuplicateVertex(IntersectedCurves);
    tri::Clean<CMesh>::RemoveUnreferencedVertex(IntersectedCurves);
    tri::Allocator<CMesh>::CompactEveryVector(IntersectedCurves);
    tri::UpdateTopology<CMesh>::VertexEdge(IntersectedCurves);
    tri::UpdateTopology<CMesh>::EdgeEdge(IntersectedCurves);
    tri::UpdateColor<CMesh>::PerVertexConstant(IntersectedCurves,Color4b::Black);
}
//try to equalize size of edges and to avoid the vertices with very sharp angles
void intersectingCurves::improveInitialLoop(std::forward_list<Point3f> &loop,bool iscircular){
    std::forward_list<Point3f>::iterator it=loop.begin();
    std::forward_list<Point3f> result;
    std::forward_list<Point3f>::iterator itr=result.before_begin();
    int n=sizefowardlist(loop);
    const float PI=3.14159265;
    while(it!=loop.end()){
        std::forward_list<Point3f>::iterator past=it;
        it++;
        if(it==loop.end()) break;
        std::forward_list<Point3f>::iterator present=it;
        it++; // next to present
        if(it==loop.end()) break;
        Point3f a=(*past)-(*present);
        Point3f b=(*it)-(*present);
        if(vcg::Angle(a,b)*180/PI<60) {
            cout<<"avoid vetices with sharp angle"<<endl;
        }
        else{
            itr=result.insert_after(itr,*present);
            it=present;
        }
    }
    if(iscircular){
        it=loop.before_begin();
        std::advance(it,n-1);
        std::forward_list<Point3f>::iterator past=it;
        it++;
        std::forward_list<Point3f>::iterator present=it;
        it=loop.begin();
        Point3f a=(*past)-(*present);
        Point3f b=(*it)-(*present);
        if(vcg::Angle(a,b)*180/PI<60) {
            cout<<"avoid vetices with sharp angle"<<endl;
        }
        else{
            itr=result.insert_after(itr,*present);
        }
    }

    loop=result;
}

bool intersectingCurves::improveloop(vector<CMesh::VertexPointer> &loop, CMesh& m){
    /*std::vector<CMesh::VertexPointer> result;
    UpdateTopology<CMesh>::FaceFace(m);
    tri::UpdateFlags<CMesh>::VertexClearS(m);
    //tri::UpdateSelection<CMesh>::VertexClear(m);

    tri::UpdateFlags<CMesh>::FaceClearS(m);
    //tri::UpdateSelection<CMesh>::FaceClear(m);
    for(int i=0;i<loop.size();i++)
        loop[i]->SetS();
    tri::UpdateSelection<CMesh>::FaceFromVertexLoose(m);
    utility::splitMesh(m,thrZero,true);
    //utility::splitMesh(m,thrZero,true);
    //tri::UpdateSelection<CMesh>::VertexFromFaceStrict(m);
    //UpdateTopology<CMesh>::VertexFace(m);
    tri::UpdateFlags<CMesh>::FaceClearS(m);
    tri::UpdateFlags<CMesh>::VertexClearS(m);
    for(int i=0;i<loop.size();i++)
        loop[i]->SetS();
    tri::UpdateSelection<CMesh>::FaceFromVertexLoose(m);
    tri::UpdateSelection<CMesh>::VertexFromFaceLoose(m);
    tri::UpdateSelection<CMesh>::FaceFromVertexStrict(m);
    tri::UpdateFlags<CMesh>::FaceBorderFromFF(m);

    for(CMesh::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
        if(!(*fi).IsD() && (*fi).IsS()){

          vcg::tri::Allocator<CMesh>::DeleteFace(m,*fi);
        }
    }
    CMesh::VertexPointer vp;
    tri::UpdateTopology<CMesh>::FaceFace(m);
    tri::UpdateTopology<CMesh>::VertexFace(m);
    tri::UpdateFlags<CMesh>::FaceBorderFromFF(m);
    tri::UpdateFlags<CMesh>::VertexBorderFromFaceAdj(m);
    for(CMesh::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
    {
        if(!(*vi).IsD() && (*vi).IsS()){
            vp=&*vi;
            face::VFIterator<CFace> vfi(vp);
            bool goodvertex=false;
            for(;!vfi.End();++vfi)
                goodvertex=true;
            if (goodvertex) break;
        }
    }
    int ind=-1;
    bool stop=false;
    face::VFIterator<CFace> vfi(vp);
    CMesh::FacePointer fv;
    for(;!vfi.End();++vfi){
        //cout<<"neighboor"<<endl;
        CMesh::FacePointer f = vfi.F();
        for(int k=0;k<3;k++){
          if(face::IsBorder(*f,k) && (f->V(k)==vp || f->V(f->Next(k))==vp)){
          //if(f->IsB(k)){
              fv=f;
              //f->SetS();
              ind=k;
              stop=true;
          }
        }
        if(stop) break;
    }
    if(ind==-1){
        cout<<"improve loop fail"<<endl;
    }
    else{
        const float PI=3.14159265;
        face::Pos<CMesh::FaceType> startPos(fv,ind);
        face::Pos<CMesh::FaceType> curPos=startPos;
        vector<face::Pos<CMesh::FaceType>> betterloop;
        result.push_back(curPos.V());
        betterloop.push_back(curPos);
        do
        {
          curPos.NextB();
          //curPos.F()->SetS();
          betterloop.push_back(curPos);
          result.push_back(curPos.V());
        }
        while(curPos!=startPos);
        result.pop_back();
        betterloop.pop_back();
        loop=result;
        int n=loop.size();
        vector<CMesh::FacePointer> todelete;
        vector<int> WrongIndexes;
        for(int k=0;k<n;k++){
             vcg::face::Pos<CMesh::FaceType> jp=betterloop[k];
             int angle=0;
             do{
                 Point3f a=(jp.VFlip()->P())-(jp.V()->P());
                 jp.FlipE();
                 Point3f b=(jp.VFlip()->P())-(jp.V()->P());
                 angle+=vcg::Angle(a,b)*180/PI;
                 jp.FlipF();
                 if(jp.VFlip()==betterloop[(k+1)%n].V())
                       break;
             }while(jp.f!=betterloop[k].f);
             if(angle<150){
                 WrongIndexes.push_back(k);
                 jp=betterloop[k];
                 do{
                     todelete.push_back(jp.F());
                     //jp.F()->SetS();
                     jp.FlipE();
                     jp.FlipF();
                     if(jp.VFlip()==betterloop[(k+1)%n].V())
                           break;
                 }while(jp.f!=betterloop[k].f);
             }

        }
        //Compute blocks to improve
        vector<int> FinalIndices;
        for(int k=0;k<WrongIndexes.size();k++){
            int sup=(WrongIndexes[k]+4)%n;
            int inf=(WrongIndexes[k]-4+n)%n;
            int next=WrongIndexes[(k+1)%WrongIndexes.size()];
            int before=WrongIndexes[(k-1+WrongIndexes.size())%WrongIndexes.size()];
            bool yes=false;
            while(inf!=WrongIndexes[k]){
                if(inf==before){
                    yes=true;
                    break;
                }
                inf=(inf==(n-1))?0:inf+1;
            }
            while(sup!=WrongIndexes[k]){
                if(sup==next){
                    yes=true;
                    break;
                }
                sup=(sup==0)?(n-1):sup-1;
            }
            if(yes)
                FinalIndices.push_back(WrongIndexes[k]);
        }
        for(int k=0;k<FinalIndices.size();k++){
            Point3f v0=betterloop[FinalIndices[k]].VFlip()->P();
            Point3f v1=betterloop[FinalIndices[k]].V()->P();
            Point3f v2=betterloop[(FinalIndices[k]+1)%n].V()->P();
            Point3f va=v1-v0;
            Point3f vb=v2-v0;
            Point3f normal=va^vb;
            Point3f med=(v0+v2)/2;
            Point3f dir=(v2-v0)^normal;
            dir.Normalize();
            float height=(v2-v0).Norm()/5;
            Point3f bestPoint=med+dir.Scale(height,height,height);
            betterloop[FinalIndices[k]].V()->P()=bestPoint;
            betterloop[FinalIndices[k]].F()->SetS();
        }
        /*tri::UpdateFlags<CMesh>::VertexClearS(m);
        for(int i=0;i<loop.size();i++)
            loop[i]->SetS();
        tri::UpdateSelection<CMesh>::FaceFromVertexLoose(m);
        //vcg::tri::UpdateColor<CMesh>::PerFaceConstant(m, vcg::Color4b::Red,true);
        //addloopToMesh(loop);
    }*/
    const float PI=3.14159265;
    int n=loop.size();
    int ind=-1;
    int i=0;
    bool stop=false;
    //tri::UpdateTopology<CMesh>::VertexFace(m);
    //tri::UpdateFlags<CMesh>::FaceClearS(m);
    tri::UpdateTopology<CMesh>::FaceFace(m);
    tri::UpdateTopology<CMesh>::VertexFace(m);
    tri::UpdateFlags<CMesh>::FaceBorderFromFF(m);
    //UpdateSelection<CMesh>::FaceFromBorderFlag(m);
    //tri::UpdateFlags<CMesh>::VertexBorderFromFaceAdj(m);
    CMesh::FacePointer fv;
    while(!stop && i<n){
        face::VFIterator<CFace> vfi(loop[i]);
        for(;!vfi.End();++vfi){
            CMesh::FacePointer f = vfi.F();
            if(f->IsD()) continue;
            for(int k=0;k<3;k++){
              if(face::IsBorder(*f,k) && (f->V(k)==loop[i]|| f->V(f->Next(k))==loop[i])){
                  fv=&*f;
                  f->SetS();
                  ind=k;
                  stop=true;
                  break;
              }
            }
            if(stop) break;
        }
        if(ind==-1)
          i++;
    }
    //vcg::tri::UpdateColor<CMesh>::PerFaceConstant(m, vcg::Color4b::Red,true);

    //It must exist some selected face with loop[0] as vertex
    //cout<<"ind "<<ind<<endl;
    //assert(ind!=-1);
    if(ind==-1){
        cout<<"Improvement loop fail"<<endl;
        return false;
    }
    std::vector<CMesh::VertexPointer> result;
    face::Pos<CMesh::FaceType> startPos(fv,ind);
    face::Pos<CMesh::FaceType> curPos=startPos;
    vector<face::Pos<CMesh::FaceType>> betterloop;
    result.push_back(curPos.V());
    betterloop.push_back(curPos);
    do
    {
      curPos.NextB();
      //curPos.F()->SetS();
      betterloop.push_back(curPos);
      result.push_back(curPos.V());
    }
    while(curPos!=startPos);
    result.pop_back();
    betterloop.pop_back();
    loop=result;
    n=loop.size();
    vector<CMesh::FacePointer> todelete;
    vector<int> WrongIndexes;
    for(int k=0;k<n;k++){
         vcg::face::Pos<CMesh::FaceType> jp=betterloop[k];
         int angle=0;
         do{
             Point3f a=(jp.VFlip()->P())-(jp.V()->P());
             jp.FlipE();
             Point3f b=(jp.VFlip()->P())-(jp.V()->P());
             angle+=vcg::Angle(a,b)*180/PI;
             jp.FlipF();
             if(jp.VFlip()==betterloop[(k+1)%n].V())
                   break;
         }while(jp.f!=betterloop[k].f);
         if(angle<150){
             WrongIndexes.push_back(k);
             jp=betterloop[k];
             do{
                 todelete.push_back(jp.F());
                 jp.F()->SetS();
                 jp.FlipE();
                 jp.FlipF();
                 if(jp.VFlip()==betterloop[(k+1)%n].V())
                       break;
             }while(jp.f!=betterloop[k].f);
         }

    }
    //Compute blocks to improve
    //This avoids smoothing corners vertices . Just stairlike sequences are smoothing
    vector<int> FinalIndices;
    for(int k=0;k<WrongIndexes.size();k++){
        int sup=(WrongIndexes[k]+4)%n;
        int inf=(WrongIndexes[k]-4+n)%n;
        int next=WrongIndexes[(k+1)%WrongIndexes.size()];
        int before=WrongIndexes[(k-1+WrongIndexes.size())%WrongIndexes.size()];
        bool yes=false;
        while(inf!=WrongIndexes[k]){
            if(inf==before){
                yes=true;
                break;
            }
            inf=(inf==(n-1))?0:inf+1;
        }
        while(sup!=WrongIndexes[k]){
            if(sup==next){
                yes=true;
                break;
            }
            sup=(sup==0)?(n-1):sup-1;
        }
        if(yes)
            FinalIndices.push_back(WrongIndexes[k]);
    }
    for(int k=0;k<FinalIndices.size();k++){
        Point3f v0=betterloop[FinalIndices[k]].VFlip()->P();
        Point3f v1=betterloop[FinalIndices[k]].V()->P();
        Point3f v2=betterloop[(FinalIndices[k]+1)%n].V()->P();
        Point3f va=v1-v0;
        Point3f vb=v2-v0;
        Point3f normal=va^vb;
        Point3f med=(v0+v2)/2;
        Point3f dir=(v2-v0)^normal;
        dir.Normalize();
        float height=(v2-v0).Norm()/10;
        Point3f bestPoint=med+dir.Scale(height,height,height);
        betterloop[FinalIndices[k]].V()->P()=bestPoint;
        betterloop[FinalIndices[k]].F()->SetS();
    }
    return true;

}
void intersectingCurves::fillloops(int index,vector<CMesh::VertexPointer> &loop, vector<CMesh::VertexPointer> &newloop,CMesh& m){
    //newloop.clear();
    //tri::UpdateFlags<CMesh>::FaceClearS(m);
    if(intersectingLoop.size()!=0){
        CMesh::VertexIterator vi=Allocator<CMesh>::AddVertices(m,intersectingLoop[index].size());
        for(int i=0;i<intersectingLoop[index].size();i++){
            vi->P()=intersectingLoop[index][i];
            newloop.push_back(&*vi);
            ++vi;
        }
        // Here we tessellate the strip
        // Very simple strategy: First choose the loop with less size
        // For each point in this find the nearest point in the another loop (with that we have edges)
        // At this point each segment in the first loop correspond to a polyline in the second loop
        // Choose the better point in this polyline(better aspect ratio)
        // For the rest of unmatched point in the second loop do the same
        if((loop.size()<newloop.size()))
            utilfillloops(loop,newloop,m,true);
        else
            utilfillloops(newloop,loop,m,false);
        //tri::Clean< CMesh >::RemoveNonManifoldFace(m);
        //tri::Clean< CMesh >::RemoveNonManifoldVertex(m);

    }
}
void intersectingCurves::utilfillloops(vector<CMesh::VertexPointer> &loop, vector<CMesh::VertexPointer> &newloop, CMesh &m,bool inmesh){
    std::unordered_map<int,int> mymap;
    int n=loop.size();
    int mn=newloop.size();
    for(int i=0;i<n;i++){
         float dist=std::numeric_limits<float>::infinity();
         int min=-1;
         int j=0;
         for(j=0;j<newloop.size();j++){
             float d=(newloop[j]->P()-loop[i]->P()).SquaredNorm();
             if(d<dist){
                 min=j;
                 dist=d;
             }
         }
         mymap.emplace(i,min);
         //utility::printPoint3f(loop[i]->P(),"point "+std::to_string(i));
         //cout<<"index after "<<tri::Index(m,loop[i])<<endl;
    }

    bool globalwalkingNext=true; // to decide what part of the circular curve we take
    vector<int> auxWalkingNext;
    auxWalkingNext.clear();
    for(int i=0;i<n;i++)
       auxWalkingNext.push_back(mymap[i]);

    globalwalkingNext=utility::moreIncrementsArray(auxWalkingNext);
    // for while we assuming that all loop is circular
    for(int i=0;i<n;i++){
        bool walkingNext=globalwalkingNext; // Even globally the curve has a circulation, some segments may have contray circulation
        int mapbegin=(mymap[i]<=mymap[(i+1)%n])?mymap[i]:mymap[(i+1)%n];
        int mapend=(mapbegin==mymap[i])?mymap[(i+1)%n]:mymap[i];
        if(mapend-mapbegin>0){
            if((mapend-mapbegin)==1){ // sure case with predominat circulation
                if(mymap[i]>mymap[(i+1)%n])
                    walkingNext=false;
                else
                    walkingNext=true;
            }
            // First attempt
            //int next=(mymap[i]==(mn-1))?0:mymap[i]+1;
            //int before=(mymap[i]==0)?mn-1:mymap[i]-1;            
            /*if(vcg::Distance(loop[(i+1)%n]->P(),newloop[next]->P())>vcg::Distance(loop[(i+1)%n]->P(),newloop[before]->P()) && vcg::Distance(loop[i]->P(),newloop[before2]->P())>vcg::Distance(loop[i]->P(),newloop[next2]->P()))
                walkingNext=false;*/

            // Second attempt
            /*int aux=mymap[i];
            int count=1;
            while(aux!=mymap[(i+1)%n]){
                 aux=(aux==(mn-1))?0:aux+1;
                 count++;
                 if(count>mn/2) break;
            }
            if(count>mn/2) walkingNext=false;*/

            //We have a segment in smaller loop. Find better third point in the bigger loop
            //float ARcurrent=utility::aspectRatioTriangle(loop[(i+1)%n]->P(),loop[i]->P(),newloop[mapbegin]->P());
            float MoreObtuseAngle=utility::cosinAngleFromPoints(loop[(i+1)%n]->P(),loop[i]->P(),newloop[mapbegin]->P(),2);
            int betterInd=mymap[i];
            int j=mymap[i];
            while(true){
                if(walkingNext) j=(j==(mn-1))?0:j+1;
                else j=(j==0)?mn-1:j-1;
                /*float AR=utility::aspectRatioTriangle(loop[(i+1)%n]->P(),loop[i]->P(),newloop[j]->P());
                if(AR<ARcurrent){
                    ARcurrent=AR;
                    betterInd=j;
                }*/
                float cosinAngle=utility::cosinAngleFromPoints(loop[(i+1)%n]->P(),loop[i]->P(),newloop[j]->P(),2);
                if(cosinAngle<MoreObtuseAngle){
                    MoreObtuseAngle=cosinAngle;
                    betterInd=j;
                }
                if(j==mymap[(i+1)%n]) break;
            }
            //Now we have probably two subpolylines in the bigger loop. Triangulate these;
            // The orientation if the faces are given by the mesh. So we take care with the orientation of the curve belonging to the mesh
            if(walkingNext){
                int k=mymap[i];
                while(true){
                    CMesh::FaceIterator f;
                    if(k==betterInd) break;
                    k=(k==(mn-1))?0:k+1;
                    //CMesh::FaceIterator f=tri::Allocator<CMesh>::AddFace(m,newloop[k],newloop[(k-1+mn)%mn],loop[i]);
                    //CMesh::FaceIterator f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k-1+mn)%mn]->P(),loop[i]->P());
                    if(inmesh){
                        f=tri::Allocator<CMesh>::AddFace(m,newloop[(k-1+mn)%mn]->P(),newloop[k]->P(),loop[i]->P());
                        //f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k-1+mn)%mn]->P(),loop[i]->P());
                    }
                    else{
                        //f=tri::Allocator<CMesh>::AddFace(m,newloop[(k-1+mn)%mn]->P(),newloop[k]->P(),loop[i]->P());
                        f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k-1+mn)%mn]->P(),loop[i]->P());
                    }
                    f->SetS();
                }
                k==betterInd;
                while(true){
                    CMesh::FaceIterator f;
                    if(k==mymap[(i+1)%n]) break;
                    k=(k==(mn-1))?0:k+1;
                    //CMesh::FaceIterator f=tri::Allocator<CMesh>::AddFace(m,newloop[k],newloop[(k-1+mn)%mn],loop[(i+1)%n]);
                    //CMesh::FaceIterator f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k-1+mn)%mn]->P(),loop[(i+1)%n]->P());
                    if(inmesh){
                        f=tri::Allocator<CMesh>::AddFace(m,newloop[(k-1+mn)%mn]->P(),newloop[k]->P(),loop[(i+1)%n]->P());
                        //f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k-1+mn)%mn]->P(),loop[(i+1)%n]->P());
                    }
                    else{
                        //f=tri::Allocator<CMesh>::AddFace(m,newloop[(k-1+mn)%mn]->P(),newloop[k]->P(),loop[(i+1)%n]->P());
                        f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k-1+mn)%mn]->P(),loop[(i+1)%n]->P());
                    }
                    f->SetS();
                }
            }
            else{
                int k=mymap[i];
                /*cout<<"SESION"<<endl;
                cout<<"k: "<<k<<endl;
                cout<<"next "<<mymap[(i+1)%n]<<endl;
                cout<<"better: "<<betterInd<<endl;
                cout<<"mn: "<<mn<<endl;*/

                while(true){
                    CMesh::FaceIterator f;
                    if(k==betterInd) break;
                    k=(k==0)?(mn-1):k-1;
                    //CMesh::FaceIterator f=tri::Allocator<CMesh>::AddFace(m,newloop[k],newloop[(k+1+mn)%mn],loop[i]);
                    //CMesh::FaceIterator f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k+1+mn)%mn]->P(),loop[i]->P());
                    if(inmesh){
                        f=tri::Allocator<CMesh>::AddFace(m,newloop[(k+1+mn)%mn]->P(),newloop[k]->P(),loop[i]->P());
                        //f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k+1+mn)%mn]->P(),loop[i]->P());
                    }
                    else{
                        f=tri::Allocator<CMesh>::AddFace(m,newloop[(k+1+mn)%mn]->P(),newloop[k]->P(),loop[i]->P());
                        //f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k+1+mn)%mn]->P(),loop[i]->P());
                    }
                    f->SetS();
                }
                k==betterInd;
                while(true){
                    CMesh::FaceIterator f;
                    if(k==mymap[(i+1)%n]) break;
                    k=(k==0)?(mn-1):k-1;
                    //CMesh::FaceIterator f=tri::Allocator<CMesh>::AddFace(m,newloop[k],newloop[(k+1+mn)% mn],loop[(i+1)%n]);
                    //CMesh::FaceIterator f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k+1+mn)% mn]->P(),loop[(i+1)%n]->P());
                    if(inmesh){
                        f=tri::Allocator<CMesh>::AddFace(m,newloop[(k+1+mn)%mn]->P(),newloop[k]->P(),loop[(i+1)%n]->P());
                        //f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k+1+mn)%mn]->P(),loop[(i+1)%n]->P());
                    }
                    else{
                        f=tri::Allocator<CMesh>::AddFace(m,newloop[(k+1+mn)%mn]->P(),newloop[k]->P(),loop[(i+1)%n]->P());
                        //f=tri::Allocator<CMesh>::AddFace(m,newloop[k]->P(),newloop[(k+1+mn)%mn]->P(),loop[(i+1)%n]->P());
                    }
                    f->SetS();
                }
            }

            CMesh::FaceIterator fi;
            if(inmesh)
                fi=tri::Allocator<CMesh>::AddFace(m,loop[(i+1)%n]->P(),loop[i]->P(),newloop[betterInd]->P());
            else{
                if(walkingNext) // if the intersection curve and the border curve have the same circulation
                   fi=tri::Allocator<CMesh>::AddFace(m,loop[i]->P(),loop[(i+1)%n]->P(),newloop[betterInd]->P());
                else
                   fi=tri::Allocator<CMesh>::AddFace(m,loop[(i+1)%n]->P(),loop[i]->P(),newloop[betterInd]->P());

            }
            //CMesh::FaceIterator fi=tri::Allocator<CMesh>::AddFace(m,loop[(i+1)%n]->P(),loop[i]->P(),newloop[betterInd]->P());
            //CMesh::FaceIterator fi=tri::Allocator<CMesh>::AddFace(m,loop[i]->P(),loop[(i+1)%n]->P(),newloop[betterInd]->P());
            fi->SetS();

        }
        else{ // if associated point is the same we already have a face
            CMesh::FaceIterator f;
            if(inmesh)
                f=tri::Allocator<CMesh>::AddFace(m,loop[(i+1)%n]->P(),loop[i]->P(),newloop[mymap[i]]->P());
            else{
               //Point3f dircurve=loop[(i+1)%n]->P()-loop[i]->P();
               /*Point3f normalthistriangle=(loop[i]->P()-newloop[mymap[i]]->P())^(loop[(i+1)%n]->P()-newloop[mymap[i]]->P());
               Point3f normalinmesh=newloop[mymap[i]]->N();
               if((normalthistriangle).dot(normalinmesh)>=0)
                   f=tri::Allocator<CMesh>::AddFace(m,loop[i]->P(),loop[(i+1)%n]->P(),newloop[mymap[i]]->P());
               else
                   f=tri::Allocator<CMesh>::AddFace(m,loop[(i+1)%n]->P(),loop[i]->P(),newloop[mymap[i]]->P());*/
               if(walkingNext) // if the intersection curve and the border curve have the same circulation
                  f=tri::Allocator<CMesh>::AddFace(m,loop[i]->P(),loop[(i+1)%n]->P(),newloop[mymap[i]]->P());
               else
                  f=tri::Allocator<CMesh>::AddFace(m,loop[(i+1)%n]->P(),loop[i]->P(),newloop[mymap[i]]->P());
            }
            f->SetS();
            //CMesh::FaceIterator f=tri::Allocator<CMesh>::AddFace(m,loop[(i+1)%n]->P(),loop[i]->P(),newloop[mymap[i]]->P());
            //CMesh::FaceIterator f=tri::Allocator<CMesh>::AddFace(m,loop[i]->P(),loop[(i+1)%n]->P(),newloop[mymap[i]]->P());

        }
    }
    bool ioriented,isorientable;
    tri::Clean< CMesh >::RemoveDuplicateVertex(m);
    UpdateTopology<CMesh>::FaceFace(m);
    //tri::Clean< CMesh >::OrientCoherentlyMesh(m,ioriented,isorientable);
    //tri::Clean< CMesh >::RemoveDegenerateVertex(m);
    //tri::Clean< CMesh >::RemoveUnreferencedVertex(m);
    //UpdateNormal<CMesh>::PerVertex(m);
    //UpdateNormal<CMesh>::PerFace(m);
}
void intersectingCurves::doIsoRemeshing(CMesh &m,bool selectedFaces){

    tri::Clean< CMesh >::RemoveDuplicateVertex(m);
    tri::Clean< CMesh >::RemoveNonManifoldFace(m);
    tri::Clean< CMesh >::RemoveUnreferencedVertex(m);
    tri::Allocator< CMesh >::CompactEveryVector(m);
    UpdateTopology<CMesh>::FaceFace(m);
    UpdateTopology<CMesh>::VertexFace(m);
    UpdateTopology<CMesh>::AllocateEdge(m);
    tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFaceNormalized(m);
    tri::UpdateBounding<CMesh>::Box(m);
    tri::MeshAssert<CMesh>::FFTwoManifoldEdge(m);
    tri::MeshAssert<CMesh>::NoUnreferencedVertex(m);
    UpdateSelection<CMesh>::FaceDilate(m);
    IsotropicRemeshing<CMesh>::Params params;
    params.SetTargetLen(thrMax/2);
    params.SetFeatureAngleDeg(90);
    params.splitFlag=false;
    params.iter=2;
    params.selectedOnly=selectedFaces;
    IsotropicRemeshing<CMesh>::Do(m,params);
}
void intersectingCurves::deleteVerticesInIntersectionCurve(int indexCurve, set<int> &indexesToDelete){
    int nsize=intersectingLoop[indexCurve].size();
    vector<CMesh::VertexPointer> newintersectingLoopPointerL;
    vector<CMesh::VertexPointer> newintersectingLoopPointerR;
    vector<Point3f> newintersectingLoop;
    for(int j=0;j<nsize;j++){
        if(indexesToDelete.find(j)==indexesToDelete.end()){
            newintersectingLoop.push_back(intersectingLoop[indexCurve][j]);
            newintersectingLoopPointerL.push_back(intersectingLoopPointers[indexCurve].first[j]);
            newintersectingLoopPointerR.push_back(intersectingLoopPointers[indexCurve].second[j]);
        }
    }
    intersectingLoop[indexCurve]=newintersectingLoop;
    intersectingLoopPointers[indexCurve]=make_pair(newintersectingLoopPointerL,newintersectingLoopPointerR);
}
