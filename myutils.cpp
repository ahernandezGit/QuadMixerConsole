#include "myutils.h"

namespace utility{

    // This predicate need the vertex border computing before
    // for example with UpdateFlags<MyMesh>::VertexBorderFromFaceAdj(mymesh)

    class EdgeSplitLenPred
    {
      public:
          float squaredlengthThr;
          EdgeSplitLenPred(){}
          EdgeSplitLenPred(float threshold) {setThr(threshold);}
          void setThr(float threshold) {squaredlengthThr = threshold*threshold; }
          bool operator()(PosType &ep)
          {
              if(ep.V()->IsB() || ep.VFlip()->IsB()){
                  if(SquaredDistance(ep.V()->P(), ep.VFlip()->P()) > squaredlengthThr)
                      return true;
                  return false;
              }
              else
                  return false;
          }
    };
    class EdgeSplitForEqualLoopsPred
    {
      public:
          float squaredlengthThr;
          EdgeSplitForEqualLoopsPred(){}
          EdgeSplitForEqualLoopsPred(float threshold) {setThr(threshold);}
          void setThr(float threshold) {squaredlengthThr = threshold*threshold; }
          bool operator()(PosType &ep)
          {
              if((ep.V()->IsB() && !ep.VFlip()->IsB())||(!ep.V()->IsB() && ep.VFlip()->IsB())){
                  if(SquaredDistance(ep.V()->P(), ep.VFlip()->P()) > squaredlengthThr)
                      return true;
                  return false;
              }
              else
                  return false;
          }
    };
    class EdgeSplitFixedEdge
    {
      public:
          PosType edge;
          EdgeSplitFixedEdge(){edge.SetNull();}
          EdgeSplitFixedEdge(PosType &iedge) {setThr(iedge);}
          void setThr(PosType &iedge) {edge=iedge;}
          bool operator()(PosType &ep)
          {
              PosType opp=edge;
              opp.VFlip();
              if(ep==edge || ep==opp)
                  return true;
              else
                  return false;
          }
    };
    void saveAsCurvePly(vector<face::Pos<CMesh::FaceType>> &curve,string  filepath){
        CMesh tcurve;
        for(size_t i=0;i<curve.size();i++)
            tri::Allocator<CMesh>::AddEdge(tcurve,curve[i].VFlip()->P(),curve[i].V()->P());
        tri::Clean<CMesh>::RemoveDuplicateVertex(tcurve);
        tri::Clean<CMesh>::RemoveUnreferencedVertex(tcurve);
        tri::Allocator<CMesh>::CompactEveryVector(tcurve);
        tri::UpdateTopology<CMesh>::VertexEdge(tcurve);
        tri::UpdateTopology<CMesh>::EdgeEdge(tcurve);
        tri::UpdateColor<CMesh>::PerVertexConstant(tcurve,Color4b::Black);
        tri::io::ExporterPLY<CMesh>::Save(tcurve,filepath.c_str(),tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY,false);
    }
    void saveAsCurvePly(vector<CMesh::CoordType> &curve,string  filepath,bool iscircular){
        CMesh tcurve;
        for(size_t i=0;i<curve.size()-1;i++)
            tri::Allocator<CMesh>::AddEdge(tcurve,curve[i],curve[i+1]);
        if(iscircular)
            tri::Allocator<CMesh>::AddEdge(tcurve,curve.back(),curve[0]);
        tri::Clean<CMesh>::RemoveDuplicateVertex(tcurve);
        tri::Clean<CMesh>::RemoveUnreferencedVertex(tcurve);
        tri::Allocator<CMesh>::CompactEveryVector(tcurve);
        tri::UpdateTopology<CMesh>::VertexEdge(tcurve);
        tri::UpdateTopology<CMesh>::EdgeEdge(tcurve);
        tri::UpdateColor<CMesh>::PerVertexConstant(tcurve,Color4b::Black);
        tri::io::ExporterPLY<CMesh>::Save(tcurve,filepath.c_str(),tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY,false);

    }
    void saveMesh(CMesh &m,string filepath){
        tri::io::ExporterPLY<CMesh>::Save(m,filepath.c_str(),tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
    }
    void ComputeRealRootsDegree2(double &root1, double &root2, double a, double b, double c,int & valid){
        double discriminant = b*b - 4*a*c;
        valid=1;
        if (discriminant > 0) {
                root1 = (-b + sqrt(discriminant)) / (2*a);
                root2 = (-b - sqrt(discriminant)) / (2*a);
                cout << "Roots are real and different." << endl;
                cout << "x1 = " << root1 << endl;
                cout << "x2 = " << root2 << endl;
        }
        else if (discriminant == 0) {
            cout << "Roots are real and same." << endl;
            root1 = (-b + sqrt(discriminant)) / (2*a);
            root2=root1;
            cout << "x1 = x2 =" << root1 << endl;
        }
        else {
            double realPart = -b/(2*a);
            double imaginaryPart =sqrt(-discriminant)/(2*a);
            cout << "Roots are complex and different."  << endl;
            root1=0.0;
            root2=0.0;
            valid=0;
        }
    }
    bool ispermuted(const vector<int>& list1,const vector<int>& list2, vector<pair<int,int>>& correspondence ){
        int n=list1.size();
        assert(list2.size()==n);
        correspondence.resize(n);
        int target=list1[0];
        for(int i=0;i<n;i++){
            if(list2[i]==target){
                bool yes=true;
                for(int j=0;j<n;j++){
                    correspondence[j]=make_pair(j,(i+j)%n);
                    yes=yes && (list2[(i+j)%n]==list1[j]);
                    if(yes==false)
                        break;
                }
                if(yes==true)
                   return true;
                yes=true;
                for(int j=0;j<n;j++){
                    correspondence[j]=make_pair(j,(i-j+n)%n);
                    yes=yes && (list2[(i-j+n)%n]==list1[j]);
                    if(yes==false)
                        break;
                }
                if(yes==true)
                   return true;
            }
        }
        return false;
    }
    void correctCornersFromCurve(vector<Point3f>& curve, bool iscircular){
       int n=curve.size();
       vector<Point3f> ncurve;
       if(n<3) return;

        //mark indexes to smooth
        set<int> indexesToSmooth;
        for(int i=0;i<n;i++){
            if(!iscircular && (i==0 || i==(n-1))){
                continue;
            }
            Point3f a=curve[(i-1+n)%n]-curve[i];
            Point3f b=curve[(i+1)%n]-curve[i];
            float angle=vcg::Angle(a,b);
            if(angle<40*3.1416/180){
                indexesToSmooth.insert(i);
                //cout<<"INDEX TO SMOOTH"<<endl;
            }
        }
        int totalindexes=indexesToSmooth.size();
        //vector<vector<int>> blocks;
        //groupIndexes(blocks,indexesToSmooth,n,iscircular);

        // for curves with insufficient amount  of vertices
        if(n-totalindexes<3) return;
        for(int i=0;i<n;i++){
            if(indexesToSmooth.find(i)==indexesToSmooth.end())
                ncurve.push_back(curve[i]);
        }
        curve=ncurve;
        /*for(std::set<int>::iterator it=indexesToSmooth.begin();it!=indexesToSmooth.end();it++)
            curve[*it]=(curve[(*it-1+n)%n]+curve[(*it+1)%n])/2;
        for(int i=0;i<blocks.size();i++){
            int nsize=blocks[i].size();
            int before=(blocks[i][0]-1+n)%n;
            int next=(blocks[i][nsize-1]+1)%n;
            Eigen::Vector3f p0(curve[before].X(),curve[before].Y(),curve[before].Z());
            Eigen::Vector3f p1(curve[next].X(),curve[next].Y(),curve[next].Z());
            Eigen::ParametrizedLine<float,3> line=Eigen::ParametrizedLine<float,3>::Through(p0,p1);
            for(int j=0;j<nsize;j++){
                Point3f tem=curve[blocks[i][j]];
                Eigen::Vector3f projection=line.projection(Eigen::Vector3f(tem.X(),tem.Y(),tem.Z()));
                curve[blocks[i][j]].X()=projection(0);
                curve[blocks[i][j]].Y()=projection(1);
                curve[blocks[i][j]].Z()=projection(2);
            }
        }*/

    }
    void groupIndexes(vector<vector<int>> &blocks, set<int> &markindexes, int nd,bool iscircular){
        std::set<int>::iterator itc=markindexes.begin();
        int count=0;
        while(count<1000){
            if(!iscircular && (*itc==0 || *itc==(nd-1))){
                itc++;
                continue;
            }
            int before=(*itc-1+nd)%nd;
            int next=(*itc+1)%nd;
            std::set<int>::iterator itbefore=markindexes.find(before);
            std::set<int>::iterator itnext=markindexes.find(next);
            vector<int> block;
            while(itbefore!=markindexes.end()){
                block.push_back(before);
                before=(before-1+nd)%nd;
                markindexes.erase(itbefore);
                itbefore=markindexes.find(before);
            }
            if(!block.empty())
                std::reverse(block.begin(),block.end());
            block.push_back(*itc);
            while(itnext!=markindexes.end()){
                block.push_back(next);
                next=(next+1)%nd;
                markindexes.erase(itnext);
                itnext=markindexes.find(next);

            }
            if(block.size()>1){
                blocks.push_back(block);
                markindexes.erase(itc);
            }
            if(!markindexes.empty()){
                /*cout<<"one pass "<<count <<endl;
                for (std::set<int>::iterator it=markindexes.begin(); it!=markindexes.end(); ++it){
                       cout<<*it<<" ";
                }
                cout<<endl;*/
                if(block.size()>1){
                    /*cout<<" encontrei grupos"<<endl;
                    for(int r=0;r<block.size();r++)
                        cout<<block[r]<<" ";
                    cout<<endl;*/
                    itc=markindexes.begin();
                }
                else{
                    //cout<<" no encontrei grupos"<<endl;
                    itc++;
                    //cout<<"novo iterator "<<*itc<<endl;
                    if(itc==markindexes.end())
                        break;
                }
            }
            else
                break;
            count++;
        }
        if(count==1000)
            cout<<"Maximun number de iteration group indexes"<<endl;
    }
    bool moreIncrementsArray(vector<int> &vec){
        int nsize=vec.size();
        int nincrements=0;
        //cout<<"Array "<<endl;
        for(int i=0;i<nsize-1;i++){
            //cout<<vec[i]<<" ";
            if(vec[i]<vec[i+1])
                nincrements++;
        }
        //cout<<vec[nsize-1]<<" ]"<<endl;
        if(nincrements>nsize/2)
            return true;
        else
            return false;
    }
    bool checkVertexIsCrease(const CMesh::VertexType* v)
    {
        vcg::edge::VEIterator<CMesh::EdgeType> vei(v);
        for (;!vei.End();++vei)
        {
            CMesh::FacePointer f=vei.E()->EFp();
            int z=vei.E()->EFi();
            if(f->IsF(z)) return true;
        }
        return false;
    }
    bool checkVertexIsCorner(const CMesh::VertexType* v,float anglethr)
    {
        assert(anglethr>=0);
        vcg::edge::VEIterator<CMesh::EdgeType> vei(v);
        for (;!vei.End();++vei)
        {
            CMesh::FacePointer f=vei.E()->EFp();
            int z=vei.E()->EFi();
            //If the edge is not border she has a adjacent face
            if(!f->IsB(z)){
                CMesh::ScalarType angle= DihedralAngleRad(*f,z);
                angle=(angle<0)?-angle:angle;
                if(angle>(anglethr-0.1) && angle<(anglethr+0.1))
                    return true;
            }
        }
        return false;
    }

    float edgeAverage(CMesh& m, bool selectedFaces)
    {
        CMesh::FaceIterator fi;
        float average=0;
        int n=0;
        for (fi = m.face.begin();  fi!=m.face.end(); fi++)
        {
            if(!(*fi).IsD() && (!selectedFaces || (*fi).IsS()))
            {
                float a=vcg::Distance(fi->P0(0),fi->P1(0));
                float b=vcg::Distance(fi->P0(1),fi->P1(1));
                float c=vcg::Distance(fi->P0(2),fi->P1(2));
                average+=((a+b)+c)/3;
                n++;
            }
        }
        if(n!=0) average=average/n;
        else{
            cout<<"non selected faces"<<endl;
        }
        return average;
    }
    bool someEdgeSmall(CMesh& m,float threshold, bool selectedFaces){
        CMesh::FaceIterator fi;
        int i=0;
        for (fi = m.face.begin();  fi!=m.face.end(); fi++,i++)
        {
            if(!(*fi).IsD() && (!selectedFaces || (*fi).IsS()))
            {
                for(int k=0;k<3;k++){
                    if(fi->IsB(k)){
                        float le=vcg::Distance((*fi).V(k)->P(),(*fi).V((*fi).Next(k))->P());
                        if(le<threshold) return true;
                    }
                }
            }
        }
        //cout<<"no edge small"<<endl;
        return false;
    }
    bool someEdgeBigger(CMesh& m,float threshold, bool selectedFaces){
        CMesh::FaceIterator fi;
        int i=0;
        for (fi = m.face.begin();  fi!=m.face.end(); fi++,i++)
        {
            if(!(*fi).IsD() && (!selectedFaces || (*fi).IsS()))
            {
                for(int k=0;k<3;k++){
                    if(fi->IsB(k)){
                        float le=vcg::Distance((*fi).V(k)->P(),(*fi).V((*fi).Next(k))->P());
                        if(le>threshold) return true;
                    }
                }
            }
        }
        //cout<<"no edge bigger"<<endl;
        return false;
    }
    vector<CMesh::CoordType> getRegularDivisionsAsBarycentricFromPos(face::Pos<CMesh::FaceType> &pos,int subdivisionsPos){
        assert(subdivisionsPos>=0);
        std::unordered_map<int,CMesh::ScalarType> cornersToCoordinate; // triangle index to barycentric coordinates
        std::unordered_map<CMesh::VertexPointer,int> pointerToCorners; //vertex to triangle index
        vector<CMesh::CoordType> result;
        result.clear();
        assert(pos.F()->VN()==3);
        int indexOut=-1;
        for(int i=0;i<3;i++){
            pointerToCorners.emplace(pos.F()->V(i),i);
            if(pos.F()->V(i)!=pos.V() && pos.F()->V(i)!=pos.VFlip()){
                indexOut=i;
            }
        }
        assert(indexOut!=-1);
        cornersToCoordinate.emplace(indexOut,0);
        cornersToCoordinate.emplace(pointerToCorners[pos.VFlip()],0);
        cornersToCoordinate.emplace(pointerToCorners[pos.V()],1);
        vector<pair<CMesh::ScalarType,CMesh::ScalarType>> coordinates;

        int numDivisions=subdivisionsPos+1;
        for(int i=1;i<numDivisions;i++){
           coordinates.push_back(make_pair(i/numDivisions,1-i/numDivisions));
        }
        for(int i=0;i<coordinates.size();i++){
           CMesh::CoordType final_coordinates;
           cornersToCoordinate[pointerToCorners[pos.VFlip()]]=coordinates[i].first;
           cornersToCoordinate[pointerToCorners[pos.V()]]=coordinates[i].second;
           final_coordinates.X()=cornersToCoordinate[0];
           final_coordinates.Y()=cornersToCoordinate[1];
           final_coordinates.Z()=cornersToCoordinate[2];
           printPoint3f(final_coordinates," barycentric coordinates");
           result.push_back(final_coordinates);
        }
        return result;

    }
    void splitSegment(const Point3f& A,const Point3f& B, vector<Point3f>& splitseg, int numberDivisions){
        vector<Point3f> result;
        splitseg.clear();
        if(numberDivisions<2)
            return;
        Point3f direction=B-A;
        direction.normalized();
        float len=vcg::Distance(A,B);
        float step=len/numberDivisions;
        //result.push_back(A);
        for(int i=1;i<numberDivisions;i++){
            result.push_back(A+direction*(i*step));
        }
        //result.push_back(B);
        splitseg=result;
    }
    void splitPosList_onPolygon(vector<Point3f> & polygon,vector<size_t> & poslist,vector<int> &timestosplit, set<int> featureOnPolygon,vector<Point3f> &result){
        assert(poslist.size()==timestosplit.size());
        size_t n=polygon.size();
        for(size_t i=0;i<poslist.size();i++){
            size_t ii=poslist[i];
            vector<Point3f> newseg;
            splitSegment(polygon[ii],polygon[(ii+1)%n],newseg,timestosplit[i]+1);
            result=polygon;
            vector<Point3f>::iterator it=polygon.begin();
            std::advance(it,(ii+1)%n);
            polygon.insert(it,newseg.begin(),newseg.end());
            //correcting feature indexes
            set<int> newfeature;
            for(const int & ind:featureOnPolygon){
                if(ind<=ii)
                    newfeature.insert(ind);
                else
                    newfeature.insert(ind+newseg.size());
            }
            featureOnPolygon=newfeature;
        }
        result=polygon;
    }
    void printPoint3f(const Point3f& p, string message)
    {
        cout<<" "<<message<<" : ["<<p.X()<<" , "<<p.Y()<<" , "<<p.Z()<<" ]"<<endl;
    }
    float aspectRatioTriangleFromFace(CMesh::FacePointer f){
        //assert(f->VN()!=3);
        if(!(*f).IsD())
        {
            float a=vcg::Distance((*f).P0(0),(*f).P1(0));
            float b=vcg::Distance((*f).P1(0),(*f).P2(0));
            float c=vcg::Distance((*f).P0(0),(*f).P2(0));
            float s=((a+b)+c)/2;
            float ar=a*b*c/(8*(s-a)*(s-b)*(s-c));
            return ar;
        }
        return 0;
    }
    float aspectRatioTriangle(const Point3f& v0,const Point3f& v1,const Point3f& v2){

            float a=vcg::Distance(v0,v1);
            float b=vcg::Distance(v0,v2);
            float c=vcg::Distance(v1,v2);
            float s=((a+b)+c)/2;
            float ar=a*b*c/(8*(s-a)*(s-b)*(s-c));
            return ar;

    }
    float signedAreaTriangle(const Point3f& v0,const Point3f& v1,const Point3f& v2){
        Point3f va=v1-v0;
        Point3f vb=v2-v1;
        Point3f result=va^vb;
        if (result.Z()>0) return 1;
        if (result.Z()==0) return 0;
        if (result.Z()<0) return -1;
    }
    void splitMeshStrip(CMesh& m, float threshold,bool selectedFaces)
    {
        //vcg::tri::UpdateTopology<CMesh>::FaceFace(m);
        vcg::tri::UpdateFlags<CMesh>::VertexBorderFromFaceAdj(m);
        EdgeSplitLenPred ep(threshold);
        tri::MidPoint<CMesh> midFunctor(&m);
        //Refine(threshold,selectedFaces);
        tri::RefineE(m,midFunctor,ep,selectedFaces);
    }
    void splitMeshStripNonBorder(CMesh& m, float threshold,bool selectedFaces)
    {
        //vcg::tri::UpdateTopology<CMesh>::FaceFace(m);
        vcg::tri::UpdateFlags<CMesh>::VertexBorderFromFaceAdj(m);
        EdgeSplitForEqualLoopsPred ep(threshold);
        tri::MidPoint<CMesh> midFunctor(&m);
        //Refine(threshold,selectedFaces);
        tri::RefineE(m,midFunctor,ep,selectedFaces);
    }

    void splitMesh(CMesh &m, float threshold, bool selectedFaces){
        tri::Refine<CMesh, tri::MidPoint<CMesh> >(m,tri::MidPoint<CMesh>(&m),threshold,selectedFaces);
    }
    void splitMeshP(PMesh &m, float threshold, bool selectedFaces){
        tri::Refine<PMesh, tri::MidPoint<PMesh> >(m,tri::MidPoint<PMesh>(&m),threshold,selectedFaces);
    }
    void easyProjectToMesh(const std::vector<CMesh::FacePointer>& faces,Point3f& p){
        std::vector<CMesh::FacePointer>::const_iterator fi;
        Point3f result;
        float mindist=std::numeric_limits<float>::infinity();
        for(fi=faces.begin();fi!=faces.end();fi++){
              Line3<float> line(p,(*fi)->N());
              float t,u,v;
              if(vcg::IntersectionLineTriangle(line,(*fi)->P0(0),(*fi)->P1(0),(*fi)->P2(0),t,u,v)){
                  Point3f dir=line.Direction();
                  Point3f projPoint=line.Origin()+dir.Scale(t,t,t);
                  float ndist=vcg::Distance(p,projPoint);
                  if(ndist<mindist){
                      result=projPoint;
                      mindist=ndist;
                  }
              }
        }
        p=result;
    }
    void easyProjectToMesh(const vector<vector<Point3f>>& faces,Point3f& p){
        vector<vector<Point3f>>::const_iterator fi;
        Point3f point=p;
        Point3f result=p;
        float mindist=std::numeric_limits<float>::infinity();
        for(fi=faces.begin();fi!=faces.end();fi++){
              Point3f normal=(*fi)[3];
              Point3f totest=point-(*fi)[0];
              if(totest.dot(normal)>=0){
                  Line3<float> line(point,normal);
                  float t,u,v;
                  if(vcg::IntersectionLineTriangle(line,(*fi)[0],(*fi)[1],(*fi)[2],t,u,v)){
                      Point3f dir=line.Direction();
                      Point3f projPoint=line.Origin()+dir.Scale(t,t,t);
                      float ndist=vcg::Distance(point,projPoint);
                      if(ndist<mindist){
                          result=projPoint;
                          mindist=ndist;
                      }
                  }
              }
        }
        p=result;
    }

    bool isCounterClockWise(vector<vcg::face::Pos<CMesh::FaceType> > &polygon){
        // if is not closed then is true
        int n=polygon.size();
        assert(polygon[n-1].V()==polygon[0].VFlip());
            float s=0;
            for(int i=0;i<n;i++){
                float x2=polygon[(i+1)%n].V()->P().X();
                float x1=polygon[i].V()->P().X();
                float y2=polygon[(i+1)%n].V()->P().Y();
                float y1=polygon[i].V()->P().Y();
                s+=(x2-x1)*(y1+y2);
            }
            if(s>0)
                return false;
            else
                return true;

    }
    bool isCounterClockWise(vector<Point3f> &polygon){
        int n=polygon.size();
        float s=0;
        for(int i=0;i<n;i++){
            float x2=polygon[(i+1)%n].X();
            float x1=polygon[i].X();
            float y2=polygon[(i+1)%n].Y();
            float y1=polygon[i].Y();
            s+=(x2-x1)*(y1+y2);
        }
        if(s>0)
            return false;
        else
            return true;
    }
    bool isCounterClockWise(vector<Point2f> &polygon){
        int n=polygon.size();
        float s=0;
        for(int i=0;i<n;i++){
            float x2=polygon[(i+1)%n].X();
            float x1=polygon[i].X();
            float y2=polygon[(i+1)%n].Y();
            float y1=polygon[i].Y();
            s+=(x2-x1)*(y1+y2);
        }
        if(s>0)
            return false;
        else
            return true;
    }
    void collapseMesh(CMesh &m, float threshold, bool selectedFaces)
    {
        // Test for learning about edge implementation and Adjacency
        /*
        CMesh::EdgeIterator ei=m.edge.begin();
        CMesh::FacePointer f=(*ei).EFp();
        int z=(*ei).EFi();
        CMesh::VertexPointer v=f->V0(z);
        CMesh::VertexPointer e0=ei->V(0);
        CMesh::VertexPointer e1=ei->V(1);
        CMesh::VertexPointer epointer=ei->V(1);
        cout<<"face V0=["<<(v->P()).X()<<" , "<<(v->P()).Y()<<" , "<<(v->P()).Z()<<" ]"<<endl;
        v=f->V1(z);
        cout<<"face V1=["<<(v->P()).X()<<" , "<<(v->P()).Y()<<" , "<<(v->P()).Z()<<" ]"<<endl;
        v=f->V2(z);
        cout<<"face V2=["<<(v->P()).X()<<" , "<<(v->P()).Y()<<" , "<<(v->P()).Z()<<" ]"<<endl;
        // e0==V0
        // e1==V1
        cout<<"e0=["<<(e0->P()).X()<<" , "<<(e0->P()).Y()<<" , "<<(e0->P()).Z()<<" ]"<<endl;
        cout<<"e1=["<<(e1->P()).X()<<" , "<<(e1->P()).Y()<<" , "<<(e1->P()).Z()<<" ]"<<endl;
        */
        //Test for learning about edge pos implementation and Adjacency
        /*vcg::edge::Pos<CMesh::EdgeType> epos(&*ei,0);
        CMesh::VertexType* v=epos.V();
        cout<<"edge vertex V=["<<(v->P()).X()<<" , "<<(v->P()).Y()<<" , "<<(v->P()).Z()<<" ]"<<endl;
        epos.NextE();
        v=epos.V();
        cout<<"edge vertex V=["<<(v->P()).X()<<" , "<<(v->P()).Y()<<" , "<<(v->P()).Z()<<" ]"<<endl;
        */

        //A collapse step
        /*vcg::tri::BasicVertexPair<CMesh::VertexType> vertexpair((*fi).V0(0),(*fi).V1(0));
        vcg::Point3f p=((*fi).P0(0)+(*fi).P1(0))/2;
        vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>> eco;
        eco.Do(m,vertexpair,p);
        */

        const double PI  =3.141592653589793238463;
        //Preprocessing
        // Mark features vertexs and edges ( borders or crease ).
        vcg::tri::UpdateFlags<CMesh>::Clear(m);
        vcg::tri::Clean<CMesh>::RemoveNonManifoldFace(m);
        vcg::tri::UpdateFlags<CMesh>::FaceBorderFromFF(m);
        vcg::tri::UpdateFlags<CMesh>::VertexBorderFromFaceAdj(m);
        //vcg::tri::UpdateFlags<CMesh>::FaceFauxCrease(m,1.74); //marks as fauxedges all the non crease edges angle<100
        //vcg::tri::UpdateFlags<CMesh>::FaceFauxBorder(m); //marks as fauxedges all the non border edge
        CMesh::EdgeIterator ei;
        for (ei = m.edge.begin();  ei!=m.edge.end(); ei++)
        {
            //Compute edge length of the edge
            float edgelen=vcg::Distance((*ei).P0(0),(*ei).P1(0));

            cout<<"el "<<edgelen<<endl;
            //Compute the opening angle at either opposing vertex

            // First face contain edge
            CMesh::FacePointer f=(*ei).EFp(); //pointer to face contain the edge
            int ez=(*ei).EFi(); // index of the edge in the found face
            CMesh::VertexPointer v0=f->V0(ez);  // endpoint of the edge
            CMesh::VertexPointer v1=f->V1(ez);  // endpoint of the edge
            CMesh::VertexPointer v2=f->V2(ez);
            CMesh::CoordType e1=v0->P()-v2->P();
            CMesh::CoordType e2=v1->P()-v2->P();
            CMesh::ScalarType angle2=0,angle1=vcg::Angle(e1,e2);

            bool angleflag=angle1<(PI/12);

            //cout<<"edge0 =["<<(ei->V(0)->P()).X()<<" , "<<(ei->V(0)->P()).Y()<<" , "<<(ei->V(0)->P()).Z()<<" ]"<<endl;
            //cout<<"edge1 =["<<(ei->V(1)->P()).X()<<" , "<<(ei->V(1)->P()).Y()<<" , "<<(ei->V(1)->P()).Z()<<" ]"<<endl;
            //cout<<"edge V0=["<<(v0->P()).X()<<" , "<<(v0->P()).Y()<<" , "<<(v0->P()).Z()<<" ]"<<endl;
            //cout<<"edge V1=["<<(v1->P()).X()<<" , "<<(v1->P()).Y()<<" , "<<(v1->P()).Z()<<" ]"<<endl;

            //If the edge is non border there is other face contain it
            if(!face::IsBorder((*f), ez)){
                CMesh::FacePointer ff=f->FFp(ez);
                int fz=ff->FFi(ez);
                CMesh::VertexPointer fv0=ff->V0(fz);  // endpoint of the edge
                CMesh::VertexPointer fv1=ff->V1(fz);  // endpoint of the edge
                CMesh::VertexPointer fv2=ff->V2(fz);
                CMesh::CoordType fe1,fe2;

                // This is awful but convention edge i=0..2 is between vertex i and (i+1)%3 didn't work
                // (fv0==v0 && fv1==v1) || (fv0==v1 && fv1==v0) is not always true
                // Find correct opposite vertex in this face
                if((fv0==v0 && fv1==v1) || (fv0==v1 && fv1==v0)){
                    fe1=fv0->P()-fv2->P();
                    fe2=fv1->P()-fv2->P();

                }
                else if((fv0==v0 && fv2==v1) || (fv0==v1 && fv2==v0)){
                    fe1=fv0->P()-fv1->P();
                    fe2=fv2->P()-fv1->P();
                }
                else if((fv1==v0 && fv2==v1) || (fv1==v1 && fv2==v0)){
                    fe1=fv1->P()-fv0->P();
                    fe2=fv2->P()-fv0->P();
                }
                else{
                    cout<< "wrong face"<<endl;
                }

                angle2=vcg::Angle(fe1,fe2);
                angleflag |=angle2<(PI/12);
            }


            // If a edge is shorther than threshold or if the opening angle at either opposing
            // vertex is less than π/12 then we try to collapse

            if(edgelen<threshold || angleflag){
                cout<<"Pode se collapse"<<endl;
                 // Test if the endpoints of the edge pertence to one crease or border edge
                 bool creaseV0=utility::checkVertexIsCrease(ei->V(0));
                 bool creaseV1=utility::checkVertexIsCrease(ei->V(1));
                 bool borderV0=ei->V(0)->IsB();
                 bool borderV1=ei->V(1)->IsB();
                 bool isfeatureV0=creaseV0 || borderV0;
                 bool isfeatureV1=creaseV1 || borderV1;
                 // flag to know if the edge is crease
                 bool creaseEdge=f->IsF(ez);
                 vcg::Point3f np;
                 if(!creaseEdge){
                     //If both vertices of an edge lie on feature edges, but the
                     //edge itself is not a feature edge then we cannot collapse
                     if(isfeatureV0 && !isfeatureV1){
                        np=ei->V(0)->P();
                        vertexPairCollapse vertexpair(ei->V(0),ei->V(1));
                        bool isFeasible=vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>>::LinkConditions(vertexpair);
                        if(isFeasible){
                             vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>> eco;
                             eco.Do(m,vertexpair,np);
                        //vcg::tri::UpdateTopology<CMesh>::AllocateEdge(m);
                        //vcg::tri::UpdateTopology<CMesh>::EdgeEdge(m);
                        //vcg::tri::UpdateTopology<CMesh>::VertexEdge(m);
                        cout<<"f0 e nao f1"<<endl;
                        //return this->collapseMesh(m,threshold,selectedFaces);
                        break;
                        }
                     }
                     else if(isfeatureV1 && !isfeatureV0){
                        np=ei->V(1)->P();
                        vertexPairCollapse vertexpair(ei->V(0),ei->V(1));
                        bool isFeasible=vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>>::LinkConditions(vertexpair);
                        if(isFeasible){
                        vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>> eco;
                        eco.Do(m,vertexpair,np);
                        //vcg::tri::UpdateTopology<CMesh>::AllocateEdge(m);
                        //vcg::tri::UpdateTopology<CMesh>::EdgeEdge(m);
                        //vcg::tri::UpdateTopology<CMesh>::VertexEdge(m);
                        cout<<"f1 e nao f0"<<endl;
                        //return this->collapseMesh(m,threshold,selectedFaces);
                        break;
                        }
                     }
                     else if(!isfeatureV1 && !isfeatureV0){
                        np=(ei->V(0)->P()+ei->V(1)->P())/2;
                        vertexPairCollapse vertexpair(ei->V(0),ei->V(1));
                        bool isFeasible=vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>>::LinkConditions(vertexpair);
                        if(isFeasible){
                        vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>> eco;
                        eco.Do(m,vertexpair,np);
                        //vcg::tri::UpdateTopology<CMesh>::AllocateEdge(m);
                        //vcg::tri::UpdateTopology<CMesh>::EdgeEdge(m);
                        //vcg::tri::UpdateTopology<CMesh>::VertexEdge(m);
                        cout<<"nao f0 e nao f1"<<endl;
                        //return this->collapseMesh(m,threshold,selectedFaces);
                        break;
                        }
                     }
                 }
                 else{
                     np=(ei->V(0)->P()+ei->V(1)->P())/2;
                     vertexPairCollapse vertexpair(ei->V(0),ei->V(1));
                     bool isFeasible=vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>>::LinkConditions(vertexpair);
                     if(isFeasible){
                     vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>> eco;
                     eco.Do(m,vertexpair,np);
                     //vcg::tri::UpdateTopology<CMesh>::AllocateEdge(m);
                     //vcg::tri::UpdateTopology<CMesh>::EdgeEdge(m);
                     //vcg::tri::UpdateTopology<CMesh>::VertexEdge(m);
                     cout<<"crease"<<endl;
                     //return this->collapseMesh(m,threshold,selectedFaces);
                     break;
                     }
                 }

            }
            else{
                cout<<"nao pode collapse"<<endl;
            }
        }


        /*CMesh::EdgeIterator ei=m.edge.begin();
        CMesh::FacePointer f=(*ei).EFp();
        vcg::BaseParameterClass param;
        vertexPairCollapse vertexpair(vertexpair((*f).V0(0),(*f).V1(0));
        TriEdgeCollapse<CMesh,vertexPairCollapse,TriEdgeCollapse> collapseSession(vertexpair,0,param);
        collapseSession.Execute(m,param);
        */
    }
    void collapseStripMesh(CMesh& m, float threshold,bool selectedFaces){
        const double PI  =3.141592653589793238463;
        //Preprocessing
        vcg::tri::UpdateTopology<CMesh>::VertexFace(m);
        vcg::tri::UpdateTopology<CMesh>::AllocateEdge(m);
        vcg::tri::UpdateTopology<CMesh>::VertexEdge(m);
        // Mark features vertexs and edges ( borders or crease ).
        vcg::tri::UpdateFlags<CMesh>::VertexClearS(m);
        //vcg::tri::Clean<CMesh>::RemoveNonManifoldFace(m);
        vcg::tri::UpdateFlags<CMesh>::FaceBorderFromFF(m);
        vcg::tri::UpdateFlags<CMesh>::VertexBorderFromFaceAdj(m); // clear Vertex border flag and after update
        float angle=PI/2;
        //vcg::tri::UpdateFlags<CMesh>::FaceFauxCrease(m,0.5); //marks as fauxedges all the non crease edges angle<100
        //vcg::tri::UpdateFlags<CMesh>::SelectVertexCornerBorder(m,angle);
        /*for(CMesh::FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
        {
          if(!(*fi).IsD()){
             for(int k=0;k<3;k++){
                 if(!fi->IsF(k)){
                     (*fi).V(k)->SetS();            //endpoints of the edge
                     (*fi).V((*fi).Next(k))->SetS();
                 }
             }
          }

        }
        */
        CMesh::EdgeIterator ei;
        for (ei = m.edge.begin();  ei!=m.edge.end(); ei++)
        {
            // First face contain edge
            CMesh::FacePointer f=(*ei).EFp(); //pointer to face contain the edge
            int ze=ei->EFi(); // index of the edge in the found face
            bool selected=f->IsS() || f->FFp(ze)->IsS();   // This edge pertence to a selected face
            if(ei->V(0)->IsB() && ei->V(1)->IsB() && (!selectedFaces || selected)){
                //cout<<" is border "<<endl;
                bool cornerV0=checkVertexIsCorner(ei->V(0),angle);
                bool cornerV1=checkVertexIsCorner(ei->V(1),angle);
                //cout<<"corner "<<cornerV0<<" "<<cornerV1<<endl;
                // Only collapse if the two endpoints are not corners
                if(!cornerV0 || !cornerV1){
                    //Compute edge length of the edge
                    float edgelen=vcg::Distance((*ei).P0(0),(*ei).P1(0));
                    if(edgelen<threshold){
                        //cout<<" border small"<<endl;
                        vcg::Point3f np;
                        if(!cornerV0 && !cornerV1) np=(ei->V(0)->P()+ei->V(1)->P())/2;
                        else if(cornerV0) np=ei->V(0)->P();
                        else if(cornerV1) np=ei->V(1)->P();

                        vertexPairCollapse vertexpair(ei->V(0),ei->V(1));
                        bool isFeasible=vcg::tri::EdgeCollapser<CMesh,vcg::tri::BasicVertexPair<CMesh::VertexType>>::LinkConditions(vertexpair);
                        if(isFeasible){
                            vcg::tri::EdgeCollapser<CMesh,vertexPairCollapse>::Do(m,vertexpair,np);
                           break;
                        }
                    }
                }
            }
        }
        //cout<<"collapse iteration finish"<<endl;
        //vcg::tri::UpdateColor<CMesh>::PerVertexConstant(m, vcg::Color4b::Red,true);
        vcg::tri::UpdateTopology<CMesh>::FaceFace(m);
        //vcg::tri::UpdateFlags<CMesh>::FaceBorderFromFF(m);
    }
    bool flipEdgesStripMesh(CMesh& m,bool selectedFaces){
        assert(vcg::tri::HasFFAdjacency(m));
        //vcg::tri::UpdateTopology<CMesh>::FaceFace(m);
        vcg::tri::UpdateTopology<CMesh>::VertexFace(m);
        CMesh::FaceIterator fi;

        for (fi = m.face.begin();  fi!=m.face.end(); fi++)
        {
            if(!(*fi).IsD() && (!selectedFaces || (*fi).IsS()))
            {
                //cout<<"selected face"<<endl;
                for(int z=0;z<3;z++){
                    if(!fi->IsB(z)){
                        //cout<<"possible edge"<<endl;
                        CMesh::VertexPointer OldDiag0 = fi->V0(z);
                        CMesh::VertexPointer OldDiag1 = fi->V1(z);
                        CMesh::VertexPointer NewDiag0 = fi->V2(z);
                        CMesh::VertexPointer NewDiag1 = fi->FFp(z)->V2(fi->FFi(z));

                        if(NewDiag0->IsS() && NewDiag1->IsS()){
                            /// Connecting opposing vertices and check if this second edge is shorter
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
                            float currentDiag=vcg::Distance(OldDiag0->P(),OldDiag1->P());
                            float newDiag=vcg::Distance(NewDiag0->P(),NewDiag1->P());
                            bool newIsShorter=newDiag<currentDiag;

                            // Check if the flip will not cause a loss of aspect ratio
                            float currentAR=aspectRatioTriangleFromFace(&*fi);
                            float newAR=aspectRatioTriangle(OldDiag0->P(),OldDiag1->P(),NewDiag0->P());
                            bool newHasBetterAR=newAR<currentAR;
                            bool isnotpossible=(NewDiag1 != NewDiag0) && (NewDiag1 != OldDiag0) && (NewDiag1 != OldDiag1);
                            //cout<<"vou checar"<<endl;
                            if(newIsShorter && newHasBetterAR && face::CheckFlipEdge(*fi, z) && !isnotpossible){// maybe can add face::CheckFlipEdgeNormal(*fi, z, math::ToRad(10.f))

                                face::FlipEdge(*fi, z);
                                //cout<<"flipped\n";
                                vcg::tri::UpdateTopology<CMesh>::FaceFace(m);
                                return true;

                            }

                        }
                    }
                }
            }
        }
        //cout<<"dont flip\n";
        return false;
    }


    float errorAB(const Point3f& A,const Point3f& B){
        const double PI = 3.141593;
        float angle=vcg::Angle(A,B);
        float a=-3/PI;
        if(angle<0) angle=-angle;
        return a*angle+3;
    }
    void correctPlanarBox(vcg::Box3<float>& aBox,float offset){
        if(aBox.DimX()==0){
           Point3f nmin(aBox.min.X()-offset,aBox.min.Y(),aBox.min.Z());
           Point3f nmax(aBox.min.X()+offset,aBox.min.Y(),aBox.min.Z());
           aBox.Add(nmin);
           aBox.Add(nmax);
        }
        if(aBox.DimY()==0){
           Point3f nmin(aBox.min.X(),aBox.min.Y()-offset,aBox.min.Z());
           Point3f nmax(aBox.min.X(),aBox.min.Y()+offset,aBox.min.Z());
           aBox.Add(nmin);
           aBox.Add(nmax);
        }
        if(aBox.DimZ()==0){
           Point3f nmin(aBox.min.X(),aBox.min.Y(),aBox.min.Z()-offset);
           Point3f nmax(aBox.min.X(),aBox.min.Y(),aBox.min.Z()+offset);
           aBox.Add(nmin);
           aBox.Add(nmax);
        }
    }
    bool pointIsInsidePoly(vector<CMesh::CoordType> &polygon,CMesh::CoordType point){
          int i, j, nvert = int(polygon.size());
          bool c = false;
          for (i = 0, j = nvert-1; i < nvert; j = i++) {
            if ( ((polygon[i].Y()>=point.Y()) != (polygon[j].Y()>=point.Y())) &&
             (point.X()<=(polygon[j].X()-polygon[i].X()) * (point.Y()-polygon[i].Y()) / (polygon[j].Y()-polygon[i].Y()) + polygon[i].X()) )
              c = !c;
          }
          return c;
    }
    void computeReflexion(CMesh::CoordType RP,CMesh::CoordType P,CMesh::CoordType A,CMesh::CoordType B){
        CMesh::CoordType direction=B-A;
        CMesh::CoordType vector=P-A;
        direction.normalized();
        CMesh::ScalarType projection_module=vector.dot(direction);
        CMesh::CoordType projection=direction*projection_module;
        CMesh::CoordType ortho=vector-projection;
        CMesh::CoordType opp=ortho*(-1);
        RP=projection+opp;
    }
    void toSingleBoundary(vector<Point2f> &outerB, vector<vector<Point2f> > &innerB, vector<Point3f> &singleB){
        singleB.clear();
        assert(isCounterClockWise(outerB));
        int closestInnerB=-1;
        vector<Point3i> diagonals; // [outerB index,inner index,innerB vector index]
        for(int i=0;i<innerB.size();i++){
            //for each hole
            assert(!isCounterClockWise(innerB[i]));
            int osize=outerB.size();
            for(int j=0;j<innerB[i].size();j++){
                for(int k=0;k<osize;k++){
                    Segment2f diagonal(outerB[k],innerB[i][j]); //diagonal
                    //check if this diagonal is valid
                    //first with the outer boundary
                    bool isvalid=true;
                    for(int r=0;r<osize;r++){
                        if(r==k || r==(k-1+osize)%osize)
                            continue;
                        Segment2f edge(outerB[r],outerB[(r+1)%osize]); //polygon's edge
                        Point2f result;
                        if(vcg::SegmentSegmentIntersection<float>(diagonal,edge,result)){
                            isvalid=false;
                            break;
                        }
                    }
                    if(!isvalid)
                        continue;

                    //second  with the other inner boundaries
                    for(int r=0;r<innerB.size();r++){
                        if(r==i)
                            continue;
                        for(int s=0;s<innerB[r].size();s++){
                            Segment2f edge(innerB[r][s],innerB[r][(s+1)%innerB[r].size()]); //hole's edge
                            Point2f result;
                            if(vcg::SegmentSegmentIntersection<float>(diagonal,edge,result)){
                                isvalid=false;
                                break;
                            }
                        }
                        if(!isvalid)
                            break;
                    }

                    if(isvalid){
                        diagonals.push_back(Point3i(k,j,i));
                    }
                    else{
                        continue;
                    }
                }
            }
        }
        assert(diagonals.size()!=0);
        Eigen::VectorXf lengths(diagonals.size());
        for(int i=0;i<diagonals.size();i++)
            lengths(i)=vcg::Distance(outerB[diagonals[i].X()],innerB[diagonals[i].Z()][diagonals[i].Y()]);
        int min=0;
        float minvalue=lengths.minCoeff(&min);
        //closest inner boundary
        closestInnerB=diagonals[min].Z();

        //constructing the new outer boundary by connecting the closest inner boundary to the outer boundary
        vector<Point2f> newouterB;
        int onOuter=diagonals[min].X();
        int start=onOuter;
        int onInner=diagonals[min].Y();
        int isize=innerB[closestInnerB].size();
        do{
           newouterB.push_back(outerB[start]);
           start=(start+1)%outerB.size();
       }while(start!=onOuter);
       newouterB.push_back(outerB[onOuter]);
       for(int k=0;k<=isize;k++)
          newouterB.push_back(innerB[closestInnerB][(k+onInner)%isize]);

       if(innerB.size()==1){
           for(int k=0;k<newouterB.size();k++)
                singleB.push_back(Point3f(newouterB[k].X(),newouterB[k].Y(),0.0));
       }
       else{
           vector<vector<Point2f>> newInnerB;
           for(int i=0;i<innerB.size();i++){
               if(i==closestInnerB)
                   continue;
               newInnerB.push_back(innerB[i]);
           }
           toSingleBoundary(newouterB,newInnerB,singleB);
       }
    }
    void toSingleBoundary(vector<Point2f> &outerB, vector<vector<Point2f> > &innerB, vector<Point3f> &singleB, set<int> &featureVertex){
        singleB.clear();
        assert(isCounterClockWise(outerB));
        int closestInnerB=-1;
        vector<Point3i> diagonals; // [outerB index,inner index,innerB vector index]
        for(int i=0;i<innerB.size();i++){
            //for each hole
            assert(!isCounterClockWise(innerB[i]));
            int osize=outerB.size();
            for(int j=0;j<innerB[i].size();j++){
                for(int k=0;k<osize;k++){
                    Segment2f diagonal(outerB[k],innerB[i][j]); //diagonal
                    //check if this diagonal is valid
                    //first with the outer boundary
                    bool isvalid=true;
                    for(int r=0;r<osize;r++){
                        if(r==k || r==(k-1+osize)%osize)
                            continue;
                        Segment2f edge(outerB[r],outerB[(r+1)%osize]); //polygon's edge
                        Point2f result;
                        if(vcg::SegmentSegmentIntersection<float>(diagonal,edge,result)){
                            isvalid=false;
                            break;
                        }
                    }
                    if(!isvalid)
                        continue;

                    //second  with the other inner boundaries
                    for(int r=0;r<innerB.size();r++){
                        if(r==i)
                            continue;
                        for(int s=0;s<innerB[r].size();s++){
                            Segment2f edge(innerB[r][s],innerB[r][(s+1)%innerB[r].size()]); //hole's edge
                            Point2f result;
                            if(vcg::SegmentSegmentIntersection<float>(diagonal,edge,result)){
                                isvalid=false;
                                break;
                            }
                        }
                        if(!isvalid)
                            break;
                    }

                    if(isvalid){
                        diagonals.push_back(Point3i(k,j,i));
                    }
                    else{
                        continue;
                    }
                }
            }
        }
        assert(diagonals.size()!=0);
        Eigen::VectorXf lengths(diagonals.size());
        for(int i=0;i<diagonals.size();i++)
            lengths(i)=vcg::Distance(outerB[diagonals[i].X()],innerB[diagonals[i].Z()][diagonals[i].Y()]);
        int min=0;
        float minvalue=lengths.minCoeff(&min);
        //closest inner boundary
        closestInnerB=diagonals[min].Z();

        //constructing the new outer boundary by connecting the closest inner boundary to the outer boundary
        vector<Point2f> newouterB;
        set<int> newfeature;
        int startNew=0;
        int onOuter=diagonals[min].X();
        int start=onOuter;
        int onInner=diagonals[min].Y();
        int isize=innerB[closestInnerB].size();
        newfeature.insert(0); // we start on the diagonal extreme on OuterB, and this is a feature
        do{
           if(featureVertex.find(start)!=featureVertex.end())
               newfeature.insert(startNew);
           newouterB.push_back(outerB[start]);
           start=(start+1)%outerB.size();
           startNew++;
       }while(start!=onOuter);
       //if(featureVertex.find(onOuter)!=featureVertex.end())
           newfeature.insert(startNew); // this position corresponds to the same diagonal extreme on outerB
       newouterB.push_back(outerB[onOuter]);
       startNew++;
       newfeature.insert(startNew);   // this position corresponds to the  diagonal extreme on closestInnerB
       for(int k=0;k<=isize;k++){
          newouterB.push_back(innerB[closestInnerB][(k+onInner)%isize]);
          startNew++;
       }
       startNew--;
       newfeature.insert(startNew);  // this position corresponds to the same diagonal extreme on closestInnerB
       if(innerB.size()==1){
           for(int k=0;k<newouterB.size();k++)
                singleB.push_back(Point3f(newouterB[k].X(),newouterB[k].Y(),0.0));
           featureVertex=newfeature;
       }
       else{
           vector<vector<Point2f>> newInnerB;
           for(int i=0;i<innerB.size();i++){
               if(i==closestInnerB)
                   continue;
               newInnerB.push_back(innerB[i]);
           }
           toSingleBoundary(newouterB,newInnerB,singleB,newfeature);
       }
    }
}
