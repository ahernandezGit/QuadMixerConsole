#include "vectorcone.h"
#include "vcg/space/intersection2.h"
#include "myutils.h"

VectorCone::VectorCone()
{
    min.SetZero();
    max.SetZero();
}
VectorCone::VectorCone(const Point3f &v){
    Point3f vec=v;
    vec.normalized();
    min=vec;
    max=vec;
}

Point3f VectorCone::bisector(){
    Point3f result=min+max;
    return result.normalized();
}
void VectorCone::set(const Point3f &v){
    Point3f vec=v;
    vec.normalized();
    min=vec;
    max=vec;
}
float VectorCone::angle(){

    float min_x_max_z=min.X()*max.Y()-min.Y()*max.X();
    float minmax=min.dot(max);
    float a;
    if(minmax!=0)
        a=std::atan(min_x_max_z/minmax);
    else
        return M_PI/2;
    if(a<0){
        a+=M_PI;
    }
    if(std::abs(a)<epsilon && minmax<0)
        a=M_PI;
    return a;
}
void VectorCone::add(const Point3f & v){
    //assert(v.X()!=0 && v.Y()!=0 && v.Z()==0);
    if(v.Z()!=0){
         utility::printPoint3f(v," point not plane");
         return;
    }
    if((v.X()==0 && v.Y()==0) && v.Z()==0){
         utility::printPoint3f(v," adding trivial vector to cone ");
         return;
    }
    if(min==max){
        if(std::abs(min.X())<epsilon && std::abs(min.Y())<epsilon){
            Point3f newvalue=v;
            newvalue.normalized();
            min=newvalue;
            max=newvalue;
            //cout<<"primeira entrada"<<endl;
            return;
        }
        else{
            float zcrosproduct=min.X()*v.Y()-min.Y()*v.X();
            if(zcrosproduct<0){
                Point3f newvalue=v;
                newvalue.normalized();
                min=newvalue;
                return;
            }
            else if(std::abs(zcrosproduct)<epsilon){
                return;
            }
            else{
                Point3f newvalue=v;
                newvalue.normalized();
                max=newvalue;
                return;
            }
        }
    }
    Point3f v_min=min^v;
    Point3f v_max=max^v;
    if(v_min.Z()>=0 && v_max.Z()<=0){  // interior of the cone.
        //cout<<"not added anything"<<endl;
        return;
    }
    else if(v_min.Z()<=0 && v_max.Z()>=0){ // in this region the cone is ambiguous. So, we set the cone flat.
        min=max;
        min.Scale(-1,-1,0);
        return;
    }
    Point3f newvalue=v;
    newvalue.normalized();
    if(v_min.Z()<0)
        min=newvalue;
    else
        max=newvalue;
}

void VectorCone::printCone(){
    cout<<"Cone "<<endl;
    cout<<"min: ["<<min.X()<<" , "<<min.Y()<<"] "<<endl;
    cout<<"max: ["<<max.X()<<" , "<<max.Y()<<"] "<<endl;
    cout<<"angle: "<<angle()<<endl;
}

vector<int> MonotonicFunctions::monotonePartition(vector<Point3f> & pts, float ang, int wherestart) {
    vector<int> result;
    ang=M_PI*ang/180;
    int start = wherestart;
    int n = pts.size();
    Point3f prev = pts[start];
    result.push_back(start);
    Point3f v = pts[(start+1)%n]-prev;
    VectorCone cone(v);
    for (int j = 1; j < n; j++){
        int i = (j+start)%n;
        Point3f p = pts[i];
        v = pts[i]-prev;
        cone.add(v);
        if (cone.angle()>ang) {
            result.push_back((i+n-1)%n);
            cone.set(v);
        }
        prev = p;
    }
    //if (result [result.length-1] != (start+n-1)%n) result.push((start+n-1)%n);
    return result;
}

// Partitions a curve in monotonic chains and returns them as an array of curves
vector<vector<int>> MonotonicFunctions::monotonicPartitionCurves(vector<Point3f> curvePoints, float ang, int wherestart) {
    vector<int> partition = monotonePartition(curvePoints,ang,wherestart);
    int n = curvePoints.size();
    int icutPrev = partition.back();
    vector<vector<int>> result;
    for(int j=0;j<partition.size();j++){
        int icut=partition[j];
        vector<int> curve;
        while (icutPrev != icut) {
            curve.push_back(icutPrev);
            icutPrev = (icutPrev+1)%n;
        }
        curve.push_back(icut);
        result.push_back(curve);
    }
    return result;
}

MonotonicFunctions::DiagonalDS MonotonicFunctions::diagonals(vector<Point3f> &curve,vector<int> &partition,float sw,float aw,float iw, int index) {
    DiagonalDS result;
    vector<DiagonalDS> candidates;

    int n = partition.size();
    int m = curve.size();


       vector<int> npartition;
       // adding most concave vertex
       float concave=0;
       int ind_convace=-1;
       for(int i = 0; i < m; i++){
           Point3f p = curve[i];
           Point3f v1= curve[(i+1)%m]-p;
           Point3f v2= curve[(i-1+m)%m]-p;
           Point3f cross=v1^v2;
           if(cross.Z()>0){
               float angle=v1.dot(v2)/v1.Norm()*v2.Norm()+M_PI;
               if(angle>concave){
                   concave=angle;
                   ind_convace=i;
               }
           }
       }
       for(int i = 0; i < n ; i++){
           if(partition[i]<=ind_convace)
               npartition.push_back(partition[i]);
           else{
               npartition.push_back(ind_convace);
               npartition.push_back(partition[i]);
           }
       }
       partition=npartition;
       n = partition.size();

    float maxSize = 0;

    float sizeWeight=sw;
    float deviationWeight=aw;
    float imbalanceWeight=iw;

    // Find bisectors for all consecutive edges
    Point3f bisector[m];
    for(int ii = 0; ii < m; ii++){
        Point3f p = curve[ii];
        VectorCone cone(curve[(ii+1)%m]-p);
        cone.add(curve[(ii+m-1)%m]-p);
        bisector[ii] = cone.bisector();
    }

    for(int i = 0; i < n; i++){
        if(index!=-1){
            i=index;
        }
        int ii = partition[i];
        Point3f p = curve[ii];
        Point3f u0 = p-curve[(ii+1)%m];
        Point3f u1 = p-curve[(ii-1+m)%m];
        u0.Normalize();
        u1.Normalize();
        for (int j = 2; j+2 < m; j++) {
            int jj = (ii+j)%m;
            Point3f q = curve[jj];
            Point3f v = q-p;
            // Ray cast and see if it is a diagonal
            bool isvaliddiagonal = isDiagonalInterior(ii,jj,curve);
            if (isvaliddiagonal) {
                // Compute diagonal energy metrics
                int k = std::abs(jj-ii);
                float imbalance = float(std::max(k,m-k))/m;
                float size = v.Norm();
                v.Normalize();
                float angularDeviation = (2-std::abs(bisector[ii].dot(v))-std::abs(bisector[jj].dot(v)))/2;

                // Discard if collinear to neighboring edges
                if (1 - std::abs(v.dot(u0)) < sin(5*M_PI/180)) continue;
                if (1 - std::abs(v.dot(u1)) < sin(5*M_PI/180)) continue;
                Point3f w0 = q-curve[(jj+1)%m];
                Point3f w1 = q-curve[(jj-1+m)%m];
                w0.Normalize();
                w1.Normalize();
                if (1 - std::abs(v.dot(w0)) < epsilon) continue;
                if (1 - std::abs(v.dot(w1)) < epsilon) continue;

                maxSize = std::max(maxSize,size);
                result.setDiag(ii,jj,size,imbalance,angularDeviation);
                candidates.push_back(result);
            }
        }
        if(index!=-1){
            break;
        }
    }

    // Getting the best diag, i.e, the diag with the least energy.
    float bestEnergy=std::numeric_limits<float>::max();
    for(int i=0;i<candidates.size();i++){
        float energy=candidates[i].size/maxSize*sizeWeight+
                     candidates[i].angularDeviation*deviationWeight+
                     candidates[i].imbalance*imbalanceWeight;
        if(energy<bestEnergy){
            bestEnergy=energy;
            result.copy(candidates[i]);
        }
    }
    return result;
}

bool MonotonicFunctions::isDiagonalInterior(int id1,int id2,vector<Point3f> &polygon){    
    int n=polygon.size();
    assert(id1>=0 && id2>=0 && id1<n && id2<n);  
    bool intersectSides=false;
    Point2f a(polygon[id1].X(),polygon[id1].Y());
    Point2f b(polygon[id2].X(),polygon[id2].Y());
    Segment2<float> diagonal(a,b);
    int nextId1=(id1+1)%n;
    int nextId2=(id2+1)%n;
    int beforeId1=(id1-1+n)%n;
    while(!intersectSides && nextId1!=beforeId1){
        int follow=(nextId1+1)%n;
        if(follow!=id2 && follow!=nextId2){
            Point2f p1(polygon[nextId1].X(),polygon[nextId1].Y());
            Point2f p2(polygon[follow].X(),polygon[follow].Y());
            Point2f intersectionPoint;
            Segment2<float> side(p1,p2);
            intersectSides=vcg::SegmentSegmentIntersection(diagonal,side,intersectionPoint);            
        }
        nextId1=follow;
    }
    if(intersectSides){
        return false;
    }
    else{
        // check if the diagonal is completely in the polygon
        // It is necessary the polygon is int counterclockwise
        vector<Point3f> copy(polygon);

        if(!utility::isCounterClockWise(polygon)){     
            std::reverse(copy.begin(),copy.end());
            id1=n-1-id1;
            id2=n-1-id2;
        }       
        nextId1=(id1+1)%n;
        beforeId1=(id1-1+n)%n;
        Point3f V1=copy[nextId1]-copy[id1];
        Point3f V2=copy[beforeId1]-copy[id1];
        Point3f V3=copy[id2]-copy[id1];        

        vector<float> crossProd={V1.X()*V2.Y()-V1.Y()*V2.X(),V1.X()*V3.Y()-V1.Y()*V3.X(),V3.X()*V2.Y()-V3.Y()*V2.X()};
        bool isInterior=crossProd[0]>=0 && crossProd[1]>=0 && crossProd[2]>=0; // the cross product V1xV2 is positive
        isInterior=isInterior || (crossProd[0]<0 && (crossProd[1]>=0 || crossProd[2]>=0));
        if(isInterior)
            return true;
        else
            return false;

    }
}

vector<vector<Point3f>> MonotonicFunctions::splitCurve(vector<Point3f> &polygon, float limitAngle, int maxPartition,float sw,float aw,float iw){
    vector<vector<Point3f>> result;
    result.clear();
    vector<int> partition = monotonePartition(polygon,limitAngle);
    if (partition.size() <= maxPartition) {
        result.push_back(polygon);
        return result;
    }
    DiagonalDS diag = diagonals(polygon,partition,sw,aw,iw);
    // Getting two curves obtained by splitting the curve along the diagonal "diag"
    int n=polygon.size();
    /*cout<<"polygon size"<<n<<endl;
    cout<<"diag i j "<<diag.i<<" "<<diag.j<<endl;
    utility::printPoint3f(polygon[diag.i],"point i");
    utility::printPoint3f(polygon[diag.j],"point j");*/
    vector<Point3f> a,b;
    for (int k = diag.i; k != diag.j; k = (k+1)%n) a.push_back(polygon[k]);
    a.push_back(polygon[diag.j]);
    for (int k = diag.j; k != diag.i; k = (k+1)%n) b.push_back(polygon[k]);
    b.push_back(polygon[diag.i]);
    result.push_back(a);
    result.push_back(b);
    return result;
}

vector<vector<Point3f>> MonotonicFunctions::recursiveSplitCurve(vector<Point3f> &polygon, float limitAngle, int maxPartition,float sw,float aw,float iw){
    vector<vector<Point3f>> firstResult=splitCurve(polygon,limitAngle,maxPartition,sw,aw,iw);
    vector<vector<Point3f>> result;
    result.clear();
    if (firstResult.size() == 1){
        return firstResult;
    }
    else{
        //cout<<"result size "<<firstResult.size()<<endl;
        for(int i=0;i<firstResult.size();i++){
            vector<vector<Point3f>> partialResult=recursiveSplitCurve(firstResult[i],limitAngle,maxPartition,sw,aw,iw);
            result.insert(result.end(),partialResult.begin(),partialResult.end());
        }
        return result;
    }
}

MonotonicFunctions::PrePatch2D::PrePatch2D(){
    precorners.clear();
    //initial_polygon.Clear();
    allpolygons.clear();
}
void MonotonicFunctions::PrePatch2D::setnull(){
    precorners.clear();
    //initial_polygon.Clear();
    allpolygons.clear();
}

MonotonicFunctions::PrePatch2D::PrePatch2D(MonotonicFunctions::Polygon& polygon){
    set(polygon);
}
void MonotonicFunctions::PrePatch2D::set(Polygon &polygon){
    precorners.clear();
    //initial_polygon.Clear();
    allpolygons.clear();
    diagonalset.clear();
   /* int n=polygon.size();
    int nsides=polygon.numbersides();
    auto vi=vcg::tri::Allocator<CMesh>::AddVertices(initial_polygon,n);
    (*vi).P()=polygon.sides[0].vertices[0];
    vi++;
    for(int i=0;i<nsides-1;i++){
        if(vi!=initial_polygon.vert.end()){
            for(j=1;j<polygon.sides[i].size;j++){
                (*vi).P()=polygon.sides[i].vertices[j];
                vi++;
            }
        }
    }*/
    processPrePatch(polygon);
}

vector<MonotonicFunctions::Polygon> MonotonicFunctions::PrePatch2D::splitCurve(MonotonicFunctions::Polygon &polygon, float limitAngle, int maxPartition,float sw,float aw,float iw){
    vector<MonotonicFunctions::Polygon> result;
    result.clear();
    if(polygon.sides.size()<= maxPartition){
        result.push_back(polygon);
        return result;
    }
    /*poly.saveAsPLY("original");
    std::set<int> ind_added_into_diagonals;
    Polygon polygon=poly.subdivideDiagonals(ind_added_into_diagonals);
    polygon.saveAsPLY("subdivided");*/

    //polygon.saveAsPLY("original"+QTime::currentTime().toString().toStdString());
    vector<int> feature;
    vector<Point3f> polygonP=polygon.get_as_points(feature);
    int n=polygon.size();
    vector<int> partition = MonotonicFunctions::monotonePartition(polygonP,limitAngle);
    DiagonalDS diag = diagonals(polygonP,partition,sw,aw,iw);
    // constructing the sides of the new polygons
    vector<Point3f> aux;
    vector<PolygonSide> sa,sb;
    sa.clear();
    sb.clear();
    int iaux;
    int countDivisions=0;
    for (int k = diag.i; k != diag.j; k = (k+1)%n){
        //if(k==diag.i || ind_added_into_diagonals.find(k)==ind_added_into_diagonals.end())
            aux.push_back(polygonP[k]);
        int ison=isOnFeatures(feature,k); // give us on what side this point is
        if(ison!=-1 && aux.size()>1){  // if aux.size==1 so the diag have a feature point as extreme
            PolygonSide tempside(aux,polygon.sides[ison].flag,polygon.sides[ison].index_3D);
            sa.push_back(tempside);
            countDivisions+=tempside.size-1;
            aux.clear();
            aux.push_back(polygonP[k]);
            iaux=ison;
        }
    }
    aux.push_back(polygonP[diag.j]);
    iaux=(iaux+1)%polygon.numbersides();
    PolygonSide tempside1(aux,polygon.sides[iaux].flag,polygon.sides[iaux].index_3D);
    sa.push_back(tempside1);
    countDivisions+=tempside1.size-1;
    vector<Point3f> diag1={polygonP[diag.j],polygonP[diag.i]};
    //split the diagonal for keep the parity and better quadrangulation
    float edgelenPoly=polygon.averageEdgeLength_fixed();
    if(edgelenPoly==0)
        edgelenPoly=polygon.averageEdgeLength();
    assert(edgelenPoly!=0);
    float edgelenDiag=vcg::Distance(polygonP[diag.j],polygonP[diag.i]);
    int numberSegments=int(edgelenDiag/edgelenPoly)+1;
    if((numberSegments+countDivisions)%2!=0)
        numberSegments++;
    if(numberSegments>1){
        vector<Point3f> newdiag;
        utility::splitSegment(polygonP[diag.j],polygonP[diag.i],newdiag,numberSegments);
        newdiag.insert(newdiag.begin(),polygonP[diag.j]);
        newdiag.push_back(polygonP[diag.i]);
        diag1=newdiag;
    }
    PolygonSide tempdiag1(diag1,S_DIAGONAL,-1);
    sa.push_back(tempdiag1);
    MonotonicFunctions::Polygon a(sa);

    aux.clear();
    countDivisions=0;
    for (int k = diag.j; k != diag.i; k = (k+1)%n){
         //if(k==diag.j || ind_added_into_diagonals.find(k)==ind_added_into_diagonals.end())
            aux.push_back(polygonP[k]);
        int ison=isOnFeatures(feature,k); // give us on what side this point is
        if(ison!=-1 && aux.size()>1){  // if aux.size==1 so the diag have a feature point as extreme
            PolygonSide tempside(aux,polygon.sides[ison].flag,polygon.sides[ison].index_3D);
            sb.push_back(tempside);
            countDivisions+=tempside.size-1;
            aux.clear();
            aux.push_back(polygonP[k]);
            iaux=ison;
        }
    }
    aux.push_back(polygonP[diag.i]);
    iaux=(iaux+1)%polygon.numbersides();
    PolygonSide tempside2(aux,polygon.sides[iaux].flag,polygon.sides[iaux].index_3D);
    sb.push_back(tempside2);
    countDivisions+=tempside2.size-1;
    vector<Point3f> diag2={polygonP[diag.i],polygonP[diag.j]};
    /*edgelenDiag=vcg::Distance(polygonP[diag.j],polygonP[diag.i]);
    numberSegments=int(edgelenDiag/edgelenPoly)+1;
    if((numberSegments+countDivisions)%2!=0)
        numberSegments++;
   */
    if(numberSegments>1){
        vector<Point3f> newdiag;
        utility::splitSegment(polygonP[diag.i],polygonP[diag.j],newdiag,numberSegments);
        newdiag.insert(newdiag.begin(),polygonP[diag.i]);
        newdiag.push_back(polygonP[diag.j]);
        diag2=newdiag;
    }
    PolygonSide tempdiag2(diag2,S_DIAGONAL,-1);
    sb.push_back(tempdiag2);
    MonotonicFunctions::Polygon b(sb);    
    //cout<<"polygon with "<<a.sides.size()<<" sides"<<endl;
    //cout<<"polygon with "<<b.sides.size()<<" sides"<<endl;
    //a.saveAsPLY("a-result"+QTime::currentTime().toString().toStdString());
    //b.saveAsPLY("b-result"+QTime::currentTime().toString().toStdString());
    result.push_back(a);
    result.push_back(b);
    return result;
}
vector<MonotonicFunctions::Polygon> MonotonicFunctions::PrePatch2D::recursiveSplitCurve(MonotonicFunctions::Polygon &polygon, float limitAngle, int maxPartition, float sw, float aw, float iw){    
    vector<MonotonicFunctions::Polygon> firstResult=splitCurve(polygon,limitAngle,maxPartition,sw,aw,iw);
    for(int i=0;i<firstResult.size();i++){
        if((firstResult[i].size()%2!=0) || (firstResult[i].size()%2!=0)){
            cout<<"in recursive"<<endl;
            cout<<"polygon size "<<polygon.size()<<endl;
            cout<<"this size "<<firstResult[i].size()<<endl;
        }
    }
    vector<MonotonicFunctions::Polygon> result;
    result.clear();
    diagonalset.clear();
    if (firstResult.size() == 1){
        return firstResult;
    }
    else{
        //cout<<"result size "<<firstResult.size()<<endl;
        for(int i=0;i<firstResult.size();i++){
            vector<MonotonicFunctions::Polygon> partialResult=recursiveSplitCurve(firstResult[i],limitAngle,maxPartition,sw,aw,iw);          
            result.insert(result.end(),partialResult.begin(),partialResult.end());
        }
        return result;
    }
}
vector<vector<vector<Point3f>>> MonotonicFunctions::PrePatch2D::getPolygonsAsPrePatch(){
    assert(!allpolygons.empty());
    vector<vector<vector<Point3f>>> result;
    for(int i=0;i<allpolygons.size();i++){
        vector<vector<Point3f>> thispoly;
        //checking coherence of the polygon sides.
        for(int j=1;j<allpolygons[i].sides.size();j++){
            Point3f b=allpolygons[i].sides[j].vertices[0];
            Point3f a=allpolygons[i].sides[j-1].vertices.back();

            //we choose the point with the least distance. This because the floating point error when compare 3d points
            if(a!=b){
                Point3f c=allpolygons[i].sides[j].vertices.back();
                float dab=vcg::Distance(a,b);
                float dac=vcg::Distance(a,c);
                if(dac<dab)
                    allpolygons[i].reverse(j);
            }
        }
        //working  with a counterclokwise polygon
        vector<Point3f> aspoly=allpolygons[i].get_as_points();
        if(!utility::isCounterClockWise(aspoly)){
            for(int j=0;j<allpolygons[i].sides.size();j++){
                    allpolygons[i].reverse(j);
            }
        }
        //we add the polygon structure
        for(int j=0;j<allpolygons[i].sides.size();j++){
            thispoly.push_back(allpolygons[i].sides[j].vertices);
        }
        result.push_back(thispoly);
    }
    return result;
}

void MonotonicFunctions::PrePatch2D::processPrePatch(Polygon& poly){
     /*int n=poly.numbersides();
    CMesh tcurve;
    for(int i=0;i<n;i++){
        CMesh tside;
        for(int j=0;j<poly.sides[i].size-1;j++)
             tri::Allocator<CMesh>::AddEdge(tside,poly.sides[i].vertices[j],poly.sides[i].vertices[j+1]);
        switch (poly.sides[i].flag) {
            case S_FIXED: UpdateColor<CMesh>::PerVertexConstant(tside,Color4b::Blue);
                          break;
            case S_DIAGONAL: UpdateColor<CMesh>::PerVertexConstant(tside,Color4b::Red);
                          break;
            case S_UNDEFINED: UpdateColor<CMesh>::PerVertexConstant(tside,Color4b::Green);
                          break;
        }
        tri::Append<CMesh,CMesh>::Mesh(tcurve,tside);
    }
    std::string patchstring="parameterization/bigpolygon.ply";
    tri::io::ExporterPLY<CMesh>::Save(tcurve,patchstring.c_str(),tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
    */
    cout<<"current components with "<<poly.size()<<" divisions "<<endl;
    if(poly.size()==172){
        cout<<"aqui"<<endl;
    }
    allpolygons=recursiveSplitCurve(poly);   
    for(int i=0;i<allpolygons.size();i++){
        if(allpolygons[i].size()==15)
            cout<<"hola"<<endl;
        cout<<"  --subdivided polygon with "<<allpolygons[i].size()<< " divisions" <<endl;
    }
    /*for(int i=0;i<allpolygons.size();i++){
        int nthis=allpolygons[i].numbersides();
        CMesh thiscurve;
        for(int j=0;j<nthis;j++){
            CMesh tside;
            for(int k=0;k<allpolygons[i].sides[j].size-1;k++)
                 tri::Allocator<CMesh>::AddEdge(tside,allpolygons[i].sides[j].vertices[k],allpolygons[i].sides[j].vertices[k+1]);
            switch (allpolygons[i].sides[j].flag) {
                case S_FIXED: UpdateColor<CMesh>::PerVertexConstant(tside,Color4b::Blue);
                              break;
                case S_DIAGONAL: UpdateColor<CMesh>::PerVertexConstant(tside,Color4b::Red);
                              break;
                case S_UNDEFINED: UpdateColor<CMesh>::PerVertexConstant(tside,Color4b::Green);
                              break;
            }
            tri::Append<CMesh,CMesh>::Mesh(thiscurve,tside);
        }
        string thispatchstring="parameterization/subpolygon"+QString::number(i).toStdString()+".ply";
        tri::io::ExporterPLY<CMesh>::Save(thiscurve,thispatchstring.c_str(),tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
    }*/

    //Here we store all the diagonals. In fact, we store two equal diagonals por each two connected polygons
    // we also store for each diagonal the corresponding polygon and side for what it belongs
    /*vector<PolygonSide> alldiagonals;
    vector<pair<int,int>> diag_to_polygons;
    for(int i=0;i<allpolygons.size();i++){
        cout<<"  --subdivided polygon with "<<allpolygons[i].size()<< " divisions" <<endl;
        for(int j=0;j<allpolygons[i].sides.size();j++){
            if(allpolygons[i].sides[j].flag==S_DIAGONAL){
                allpolygons[i].sides[j].id=i;
                alldiagonals.push_back(allpolygons[i].sides[j]);
                diag_to_polygons.push_back(make_pair(i,j));
            }
        }
    }
    // We construct a correspondence between pair of digonals that are the same
    vector<int> correspondencediagonals(alldiagonals.size());
    for(int i=0;i<alldiagonals.size();i++){
        assert(alldiagonals[i].size==2);
        for(int j=0;j<alldiagonals.size();j++){
            if(j==i) continue;
            assert(alldiagonals[j].size==2);
            Point3f a=alldiagonals[i].vertices[0];
            Point3f b=alldiagonals[i].vertices[1];
            Point3f c=alldiagonals[j].vertices[0];
            Point3f d=alldiagonals[j].vertices[1];
            if(a==c && b==d){
                  correspondencediagonals[i]=j;
                  break;
            }
            if(a==d && b==c){
                  correspondencediagonals[i]=j;
                  break;
            }
        }

    }
    //Construct graph
    //cout<<"construct graph"<<endl;
    vector<myutils::Edge> edgesg;
    vector<bool> valuesg,constantsg;
    for(int i=0;i<correspondencediagonals.size();i++){
        MonotonicFunctions::PolygonSide thisdiag=alldiagonals[i];
        int src=thisdiag.id;
        int dest=alldiagonals[correspondencediagonals[i]].id;
        myutils::Edge temp={src,dest,1};
        edgesg.push_back(temp);
    }
    //resampliong diagonals to better quadrangulations
    std::set<int> indexes_diagonals;
    for(int i=0;i<alldiagonals.size();i++)
        indexes_diagonals.insert(i);
    int k=0;
    while(!indexes_diagonals.empty() && k<alldiagonals.size()){
        if(indexes_diagonals.find(k)!=indexes_diagonals.end()){
            int pa=diag_to_polygons[k].first;  // corresponding polygon to this diagonal
            int pa_s=diag_to_polygons[k].second;  // side on the corresponding polygon that corresponds to this diagonal
            int pb=diag_to_polygons[correspondencediagonals[k]].first;  // corresponding polygon to this diagonal (there are always two repeated diagonals)
            int pb_s=diag_to_polygons[correspondencediagonals[k]].second;  // side on the corresponding polygon that corresponds to this diagonal
            float edgelenPoly=(allpolygons[pa].averageEdgeLength_fixed()+allpolygons[pb].averageEdgeLength_fixed())/2;
            float edgelenDiag=alldiagonals[k].length();
            int numberSegments=int(edgelenDiag/edgelenPoly)+1;
            //int numberSegments=5;
            if(numberSegments>1){
                vector<Point3f> newdiag;
                Point3f direction=alldiagonals[k].vertices.back()-alldiagonals[k].vertices[0];
                direction.normalized();
                float step=edgelenDiag/numberSegments;
                for(int j=0;j<=numberSegments;j++){
                    newdiag.push_back(alldiagonals[k].vertices[0]+direction*(j*step));
                }
                allpolygons[pa].sides[pa_s].vertices=newdiag;
                allpolygons[pb].sides[pb_s].vertices=newdiag;
            }
            auto element=indexes_diagonals.find(k);
            indexes_diagonals.erase(element);
            element=indexes_diagonals.find(correspondencediagonals[k]);
            indexes_diagonals.erase(element);
        }
        k++;
    }
    for(int i=0;i<allpolygons.size();i++){
        valuesg.push_back(allpolygons[i].size()%2);
        constantsg.push_back(allpolygons[i].getConstant());
    }
    graph.set(edgesg,valuesg,constantsg);
    Eigen::MatrixXi matrixgraph=graph.toAdjacencyMatrix();*/
    //graph.export_as_dot(matrixgraph,graph.values,"graph");
    //graph.GurobiSolver();
}
int MonotonicFunctions::PrePatch2D::isOnFeatures(const vector<int> &feature, int ind){
    int s=feature.size();
    for(int k=0;k<s;k++){
        if(ind==feature[k]){
            return (k-1+s)%s;
        }
    }
    return -1;
}
