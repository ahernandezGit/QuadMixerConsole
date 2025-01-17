#ifndef QUAD_TRACER
#define QUAD_TRACER
#include <unordered_map>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/append.h>

template <class PolyMeshType>
class QuadMeshTracer
{
    typedef typename PolyMeshType::FaceType PFace;

    PolyMeshType *PolyM;

    std::deque<vcg::face::Pos<PFace> > TracingStack;
    std::set<size_t> VisistedVertices;
    std::set<std::pair<size_t,size_t> > TracedEdges;

public:

    bool MotorCycle=true;
    std::vector<int> FacePatch;
    bool DebugMessages=false;
private:

    void InitSingularities()
    {
        TracedEdges.clear();
        TracingStack.clear();
        vcg::tri::UpdateQuality<PolyMeshType>::VertexValence(*PolyM);
        for (size_t i=0;i<(*PolyM).face.size();i++)
            for (int j=0;j<(*PolyM).face[i].VN();j++)
            {
                if ((*PolyM).face[i].V(j)->IsB())continue;
                if ((*PolyM).face[i].V(j)->Q()==4)continue;
                vcg::face::Pos<PFace> CurrP(&(*PolyM).face[i],j);

                size_t IndexV=vcg::tri::Index(*PolyM,CurrP.V());
                VisistedVertices.insert(IndexV);

                CurrP.FlipV();
                TracingStack.push_back(CurrP);

            }
        if (DebugMessages)
        std::cout<<"Inserted "<<TracingStack.size()<<" portals"<<std::endl;
    }

    int SplitIntoPatches()
    {
        FacePatch.resize((*PolyM).face.size(),-1);
        std::vector<size_t> StackMoves;
        vcg::tri::UpdateFlags<PolyMeshType>::FaceClearS(*PolyM);
        int currPartition=0;
        while (true)
        {
            //find next seed
            for (size_t i=0;i<(*PolyM).face.size();i++)
            {
                if ((*PolyM).face[i].IsS())continue;
                StackMoves.push_back(i);break;
            }
            if (StackMoves.empty())return currPartition;
            if (DebugMessages)
                std::cout<<"Partition "<<currPartition<<std::endl;
            while (!StackMoves.empty())
            {
                size_t currF=StackMoves.back();
                (*PolyM).face[currF].SetS();
                StackMoves.pop_back();
                FacePatch[currF]=currPartition;
                for (int j=0;j<(*PolyM).face[currF].VN();j++)
                {
                    size_t IndexV0=vcg::tri::Index(*PolyM,(*PolyM).face[currF].V0(j));
                    size_t IndexV1=vcg::tri::Index(*PolyM,(*PolyM).face[currF].V1(j));
                    std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                    //found a border edge
                    if (TracedEdges.count(key)>0)continue;
                    // go on the other size
                    size_t nextF=vcg::tri::Index(*PolyM,(*PolyM).face[currF].FFp(j));
                    if ((*PolyM).face[nextF].IsS())continue;
                    if (nextF==currF)continue;
                    StackMoves.push_back(nextF);
                }
            }
            currPartition++;
        }
        return -1;
    }

    void DoTrace()
    {
        if (DebugMessages)
            std::cout<<TracingStack.size() << " need to be traced! "<<std::endl;

        int step=0;
        do{
            vcg::face::Pos<PFace> CurrP=TracingStack.front();
            TracingStack.pop_front();
            size_t IndexV0=vcg::tri::Index(*PolyM,CurrP.V());
            size_t IndexV1=vcg::tri::Index(*PolyM,CurrP.VFlip());
            assert(IndexV1!=IndexV0);

            // add to traced paths
            std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            TracedEdges.insert(key);

            //already visited
            if (VisistedVertices.count(IndexV0)>0)continue;

            //border reached
            if ((*PolyM).vert[IndexV0].IsB())continue;

            //must been visited already
            if (MotorCycle)
            {
                assert(VisistedVertices.count(IndexV1)>0);
            }

            if (MotorCycle)
                VisistedVertices.insert(IndexV0);

            //then go to next one
            CurrP.FlipE();
            CurrP.FlipF();
            CurrP.FlipE();
            CurrP.FlipV();

            TracingStack.push_front(CurrP);

            step++;
            if (((step%100)==0)&&(DebugMessages))
                std::cout<<"Traced "<<step<<" separatrices"<<std::endl;
        }while (!TracingStack.empty());
    }



public:

    void SaveColoredMesh()
    {
        int NumPartitions=0;
        for (size_t i=0;i<FacePatch.size();i++)
            NumPartitions=std::max(NumPartitions,FacePatch[i]);

        for (size_t i=0;i<(*PolyM).face.size();i++)
        {
            assert(FacePatch[i]>=0);
            vcg::Color4b currC=vcg::Color4b::Scatter(NumPartitions,FacePatch[i]);
            (*PolyM).face[i].C()=currC;
        }

        //vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(PolyM,PatchPath.c_str(),vcg::tri::io::Mask::IOM_FACECOLOR);
    }

    void SavePatchDecomposition(std::string &PatchFile)
    {
        FILE *f=fopen(PatchFile.c_str(),"wt");
        for (size_t i=0;i<FacePatch.size();i++)
            fprintf(f,"%d\n",FacePatch[i]);
        fclose(f);
    }

    void SaveTracedEdgeMesh(std::string &EdgePath)
    {
        PolyMeshType TestEdge;
        std::set<std::pair<size_t,size_t> >::iterator iteE;
        for (iteE=TracedEdges.begin();iteE!=TracedEdges.end();iteE++)
        {
            typename PolyMeshType::CoordType pos0=(*PolyM).vert[(*iteE).first].P();
            typename PolyMeshType::CoordType pos1=(*PolyM).vert[(*iteE).second].P();
            vcg::tri::Allocator<PolyMeshType>::AddEdge(TestEdge,pos0,pos1);
        }
        //vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(TestEdge,EdgePath.c_str(),vcg::tri::io::Mask::IOM_EDGEINDEX);
    }

    std::set<std::pair<size_t,size_t>> getTraceEdges(){
        return TracedEdges;
    }

    void TracePartitions()
    {
        InitSingularities();

        DoTrace();

        SplitIntoPatches();
    }

    //QuadMeshTracer(PolyMeshType &_PolyM):PolyM(_PolyM){}
    QuadMeshTracer(){
        PolyM=NULL;
    }
    void setModel(PolyMeshType & _PolyM){
        PolyM=&_PolyM;
        VisistedVertices.clear();
    }
};



#endif
