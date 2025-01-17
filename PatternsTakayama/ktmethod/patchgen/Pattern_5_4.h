#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <sstream>
#include <vcg/complex/algorithms/create/platonic.h>
#include "edgeloop.h"
/*
equation for pattern 4: Fig 9(b) of Yassen et al., Computer-Aided Design 45 (2013)
  |0|     |1|     |0|     |0|     |1|     |0|     |1|   |3|    |1|   |l0|
  |1|     |0|     |1|     |0|     |0|     |1|     |0|   |3|    |1|   |l1|
p0|0| + p1|1| + p2|0| + p3|1| + p4|0| + q0|0| + q1|1| + |3| + x|0| = |l2|
  |0|     |0|     |1|     |0|     |1|     |0|     |0|   |1|    |0|   |l3|
  |1|     |0|     |0|     |1|     |0|     |1|     |0|   |2|    |0|   |l4|
*/
namespace patchgen {
    template <>
    struct Pattern<5, 4> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(5, 8);
                constraint_matrix << 0, 1, 0, 0, 1, 0, 1, 1,
                                     1, 0, 1, 0, 0, 1, 0, 1,
                                     0, 1, 0, 1, 0, 0, 1, 0,
                                     0, 0, 1, 0, 1, 0, 0, 0,
                                     1, 0, 0, 1, 0, 1, 0, 0;
            }
            return constraint_matrix;
        }
    
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return kt84::make_Vector5d(l[0] - 3,
                                       l[1] - 3,
                                       l[2] - 3,
                                       l[3] - 1,
                                       l[4] - 2);
        }
    
        static int& get_variable(PatchParam& param, int index) {
            if (index < 5) return param.p[index];
            if (index == 5) return param.q[0];
            if (index == 6) return param.q[1];
            if (index == 7) return param.x;
            return param.y;
        }

        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param) {
            auto& constraint_matrix = get_constraint_matrix();
            const int num_variables = constraint_matrix.cols();
            ILP ilp(num_variables);
            ilp.add_constraint(constraint_matrix, EQ, get_constraint_rhs(l));
        
            // arbitrary constraints and objective
            // xmin
            ilp.set_objective(kt84::make_Vector8d(0, 0, 0, 0, 0, 0, 0, 1), false);
            if (!ilp.solve()) return false;
            int xmin = ilp.get_variables()[7];
            // xmax
            ilp.refresh();
            ilp.set_objective(kt84::make_Vector8d(0, 0, 0, 0, 0, 0, 0, 1), true);
            if (!ilp.solve()) return false;
            int xmax = ilp.get_variables()[7];
            // xmid-1<=x<=xmid+1
            ilp.refresh();
            int xmid = (xmin + xmax) / 2;
            ilp.add_constraint(kt84::make_Vector8d(0, 0, 0, 0, 0, 0, 0, 1), LE, xmid + 1);
            ilp.add_constraint(kt84::make_Vector8d(0, 0, 0, 0, 0, 0, 0, 1), GE, xmid - 1);
            // p1<=q1
            ilp.add_constraint(kt84::make_Vector8d(0, 1, 0, 0, 0, -1, 0, 0), LE, 0);
            // p4<=q4
            ilp.add_constraint(kt84::make_Vector8d(0, 0, 0, 0, 1, 0, -1, 0), LE, 0);
            // maximize p0+p1+p2+p3+p4
            ilp.set_objective(kt84::make_Vector8d(1, 1, 1, 1, 1, 0, 0, 0), true);
            if (!ilp.solve()) return false;
        
            auto variables = ilp.get_variables();
            for (int i = 0; i < num_variables; ++i)
                get_variable(param, i) = variables[i];
        
            param.pattern_id = 4;
            return true;
        }
        template <typename PatchT>
        static void generate_subtopology(const PatchParam& param, PatchT& patch) {
            /*
            |             __C3--__    v--q1
            |         __--      __V5--__
            |     __--       __/        V4___   
            | C4--        __/          /     ----___
            |  \        _/            /             C2
            |   \    __/             /               /
            |    \ _/               /               /
            |    V6---------------V7--------------V3
            |      \             /                /
            |       \            |               /<--q0
            |   q0-->\          /         ___--V2
            |         \         |   ___---     / <--x
            |          \       /__--          /
            |          C0-----V0------V1-----C1
            |          q1--^    x--^ 
            */
         /*   patch.clear();
            typename PatchT::VHandle C[5];
            typename PatchT::VHandle V[8];
            for (int i = 0; i < 5; ++i) C[i] = add_tagged_vertex(patch, i, true );
            for (int i = 0; i < 8; ++i) V[i] = add_tagged_vertex(patch, i, false);
        
            patch.add_face(C[0], V[0], V[7], V[6]);
            patch.add_face(V[0], V[1], C[1], V[2]);
            patch.add_face(V[2], V[3], V[7], V[0]);
            patch.add_face(V[3], C[2], V[4], V[7]);
            patch.add_face(V[4], V[5], V[6], V[7]);
            patch.add_face(V[5], C[3], C[4], V[6]);
        
            auto h_insert_x = patch.halfedge_handle(V[1]);  // corresponds to V1-V0
            for (int i = 0; i < param.x; ++i)
                insert_edgeloop(patch, h_insert_x);
            auto h_insert_q0 = patch.halfedge_handle(V[3]);  // corresponds to V3-V2
            for (int i = 0; i < param.q[0]; ++i)
                insert_edgeloop(patch, h_insert_q0);
            auto h_insert_q1 = patch.halfedge_handle(V[0]);  // corresponds to V0-C0
            for (int i = 0; i < param.q[1]; ++i)
                insert_edgeloop(patch, h_insert_q1);
                */
            patch.Clear();
            bool hasAttributeL = vcg::tri::HasPerVertexAttribute(patch,"LeftSide");
            bool hasAttributeR = vcg::tri::HasPerVertexAttribute(patch,"RightSide");
            if(hasAttributeL)
                vcg::tri::Allocator<PatchT>::DeletePerVertexAttribute(patch,"LeftSide");
            if(hasAttributeR)
                vcg::tri::Allocator<PatchT>::DeletePerVertexAttribute(patch,"RightSide");

            auto side_indexL=vcg::tri::Allocator<PatchT>::template GetPerVertexAttribute<int>(patch,std::string("LeftSide"));
            auto side_indexR=vcg::tri::Allocator<PatchT>::template GetPerVertexAttribute<int>(patch,std::string("RightSide"));
            vcg::tri::Allocator<PatchT>::AddVertices(patch,13);
            patch.vert[0].SetS(); //C0
            patch.vert[1].SetS(); //C1
            patch.vert[2].SetS(); //C2
            patch.vert[3].SetS(); //C3
            patch.vert[4].SetS(); //C4
            // Populate index sides
            side_indexL[0]=1;
            side_indexR[0]=0;
            side_indexL[1]=0;
            side_indexR[1]=4;
            side_indexL[2]=4;
            side_indexR[2]=3;
            side_indexL[3]=3;
            side_indexR[3]=2;
            side_indexL[4]=2;
            side_indexR[4]=1;
            side_indexL[5]=0;   //V0
            side_indexR[5]=0;
            side_indexL[6]=0;   //V1
            side_indexR[6]=0;
            side_indexL[7]=4;   //V2
            side_indexR[7]=4;
            side_indexL[8]=4;   //V3
            side_indexR[8]=4;
            side_indexL[9]=3;   //V4
            side_indexR[9]=3;
            side_indexL[10]=3;   //V5
            side_indexR[10]=3;
            side_indexL[11]=1;   //V6
            side_indexR[11]=1;
            side_indexL[12]=-1;   //V7
            side_indexR[12]=-1;

            typename PatchT::FaceIterator pfi=vcg::tri::Allocator<PatchT>::AddFaces(patch,6);
            (*pfi).Alloc(4); //C[0], V[0], V[7], V[6]
            (*pfi).V(0)=&patch.vert[0];(*pfi).V(1)=&patch.vert[5];(*pfi).V(2)=&patch.vert[12];(*pfi).V(3)=&patch.vert[11];
            pfi++;
            (*pfi).Alloc(4);  //V[0], V[1], C[1], V[2]
            (*pfi).V(0)=&patch.vert[5];(*pfi).V(1)=&patch.vert[6];(*pfi).V(2)=&patch.vert[1];(*pfi).V(3)=&patch.vert[7];
            pfi++;
            (*pfi).Alloc(4);  //V[2], V[3], V[7], V[0]
            (*pfi).V(0)=&patch.vert[7];(*pfi).V(1)=&patch.vert[8];(*pfi).V(2)=&patch.vert[12];(*pfi).V(3)=&patch.vert[5];
            pfi++;
            (*pfi).Alloc(4);  //V[3], C[2], V[4], V[7]
            (*pfi).V(0)=&patch.vert[8];(*pfi).V(1)=&patch.vert[2];(*pfi).V(2)=&patch.vert[9];(*pfi).V(3)=&patch.vert[12];
            pfi++;
            (*pfi).Alloc(4);  //V[4], V[5], V[6], V[7]
            (*pfi).V(0)=&patch.vert[9];(*pfi).V(1)=&patch.vert[10];(*pfi).V(2)=&patch.vert[11];(*pfi).V(3)=&patch.vert[12];
            pfi++;
            (*pfi).Alloc(4);  //V[5], C[3], C[4], V[6]
            (*pfi).V(0)=&patch.vert[10];(*pfi).V(1)=&patch.vert[3];(*pfi).V(2)=&patch.vert[4];(*pfi).V(3)=&patch.vert[11];

            vcg::face::Pos<typename PatchT::FaceType> startPos=getPosFromIndex(patch,1,7);
            for (int i = 0; i < param.x; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

            startPos=getPosFromIndex(patch,7,8);
            for (int i = 0; i < param.q[0]; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

            startPos=getPosFromIndex(patch,0,5);
            for (int i = 0; i < param.q[1]; ++i)
               insert_edgeloop<PatchT>(patch,startPos);


            cout<<"patch 5 -- 4"<<endl;
        }
        /*static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            if (variable_indicators.empty()) {
                variable_indicators.resize(3);
                // for q0
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::V2, PatchVertexTag::V3));
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::V0, PatchVertexTag::V7));
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V6));
                // for q1
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V0));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::V6, PatchVertexTag::V7));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::V5, PatchVertexTag::V4));
                // for x
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::V0, PatchVertexTag::V1));
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::V2, PatchVertexTag::C1));
            }
            return variable_indicators;
        }*/
        static std::string get_param_str(const PatchParam& param) {
            std::stringstream ss;
            ss << "q0=" << param.q[0]
               << "_q1=" << param.q[1]
               << "_x=" << param.x;
            return ss.str();
        }
    };
}
