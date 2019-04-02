#ifndef ILPGUROBI_H
#define ILPGUROBI_H
#include "gurobi_c++.h"
#include <Eigen/Core>

const int EQ=3;
const int LE=1;
const int GE=2;
struct ILPGu{
     int num_variables;
     GRBEnv *env=NULL;
     GRBModel *model=NULL;
     GRBVar* vars=NULL;
     ILPGu(int inum_variables): num_variables(inum_variables){
         env = new GRBEnv ();
         env->set(GRB_IntParam_OutputFlag, 0);
         model = new GRBModel(*env);
         // set all variables to be integer
         char intergertypes[num_variables];
         for(int i=0;i<num_variables;i++){
              intergertypes[i]=GRB_INTEGER;
         }
         vars= model->addVars(NULL,NULL,NULL,intergertypes,NULL,num_variables);
         model->setObjective(GRBLinExpr(0));
         // default ojective
         //set_objective(Eigen::VectorXd::Ones(num_variables), true);
     }
     ~ILPGu(){
         delete [] vars ;
         delete model;
         delete env ;
     }
     void add_constraint(const Eigen::VectorXd& row, int constr_type, double rhs) const {
         assert(row.size() == num_variables);
         char sense;
         switch (constr_type) {
            case 1: sense=GRB_LESS_EQUAL; break;
            case 3 : sense=GRB_EQUAL; break;
            case 2 : sense=GRB_GREATER_EQUAL; break;
            default: assert(false);
         }
         //construnctin the expresion
         GRBLinExpr lhs=0;
         for(int j = 0; j<row.size();j++)
             lhs+=row(j)*vars[j];
         model->addConstr(lhs,sense,rhs);
     }
     void add_constraint(const Eigen::MatrixXd& rows, int constr_type, const Eigen::VectorXd& rhs) const {
         for (int i = 0; i < rows.rows(); ++i)
             add_constraint(rows.row(i), constr_type, rhs[i]);
     }
     void set_objective(const Eigen::VectorXd& row, bool is_maxim) const {
            assert(row.size() == num_variables);
            // Set objective
            GRBLinExpr funobj=0;
            for(int i=0;i<row.size();i++){
                funobj+=row[i]*vars[i];
            }
            if(is_maxim)
                model->setObjective(funobj,GRB_MAXIMIZE);
            else
                model->setObjective(funobj,GRB_MINIMIZE);
     }
     bool solve() const {
         model->optimize();
         // Status checking
         auto status = model->get(GRB_IntAttr_Status);
         if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED) {
            std::cout << " The model cannot be solved " <<
                 " because it is infeasible or unbounded " << std::endl ;
         }
         return model->get(GRB_IntAttr_Status)== GRB_OPTIMAL;
     }
        void refresh() {            // workaround for weird error when re-solving
            delete [] vars ;
            delete model;

            model = new GRBModel(*env);
            // set all variables to be integer
            char intergertypes[num_variables];
            for(int i=0;i<num_variables;i++){
                 intergertypes[i]=GRB_INTEGER;
            }
            vars= model->addVars(NULL,NULL,NULL,intergertypes,NULL,num_variables);
        }
        Eigen::VectorXi get_variables() const {
            Eigen::VectorXd variables = Eigen::VectorXd::Zero(num_variables);
            for ( int i = 0; i <variables.size() ; i++)
                variables(i) = vars[i].get(GRB_DoubleAttr_X);
            return variables.cast<int>();
        }
};

#endif // ILPGUROBI_H
