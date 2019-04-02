#ifndef GRAPH_H
#define GRAPH_H
//#include "gurobi_c++.h"
#include <iostream>
#include <vector>
#include <QString>
#include <QTime>
#include <QDate>
#include <fstream>
#include <assert.h>
#include <Eigen/Core>


using namespace std;

namespace  myutils{

// data structure to store graph edges
struct Edge {
    int src;
    int dest;
    int weight;
};

// class to represent a graph object
class Graph
{
    public:
    // An array of vectors to represent adjacency list
    vector<vector<pair<int, int>>> adjList;
    vector<bool> values;
    vector<bool> constants;
    // Default Constructor
    Graph(){
        adjList.clear();
        values.clear();
        constants.clear();
    }
    // Constructor
    Graph(vector<Edge> &edges, vector<bool> &ival,vector<bool> &iconst)
    {
       values.clear();
       constants.clear();
       adjList.clear();
       set(edges,ival,iconst);
    }
    //Initializator
    void set(vector<Edge> &edges, vector<bool> &ival,vector<bool> &iconst){

        assert(ival.size()==iconst.size());
        assert(ival.size()!=0);

        int N=ival.size();
        adjList.resize(N);
        values=ival;
        constants=iconst;
        // add edges to the directed graph
        for (unsigned i = 0; i < edges.size(); i++)
        {
            int src = edges[i].src;
            int dest = edges[i].dest;
            int weight = edges[i].weight;

            // insert at the end
            adjList[src].push_back(make_pair(dest, weight));

            // Because it is an undirected graph
            adjList[dest].push_back(make_pair(src, weight));
        }
    }
    void addGraph(Graph & added){
        //assert(added.values.size()==added.values.size()==)
        int na=added.values.size();
        int nn=values.size();
        // updating values and constants
        for(int i=0;i<na;i++){
            values.push_back(added.values[i]);
            constants.push_back(added.constants[i]);
        }
        vector<vector<pair<int, int>>> newadjList=added.adjList;
        for(int i=0;i<na;i++){
            vector<pair<int,int>> toadd;
            for (pair<int, int> v : newadjList[i]){
                toadd.push_back(make_pair(v.first+nn,v.second));
            }
            adjList.push_back(toadd);
        }
        //adjList.insert(adjList.end(),newadjList.begin(),newadjList.end());
    }
    void save(){
        QDate cd = QDate::currentDate();
        QTime ct = QTime::currentTime();
        QString filename ="../exported/"+ct.toString("hh:mm:ss")+cd.toString("'graph_'MM_dd_yyyy'.txt'");
        std::ofstream ofile(filename.toStdString(), ios::out);
        if(ofile.is_open()){
            ofile<<"vertices "<<values.size()<<"\n";
            for(int i=0;i<values.size();i++){
               for (pair<int, int> v : adjList[i])
                   ofile<<i<<" "<<v.first<<" "<<v.second<<" "<<values[i]<<
                          " "<<values[v.first]<<" "<<constants[i]<<" "<<constants[v.first]<<"\n";
            }
            ofile.close();
        }
    }
    Eigen::MatrixXi toAdjacencyMatrix(){
        Eigen::MatrixXi result=Eigen::MatrixXi::Constant(values.size(),values.size(),-1);
        for (int i = 0; i < values.size(); i++)
        {
            for (pair<int, int> v : adjList[i]){
                result(i,v.first)=v.second;
            }
        }
        return result;
    }
    // Destructor
    ~Graph() {
    }
    void export_as_dot(Eigen::MatrixXi &graph,vector<bool> &ivalues, QString name){
        int n=values.size();
        assert(n==adjList.size());
        cout<<"vertices this graph "<<n<<endl;

        string filenamefinal=name.toStdString();
        filenamefinal+=QTime::currentTime().toString().toStdString() + QString::number(n).toStdString()+".txt";
        std::ofstream fileGraph(filenamefinal);
        if (fileGraph.is_open()){
            fileGraph<<"graph { \n";
            for(int j=0;j<graph.rows();j++){
                for(int k=j;k<graph.rows();k++){
                    if(graph(j,k)!=-1){
                        if(graph(j,k)==0)
                            fileGraph<<"v"<<j<<" -- "<<"v"<<k<<" [color=green,penwidth=3.0] \n";
                        else if(graph(j,k)==1)
                            fileGraph<<"v"<<j<<" -- "<<"v"<<k<<" [color=red,penwidth=3.0] \n";
                    }
                }
            }
            fileGraph<<"{\n";
            for(int j=0;j<ivalues.size();j++){
                if(ivalues[j])
                    fileGraph<<"v"<<j<<" [style=filled,fillcolor=red,fontcolor=white]\n";
                else
                    fileGraph<<"v"<<j<<" [style=filled,fillcolor=green]\n";
            }
            fileGraph<<"}\n";
            fileGraph<<"} \n";
            fileGraph.close();
        }
    }
    void export_as_dot(QString name){
       Eigen::MatrixXi thismatrix=toAdjacencyMatrix();
       vector<bool> thisvalues=values;
       export_as_dot(thismatrix,thisvalues,name);
    }
    // print adjacency list representation of graph
    void printGraph()
    {
        for (int i = 0; i < values.size(); i++)
        {
            // print all neighboring vertices of given vertex
            for (pair<int, int> v : adjList[i])
                cout << "(" << i << ", " << v.first <<
                        ", " << v.second << ") ";
            cout << endl;
        }
    }
    void GurobiSolver(){
        /*GRBEnv env = GRBEnv();
        env.set(GRB_IntParam_OutputFlag, 0); // to avoid information printing
        GRBModel model = GRBModel(env);
        Eigen::MatrixXi adjacency=toAdjacencyMatrix();
        // Create variables
        int numberEdges=0;
        int numberVertices=values.size();
        for(int j=0;j<adjacency.rows();j++){
            for(int k=0;k<adjacency.rows();k++){
                if(adjacency(j,k)!=-1){
                    numberEdges++;
                }
            }
        }
        assert(numberEdges%2==0);
        numberEdges=numberEdges/2; // the adjacency matrix is symmetric
        double elb[numberEdges];
        double clb[numberVertices];
        double erb[numberEdges];
        double crb[numberVertices];
        char binarytypes[numberEdges];
        char intergertypes[numberVertices];
        for(int i=0;i<numberEdges;i++){
            clb[i]=0.0;
            elb[i]=0.0;
        }
        for(int i=0;i<numberEdges;i++){
           binarytypes[i]=GRB_BINARY;
           erb[i]=1.0;
        }
        for(int i=0;i<numberVertices;i++){
           intergertypes[i]=GRB_INTEGER;
           crb[i]=3.0;
        }

        GRBVar* e= model.addVars(elb,erb,NULL,binarytypes,NULL,numberEdges);
        GRBVar* c= model.addVars(clb,crb,NULL,intergertypes,NULL,numberVertices);

        //Enumering the e_i variables
        vector<std::pair<int,int>> ei;
        for(int j=0;j<adjacency.rows();j++){
            for(int k=j;k<adjacency.rows();k++){
                if(adjacency(j,k)!=-1){
                   ei.push_back(make_pair(j,k));
                }
            }
        }
        assert(numberEdges==ei.size());
        //Adding constraints
        for(int i=0;i<numberVertices;i++){
            // finding edges from this node
            vector<int> associated;
            for(int j=0;j<ei.size();j++){
                if(ei[j].first==i || ei[j].second==i){
                    associated.push_back(j);
                }
            }
            if(!associated.empty()){
                GRBLinExpr lhs=0;
                for(int j = 0; j < associated.size();j++)
                    lhs+=e[associated[j]];
                lhs+=-2*c[i];
                model.addConstr(lhs,GRB_EQUAL,double(-constants[i]));  //ej+ek+...ep-2pi=-constant
            }

        }
        // Set objective
        GRBLinExpr funobj=0;
        for(int i=0;i<ei.size();i++){
            funobj+=e[i];
        }
        model.setObjective(funobj,GRB_MAXIMIZE);
        // Optimize model
        model.optimize ();

        if(model.get(GRB_IntAttr_Status)==GRB_OPTIMAL){
            for(int i=0;i<ei.size();i++){
                cout<<"e"<<i<<" "<<e[i].get(GRB_DoubleAttr_X)<<endl;
            }

            //Computing the new configuration
            for(int i=0;i<ei.size();i++){
                int newweight=int(e[i].get(GRB_DoubleAttr_X));
                assert(newweight==0 || newweight==1 );
                adjacency(ei[i].first,ei[i].second)=newweight;
                adjacency(ei[i].second,ei[i].first)=newweight;
            }
            for(int i=0;i<numberVertices;i++){
                // finding edges from this node
                vector<int> associated;
                for(int j=0;j<ei.size();j++){
                    if(ei[j].first==i || ei[j].second==i){
                        associated.push_back(j);
                    }
                }
                if(!associated.empty()){
                    int cyclesum=constants[i]?1:0;
                    for(int j = 0; j < associated.size();j++){
                        int newweight=int(e[associated[j]].get(GRB_DoubleAttr_X));
                        cyclesum+=newweight;
                    }
                    cyclesum=cyclesum%2;
                    values[i]=cyclesum;
                }
            }
            export_as_dot(adjacency,values,"optimizedSolution");
        }
        else{
            cout<<"MODEL IMPOSIBLE TO OPTIMIZE"<<endl;
        }
        delete [] e;
        delete [] c;*/
    }
};

}
#endif // GRAPH_H
