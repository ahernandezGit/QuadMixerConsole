#ifndef AABB_H
#define AABB_H

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#include <stdlib.h>
//#include <ctime>
#include <vector>
#include "meshtypes.h"
using namespace std;
//
// Axis aligned bounding boxes
//
namespace AABB {

    // Number of dimensions (typically 3)
    const int dim = 3;

    //
    // Class that represents a point in dim dimensions
    //
    struct Point {
        // Where the coordinates are stored
        double coord[dim];

        // Indexing operator
        double& operator[] (const int index) {
            return coord[index];
        }

        // Indexing operator
        const double& operator[] (const int index) const {
            return coord[index];
        }

        // Default constructor
        Point () {
            for (int i = 0; i < dim; i++) coord[i] = 0;
        }

        // Copy constructor
        Point (const Point& other) {
            for (int i = 0; i < dim; i++) coord[i] = other[i];
        }

        // Constructor from array
        Point (double * array) {
            for (int i = 0; i < dim; i++) coord[i] = array[i];
        }

        // Assignment operator
        Point& operator= (const Point& other) {
            for (int i = 0; i < dim; i++) coord[i] = other[i];
            return *this;
        }

    };

    // Prints a point
    ostream& operator<<(ostream& out, const Point& p) {
        out << "[" << p[0];
        for (int i = 1; i < dim; i++) {
            out << "," << p[i];
        }
        out << "]";
        return out;
    }

    //
    // Class that represents a single axis-aligned bounding box
    //
    //template <class FaceTypeMesh>
    struct Box {

        Point min, max; // Limits
        float x1,x2,x3,y1,y2,y3,z1,z2,z3;
        CMesh::FaceType* faceobject;
        PMesh::FaceType* faceobjectPolygon;
        //FaceTypeMesh faceobject;

        // Default Constructor
        Box() {
            x1=0;
            x2=0;
            x3=0;
            y1=0;
            y2=0;
            y3=0;
            z1=0;
            z2=0;
            z3=0;
        }

        // Constructor from one point
        Box (const Point& p) { min = max = p; }

        // Constructor from two points
        Box (const Point& p, const Point& q) {
            for (int i = 0; i < dim; i++) {
                if (p[i] < q[i]) {
                    min [i] = p[i];
                    max [i] = q[i];
                }
                else  {
                    min [i] = q[i];
                    max [i] = p[i];
                }
            }
        }

        // Enlarge this box to include point p
        Box& add (const Point& p) {
            for (int i = 0; i <  dim; i++) {
                if (p[i] < min[i]) {
                    min [i] = p[i];
                }
                else if (p[i] > max[i]) {
                    max[i] = p[i];
                }
            }
            return *this;
        };

        // Enlarge this box so that it also contains box other
        Box& add (const Box& other) {
            add (other.min);
            add (other.max);
            return *this;
        }

        // Returns true iff this box has a non-empty intersection with other
        bool intersects (const Box& other) {
            for (int i = 0; i < dim; i++) {
                if (MAX(min[i] , other.min[i]) > MIN(max[i],other.max[i])) return false;
            }
            return true;
        }

        // Possibly shrink this box so that it contains only points common with other box
        // Requirement: the boxes must intersect!
        Box& intersection (const Box& other) {
            assert (intersects(other));
            for (int i = 0; i < dim; i++) {
                min[i]=MAX(min[i], other.min[i]);
                max[i]=MIN(max[i], other.max[i]);
            }
            return *this;
        }

        // Returns the length of a box in the i'th dimension
        double length(const int i) const {
            return max[i]-min[i];
        }

        // Returns the axis with largest length
        int largestAxis () const {

            int largest = 0;
            double lmax = length(0);
            for (int i = 1; i < dim; i++) {
                double len = length(i);
                if (len > lmax) {
                    lmax = len;
                    largest = i;
                }
            }
            return largest;
        }

        void setTriangleFace(float a,float b,float c,float d,float e,float f,float g,float h,float k){
            x1=a;
            x2=b;
            x3=c;
            y1=d;
            y2=e;
            y3=f;
            z1=g;
            z2=h;
            z3=k;
        }
        void setTriangle(CMesh::FaceType *tface){
            faceobject=tface;
        }
        void setPolygonFace(PMesh::FaceType *tface){
            faceobjectPolygon=tface;
        }
    };

    // Prints a box
    ostream& operator<<(ostream& out, const Box& box) {
        out << "box(" << box.min << " , " << box.max << ")";
        return out;
    }

    // This comparison object encapsulates a function for ordering
    // Boxes according to the given axis
    struct BoxComparisonClass {
        int axis; // Which axis is used to compare min points of boxes
        bool operator() (const Box& a, const Box& b) {
            return a.min[axis] < b.min[axis];
        }
    } boxcompare;

    //
    // An AABB tree node
    //
    struct Tree {

        Box box;        // A bounding box for the whole tree
        Tree *son[2];	// The two sons (null if a leaf node)

        // Creates a tree for storing a collection of boxes
        Tree (vector<Box>& data) {

            // Find a bounding box for all elements
            int n = data.size();
            //cout<<"n "<<n<<endl;
            assert (n>0);
            box = data[0];
            if (n == 1) {
                // Leaf node : sons are null
                son[0] = son[1] = 0;
                return;
            }
            for (int i = 1; i < n; i++) {
                box.add(data[i]);
            }

            // Internal node:
            // Split the data collection in two by sorting them according
            // to the largest axis
            boxcompare.axis = box.largestAxis();
            sort (data.begin(),data.end(),boxcompare);

            vector<Box> sonData[2];
            int m = n/2;
            for (int i = 0; i < m; i++) sonData[0].push_back(data[i]);
            for (int i = m; i < n; i++) sonData[1].push_back(data[i]);

            son[0] = new Tree(sonData[0]);
            son[1] = new Tree(sonData[1]);
        }

        // Tells if this node is a leaf node
        bool isLeaf () const {
            assert (son[0] == 0 && son[1] == 0 || son[0] != 0 && son[1] != 0);
            return son[0] == 0;
        }

        // recursively collects all leaf boxes that intersect the query box and adds them
        // to result
        void intersectedBoxes (const Box& query, vector<Box*>& result) {
            if (isLeaf()) {
                if (box.intersects(query)) result.push_back (&box);
            }
            else {
                if (son[0]->box.intersects(query)) son[0]->intersectedBoxes(query, result);
                if (son[1]->box.intersects(query)) son[1]->intersectedBoxes(query, result);
            }
        }

    };

    // Prints a tree
    ostream& operator<<(ostream& out, const Tree& tree) {
        if (tree.isLeaf()) {
            out << "Leaf(" << tree.box << ")";
        }
        else {
            out << "Node(" << tree.box << ", left " << *tree.son[0] << ", right " << *tree.son[1] << ")";
        }
        return out;
    }

}

const int n = 1000000; // Number of boxes to generate
const int max_world_size = 10000; // Size of world
const int max_box_size = 100; // maximum box size

// Generates a random box of size up to max_box_size in each dimension
// somewhere in a virtual world of size max_world_size
AABB::Box randombox () {
    AABB::Point pos;
    for (int i = 0; i < AABB::dim; i++) {
        pos[i] = rand()%max_world_size;
    }
    AABB::Box box(pos);
    AABB::Point p;
    for (int i = 0; i < AABB::dim; i++) {
        p[i] = pos[i]+rand()%max_box_size;
    }
    box.add(p);
    return box;
}
// Creates a point from explicit coordinates
AABB::Point make_point (double x, double y, double z) {
    double v[] = { x, y, z };
    return AABB::Point (v);
}



#endif // AABB_H
