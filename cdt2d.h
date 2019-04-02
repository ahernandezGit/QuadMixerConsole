#ifndef CDT2D_H
#define CDT2D_H
#include <list>
#include <cmath>
//#include <unordered_map>
//#include <iterator>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Polygon_2.h>
//#include <CGAL/number_utils.h>
//#include <CGAL/Segment_2.h>
//#include <CGAL/intersections.h>


template<class Kernel, class Container>
void print_polygon (const CGAL::Polygon_2<Kernel, Container>& P)
{
  typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator vit;
  std::cout << "[ " << P.size() << " vertices:";
  for (vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    std::cout << " (" << *vit << ')';
  std::cout << " ]" << std::endl;
}
template<class PointType>
class cdt2d
{
public:
    //typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
    typedef CGAL::Exact_predicates_exact_constructions_kernel                   KernelE;
    typedef KernelE::Point_2                                                    PointE;
    //typedef KernelE::Segment_2                                               Segment_2;
    //typedef Kernel::Point_2                                                      Point;
    typedef CGAL::Polygon_2<KernelE>                                         Polygon_2;
    typedef Polygon_2::Edge_const_iterator                                EdgeIterator;
    //typedef CGAL::Polygon_with_holes_2<KernelE>                   Polygon_with_holes_2;
    //typedef std::list<Polygon_with_holes_2>                                 Pwh_list_2;
    typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int,KernelE>      Vb;
    typedef CGAL::Triangulation_face_base_2<KernelE>                               Fbb;
    typedef CGAL::Constrained_triangulation_face_base_2<KernelE,Fbb>                Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                            Tds;
    typedef CGAL::Constrained_Delaunay_triangulation_2<KernelE,Tds,CGAL::Exact_predicates_tag>                CDT;
    typedef CDT::Vertex_handle                                           Vertex_handle;
    //typedef CGAL::cpp11::result_of<KernelE::Intersect_2(Segment_2,Segment_2)>::type intersection_type;

public:
    cdt2d(){
        inner_points.clear();
        outer_poly.clear();
    }
    cdt2d(std::vector<PointType> &polygon,std::vector<PointType> &innerPoints){
        inner_points.clear();
        outer_poly.clear();
        //storing data
        for(size_t i=0;i<polygon.size();i++)
            outer_poly.push_back(PointE(polygon[i].X(),polygon[i].Y()));

        if(outer_poly.orientation() == CGAL::CLOCKWISE){
            outer_poly.reverse_orientation();
        }
        for(size_t i=0;i<innerPoints.size();i++)
            inner_points.push_back(PointE(innerPoints[i].X(),innerPoints[i].Y()));
    }
    void get_triangulation(std::vector<PointType> &vertices, std::vector<std::vector<unsigned int>> &faces){
        CDT cdt;
        cdt.insert_constraint(outer_poly.vertices_begin(),outer_poly.vertices_end(),true);
        for(size_t i=0;i<inner_points.size();i++)
            cdt.insert(inner_points[i]);

        assert(cdt.is_valid());
        std::cout<<"triangulation done "<<std::endl;
        CDT::Finite_vertices_iterator vit;
        unsigned int index=0;
        for (vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); vit++,index++){
             vit->info()=index;
             double cx=CGAL::to_double(vit->point().x());
             double cy=CGAL::to_double(vit->point().y());
             vertices.push_back(PointType(cx,cy,0));
        }
        CDT::Finite_faces_iterator fit;
        for (fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++){
            //To discard wrong faces
            PointE centroid=CGAL::centroid(cdt.triangle(fit));
            if(outer_poly.bounded_side(centroid)==CGAL::ON_BOUNDED_SIDE){
                std::vector<unsigned int> face={fit->vertex(0)->info(),fit->vertex(1)->info(),fit->vertex(2)->info()};
                faces.push_back(face);
            }
        }
    }
private:
    Polygon_2 outer_poly;
    std::vector<PointE> inner_points;
};

#endif // CDT2D_H
