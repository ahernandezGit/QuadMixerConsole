#ifndef PRECOMPILED_LIBRARIES_H
#define PRECOMPILED_LIBRARIES_H

#include <QTemporaryFile>
#include <QDir>
#include <QLibraryInfo>
#include <ctime>
#include <cstring>
#include <fstream>
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <list>
#include <memory>



#include <CGAL/config.h>
#include <igl/copyleft/cgal/mesh_boolean.h>

#include "vcg/space/intersection2.h"
#include <vcg/complex/algorithms/intersection.h>
#include <vcg/complex/algorithms/curve_on_manifold.h>
#include <wrap/io_trimesh/export.h>
#include<vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/simplex/face/jumping_pos.h>
#include<vcg/complex/algorithms/isotropic_remeshing.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/append.h>
#include <vcg/space/intersection2.h> // for IsInsideTrianglePoint()
#include <vcg/space/segment2.h>      // for using segment2
#include <vcg/space/triangle2.h>     // for using triangle2
 #include <vcg/space/triangle3.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/refine_loop.h>
#include <vcg/complex/algorithms/edge_collapse.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/inside.h>
#include <vcg/space/intersection2.h>
#include <vcg/space/triangle3.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/refine_loop.h>
#include <vcg/complex/algorithms/update/selection.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/complex/algorithms/parametrization/poisson_solver.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>
#include<vcg/complex/algorithms/isotropic_remeshing.h>
#include<vcg/complex/algorithms/smooth.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/intersections.h>

#endif // PRECOMPILED_LIBRARIES_H
