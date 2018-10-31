#ifdef __cplusplus
extern "C" {
#endif
void libsupermesh_sort_intersection_finder_set_input(long* nnodes_a, int* dim_a, long* nelements_a, int* loc_a, long* nnodes_b, int* dim_b, long* nelements_b, int* loc_b, double* positions_a, long* enlist_a, double* positions_b, long* enlist_b);
void libsupermesh_sort_intersection_finder_query_output(long* nindices);
void libsupermesh_sort_intersection_finder_get_output(long* nelements, long* nindices, long* indices, long* ind_ptr);
void libsupermesh_tree_intersection_finder_set_input(long* nnodes_a, int* dim_a, long* nelements_a, int* loc_a, long* nnodes_b, int* dim_b, long* nelements_b, int* loc_b, double* positions_a, long* enlist_a, double* positions_b, long* enlist_b);
void libsupermesh_tree_intersection_finder_query_output(long* nindices);
void libsupermesh_tree_intersection_finder_get_output(long* nelements, long* nindices, long* indices, long* ind_ptr);
void libsupermesh_intersect_tris_real(double* tri_a, double* tri_b, double* tris_c, int* n_tris_c);
void libsupermesh_intersect_tets_real(double* tet_a, double* tet_b, double* tets_c, int* n_tets_c);
void libsupermesh_triangle_area(double* tri, double* area);
void libsupermesh_tetrahedron_volume(double* tet, double* volume);
void libsupermesh_interval_size(double* interval, double* size);
void libsupermesh_intersect_intervals(double* interval_a, double* interval_b, double* intervals_c, int* n_intervals_c);
void libsupermesh_intersect_tet_hex(double* tet_a, double* hex_b, double* tets_c, int* n_tets_c);
void libsupermesh_intersect_tri_quad(double* tri_a, double* quad_b, double* tris_c, int* n_tris_c);

#ifdef __cplusplus
}

int libsupermesh_max_n_tris_c(int n_lines_b, int n_lines_a = 3){
   int out = n_lines_a*std::pow(2, n_lines_b) - 2;
   return out;
}

int libsupermesh_max_n_tets_c(int n_planes_b){
   int out = std::pow(3, n_planes_b);
   return out;
}

#endif
