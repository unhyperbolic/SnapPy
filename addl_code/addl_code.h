/* Prototypes for functions defined in gluing_equations.c */

extern int** get_gluing_equations(Triangulation* manifold, 
				  int* num_rows, 
				  int* num_cols);

extern void free_gluing_equations(int** equations, 
				  int num_rows);

extern int* get_cusp_equation(Triangulation* manifold,
			      int cusp_num,
			      int m,
			      int l,
			      int* num_rows);

extern void free_cusp_equation(int* equation);

/* Prototype for the function defined in load_link_proj.cc */

extern Triangulation*    triangulate_link_complement_from_file(char* file_name,
							       char* path);

/* Prototype for the function defined in braid.cc */

extern Triangulation* fibered_manifold_associated_to_braid(int numStrands,
						           int braidLength,
							   int* word);
/* Prototypes for functions defined in set_tet_shapes.c */

extern void set_tet_shapes(Triangulation* manifold, Complex *shapes);
extern void set_target_holonomy(Triangulation* manifold,
                                int            theCuspIndex,
                                Complex        theTarget,
                                int            theRecomputeFlag);

/* Prototype for the function defined in dt2snap.cc */

extern Triangulation* DT2Triangulation(char* c_link_record);

/* Prototype for the functions defined in choose_gen_helper.c */

extern void choose_gen_tetrahedron_info(Triangulation    *manifold, 
				 int tet_index, 
				 int *generator_path,
				 int *face0_gen,
				 int *face1_gen,
				 int *face2_gen,
				 int *face3_gen,
				 Complex *corner0, 
				 Complex *corner1, 
				 Complex *corner2, 
					Complex *corner3);
