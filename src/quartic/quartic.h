/*----------------------------------------------------------------------------*/
/*		      Excessive Plane Quartic Search Code                     */
/*----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

This library of functions is used to search for excessive curves of genus 3 over
GF(q). See

  [1] Faber, Xander and Grantham, Jon and Howe, Everett. "On the Maximum
      Gonality of a Curve over a Finite Field". Preprint, 2022. 

for descriptions of the algorithms used.  

TABLE OF CONTENTS

- MACROS

  - MEMORY_ERROR
  - DEBUG
  - QUADBOUND
  - SETBITS

- UTILITIES (utilities.c)

  - format_time
  - start_and_stop_work
  - loop_printer

  - next_vector
  - start_vector
  - next_projective_vector
  - test_next_projective_vector

  - key

- FIELD ARITHMETIC AND VECTOR FUNCTIONS (fields.c)

  - field_cardinality_to_int
  - field_elements_array
  - construct_field
  - clear_field
  - irreducible_quadratic

  - ops_store_t
  - field_element_to_index
  - construct_ops_store
  - clear_ops_store

  - fq_idx_add
  - fq_idx_sub
  - fq_idx_mul
  - fq_idx_inv
  - fq_idx_square
  - fq_idx_cube
  - fq_idx_quad
  
  - fq_idx_vec_t
  - fq_idx_vec_init
  - fq_idx_vec_clear
  - fq_idx_vec_dot
  - write_index_vector_to_file
  - write_index_vector_to_stdout


- POINTS (points.c)

  - point_t
  - point_init
  - point_clear
  - point_pretty_print
  - point_eq_raw
  - point_quadratic_frobenius

  - construct_rational_points
  - construct_rational_point_arrays

- CURVES (curves.c)

  - quartic_basis_vec_t
  - quartic_basis_vec_init
  - quartic_basis_vec_clear
  - quartic_basis_vec_pretty_print
  - quartic_basis_vec_pretty_print_file

  - quartic_basis_pol_evaluations

  - point_data_test

  - quartic_has_rational_point
  - quadratic_places_on_quartic

- PGL_3 ACTION (pgl_action.c)

  - zzzz_coef
  - xxxx_coef
  - yyyy_coef
  - xxxy_coef
  - xxyy_coef
  - xyyy_coef
  - xxyz_coef
  - xyyz_coef
  - xyzz_coef
  - xxxz_coef
  - yyyz_coef
  - xxzz_coef
  - yyzz_coef
  - xzzz_coef
  - yzzz_coef

  - mark_good_translates_naive
  - mark_good_translates_smart

------------------------------------------------------------------------------*/

#ifndef QUARTIC_H
#define QUARTIC_H

#include <fmpz.h>
#include <fq.h>
#include <fq_vec.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*----------------------------------------------------------------------------*/
/*---------------------------------  MACROS  ---------------------------------*/
/*----------------------------------------------------------------------------*/

#define MEMORY_ERROR {{fprintf(stdout,"ERROR: Cannot allocate memory at line %d\n",__LINE__); exit(1);}}

#define DEBUG 0  // Set to 1 for debug printing and file output

#define QUADBOUND 530 // Max number of quadratic places for q \le 23:
                       // (q^2 + 1 + 2*3*q)/2 \le 334 in the smooth case
                       // q^2 + 1 \le 530 in the case of a product of conjugate conics

#define SETBITS 28 // log_2(Number of entries) for set/hash table construction
                   // Setting to 28 gives a table size of 2GB

/*----------------------------------------------------------------------------*/
/*-------------------------------  UTILITIES  --------------------------------*/
/*----------------------------------------------------------------------------*/

// Format time t as hrs,mins,secs and write to character array time_str
void format_time(char *time_str, time_t t);

// Assign start and stop indices for the work this job should do.  Write
// total_work = q*num_jobs + r. Assign q+1 tasks to the first r jobs and q tasks
// to each of the remaining jobs.
void start_and_stop_work(uint64_t *start, uint64_t *stop,
			 uint64_t total_work, int num_jobs, int job);


// Print a percentage if the loop index is divisible by steps_until_print. If
// num_prints is divisible by 10, also print some timing info. This updates the
// value of num_prints and loop_timer.
// Set num_prints=0 and loop_timer = time(NULL) before using this function.
void loop_printer(uint64_t loop_idx, uint64_t start, uint64_t stop,
		  uint64_t steps_until_print, int *num_prints, time_t *loop_timer);

/*----------------------------------------------------------------------------*/

// Write start_idx in base q and assign the digits to the entries of vv in
// big-endian form. No checking is done to insure start_idx is in
// [0,q)^len_vv.
void start_vector(uint16_t *vv, int len_vv, int q, int start_idx);

// View the array *vv as a vector of integers in [0,q) of length len_vv.
// Modify the vector to be the next in the lex ordering if such exists and
// return 1. The first vector in this ordering is [0,...,0]. If
// *vv = [q-1,...,q-1], then there is no next vector; return 0.
// No validity testing is done on the array.
int next_vector(uint16_t *vv, int len_vv, int bound);

// View the array *vv as a vector of integers in [0,q) of length len_vv.
// Modify the vector to be the next in the lex ordering subject to the
// additional restriction that the leftmost nonzero entry must be 1, if such
// exists, and return 1. The first vector in this ordering is [0,...,0,1].
// If *vv = [1,q-1,...,q-1], then there is no next vector; return 0.
// No validity testing is done on the array.
// WARNING: We treat 0 and 1 like they are field elements, and in order to apply
//   this to actual fields, we will need to ensure that our lists of field
//   elements start with 0,1. 
int next_projective_vector(uint16_t *vv, int len_vv, int q);

// Simple code for testing 
void test_next_projective_vector(int m, int q);

/*----------------------------------------------------------------------------*/

// A simple hash table to implement a set for storing quartics.
// Use linear probing to deal with collisions.
typedef struct {
  int64_t*   members;
  size_t num_members; // 2^SETBITS
  int64_t       mask; // 2^SETBITS - 1, for easy modular reduction
  int              q; // For base-q arithmetic
} qset_t;

// Initialize/clear a qset structure
void qset_initialize(qset_t *S, int q);
#define qset_clear(S) free(S.members)

// Add Q to the set
void qset_add_item(qset_t *S, uint16_t *Q);

// Return 1 if Q is a member of the set, 0 otherwise
int qset_is_member(qset_t *S, uint16_t *Q);

/*----------------------------------------------------------------------------*/
/*-----------------  FIELD ARITHMETIC AND VECTOR FUNCTIONS  ------------------*/
/*----------------------------------------------------------------------------*/

// Utility for grabbing the cardinality of a flint field
int field_cardinality_to_int(fq_ctx_t FF);

// Utility for grabbing the cardinality of flint field representing GF(q^2)
int get_q(fq_ctx_t GFq_squared);

// Set rop = op^q inside GF(q^2)
void quadratic_frobenius(fq_t rop, fq_t op, fq_ctx_t GFq_squared);

// Element op is in GF(q^2); return 1 if op is in GF(q) and 0 otherwise.
int fq_is_rational(fq_t op, fq_ctx_t GFq_squared);

// Allocate memory for and populate an array of the field elements in
// GF(q^2). We order the elements of GF(q^2) as follows. If p is the
// characteristic, then we write an integer i in base p as
//   i = j_0 + j_1*p + j_2*p^2 + ...
// The i-th element of GF(q^2) is then given by
//   j_0 + j_1*a + j_2*a^2 + ...,  where GF(q^2)= GF(p)(a).
// We record these elements in the array *elts in order, subject to the
// additional restriction that we put all of the elements of GF(q) first.
// In particular, the i-th element is i for i < p.
void field_elements_array(fq_t **elts, fq_ctx_t GFq_squared);

// Initialize the field GF(q^2), allocate memory for the array of its elements,
// and populate the array. This errors out if q is not a prime power <= 32
// WARNING/FIXME: This messes up many calls to field_cardinality_to_int
void construct_quadratic_extension(fq_ctx_t GFq_squared, fq_t **elts, int q);

// Free up all memory involved in elts and FF
void clear_field(fq_ctx_t FF, fq_t *elts);

// Set a to be the first elements of elts such that a \in GF(q) and
// x^2 + ax + 1 is irreducible over GF(q). 
void irreducible_quadratic(fq_t a, fq_t *elts, fq_ctx_t GFq_squared);

/*----------------------------------------------------------------------------*/

// Arithmetic operations tables for field indices for GF(q^2)
typedef struct {
  int               q; // Cardinality of field GF(q)
  int              qq; // Cardinality of GF(q^2)
  int               p; // Field characteristic
  uint16_t          a; // Index of a where T^2 + a*T + 1 is irreducible over GF(q)
  uint16_t       root; // Root of T^2 + a*T + 1
  uint16_t *add_table; // q^2 entries for addition table
  uint16_t *sub_table; // q^2 entries for addition table  
  uint16_t *mul_table; // q^2 entries for subtraction table
  uint16_t     *conjs; // q entries for conjugates array    
  uint16_t      *invs; // q entries for inverse array  
  uint16_t   *squares; // q entries for squares array
  uint16_t     *cubes; // q entries for cubes array
  uint16_t     *quads; // q entries for quads array
} ops_store_t;

// The element elt appears at position j in the array fld_elts; return j.
// Return cardinality of FF on error. 
uint16_t field_element_to_index(fq_t elt, fq_t *fld_elts, fq_ctx_t FF);

// Construct arithmetic operations tables for field GF(q^2), whose elements
// occupy the array fld_elts. Space will be allocated for the entries in the
// tables. The sum/difference/product of fld_elts[i] and fld_elts[j] is stored
// at position i + qq*j, where qq = q^2. The inverse of fld_elts[i] is stored at
// position i in the array negs, and similarly for squares, cubes, and
// quads. The element 0 is stored at position 0 of the inverse array. We also
// pass in a special element a whose index is stored, and a root of
// T^2 + aT + 1 is computed and stored. 
void construct_ops_store(ops_store_t *tables, fq_t a,
			 fq_t *fld_elts, fq_ctx_t GFq_squared);

// Clear all memory associated with tables
void clear_ops_store(ops_store_t *tables);

// Print tables and arrays in a useful manner for debugging
void print_ops_store(ops_store_t *tables);

// Binary operations using indices; return the index given by combining
// op1 and op2 using the specified operation. 
uint16_t fq_idx_add(uint16_t *op1, uint16_t *op2, ops_store_t *tables);
uint16_t fq_idx_sub(uint16_t *op1, uint16_t *op2, ops_store_t *tables);
uint16_t fq_idx_mul(uint16_t *op1, uint16_t *op2, ops_store_t *tables);
uint16_t fq_idx_mul_ui(uint16_t *op1, int16_t op2, ops_store_t *tables);

// Unary operations using indices; return the index given by applying the
// specified operator to op. (Inverse returns 0 on input 0).
uint16_t fq_idx_conj(uint16_t *op, ops_store_t *tables);   // q-power frobenius
uint16_t fq_idx_inv(uint16_t *op, ops_store_t *tables);    // inverse
uint16_t fq_idx_square(uint16_t *op, ops_store_t *tables); // square
uint16_t fq_idx_cube(uint16_t *op, ops_store_t *tables);   // cube
uint16_t fq_idx_quad(uint16_t *op, ops_store_t *tables);   // 4-th power

/*----------------------------------------------------------------------------*/

// Vector of field indices
typedef uint16_t* fq_idx_vec_t;

// Allocate memory for a table vector with num_entries entries
void fq_idx_vec_init(fq_idx_vec_t *vec, int num_entries);

// Free memory from a coefficient vector 
#define fq_idx_vec_clear(vec) free(vec)

// Dot product of vectors op1 and op2, with all arithmetic performed using field
// indices
uint16_t fq_idx_vec_dot(fq_idx_vec_t op1, fq_idx_vec_t op2, int len,
			ops_store_t *tables);

// Each vector is written as a space-separated sequence of indices. If writing
// to file, a newline is added. 
void write_index_vector_to_file(FILE *fp, fq_idx_vec_t vv, int len);
void write_index_vector_to_stdout(fq_idx_vec_t vv, int len);

/*----------------------------------------------------------------------------*/
/*---------------------------------  POINTS  ---------------------------------*/
/*----------------------------------------------------------------------------*/

// Point of PP^2(GF(q))
typedef struct {
  fq_t x;
  fq_t y;
  fq_t z;
} point_t;

// Initialize and clear P in PP^2(FF)
void point_init(point_t *P, fq_ctx_t FF);
void point_clear(point_t *P, fq_ctx_t FF);

// Print P to stdout in (*, *, *) nice format (with no line break)
void point_pretty_print(point_t *P, fq_ctx_t FF);

// Return 1 if P = Q as points in AA^3(FF), 0 otherwise.  (This can be used to
// test equality of projective points if we already know they are normalized.)
int point_eq_raw(const point_t *P,const point_t *Q, fq_ctx_t FF);

// Set Q = f(P), where f is the q-Frobenius map on FF = GF(q^2).
void point_quadratic_frobenius(point_t *Q, const point_t *P, fq_ctx_t FF);

/*----------------------------------------------------------------------------*/

// Construct four arrays of rational points for the different phases of testing
// xpts are the   q points on the line x=0 not equal to (0:0:1)
// ypts are the   q points on the line y=0 not equal to (0:0:1)
// zpts are the q-1 points on the line z=0 not equal to (1:0:0) or (0:1:0)
// genpts are the remaining (q-1)^2 points not on a coordinate line
void construct_rational_point_arrays(point_t **xpts, point_t **ypts, point_t **zpts,
				     point_t **genpts, fq_t *elts,
				     fq_ctx_t GFq_squared);

// Allocate memory for and store one representative for each Galois orbit of
// points in P^2(GF(q^2)) - P^2(GF(q)). 
void construct_quadratic_points(point_t **pts, fq_t *elts, fq_ctx_t GFq_squared);

// Allocate space and set array of points in index form. We allocate QUADBOUND
// slots, which is enough for q \le 23 by the Weil bound plus some analysis of
// geometrically reducible quartics.
void quad_point_indices(fq_idx_vec_t **point_inds, point_t *quad_pts,
			fq_t *fld_elts, fq_ctx_t GFq_squared);

// Return 0 if the distinct, non-conjugate quadratic points P and Q lie on a
// line over GF(q); return 1 otherwise. (In the latter case, the corresponding
// quadratic places are in general position.)
int quad_places_in_general_position(fq_idx_vec_t P, fq_idx_vec_t Q,
				    ops_store_t *tables);

/*----------------------------------------------------------------------------*/
/*---------------------------------  CURVES  ---------------------------------*/
/*----------------------------------------------------------------------------*/

// Special Quartic Polynomial Basis (ordered):
// [z^4 + x^2*z^2 + a*x*z^3 + y^2*z^2 + a*y*z^3,
//  x^4 + a*x^3*z + x^2*z^2,
//  x^3*z + a*x^2*z^2 + x*z^3,
//  y^4 + a*y^3*z + y^2*z^2,
//  y^3*z + a*y^2*z^2 + y*z^3,
//  x^3*y,
//  x^2*y^2,
//  x*y^3,
//  x^2*y*z,
//  x*y^2*z,
//  x*y*z^2]
//  where the polynomial T^2 + aT + 1 is irreducible over GF(q)

// Struct for holding vector of coefficients or evaluations of 
// this special basis.
typedef fq_struct* quartic_basis_vec_t;

void quartic_basis_vec_init(quartic_basis_vec_t *vec, fq_ctx_t FF);
void quartic_basis_vec_clear(quartic_basis_vec_t *vec, fq_ctx_t FF);

// Print vv to standard out
void quartic_basis_vec_pretty_print(quartic_basis_vec_t vv, fq_ctx_t FF);

// Print vv to file
void quartic_basis_vec_pretty_print_file(FILE* fp, quartic_basis_vec_t vv,
					 fq_ctx_t FF);

/*----------------------------------------------------------------------------*/

// Allocate memory and store the vector of indices of the quartic monomial
// evaluations v_P = (m(P) : m \in quartic basis polynomials), for P in pts
void quartic_basis_pol_evaluations(fq_idx_vec_t **vecs, point_t *pts, int num_pts,
				   fq_t *elts, fq_t a, fq_ctx_t FF);

/*----------------------------------------------------------------------------*/

// Write debug data to file for use in Sage script debug_tests.sage
void point_data_test(FILE *fp, fq_t *elts,
		     point_t *xpts, point_t *ypts, point_t *zpts, point_t *genpts,
		     fq_idx_vec_t *xpt_data, fq_idx_vec_t *ypt_data,
		     fq_idx_vec_t *zpt_data, fq_idx_vec_t *genpt_data,
		     point_t *quadpts, fq_idx_vec_t *quadpt_data, 
		     fq_t a, fq_ctx_t GFq_squared);

/*----------------------------------------------------------------------------*/

// Return 1 if the quartic with coefficients given by Q, and 0 otherwise. For
// speed, all arithmetic is done with field indices rather than elements.
int quartic_has_rational_point(uint16_t *Q, fq_idx_vec_t *pt_data, int num_pts,
			       ops_store_t *tables);


// Compute all quadratic places on the quartic Q, stored as a list of indices
// into the full list of quadratic places of the plane. The array Qquad_inds
// holds these indices, and the value of num_places holds how many quadratic
// places there are on Q. We assume that 80 slots have been allocated for
// Qplace_inds, which is enough for q=23 by the Weil bound.
void quadratic_places_on_quartic(int *Qquad_inds, int *num_places,
				 uint16_t * Q, fq_idx_vec_t *quadpt_data,
				 fq_idx_vec_t *quadpt_inds, ops_store_t *tables);


/*----------------------------------------------------------------------------*/
/*				 MATRIX ENTRIES                               */
/*----------------------------------------------------------------------------*/

// Compute the coefficients of the matrix M in PGL(3,q) that maps the quadratic
// points pt1, pt2 to the points (c:0:1) and (0:c:1). The code for these
// functions was automatically generated by the Sage script
// quad_point_helper.sage.  The formulas for these entries are only guaranteed
// to be correct if pt1 and pt2 are in general position. 
uint16_t quad_transmat_entry00(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables);
uint16_t quad_transmat_entry01(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables);
uint16_t quad_transmat_entry02(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables);
uint16_t quad_transmat_entry10(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables);
uint16_t quad_transmat_entry11(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables);
uint16_t quad_transmat_entry12(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables);
uint16_t quad_transmat_entry20(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables);
uint16_t quad_transmat_entry21(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables);
uint16_t quad_transmat_entry22(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables);

/*----------------------------------------------------------------------------*/
/*				  PGL_3 ACTION                                */
/*----------------------------------------------------------------------------*/

// Compute the coefficients on the various monomials for F(M(x,y,z)), for F =
// cc(x,y,z) a quartic polynomial and M = matrix(mm). The code for these
// functions was automatically generated by the Sage script pgl3_helper.sage
// using explicit formulas for these coefficients (quartic in coefficients mm
// and linear in cc).
uint16_t zzzz_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t xxxx_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t yyyy_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t xxxy_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t xxyy_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t xyyy_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t xxyz_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t xyyz_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t xyzz_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t xxxz_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t yyyz_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t xxzz_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t yyzz_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t xzzz_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);
uint16_t yzzz_coef(fq_idx_vec_t cc, fq_idx_vec_t mm, ops_store_t *tables);

/*----------------------------------------------------------------------------*/

// Simple test whether the above coefficient functions are computing the correct
// thing. We use a small number of simple quartices and the identity map from
// PGL(3,q)
void pgl_tests(ops_store_t *tables);

/*----------------------------------------------------------------------------*/

// Act on a quartic curve F by all elements of PGL(3,q) and record which of
// these translates are of the desired form in the seen array. This looks at
// O(q^8) elements of PGL(3,q). 
void mark_good_translates_naive(qset_t *seen, uint16_t *Q, ops_store_t *tables);  

// Act on a quartic curve F by all elements of PGL(3,q) that preserve the
// property "passes through (c:0:1) and (0:c:1)". Record which of these
// translates are of the desired form in the seen array. Return 1 if the number
// of quadratic places satisfies the Weil bounds for a genus-3 curve and 0
// otherwise. The information need to rapidly determine the quadratic points on
// the curve is given in quadpt_data (monomial evaluations at quadratic points)
// and quadpt_inds (the full list of quadratic points of the plane, up to Galois
// action). This looks at O(q^4) elements of PGL(3,q)
int mark_good_translates_smart(qset_t *seen, uint16_t *Q, fq_idx_vec_t *quadpt_data,
			       fq_idx_vec_t *quadpt_inds, ops_store_t *tables);

#endif /* QUARTIC_H */
