/*----------------------------------------------------------------------------*/
/*			  Excessive Curve Search Code                         */
/*----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

An excessive curve is a curve over a finite field of genus g and gonality g+1.
This library of functions is used to search for excessive curves of genus 6 or 7
over GF(2) that are represented as plane curves of specific degrees. See 

  [1] Faber, Xander and Grantham, Jon and Howe, Everett. "On the Maximum
      Gonality of a Curve over a Finite Field". Preprint, 2021. 

for descriptions of the algorithms used.  


TABLE OF CONTENTS

- MACROS

  - MEMORY_ERROR
  - DEBUG
  - QUEUE_VECS

- UTILITIES (utilities.c)

  - format_time
  - start_and_stop_work

- FIELD ARITHMETIC AND VECTOR FUNCTIONS (fields.c)

  - field_cardinality_to_long
  - field_elements_array
  - construct_field
  - clear_field

  - table_t
  - construct_op_table
  - fq_print_table
  - fq_op_table
  - fq_sum_terms_table

- POINTS (points.c)

  - point_t
  - point_init
  - point_clear
  - point_set
  - point_pretty_print
  - point_normalize
  - point_eq_raw
  - point_eq
  - point_frobenius
  - galois_orbits

- BIT VECTOR ITERATORS (next_vec.c)

  - next_lex_vector_of_same_weight_with_offset
  - next_odd_weight_vector_with_offset
  - print_bits
  - set_first_odd_weight_vector_with_offset
  - test_next_odd_weight_vector_with_offset
  - odd_weight_vector_with_offset_and_index

- PLANE CURVES (curves.c)

  - monomial_t
  - construct_monomials
  - coef_vec_t

  - table_vec_t
  - table_vec_init
  - table_vec_clear

  - point_data_t
  - point_data_init
  - point_data_clear
  - point_data_populate_all
  - point_data_test_short
  - point_data_test_full
  - point_is_smooth_on_curve

- PRINT QUEUE  (queue.c)

  - print_queue_t;
  - print_queue_init
  - print_queue_clear
  - write_vector_to_string
  - print_queue_purge
  - print_queue_write  


------------------------------------------------------------------------------*/

#ifndef PLANE_H
#define PLANE_H

#include <fmpz.h>
#include <fq.h>
#include <fq_vec.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*----------------------------------------------------------------------------*/
/*---------------------------------  MACROS  ---------------------------------*/
/*----------------------------------------------------------------------------*/

#define MEMORY_ERROR {{fprintf(stdout,"ERROR: Cannot allocate memory at line %d\n",__LINE__); exit(1);}}

#define DEBUG 0  // Set to 1 for debug printing

#define QUEUE_VECS 1024 // Maximum number of vectors in the print queue

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

/*----------------------------------------------------------------------------*/
/*-----------------  FIELD ARITHMETIC AND VECTOR FUNCTIONS  ------------------*/
/*----------------------------------------------------------------------------*/

// Utility for grabbing the cardinality of a flint field
int field_cardinality_to_int(fq_ctx_t FF);
  
// Allocate memory for an array of the field elements in FF and populate the
// array. If FF has characteristic p, then the i-th element of the array is
// given by writing
//   i = j_0 + j_1*p + j_2*p^2 + ...
// and setting
//   *elts[i] = j_0 + j_1*a + j_2*a^2 + ...
// where FF = GF(p)(a). 
void field_elements_array(fq_t **elts, fq_ctx_t FF);

// Construct the field GF(2^r) and make an array of its elements, with the first
// two elements being 0 and 1.  We use Conway polynomials to construct the
// extension. The ordering of the field elements is given in the description of
// field_elements_array.
void construct_field(fq_ctx_t FF, fq_t **elts, int r);

// Free up all memory involved in elts and FF
void clear_field(fq_ctx_t FF, fq_t *elts);

/*----------------------------------------------------------------------------*/

typedef struct {
  int          width; // table width
  uint16_t  *entries; // width^2 entries in table
} table_t;

// Construct addition/multiplication tables for field FF, whose elements occupy
// the array fld_elts.  Space will be allocated for the entries in the
// table. The sum/product of fld_elts[i] and fld_elts[j] is stored at position
// i + q*j, where q is the field cardinality. This constructs an addition table
// if is_add is nonzero, and otherwise it constructs a multiplication table. 
void construct_op_table(table_t *add_table, fq_t *fld_elts, fq_ctx_t FF, int is_add);

// Very simple table printing, for debug purposes
void fq_print_table(table_t *table);

// Combine field indices for GF(q) using a table constructed with
// construct_op_table. Working with indices and a lookup table is far more
// efficient than performing field arithmetic on small fields. 
uint16_t fq_op_table(uint16_t *op1, uint16_t *op2, table_t *table);

// Return a specified sum of field elements. If the i-th bit of mask is set, we
// add the field element whose table index is inds[i]. For maximum efficiency,
// this function DOES NOT check that mask and inds have compatible lengths.
uint16_t fq_sum_terms_table(uint64_t mask, uint16_t *inds, int len, table_t *add_table);

/*----------------------------------------------------------------------------*/

// Vector of table coefficients
typedef uint16_t* table_vec_t;

// Allocate memory for a table vector with num_entries entries
void table_vec_init(table_vec_t *vec, int num_entries);

// Free memory from a coefficient vector 
#define table_vec_clear(vec) free(vec)

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

// Set P = Q
void point_set(point_t *P, point_t *Q, fq_ctx_t FF);

// Print P to stdout in (*, *, *) nice format (with no line break)
void point_pretty_print(point_t *P, fq_ctx_t FF);

// Rescale P in place so that left-most nonzero coordinate is 1
void point_normalize(point_t *P, fq_ctx_t FF);

// Return 1 if P = Q as points in AA^3(FF), 0 otherwise.
// (This can be used to test equality of projective points
// if we already know they are normalized.)
int point_eq_raw(const point_t *P, const point_t *Q, fq_ctx_t FF);

// Return 1 if P = Q, 0 otherwise.
// This normalizes P and Q in the process.
int point_eq(point_t *P, point_t *Q, fq_ctx_t FF);

// Set Q = f(P), where f is the p-power Frobenius map on FF = GF(p^*)
void point_frobenius(point_t *Q, const point_t *P, fq_ctx_t FF);

// Initialize and construct an array of representatives of Galois orbits of
// points in PP^2(FF). Return the number of orbits found. If rat_pts=0,
// do not include Frobenius-invariant points among the representatives.
int galois_orbits(point_t **orbits, fq_ctx_t FF, fq_t *elts, int rat_pts);


/*----------------------------------------------------------------------------*/
/*--------------------------  BIT VECTOR ITERATORS  --------------------------*/
/*----------------------------------------------------------------------------*/

// Modify the substring of *vv of length len with given offset (from right)
// to be the next lex vector of the same Haming weight. Return 1 on success
// and 0 if there is no such substring (or fail). Example:
//    *vv = 00101011, len = 4, offset = 2.
// The updated value is 
//    *vv = 00110011.
int next_lex_vector_of_same_weight_with_offset(uint64_t *vv, int len, int offset);

// Modify the substring of *vv of length len with given offset (from right) to
// be the next lex vector of odd Hamming weight. Here v > w if wt(v) > wt(w) or
// if wt(v) = wt(w) and v > w for the lex ordering. No validity testing is done
// on the input word.
int next_odd_weight_vector_with_offset(uint64_t *vv, int len, int offset);

// Print the bit representation of n, using at least min_bits to do so
void print_bits(uint64_t n, int min_bits);

// Set field of length len at given offset to be 0...01.
void set_first_odd_weight_vector_with_offset(uint64_t *vv, int len,
					     int offset);

// Code to test next_odd_weight_vector_with_offset: print all odd-weight binary
// vectors of length total_len with active field of length len and given offset.
// Remaining entries are populated with random bits. 
void test_next_odd_weight_vector_with_offset(int total_len, int len,
					     int offset);

// Set the active field of *vv to be the j-th odd-weight vector of length len in
// the lex ordering, where j = index, and return 1. If no such vector exists,
// return 0. (This implementation just advances through the odd_weight vectors
// until it gets to the one with given index.)
int odd_weight_vector_with_offset_and_index(uint64_t *vv, int len,
					    int offset, uint64_t index);

/*----------------------------------------------------------------------------*/
/*------------------------------  PLANE CURVES  ------------------------------*/
/*----------------------------------------------------------------------------*/

// Exponents of a monomial x^a * y^b * z^c
typedef struct {
  int a;
  int b;
  int c;
} monomial_t;

// Allocate memory for the array of monomials of degree d=deg and construct the
// list in a particular order:
// x^d, y^d, z^d,
// x^{d-1}y, ..., xy^{d-1},
// x^{d-1}z, ..., xz^{d-1},
// y^{d-1}z, ..., yz^{d-1},
// xyz^(d-2), xy^2z^{d-3}, ..., xy^{d-2}z, x^2yz^{d-3}, ..., x^{d-2}yz
// In the last line, we have terms ^a y^b z^c with a+b+c=d and min(a,b,c) > 0. 
// Return the number of such monomials, (d+2)*(d+1)/2.
int construct_monomials(monomial_t **terms, int deg);

// Coefficient vector for a polynomial over GF(2) with at most 64 terms, in an
// implicit basis. (Make the basis explicit when working with these!)
typedef uint64_t coef_vec_t;

// Write a GF(2)-polynomial in the standard monomial basis to string.  Monomials
// are given in terms, and coefficients are given as the bits of coefs. 
void coef_vec_to_poly_str(char *str, monomial_t *terms, coef_vec_t coefs);

/*----------------------------------------------------------------------------*/

// Point of PP^2(FF), augmented with its table representation, a list of values
// of polynomials evaluated at the point, and a list of the x-,y-,z-derivatives
// of those same polynomials evaluated at the point. The list of polynomials is
// unspecified in this data structure. (In practice, it is either the list of
// monomials in standard order, or a special basis of polynomials.)
typedef struct
{
  point_t *pt;         // pointer to the actual point
  uint16_t  x;         // table index of x-coordinate
  uint16_t  y;         // table index of y-coordinate
  uint16_t  z;         // table index of z-coordinate
  int len;             // vector length
  table_vec_t pt_vals; // vector of table indices of polynomial values
  table_vec_t x_deriv; // vector of table indices of values of x-derivatives
  table_vec_t y_deriv; // vector of table indices of values of y-derivatives
  table_vec_t z_deriv; // vector of table indices of values of z-derivatives
} point_data_t;

// Allocate/clear space for coefficient vector fields in point_data_t object
void point_data_init(point_data_t *data, int num_terms);
void point_data_clear(point_data_t *data);

// Allocate space for an array of point_data_t objects and populate the
// data for the *standard* basis of monomials (see above). 
void point_data_populate_all(point_data_t **data, point_t *pts, int num_pts,
			     fq_t *elts, fq_ctx_t FF, table_t *mul_table,
			     monomial_t *terms, int num_terms);

// Write polynomial evaluation data for each point to file fp for testing with
// the Sage script point_test.sage
void point_data_test_short(FILE *fp, point_data_t *data, int num_pts);

// Write additional data to file fp for testing with the Sage script
// point_test.sage (This function does more testing than the short test.)
void point_data_test_full(FILE *fp, point_data_t *data, int num_pts,
			  fq_ctx_t FF, fq_t *elts);

// Return 1 if P is a smooth point on the plane curve defined by pol, and 0
// otherwise.
int point_is_smooth_on_curve(point_data_t *P, coef_vec_t pol, table_t *add_table);

/*----------------------------------------------------------------------------*/
/*------------------------------  PRINT QUEUE  -------------------------------*/
/*----------------------------------------------------------------------------*/

typedef struct
{
  int num_vecs;                // Number of vectors currently in the queue
  coef_vec_t vecs[QUEUE_VECS]; // Array of vectors
} print_queue_t;

// Initialize/clear all of the fields in the print queue. 
void print_queue_init(print_queue_t *queue);
void print_queue_clear(print_queue_t *queue);

// Write all vector coefficients in the queue to the file object fp
// as a string of entries, each vector on its own line. 
void print_queue_purge(print_queue_t *queue, FILE *fp);

// Write the vector v to the print queue. If the print queue is full,
// purge the queue to the file fp. 
void print_queue_write(print_queue_t *queue, coef_vec_t v, FILE *fp);


#endif /* PLANE_H */
