/*----------------------------------------------------------------------------*/
/*				 MATRIX ENTRIES                               */
/*----------------------------------------------------------------------------*/

#include "quartic.h"

/*----------------------------------------------------------------------------*/

uint16_t quad_transmat_entry00(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables)
{
  int p = tables->p;
  uint16_t c = tables->root;
  uint16_t rop, s, t;
  uint16_t pts[12];
   pts[0] = pt1[0];
   pts[1] = pt1[1];
   pts[2] = pt1[2];
   pts[3] = tables->conjs[pt1[0]];
   pts[4] = tables->conjs[pt1[1]];
   pts[5] = tables->conjs[pt1[2]];
   pts[6] = pt2[0];
   pts[7] = pt2[1];
   pts[8] = pt2[2];
   pts[9] = tables->conjs[pt2[0]];
  pts[10] = tables->conjs[pt2[1]];
  pts[11] = tables->conjs[pt2[2]];
  rop = 0;
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  return rop;
}

/*----------------------------------------------------------------------------*/

uint16_t quad_transmat_entry01(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables)
{
  int p = tables->p;
  uint16_t c = tables->root;
  uint16_t rop, s, t;
  uint16_t pts[12];
   pts[0] = pt1[0];
   pts[1] = pt1[1];
   pts[2] = pt1[2];
   pts[3] = tables->conjs[pt1[0]];
   pts[4] = tables->conjs[pt1[1]];
   pts[5] = tables->conjs[pt1[2]];
   pts[6] = pt2[0];
   pts[7] = pt2[1];
   pts[8] = pt2[2];
   pts[9] = tables->conjs[pt2[0]];
  pts[10] = tables->conjs[pt2[1]];
  pts[11] = tables->conjs[pt2[2]];
  rop = 0;
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  return rop;
}

/*----------------------------------------------------------------------------*/

uint16_t quad_transmat_entry02(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables)
{
  int p = tables->p;
  uint16_t c = tables->root;
  uint16_t rop, s, t;
  uint16_t pts[12];
   pts[0] = pt1[0];
   pts[1] = pt1[1];
   pts[2] = pt1[2];
   pts[3] = tables->conjs[pt1[0]];
   pts[4] = tables->conjs[pt1[1]];
   pts[5] = tables->conjs[pt1[2]];
   pts[6] = pt2[0];
   pts[7] = pt2[1];
   pts[8] = pt2[2];
   pts[9] = tables->conjs[pt2[0]];
  pts[10] = tables->conjs[pt2[1]];
  pts[11] = tables->conjs[pt2[2]];
  rop = 0;
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  return rop;
}

/*----------------------------------------------------------------------------*/

uint16_t quad_transmat_entry10(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables)
{
  int p = tables->p;
  uint16_t c = tables->root;
  uint16_t rop, s, t;
  uint16_t pts[12];
   pts[0] = pt1[0];
   pts[1] = pt1[1];
   pts[2] = pt1[2];
   pts[3] = tables->conjs[pt1[0]];
   pts[4] = tables->conjs[pt1[1]];
   pts[5] = tables->conjs[pt1[2]];
   pts[6] = pt2[0];
   pts[7] = pt2[1];
   pts[8] = pt2[2];
   pts[9] = tables->conjs[pt2[0]];
  pts[10] = tables->conjs[pt2[1]];
  pts[11] = tables->conjs[pt2[2]];
  rop = 0;
  s = (2*p +(1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  return rop;
}

/*----------------------------------------------------------------------------*/

uint16_t quad_transmat_entry11(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables)
{
  int p = tables->p;
  uint16_t c = tables->root;
  uint16_t rop, s, t;
  uint16_t pts[12];
   pts[0] = pt1[0];
   pts[1] = pt1[1];
   pts[2] = pt1[2];
   pts[3] = tables->conjs[pt1[0]];
   pts[4] = tables->conjs[pt1[1]];
   pts[5] = tables->conjs[pt1[2]];
   pts[6] = pt2[0];
   pts[7] = pt2[1];
   pts[8] = pt2[2];
   pts[9] = tables->conjs[pt2[0]];
  pts[10] = tables->conjs[pt2[1]];
  pts[11] = tables->conjs[pt2[2]];
  rop = 0;
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  return rop;
}

/*----------------------------------------------------------------------------*/

uint16_t quad_transmat_entry12(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables)
{
  int p = tables->p;
  uint16_t c = tables->root;
  uint16_t rop, s, t;
  uint16_t pts[12];
   pts[0] = pt1[0];
   pts[1] = pt1[1];
   pts[2] = pt1[2];
   pts[3] = tables->conjs[pt1[0]];
   pts[4] = tables->conjs[pt1[1]];
   pts[5] = tables->conjs[pt1[2]];
   pts[6] = pt2[0];
   pts[7] = pt2[1];
   pts[8] = pt2[2];
   pts[9] = tables->conjs[pt2[0]];
  pts[10] = tables->conjs[pt2[1]];
  pts[11] = tables->conjs[pt2[2]];
  rop = 0;
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  return rop;
}

/*----------------------------------------------------------------------------*/

uint16_t quad_transmat_entry20(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables)
{
  int p = tables->p;
  uint16_t c = tables->root;
  uint16_t rop, s, t;
  uint16_t pts[12];
   pts[0] = pt1[0];
   pts[1] = pt1[1];
   pts[2] = pt1[2];
   pts[3] = tables->conjs[pt1[0]];
   pts[4] = tables->conjs[pt1[1]];
   pts[5] = tables->conjs[pt1[2]];
   pts[6] = pt2[0];
   pts[7] = pt2[1];
   pts[8] = pt2[2];
   pts[9] = tables->conjs[pt2[0]];
  pts[10] = tables->conjs[pt2[1]];
  pts[11] = tables->conjs[pt2[2]];
  rop = 0;
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  return rop;
}

/*----------------------------------------------------------------------------*/

uint16_t quad_transmat_entry21(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables)
{
  int p = tables->p;
  uint16_t c = tables->root;
  uint16_t rop, s, t;
  uint16_t pts[12];
   pts[0] = pt1[0];
   pts[1] = pt1[1];
   pts[2] = pt1[2];
   pts[3] = tables->conjs[pt1[0]];
   pts[4] = tables->conjs[pt1[1]];
   pts[5] = tables->conjs[pt1[2]];
   pts[6] = pt2[0];
   pts[7] = pt2[1];
   pts[8] = pt2[2];
   pts[9] = tables->conjs[pt2[0]];
  pts[10] = tables->conjs[pt2[1]];
  pts[11] = tables->conjs[pt2[2]];
  rop = 0;
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(-1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = (2*p +(1))%p;
  t = fq_idx_mul_ui(tables->squares+c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  return rop;
}

/*----------------------------------------------------------------------------*/

uint16_t quad_transmat_entry22(fq_idx_vec_t pt1, fq_idx_vec_t pt2, ops_store_t *tables)
{
  int p = tables->p;
  uint16_t c = tables->root;
  uint16_t rop, s, t;
  uint16_t pts[12];
   pts[0] = pt1[0];
   pts[1] = pt1[1];
   pts[2] = pt1[2];
   pts[3] = tables->conjs[pt1[0]];
   pts[4] = tables->conjs[pt1[1]];
   pts[5] = tables->conjs[pt1[2]];
   pts[6] = pt2[0];
   pts[7] = pt2[1];
   pts[8] = pt2[2];
   pts[9] = tables->conjs[pt2[0]];
  pts[10] = tables->conjs[pt2[1]];
  pts[11] = tables->conjs[pt2[2]];
  rop = 0;
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+3,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+2,tables);
  s = fq_idx_mul(&s,pts+4,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+10,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+0,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+7,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(-1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+8,tables);
  s = fq_idx_mul(&s,pts+9,tables);
  rop = fq_idx_add(&s,&rop,tables);
  s = 0;
  t = fq_idx_mul_ui(&c,(2*p +(1))%p,tables);
  s = fq_idx_add(&s,&t,tables);
  s = fq_idx_mul(&s,pts+1,tables);
  s = fq_idx_mul(&s,pts+5,tables);
  s = fq_idx_mul(&s,pts+6,tables);
  s = fq_idx_mul(&s,pts+11,tables);
  rop = fq_idx_add(&s,&rop,tables);
  return rop;
}

/*----------------------------------------------------------------------------*/

