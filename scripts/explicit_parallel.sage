r""".

Sage script for constructing a complete list of real Weil polynomials of
Jacobians for hypothetical excessive curves. It is a version of
explicit.find_real_weil_pols that parallelizes the work done in the final step
of the computation (corresponding to looping over a_g).

This can be used in tandem with sage_launcher.py.

"""

# Best bounds for a_g. These were located by repeated hand-use of
# explicit.a_bounds. Keys are (g,q),
# values are (a_g bound, u-values that produced that bound)

best_bds = {
  (4,5): (213, [0.006360184576491568, 9.951900809297541, 0.0014836763372916284]),
  (4,7): (724, [0.07807611789308222, 8.617467309593438, 0]), 
  (5,4): (512,[0.02820200409635598, 9.939652506364535, 6.068877586817745]),
  (6,2): (28, [0.07195814264202349, 0.1485826850743588, 7.8632809023238455, 0.04107214679586302]),
  (6,3): (194, [0.057546657809895985, 0.3802090759313004, 8.643303272541122, 0.09671321698647373]),
  (7,2): (65, [0.06842676668710079, 0.08034697808597668, 7.473955470453576, 6.431840827679441]),
  (8,2): (69, [0.2316352795058063, 0.07377849871709419, 0.31852364129229316, 9.141930504741483, 0.041257465957017514]),
  (9,2): (171, [0.05521085867995157, 0.04821115255780217, 0.10215406597031351, 9.70274669074866, 9.157094031020181]),  
  (10,2): (185, [0.05739818168041255, 0.19920229761902664, 0.12463774343493772, 0.227410982013998, 8.640187084738141, 0.1384629059431508])
}
  

################################################################################

import argparse
import explicit
import progress
import search_tools

parser = argparse.ArgumentParser()
parser.add_argument('num_jobs',type=int,help='break up work into this many jobs')
parser.add_argument('job',type=int,help='which job is the current one?')
parser.add_argument('filename',type=str,help='where to write output')
parser.add_argument('genus',type=int,help='genus of curves')
parser.add_argument('field_cardinality',type=int,help='prime power')
args = parser.parse_args()

num_jobs = args.num_jobs
job = args.job
filename = args.filename
g = args.genus
q = args.field_cardinality

################################################################################

fp = open(filename,'w')

if g <= 2:
  print('There is no excessive curve of genus g = {}'.format(g))
  fp.close()
  exit(int(0))

# Get pointless conditions for excessive curve  
a_vals = {}
conds = explicit.excessive_conditions(g)
for d in range(1,g-1):
  if (d,) not in conds: continue
  a_vals[d] = 0
  conds.remove((d,))

# Set up checking for excessive curve pointless properties
# A condition (i1, ..., ir) only needs to be checked for a_vector
# if the max_index equals ir (since the condition is not enforceable if
# ir > max_index and it is already forced to hold if ir < max_index).
if len(conds) == 0:
  seems_excessive = lambda a_vector,max_index : True
else:
  def seems_excessive(a_vector, max_index):
    for cond in conds:
      if max_index != max(cond): continue
      # a_vector has 0-based indexing
      if all(a_vector[d-1] != 0 for d in cond): return False
    return True

# Print the requirements imposed on the a_d's
msg = ''
num_vals = len(a_vals)
for i,(key,val) in enumerate(a_vals.items(),1):
  msg += 'a_{} = {}'.format(key,val)
  if i < num_vals - 1:
    msg += ', '
  elif i == num_vals - 1:
    if num_vals == 1:
      msg += '.'
    elif num_vals == 2:
      msg += ' and '
    else:
      msg += ', and '
  else:
    msg += '.'      
print('\nComputing real Weil polynomials for g = {}, q = {},\n  and a-invariants {}'.format(g,q,msg))
if len(conds) > 0:
  msg = 'Also require the following conditions:\n'
  for cond in conds:
    msg += '  at least one of ('
    for i, y in enumerate(cond):
      msg += 'a_{}'.format(y)
      if i < len(cond) - 1:
        msg += ', '
      else:
        msg += ') vanishes\n'
  print(msg)
print('')

# Require that we have a fixed bound for a_g for this choice of (g,q).
# We don't want different bounds on different Sage instances. 
initial_bounds,_ = explicit.a_bounds(g,q,0,num_tries=2**15)
if (g,q) not in best_bds:
  raise RuntimeError('No bound for a_g set in the dict best_bds')
initial_bounds[g],_ = best_bds[(g,q)]


# We know next a_d lies in [lower,upper], but we refine upper bound
# further using previous a_i's
def better_upper_bound(lower, upper, g, q, aa):
  if lower == upper:
    return upper
  new_upper = explicit.random_ad_bound(g,q,aa,num_tries=64)
  return min(upper,new_upper)

# Rule out various choices of a-values using derivative information
F = explicit.real_weil_pol(g,q)  # Indeterminate coefficients
real_pols = PolynomialRing(QQ,names=('T',))
good_as = [()]
msg = 'Determining information about a_{}:\n'
msg += '  Using range [{}, {})\n'
msg += '  At most {} polynomials to check'
for j in range(g-1,0,-1):
  Fj = F.derivative(j)
  Fjcoefs = Fj.list()
  next_good_as = []
  index = g - j # Determining a_{g-j}
  # Set upper and lower bounds
  if index in a_vals:
    lower = a_vals[index]
    upper = a_vals[index]
  else:
    lower = 0
    upper = initial_bounds[index]
  work = len(good_as)*(upper-lower+1)
  print(msg.format(index,lower,upper+1,work))
  prog = progress.Progress(len(good_as))
  for aa in good_as:
    # Try to improve the bound for the next a_d
    new_upper = better_upper_bound(lower,upper,g,q,aa)
    next_coef_range = (lower,new_upper+1)
    for b in range(*next_coef_range): 
      subs = aa + (b,) + (0,)*(j) # pad to get length g
      if not seems_excessive(subs,g-j): continue
      # Specialize F_j using these choices of a_d      
      Fj_spec_coefs = [search_tools.evaluate_poly(c,subs) for c in Fjcoefs]
      Fj_spec = real_pols(Fj_spec_coefs)
      if explicit.real_roots_check(Fj_spec,q):
        next_good_as.append(aa + (b,))
    prog()
  good_as = next_good_as    
  prog.finalize()


# Finally, check the full polynomial F
Fcoefs = F.list()
start,stop = search_tools.start_and_stop_work(1+initial_bounds[g],num_jobs,job)
lower = start
upper = stop-1 # start/stop uses pythonic indexing
work = len(good_as)*(stop-start)
print(msg.format(g,start,stop,work))    
prog = progress.Progress(len(good_as))
for aa in good_as:
  # Try to improve the bound for the next a_d
  new_upper = better_upper_bound(lower,upper,g,q,aa)
  next_coef_range = (lower,new_upper+1)
  for b in range(*next_coef_range):
    subs = aa + (b,)
    Fspec_coefs = [search_tools.evaluate_poly(c,subs) for c in Fcoefs]
    Fspec = real_pols(Fspec_coefs)
    if explicit.real_roots_check(Fspec,q):
      fp.write('{}\n'.format(Fspec))
  prog()
prog.finalize()
