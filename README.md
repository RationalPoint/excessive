### On the Maximum Gonality of a Curve over a Finite Field

This README gives an overview of the software used to perform the computations
in

[1] Xander Faber, Jon Grantham, and Everett Howe. "On the Maximum Gonality of
    a Curve Over a Finite Field." Preprint, 2022. arXiv: ??

CODE AUTHORS: Xander Faber (awfaber@super.org) and Everett Howe(however@alumni.caltech.edu)

LICENSE: [GPL-3.0 License](LICENSE)

#### Brief software description

* **Lauter's Algorithm**: Sage code for finding real Weil polynomials for putative
  excessive curves. This is an implementation of an algorithm of Lauter, based
  on Serre's explicit method.

  * `explicit.py`: A Sage module containing the primary functions for finding
    real Weil polynomials with specified numerical constraints.

  * `explicit_parallel.sage`: A Sage script giving an alternate implementation
     of the function `explicit.find_real_weil_pols` that parallelizes the final
     loop. While it can be used for any genus, we only needed it to break up the
     computation for curves of genus 9 and 10.

* **Quartic Curve Search**: Code for enumerating isomorphism classes of smooth plane
  quartic curves over finite fields with no rational point.

  * `PointlessQuartics.magma`: A Magma script for finding all excessive curves of
    genus 3 over GF(q) for q < 17.

  * `AllGenus3.magma`: A Magma script containing information about the 215
  excessive curves of genus 3 over finite fields, in a form that we hope is both
  human- and computer-readable. 

  * `PointlessDoubleCovers.magma`: A Magma script for computing pointless double
  covers of an elliptic curve over a finite field. It is used in the paper to
  enumerate all excessive curves of genus 3 over GF(29). (There is only one, up
  to isomorphism.)

  * `quartic_search.c`: C code for enumerating excessive curves of genus 3 over
    GF(q) for 5 <= q <= 23. It will also find non-smooth curves and multiple
    representatives of the same isomorphism class. This is housed in the `quartic`
    subdirectory.

  * `quartic_dedupe.py`: A Sage script for determining the smooth curves in the
    output of `quartic_search.c`, and for removing duplicate representatives of
    the same isomorphism class.

* **Plane Curve Search**: Code for finding isomorphism classes of excessive curves of
  genus 6 and 7 with certain extra conditions, as detailed in sections 5 and 6
  of [1].

  * `plane_search`: C code for finding all (possibly singular) plane curves of
     degree g+1 over GF(2) with certain numeric constraints, all of which must
     be satisfied by an excessive curve of genus g possessing a divisor of
     degree g-3.

  * `xp_general.magma`: A Magma script that determines which (if any) of the
     curves found by `plane_search` are excessive.

  * `genus7_plane_search`: C code for finding all (possibly singular) degree-9
     plane curves of genus 7 over GF(2) with at least three distinct cubic
     points.

  * `xp_genus7_special.magma`: A Magma script that determines which (if any) of
     the curves found by `genus7_plane_search` are excessive.

  * `ps_launcher.py`: A Python script for launching multiple instances of
     `plane_search` and `genus7_plane_search` and tracking their
     completion. This is used to split the search space up.
  
Some of the Sage scripts depend on `search_tools.py` and `progress.py`, which
can be found in the [Gonality](https://github.com/RationalPoint/gonality)
repo. In addition, the `sage_launcher.py` script is valuable for tiling a
search space across multiple instances of Sage; it can also be found in the
[Gonality](https://github.com/RationalPoint/gonality) repo.

Additional detail on how to run this code is given below. Code blocks that begin
with `$` should be run on the command line; those that begin with `sage:` should
be run in an interactive Sage session; those that begin with `>` should be run
in an interactive Magma session. For documentation, see the headers of the
appropriate `.c`, `.magma`, `.py`, or `.sage` files.

#### Enumerating Real Weil Polynomials

To find all real Weil polynomials associated with excessive curves as in Section
3 of [1], do the following.

1. In an interactive Sage session, compute the real Weil polynomials for curves
   of genus g <= 8 and q as in Proposition 2.2. For example, to compute the real
   Weil polynomials for genus 7 over GF(2), use the command
   
   ```
   sage: import explicit
   sage: explicit.find_real_weil_pols(7,2,0,excessive=True)
   ```
   Grabbing the `progress` module from the
   [Gonality](https://github.com/RationalPoint/gonality) repository will make it
   a little more pleasant to wait for this computation to complete.

2. The function in the previous step will work for genus 9 and 10, but it will
   take far too long in a single Sage instance. Instead, we recommend running 32
   instances of the script `explicit_parallel.sage` using `sage_launcher.py` (on
   a machine with at least 32 CPUs). For example, this can be run for genus 9
   over GF(2) with the following command:
   ```
   $ sage_launcher.py explicit_parallel.sage data92 32 -e '9 2' -c
   ```
   The log files and the output will appear in the subdirectory `data92`. In
   particular, the unique real Weil polynomial it finds will be in the file `out`.
   
#### Genus 3

Theorem 4.1 of [1] asserts that there exists an excessive curve of genus 3 over
GF(q) if and only if q <= 23 or q=29 or q = 32. Any such curve can be written as a
pointless plane quartic curve.

* For q <= 5, a brute force search over plane quartics can be executed in order to
  find all of them up to isomorphism. For example, to compute the pointless
  quartic curves over GF(2), do the following in an interactive Magma session:
  ```
  > load 'PointlessQuartics.magma';
  > brute_force(2);
  ```

* For 7 <= q <= 23, we use a more efficient search as described in section 8 of
  [1]. For q <= 16, this can be performed in Magma in a reasonable amount
  of time. For example, if q = 13, the commands are:
  ```
  > load 'PointlessQuartics.magma';
  > pointless_quartics(13);
  ```

  For q > 16, we instead use C code that performs the same algorithm. (In fact,
  the C implementation works for all q >= 5.) First, modify the Makefile template
  as appropriate in the directory called `quartic` and run `make`. Running `make
  test` will do the case q=3, though the results will be incomplete.  (There is a
  theoretical obstruction to finding all excessive curves over GF(3) with this
  algorithm.)  To run the search for q = 17, do the following:
  ```
  $ quartic_search 17 q17.data 1 0
  ```
  This will take between 10 and 11 hours to complete on a modestly new CPU. To
  parallelize the search across 32 CPUs and put the data in the file `data17/out`,
  use the `ps_launcher.py` utility as follows:
  ```
  $ ps_launcher.py quartic_search data17 32 17 -c
  ```
  (Running ps_launcher.py with no arguments gives the syntax.)
  The output of the `quartic_search` executable will contain singular curves, as
  well as duplicate representatives of the same isomorphism class. We have a Sage
  script that selects representatives for the smooth curves in the output:
  ```
  $ sage quartic_dedupe.py 17 data17/out
  ```
  The results will be written to the file `data17/out.excessive`.

* In section 8 of [1], it is shown that any excessive curve of genus 3 over GF(29)
  must be a double cover of the elliptic curve y^2 = x^3 + x. One can enumerate
  all such covers up to isomorphism (there is only one) using the Magma commands
  ```
  > load 'PointlessDoubleCovers.magma';
  > pointless_quartics_over_E(EllipticCurve([GF(29)| 0,0,0,1,0]));

* Finally, the case q = 32 is handled by purely theoretical means in [1].

To see a summary of the 215 isomorphism classes of excessive curves of genus 3,
as well as to obtain their real Weil polynomials and automorphism groups, look
in `AllGenus3.magma`.

#### Genus 4 and 5

There is nothing new to do here. See Section 4 of [1].

#### Genus 6 and 7

We begin by searching for curves of genus g with a divisor of degree g-3, as in
sections 5 and 6.1 of [1].  TO that end, we use the `plane_search` app:

1. Modify the Makefile template as appropriate in the directory called
   `plane` and run `make`.

2. For g < 7, this app can be used on a single CPU. For example, for genus 6, it
   is called as follows:
   ```
   $ plane_search 6 genus6.data 1 0
   ```
   For genus 7, the computation will take around 2.5 days on a single CPU.
   On a machine with multiple CPUs, we recommend splitting it into around 16 jobs.
   For example, to break it into 16 jobs and run the third tile of the search space,
   use the command
   ```
   $ plane_search 7 genus7.data.3 16 3
   ```
   To automate the running of all jobs, use `ps_launcher.py`. Type
   `ps_launcher.py -h` for syntax.

3. Run the Magma script `xp_general.magma` to determine which curves in the
   output of the previous step are smooth of the correct genus and gonality. For
   example, to run it on the file generated from the third tile of the search
   space, we would use the command
   ```
   $ magma g:=7 infile:=genus7.data.3 outfile:=genus7.data.3.magma xp_general.magma
   ```
---

Any remaining excessive curve of genus 7 must have at least three distinct cubic
points. To search for such curves, we use the `genus7_plane_search` app.

1. Modify the Makefile template as appropriate in the directory called
   `plane` and run `make`. (This may have already been done above.)

2. Run `genus7_plane_search` three times: once for each case, among 1, 2, 3. Cases
   2 and 3 should be broken into several computation and run in parallel.  For
   example, to divide case 2 into 64 jobs and examine tile 13 of the search space,
   use the command
   ```
   $ genus7_plane_search 2 case2.data.13 64 13
   ```
   To automate the running of all jobs, use `ps_launcher.py`. Type
   `ps_launcher.py -h` for syntax.

3. Run the Magma script `xp_genus7_special.magma` to determine which curves in
   the output of the previous step are smooth of the correct genus and
   gonality. For example, to run it on the file generated from the thirteenth
   tile of the search space, we would use the command
   ```
   $ magma infile:=case2.data.13 outfile:=case2.data.13.magma xp_genus7_special.magma
   ```
