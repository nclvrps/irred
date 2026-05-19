/* irred.c	 					irred version 3.15
 
==========================================================================
|                                                                        |
|  Copyright (C) 2003 R. P. Brent.                                       |
|                                                                        | 
|  This program is free software; you can redistribute it and/or         |
|  modify it under the terms of the GNU General Public License,          |
|  version 2, June 1991, as published by the Free Software Foundation.   |
|  For details see http://www.gnu.org/copyleft/gpl.html .                |    
|                                                                        |
|  External assembler routines (see ASM and comments below) are covered  |
|  by the same copyright conditions as this C code.                      |
|                                                                        |
==========================================================================
 
Author: Richard P. Brent (rpb@comlab.ox.ac.uk)
      	with assistance and contributions from
   	Samuli Larvala and Paul Zimmermann

Dates:  20000528..20030328

Summary:

   An efficient test for irreducibility of trinomials
   x^r + x^s +1 (mod 2) where r is a prime, 
   or primitivity in the case that 2^r-1 is a Mersenne prime,
   using (32-bit or 64-bit) word-length bit operations.
   Written in C (with optional assembler routines for some machines).

Restrictions, Assumptions and Limitations:

   Assumes 0 < s < r, s >= LIM, r-s >= LIM,
   where LIM is usually twice the wordlength in bits,
   or six times the wordlength for ULTRA option,
   or ten times the wordlength for UNROLL option,
   or 256 for IBMPC ASM option
   (other cases can be handled by less efficient programs,
   e.g. irred version 1.35 or earlier versions,
     or irred version 1.50 with ODDS false).

   It is assumed that r is prime - if not true then a warning is printed.
   There is no check on whether or not 2^r-1 is prime.
   
   When output says A(x) is "irreducible" this actually means that
   x^(2^r) = x mod A(x). If r is not prime, we should also check
   that GCD(x^(2^d) - x, A(x)) = 1 for divisors d of r, but the program
   does not do this (it would be easy enough to add code to perform this
   check if desired).

   Similarly, if r is prime but 2^r-1 is not prime, then an
   irreducible A(x) might not be primitive, because x^((2^r-1)/p) mod A(x) 
   is not checked for any prime divisors p of 2^r-1. (In this case we
   can not complete the primitivity test without factoring 2^r-1.)

Usage:
   	  irred restart-file log-file time-limit skip-file
   	  					(all arguments optional)
Arguments:

          A restart file can be specified as first command line argument.
          If no file is specified, input data is read from stdin.

            Format of restart file or stdin is one line
              r s1 s2
            meaning to start next run at this (r,s1) and go to (r,s2-1)
            [or, if s1 > s2, go from (r,s1) to (r,s2+1)].
            If s1 = s2 then run is finished.

	  A log file can be specified as second command line argument.
	  If no file is specified, then log output -> stdout with an
	  identifier "LOG" (can be extracted using grep).

 	  Max minutes to run (m) can be specified as 3rd command line argument.
	
	    Maximum runtime is m minutes (infinite if m = 0).
	
	    The parameters m and s2 may be omitted (in which case they
	    assumed to be zero).
	
	  A skip file can be specified as the fourth command line argument.
	  This file is a list of pairs (u, v) (one per line) indicating that
	  s in the interval [u, v] should not be tested. This is useful if
	  certain intervals are already done or reserved for another machine.
	  If a skip file is desired, the parameter m must be specified
	  (possibly zero).
	  
	  A maximum number of interactive users can be specified as the
	  fifth command line argument to override the default MAXUSERS.

	  The sixth and seventh command-line arguments can be specified
	  to override the defaults LOADTOLR and LOADTOLS.

	  For a "parallel" version using MPI, see irredpar.c

Comments on timing and memory requirements:

	Sieving quickly discards about 93% of polynomials.
	Each polynomial not discarded by sieving takes time O(r^2)
	but there may be a slow-down as r increases (this is a cache effect).

	For a 300 MHz P-II with 512KB L2 cache, times range from
	1.07r^2 nsec (for r = 44497) to 1.72r^2 nsec (for r = 3021377).

	For a 500 MHz P-III with 512KB L2 cache, times range from
	0.67r^2 nsec (for r = 44497) to 0.80r^2 nsec (for r = 3021377).

	If we measure time in units of cycles, then the time on IBM PCs
	is about 0.33r^2 if the problem is small enough for the L1 cache,
	degrading to about 0.6r^2 if the L2 cache is half as large as
	necessary (e.g. 512KB for r = 6972593).

	Detailed timing results for different versions and (r,s)
	are in the separate file times.dat.

	Space requirement is 7r/16 + (constant) bytes
	(temporarily increasing to 3r.FASTTRY/16 if FASTTRY > 2, see below).
	
	The working set after sieving is 3r/16 + (constant) bytes.

References:

        R. P. Brent, Search for primitive trinomials (mod 2),
        http://www.comlab.ox.ac.uk/oucl/work/richard.brent/trinom.html

        R. P. Brent, S. Larvala and P. Zimmermann, A fast algorithm for
        testing irreducibility of trinomials mod 2 (preliminary report),
        Report PRG-TR-13-00, Oxford University Computing Laboratory,
        30 Dec 2000. Revision to appear in Mathematics of Computation
        (posted electronically 18 Dec 2002).     
        See http://www.comlab.ox.ac.uk/oucl/work/richard.brent/pub/pub199.html  

        S. W. Golomb, Shift register sequences, Holden-Day, San Francisco,
        1967. Revised edition, Aegean Park Press, 1982.

   	J. W. Heringa, H. W. J. Bl\"ote and A. Compagner, New primitive 
        trinomials of Mersenne-exponent degrees for random-number generation,
        Int. Journal of Modern Physics C 3 (1992), 561-564.

	D. E. Knuth, The Art of Computer Programming, Vol. 2, third ed., 
	Sec. 3.2.2.

	T. Kumada, H. Leeb, Y. Kurita and M. Matsumoto, New primitive
        t-nomials (t = 3, 5) over GF(2) whose degree is a Mersenne exponent,
	Mathematics of Computation 69 (2000), 811-814.
	Corrigenda: ibid 71 (2002), 1337-1338.

	N. Zierler, Primitive trinomials whose degree is a Mersenne exponent,
        Information and Control 15 (1969), 67-69.

*/

#include <stdio.h>
#include <stdlib.h>		/* For malloc */
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>     /* For sleep */
#include <time.h>   		/* For clock */

#define true 1
#define false 0

/* VERBOSE, CONTINUE, GNU determine program behaviour */

#define VERBOSE true		/* If true give more informative output */

#define CONTINUE false		/* If true, continue after finding a
                                   primitive trinomial (useful for small r) */

#define GNU true		/* Determines if copyright notice is printed
				   (but the program is copyright anyway) */

/* ===>>> Set exactly one of the following machine types to true <<<=== */

  #define IBMPC true /* 32-bit CISC, little-endian, e.g. Pentium.

		      With gcc 2.95.2, compile irred.c with
		      gcc -O1 -fomit-frame-pointer -funroll-loops
		      	-o irred irred.c [irred.o]

		      With other versions of gcc, e.g. 2.7.2.3, try
		      gcc -O1 -fomit-frame-pointer
		        -o irred irred.c [irred.o]

		      where the assembler routines (if any) have been
		      assembled to give relocatables irred.o
	 	      */

  #define SPARC false  /* 32-bit RISC, e.g. Old (32-bit) Sun Sparc. 
  			 Compile with gcc -O6 or (better) cc -Ofast */

  #define ALPHA false /* 64-bit RISC, e.g. DEC Alpha. Compile with
			   gcc -O4 -mcpu=ev6 -Wa,-arch -Wa,ev6 -funroll-loops 
			   or cc -fast (on Compaq alphaserver).
			   Also recommended on recent Sparc processors 
			   with 64-bit options: 
			   cc -Ofast -xtarget=ultra -xarch=v9 */

  #define SGI false   /* 64-bit SGI R10000. Compile with cc -64 -Ofast */


#define ULTRA (SPARC || SGI)	/* Try on Sparc Ultra and MIPS R12000 */  	
#define UNROLL true		/* Best to set false on Sun-Blade 1000,
				   but true on Ultra-80 */

/* Following are only relevant on IBM PCs */

  #define LINUX true		/* Set true if running under Linux,
				   false otherwise. Relevant to choice
				   of FASTTRY (see below). */

/*

Comments on different versions:

   Version 3.15 has
   
   1) Extra 4 words allocated before arrays p and q in main
      and a0 in main/fastmem to avoid (probably harmless) valgrind 
      diagnostics - thanks to Julian Seward and Andrew Tridgell for 
      pointing out the problem.

   Version 3.14 has
   
   1) UNROLL flag to permit 4-way unrolled version of reducer
      (good on Sparc Ultra-80).

   Version 3.13 has
   
   1) ULTRA flag to permit another version of reducer
      (good on Sparc Ultra-80 and MIPS R12000).

   Version 3.12 has
   
   1) Argument of prime(n) is type int instead of uint64_t
      to avoid diagnostic on 64-bit (e.g. SGI) machines.

   Version 3.11 has
   
   1) Better documentation, no significant changes to the code.

   Version 3.10 has

   1) Explicit reference to the GNU General Public License (previously
      this information was communicated separately). 

   2) Test for primality of r (a warning is printed if r is composite).
      The irreducibility test assumes that r is prime, and is incomplete
      if r is composite (see comments above). Earlier versions assumed that
      r was prime without any check.

   Version 3.07 has
   
   1) More accurate output (only claims irreducible rather than primitive).
      The earlier versions assume r is a Mersenne exponent (so 2^r-1 prime),
      this version only assumes that r is prime and tests if the trinomial
      is irreducible. In general, if 2^r-1 is composite, further tests 
      depending on the factors of 2^r-1 are necessary to check primitivity.
      If 2^r-1 is prime, no further checks are necessary.
      
   Version 3.06 has
   
   1) Skip list handling fixed so can skip a single item 
   (previously low and high bounds equal terminated list processing).

   Versions 3.04 and 3.05 have

   1) New CPU info for booth (1-35, 51-55) and tosca.

   Version 3.03 has

   1) Option to modify load tolerances using command-line arguments.

   Version 3.02 has

   1) Improved version of sysload, no longer limited to range [0, 2.55].

   2) SLEEPTIME becomes two different constants SLEEPU, SLEEPL as the load
      is likely to change more often than the number of users.

   3) Option -2 for MAXUSERS (checks both userkt and sysload).

   Version 3.01 has

   1) More efficient versions of userkt and sysload, avoiding use
      of a temporary file.

   Version 3.00 has

   1) uptime instead of who|wc to count interactive users (this can
      give different results if one user has several sessions).

   2) Counters for clock and system calls.

   3) Better error handling for system calls (e.g. will not bomb if
      the file /tmp/irred.tmp is not writeable).

   4) Better error handling for command-line arguments, e.g.  skip file
      (so can use 0 or /dev/null if don't desire skip file (4th argument)
      but need 5th argument). All sscanf calls are now checked to see
      if the expected number of fields was found.

   5) Name changed: irredg -> irred (as it used to be).

   Version 2.93 has

   1) Flag BOOTH as shortcut for options on boothxx machines.

   2) option to sleep while load average is high (set MAXUSERS = -1).
      
   Versions 2.92 and 2.91 are experimental variants on the "sleeping" theme.   

   Version 2.90 has 
   
   1) option to go to sleep while the number of interactive users exceeds
      a tolerance (set MAXUSERS to a nonnegative value).

   Version 2.82 has
   
   1) Changed format of log file - hex numbers now have
      leading zeros instead of blanks. For example, 
      "x   12345" -> "x00012345".
      
      Old log files can be checked and converted to the
      format used by version 2.82 with the separate program "fixlog".

   Version 2.81 has
   
   1) Changed format of log file - if sodd = r-s > r/2 use
      r s yabcdefgh instead of r sodd xabcdefgh

   Version 2.80 has
   
   1) Dynamic choice of sieving cutoff, only use NEXTRA for upper bound.
      This should improve overall efficiency because the cutoff can 
      adapt to different machines/compilers and also to exceptional
      (r, s) values, e.g. if s is small or close to r/2.
   
   2) Printing of skip list commented out (as may be too long).

   Version 2.76 has

   1) Constant NEXTRA reduced from 4 to 3 in some cases, including
      (ASM AND IBMPC). This should save about 3% on average for time 
      per trinomial on IBM P-III such as booth, at the expense of a 
      small inconsistency in the log files produced on different machines.

   2) Introduction of LINUX flag to distinguish between Linux and Solaris
      (only use FASTTRY > 1 under Linux).

   Version 2.75 has

   1) Optional empirical choice of a "good" memory location for data arrays
      by timing of some dummy iterations. Perhaps surprisingly, this helps
      on Portland (P-II) when r is large. May depend on the operating system
      (seems less effective under Solaris than Linux).
      See FASTTRY and fastmem below.

   Version 2.72 has

   1) Memory allocation on stack instead of heap (this turned out to be
      no faster, so reverted to heap in later versions).

   Version 2.71 has
   
   1) Option CLEAR to clear space allocated by malloc before use
      (may also have effect of clearing L1 and L2 caches).

   2) Rounding up of size of data arrays to multiple of cache line
      size (32 bytes) for IBM PC (this did not help, so omitted in
      later versions).

   Version 2.70 has
   
   1) New MMX routine reducemx, similar to reducemr but making
      different alignment assumptions.
      
   2) Interface routine reduceml to call reducemr/reducemx as
      appropriate, ensuring the correct data alignment.
      
   3) prefetch instructions removed in reducer to speed small cases
      on P-II machines.  [Restored but commented out in version 2.76]

   Version 2.65 has

   1) MMX routine reducemr used instead of reducer in some case.

   Version 2.61 has
   
   1) Call to reducer if shift zero but (IBMPC AND ASM)
      (saves about 5% in the special case of shift zero,
       which should occur about 1 in 32 cases).

   Version 2.60 has 
   
   1) interlvr restored as versions without it were 
      sometimes slower than versions with it (at least this is
      true for r = 3021377 on booth, where 3r/16 bytes is very 
      close to the cache size - see comments re version 2.40 below).

   2) Space requirement decreased to 7r/16 bytes and
      working set after sieving decreased to 3r/16 bytes.

   Versions 2.51-2.54 have minor variations in interlvf MMX code.

   Version 2.50 has

   1) interlvr abolished, instead we use a cycle of length three
      involving interlvr to get better cache performance in reducer.

   2) Space requirement increased to r/2 bytes (formerly 7r/16 bytes),
      working set after sieving increased to r/4 bytes (formerly 3r/16 bytes).

      Version 2.50 seems about 10% slower than version 2.40 on booth, 
      for r = 3021377.  Slightly faster for r = 6972593 so may revert 
      to version 2.50 in the future.

   Version 2.40 has

   1) Call to MMX routine reducemr removed as non-MMX reducer is faster.
      [This may depend on data alignment. The problem is that we can't ensure 
       that data is aligned on 8-byte boundaries in calls to reducemr.]

   Version 2.34 has

   1) Improved reducer loop in irred.s (non-MMX code).

   2) Call to reducer instead of reducemr.

   Version 2.33 has

   1) Unrolled reducer restored from version 2.00 because faster
      on SPARC (no change in MMX version).

   Version 2.32 has

   1) Reordered loops in MMX code (version 2.31 is slightly different)

   Version 2.30 has

   1) Combined align and interleave, alternating forward and back.
      This gives improved performance on IBM PC, not much difference
      on other machines.

   2) Space requirement increased to 7r/16 bytes (formerly 3r/8 bytes)
      but working set after sieving is the same, 3r/16 bytes.

   3) Options FORWARD, MMX, TABLE, UNROLL eliminated as now redundant.

   Version 2.20 has
   
   1) FORWARD parameter for optional forward loop in alignmr[f] (MMX version)
      and aligner (C non-unrolled version).
      This improves cache performance in some cases.

   Version 2.12 has
   
   1) Faster (at least on Portland PC) interlvm with shorter inner loop
      (may be slower on booth)
	
   Version 2.11 has
   
   1) Replacement of some shifts and ors by adds in interlvm

   Version 2.10 has
   
   1) Option to use ASM without MMX, and new non-MMX IBM assembler
      routines interlv2, interlv3

   Version 2.01 has

   1) Improved comments.

   2) Improved loop control in reducemr, alignmr (now OK to call with
      kt negative).

   Version 2.00 has
   
   1) Improved algorithm with reduced memory references (C only so far).
      This algorithm is faster overall than the one used in earlier versions
      although sieving is slightly slower.
      
   2) Restriction that min(s, r-s) can not be too small (previously
      only r-s was restricted) and that r is odd. Use version 1.35 or 
      earlier if necessary to circumvent these restrictions.
      
   3) Because even s is replaced by r-s, log files produced by version 2.00
      are slightly different from those produced by version 1.35 and earlier.
      
   4) PRELOAD and TABLEZ flags abolished.   
   
   Version 1.50
     
   1) Optionally uses odd s (useful for checking version 2.00 and later),
      otherwise the same as version 1.35.

   Version 1.35 has
      
   1) Better detection of illegal parameters (if ASM AND IBMPC)
         
   Version 1.34 has
   
   1) More informative output re CPU type (if IBM PC).

   Version 1.33 has
   
   1) Faster code in case that r-s is divisible by 32 and IBM PC assembler
      is available (previously this case did not call assembler).

   2) ASM flag reinstated to replace MMX and MMX2 flags.

   Version 1.32 has
   
   1) More informative output - times for last sieve iteration
      and fraction of time spent sieving.

   2) Bug in versions 1.30 & 1.31 (which gave wrong "n" for log on stdout
      if printed in phase 2) fixed.
      
   3) Storage allocation bug (which caused segmentation fault if only one
      command-line argument on some machines) fixed.
        
   Version 1.31 has
   
   1) Unrolled option in reducer (may help on Sun etc).

   Version 1.30 has 
   
   1) Faster relprime calling reducer for inner loop
      (with assembler version for IBM PC).
   2) Introduction of NEXTRA and abolition of RATIO.
   3) Storage requirement reduced to 3r/8 + O(1) bytes.

   Version 1.25 has
   
   1) More CPU information in case of IBM PC.
   2) Two assembler versions - 
      one (irrednas.s) uses prefetching (works on PIII and above)
      the other (irredasm.s) without prefetching.
      (Note: later irredasm.s was abolished and irrednas.s renamed irred.s)

   Version 1.24 has
   
   1) MMX2 option to allow reduce2.s in case q1 and q1-q2 both odd.

   Version 1.23 has
   
   1) Alignment of most 8-byte accesses on 8-byte boundaries.
   2) Abolition of PORTLAND, ASM, OLDGCC, OLDR options.
   3) Merging of MMXR and MMXS options to MMX.

   Version 1.22 has
   
   1) Larvala's new MMX routines (renamed squarem and reducem)

   Version 1.21 has
   
   1) Alternative version of MMX (reduce2.s) modified RPB  

   Version 1.20 has

   1) MMXR and MMXS options for use of Multi Media extension (mmx) 
      instructions on Pentium etc (written by Samuli Larvala).

   Version 1.12 has
   
   1) Slightly faster reducea (though still some room for improvement ?)
      Set OLDR to use the old version of reducea (may be faster on machines
      other than IBM PC).
      
   2) Comments and output improved (2^p is never prime for p > 1 !).   

   Version 1.11 has
       
   1) ASM flag added and square.s implemented (IBM PC version only)
   2) comments improved 

   Compared to the previous versions 1.0 and earlier (unnumbered),
   version 1.1 has

   1) a fix to avoid timer overflow for large r.
   2) memory requirements for sieving reduced by about 46%.
   3) addition of a "skip file" as optional 4th command-line argument.
   4) various small enhancements and improvements in the comments.

Compilation flags:

   Set following flags to optimise performance depending on machine type
   and compiler.
   
   On an IBM PC which can execute MMX instructions it is best 
   to set ASM = true and use the assembler routines in irred.s
   (version 3.10 or higher). 
   
*/

#if UNROLL
  #define LIM (10*WLEN)
#elif ULTRA
 #define LIM (6*WLEN)
#else			
 #define LIM (2*WLEN) 
#endif

#define CLEAR true			   /* To clear allocated space */

/* Some constants */

#define TIMEOUT 10		/* Seconds to timeout when opening files */
#define CPUTOL 0.01		/* Times less than CPUTOL sec are negligible */
#if (IBMPC && LINUX)				

  #define FASTTRY 16		/* Number of a0 to try in fastmem.
				   Set to 1 if trivial selection desired. */
#else

  #define FASTTRY 1		/* Selection of a0 only seems useful on IBM PC
  				   running under Linux (not Solaris),
  				   but in any case do one timing run to 
  				   estimate non-sieving time and thus 
  				   determine a good cutoff for sieving. */
#endif

#define SMALLR  200000		/* Trivial selection if r < SMALLR */
#define CPUTEST 1		/* Seconds for each trial in fastmem.
				   Assume that CPUTOL is much less than 
				   CPUTEST. */

#define ERROR 1
#define OK 0

/* Useful macro definitions */

#define ODD(i)		((i)&1)

#define WLEN  (64)	/* Bits in a long word, always assume 64-bit */
#define WLENM (WLEN-1)		/* Ditto less 1, i.e. 31 or 63 */
#define WD    (6)		/* 5 or 6, so WLEN = 2^WD */

	
#define	NEXTRA 0		/* Do up to NEXTRA GCD computations for
				   2^(degree) >= r. Usually 3 or 4 is optimal.
				   Since version 2.80 the cutoff has been 
				   estimated at runtime (see DYNAMIC below)
				   and NEXTRA is only used as an upper bound.
				   This can result in minor differences in 
				   the log files produced in different runs
				   or on different machines. */

#define NMAX (WLENM + NEXTRA)	/* Sieve over x^kn - 1 for kn = 2^n-1,
				   2^n < r, n <= NMAX (but see also NEXTRA) */

#define MC 256			/* Max string length */

/* Struct for the skip list (a singly-linked list of pairs) */

struct skip {
  int low, high;
  struct skip *next;
  };

/* Global variables */

int r, q1, sodd;			/* Main parameter is r */

double CPUtime;				/* For timer */
clock_t cstart;

	/* 32 bits should suffice for these counters */

unsigned int syscalls = 0;			/* Counts calls to system/sleep */
unsigned int clockcalls = 0;			/* Counts calls to clock */
unsigned int space = 0;				/* Counts bytes allocated by malloc */

/* Routines start here */

double clockd(starta, first)
clock_t *starta;
bool first;

/* If first = true, initialises *starta and returns zero.
   If first = false, returns time in sec since last call and updates starta. 
   Note that this routine should be called at least each 4290 sec on
   IBM PC (under Linux) since CLOCKS_PER_SEC = 1000000 and 32-bit integer
   overflow may occur if called infrequently */
  
  {
  clock_t stop;
  uint64_t diff;
  stop = (clock_t) clock();
  clockcalls++;
  if (first) *starta = stop;
  diff = (uint64_t)(stop - *starta);	/* May wrap around in 2^32/10^6 sec */
  *starta = stop;
  return ((double)diff*((double)1/(double)CLOCKS_PER_SEC));
  }


char *mymalloc(size)
int size;

/* Same as malloc but gives error exit if can not allocate space.
   Optionally clears space allocated. */

  {
  char *ptr;
#if CLEAR
  int j;
#endif
  space += size;
  ptr = (char *)malloc(size);
  if (ptr == NULL) {
    printf("malloc unable to allocate %d bytes\n", space);
    exit(ERROR);
    }
#if CLEAR
  for (j = 0; j < size; j++) ptr[j] = (char)0;
#endif
  return(ptr);
  }

bool comparex(a)

/* Returns true if poly a of degree r-1 is x */

uint64_t *a;
  {
  uint64_t mask1;
  int j;
  if (a[0] != 2L) return(false);
  for (j = q1-1; j > 0; j--) {
    if (a[j] != 0) return(false);
    }
  mask1 = (uint64_t)(~0L) >> (WLENM - ((r-1) & WLENM));  
  if ((a[q1] & mask1) != 0) return(false);
  return(true);
  }
  
#if UNROLL

void reducer(a, b, kt, shift, prev)

/* 4-way unrolled version of reducer, optimised for Sparc Ultra-80.
   Note that LIM must be at least 10*WLEN. 
   Called by relprime and reducep. Assumes 0 < shift < WLEN. */

uint64_t *a, *b;
int kt, shift;
uint64_t *prev;		/* Previous -> last value of new */

  {
  int j, shiftc;
  uint64_t bj, bj2, bj3, bj4, old, new, old2, new2, old4, new4;
  new = *prev;
  shiftc = WLEN - shift;
  for (j = kt; (j >= 0) && ((j & 3) != 3); j--) { /* Up to 3 iterations */
    old = new;
    new = a[j];
    b[j] ^= (new >> shift) | (old << shiftc);
    }

  /* Now j+1 is divisible by 4 */
  
  old2 = a[j];			/* Preload some registers */
  new2 = a[j-1];  
  old4 = a[j-2];
  new4 = a[j-3];
  bj   = b[j];
  bj2  = b[j-1];
  bj3  = b[j-2];
  bj4  = b[j-3];
  				/* Main loop.  Try to load ahead as far
  				   as possible. This will involve some
  				   harmless "out of bound" loads */
  for (; j >= 0; j -= 4) {
    old  = old2;
    old2 = a[j-4];
    b[j] = bj ^ ((old >> shift) | (new << shiftc));
    bj   = b[j-4];
    new  = new2;
    new2 = a[j-5];
    b[j-1] = bj2 ^ ((new >> shift) | (old << shiftc));
    bj2  = b[j-5];
    old  = old4;
    old4 = a[j-6];
    b[j-2] = bj3 ^ ((old >> shift) | (new << shiftc));
    bj3  = b[j-6];
    new  = new4;
    new4 = a[j-7];
    b[j-3] = bj4 ^ ((new >> shift) | (old << shiftc));
    bj4  = b[j-7];
    }
  *prev = new;
  return;
  }

#elif ULTRA

void reducer(a, b, kt, shift, prev)

/* Unrolled version of reducer, good on Sparc Ultra-80 and R12000.
   Called by relprime and reducep. Assumes 0 < shift < WLEN. */

uint64_t *a, *b;
int kt, shift;
uint64_t *prev;		/* Previous -> last value of new */

  {
  int j, shiftc;
  uint64_t bj, bj2, old, new, old2, new2;
  new = *prev;
  shiftc = WLEN - shift;
  for (j = kt; (j >= 0) && ((j & 1) == 0); j--) { /* One iteration if
  						      kt even */
    old = new;
    new = a[j];
    b[j] ^= (new >> shift) | (old << shiftc);
    }

  old2 = a[j];			/* Preload some registers */
  new2 = a[j-1];  
  bj   = b[j];
  bj2  = b[j-1];
  				/* Main loop.  Try to load ahead as far
  				   as possible. This will involve some
  				   harmless "out of bound" loads */
  for (; j >= 0; j -= 2) {
    old = old2;
    old2 = a[j-2];
    b[j] = bj ^ ((old >> shift) | (new << shiftc));
    bj = b[j-2];
    new = new2;
    new2 = a[j-3];
    b[j-1] = bj2 ^ ((new >> shift) | (old << shiftc));
    bj2 = b[j-3];
    }
  *prev = new;
  return;
  }

#else			

void reducer(a, b, kt, shift, prev)

/* Unrolled version.
   Simplified version of reducea, called by relprime and reducep.
   Assumes 0 < shift < WLEN. */

uint64_t *a, *b;
int kt, shift;
uint64_t *prev;		/* Previous -> last value of new */

  {
  int j, shiftc;
  uint64_t old, new;
  new = *prev;
  shiftc = WLEN - shift;
  for (j = kt; (j >= 0) && ((j & 1) == 0); j--) { /* One iteration if
  						      kt even */
    old = new;
    new = a[j];
    b[j] ^= (new >> shift) | (old << shiftc);
    }
  for (; j >= 0; j -= 2) {
    old = a[j];					/* old <-> new here ! */
    b[j] ^= (old >> shift) | (new << shiftc);
    new = a[j-1];
    b[j-1] ^= (new >> shift) | (old << shiftc);
    }
  *prev = new;
  return;
  }

#endif

void reducep(a)	

/* Reduces implicit square of polynomial a of degree < r "in place"
   given x^r = x^sodd + 1 (mod 2), result left in scrambled order
   0 2 4 ... 1 3 5 ...
   
   Assumes 0 < sodd < r && r-sodd >= 2*WLEN && ODD(r) && ODD(sodd).
    
   Uses about (2.load + 1.store)r/(2.WLEN) memory ops.
   
   RPB, 20000815. */

uint64_t *a;

  {
  uint64_t new;
  uint64_t temp;
  int j;
  int alpha, delta;		/* Could be global */
  int deltaw, deltaq, deltaqc;	/* ditto */
  int q1, q4;			/* ditto */
  uint64_t mask1, mask2;		/* ditto */
  
  alpha = r >> 1;		/* alpha = (r-1)/2 */
  delta = (r - sodd) >> 1;	/* delta = (r - sodd)/2 */
  deltaw = delta >> WD;		/* deltaw = delta div WLEN */
  deltaq = delta & WLENM;	/* deltaq = delta mod WLEN */
  deltaqc = WLEN - deltaq;	/* Special case if deltaq is zero */
  q1 = (r-1) >> WD;		/* q1 = (r-1) div WLEN */
  q4 = alpha >> WD;		/* q4 = alpha div WLEN */
  
  mask1 = (uint64_t)(~0L) >> (WLENM - ((r-1) & WLENM));
  				/* mask1 has WLEN-1 - ((r-1) mod WLEN)
  				   zero bits in high positions */
  mask2 = (~1L) << (alpha & WLENM);
  				/* mask2 has (alpha mod WLEN) + 1 zero bits
  				   in low positions */

  /* To be safe, mask any high bits of a which are irrelevant */

  a[q1] &= mask1;
  a[q1+1] = 0;			/* In case reducemx called */
  a[q1+2] = 0;			/* Ditto */
  new = 0;			/* Assumed by reducemx */
  
  if (deltaq == 0) {		/* Special case deltaqc == WLEN */
    
    for (j = q1; j > q4; j--)	/* C reducer does not work here */
      a[j-deltaw] ^= a[j];
    a[q4-deltaw] ^= a[q4] & mask2;
    }

  else {			/* Usual case, deltaqc < WLEN */

    reducer (a+q4+1, a+q4+1-deltaw, q1-q4-1, deltaq, &new);

    /* The last two iterations are special as need to mask some bits */

    temp = new;
    new = a[q4] & mask2;
    a[q4-deltaw]   ^= (new >> deltaq) | (temp << deltaqc);
    a[q4-deltaw-1] ^=  new << deltaqc;
    }
  }
  
void setupx(a)

/* Sets up a = x (polynomial of degree r-1, all mod 2) */

uint64_t *a;
  {
  int j;
  for (j = q1; j > 0; j--)
    a[j] = 0UL;
  a[0] = 2UL;			/* 0...010 represents x */
  }

void interlvf(a, b, r)

/* 64-bit version of interleave. Loop index runs up.

   If bits 0, 2, 4, ... , r-3, r-1, 1, 3, 5, ..., r-4, r-2 in a,
   moves them to b in correct order.

   Work is about (1.load + 1.store + 32.ops)r/WLEN.

   RPB, 20000907 */

uint64_t *a, *b;
int r;

  {
  int j, s1, s2, q4;
  uint64_t t, u, v, w, next1, next2, old, new;
  uint64_t c0, c1, c2, c3, c4, c5;
  int alpha;

  c0 = 0x00000000FFFFFFFFUL;		/* Some 64-bit constants */
  c1 = 0x0000FFFF0000FFFFUL;
  c2 = 0x00FF00FF00FF00FFUL;
  c3 = 0x0F0F0F0F0F0F0F0FUL;
  c4 = 0x3333333333333333UL;
  c5 = 0x2222222222222222UL;

  alpha = r >> 1;			/* alpha = (r-1)/2 */
  q4 = (alpha+1) >> WD;			/* q4 = (alpha+1) div WLEN */
  s1 = (alpha+1) & WLENM;		/* s1 = (alpha+1) mod WLEN */
  s2 = WLENM - s1;			/* In [0, WLEN) */
  
  next1 = a[0];
  old = a[q4];
  new = a[q4+1];

  for (j = 0; j <= q4; j++) {

    next2 = (old >> s1) | ((new << 1) << s2);	/* Beware case s2 == 63 */

    u = next1>>32;			/* High order 32 bits low part of a */
    w = next2>>32;			/* Ditto high part of a */
    t = next1 & c0;			/* Low order 32 bits low part of a */
    v = next2 & c0;			/* Ditto high part of a */

    next1 = a[j+1];
    old   = new;
    new   = a[j+q4+2];

    u = (u | u<<16) & c1;		 /* Operations on t,u,v,w   */
    t = (t | t<<16) & c1;		 /* can be done in parallel */
    w = (w | w<<16) & c1;
    v = (v | v<<16) & c1;

    u = (u | u<<8) & c2;
    t = (t | t<<8) & c2;
    w = (w | w<<8) & c2;
    v = (v | v<<8) & c2;

    u = (u | u<<4) & c3;
    t = (t | t<<4) & c3;
    w = (w | w<<4) & c3;
    v = (v | v<<4) & c3;

    u = (u | u<<2) & c4;
    t = (t | t<<2) & c4;
    w = (w | w<<2) & c4;
    v = (v | v<<2) & c4;

    u += u & c5;
    t += t & c5;
    w += w & c5;
    v += v & c5;

    b[2*j+1] = u | (w << 1);
    b[2*j]   = t | (v << 1);
    }
  }
  
void interlvr(a, b, r)

/* 64-bit version of interleave. Loop index runs down (compare interlvf).

   If bits 0, 2, 4, ... , r-3, r-1, 1, 3, 5, ..., r-4, r-2 in a,
   moves them to b in correct order.

   Work is about (1.load + 1.store + 32.ops)r/WLEN.

   Because of a loop optimisation, it should be valid to access a[-1].

   RPB, 20000907 */

uint64_t *a, *b;
int r;

  {
  int j, s1, s2;
  uint64_t t, u, v, w, next1, next2, old, new;
  uint64_t c0, c1, c2, c3, c4, c5;
  int q4, alpha;

  c0 = 0x00000000FFFFFFFFUL;		/* Some 64-bit constants */
  c1 = 0x0000FFFF0000FFFFUL;
  c2 = 0x00FF00FF00FF00FFUL;
  c3 = 0x0F0F0F0F0F0F0F0FUL;
  c4 = 0x3333333333333333UL;
  c5 = 0x2222222222222222UL;

  alpha = r >> 1;			/* alpha = (r-1)/2 */
  q4 = (alpha+1) >> WD;			/* q4 = (alpha+1) div WLEN */
  s1 = (alpha+1) & WLENM;		/* s1 = (alpha+1) mod WLEN */
  s2 = WLENM - s1;			/* In [0, WLEN) */
  
  next1 = a[q4];
  old = a[2*q4+1];
  new = a[2*q4];

  for (j = q4; j >= 0; j--) {

    next2 = (new >> s1) | ((old << 1) << s2);	/* Beware case s2 == 63 */

    u = next1>>32;			/* High order 32 bits low part of a */
    w = next2>>32;			/* Ditto high part of a */
    t = next1 & c0;			/* Low order 32 bits low part of a */
    v = next2 & c0;			/* Ditto high part of a */

    next1 = a[j-1];
    old   = new;
    new   = a[j+q4-1];

    u = (u | u<<16) & c1;		 /* Operations on t,u,v,w   */
    t = (t | t<<16) & c1;		 /* can be done in parallel */
    w = (w | w<<16) & c1;
    v = (v | v<<16) & c1;

    u = (u | u<<8) & c2;
    t = (t | t<<8) & c2;
    w = (w | w<<8) & c2;
    v = (v | v<<8) & c2;

    u = (u | u<<4) & c3;
    t = (t | t<<4) & c3;
    w = (w | w<<4) & c3;
    v = (v | v<<4) & c3;

    u = (u | u<<2) & c4;
    t = (t | t<<2) & c4;
    w = (w | w<<2) & c4;
    v = (v | v<<2) & c4;

    u += u & c5;
    t += t & c5;
    w += w & c5;
    v += v & c5;

    b[2*j+1] = u | (w << 1);
    b[2*j]   = t | (v << 1);
    }
  }
  
FILE *myfopen(fname, flag)

/* Attempts to open file fname. If not successful at first, keeps trying
   for TIMEOUT seconds, then error exits.
*/
   
char *fname, *flag;

  {
  FILE *fp;
  clock_t openstart;
  double time;
  time = clockd(&openstart, true);	/* Independent timer here */
  
  for (;;) {
    fp = fopen(fname, flag);
    if (fp != NULL) break;
    if ((time += clockd(&openstart, false)) > TIMEOUT) {
      printf("Could not open %s after trying for %d seconds\n", 
	fname, TIMEOUT);
      exit(ERROR);
      }
    }
  return(fp);
  }

bool skips(skiplist, s)
struct skip *skiplist;
int s;

/* Returns true if s is in the skip list */

  {
  struct skip *skiprec = skiplist;
  while (skiprec != NULL) {
    if ((skiprec->low <= s) && (s <= skiprec->high))
      return(true);
    skiprec = skiprec->next;
    }  
  return(false);
  }

uint64_t *fastmem(r, sodd, sizeah, CPUest)

int r, sodd, sizeah;
double *CPUest;

/* Selects and returns a "good" pointer a0 by performing some timing runs  -
   it is not obvious why this works (perhaps due to real <-> virtual page
   mapping, perhaps due to conflict with system tables in L2 cache), but it
   does work on Portland/bsp0/tosca for r = 3021377. Seems to have little
   effect on booth (P-III running Solaris) or on non-IBM machines. 
   
   Also (for use in computing sieve cutoff) returns estimated overall time
   in *CPUest. 
*/

  {
  clock_t cstart;
  double CPUbest, CPUworst, CPUtime;
  int nkt, kt, bestkt, nits;
  uint64_t *a0, *a1, *a;
  uint64_t *savea0[FASTTRY+1];

  nits = ((FASTTRY <= 0) || (r < SMALLR)) ? 1 : FASTTRY;

  bestkt = 1;		/* Index of best a0 yet */
  CPUworst = 0;
  CPUbest = 1;
  for (kt = nits; kt > 0; kt--) {
    CPUtime = clockd(&cstart, true);
    a0 = (uint64_t *)mymalloc(3*sizeah*(int)sizeof(uint64_t));
    savea0[kt] = a0; 	/* Save for later free or return */
    /* Changed 2 to 6 in version 3.15 to avoid out of bounds problem 
       (sizeah also increased by 4 before call to fastmem) */
    a0 += 6;		/* For optimisation in interlvr, may access a[-1] */	
    a1 = a0 + sizeah;
    a = a1;		/* Will cycle through a1, a0 */

    setupx (a);	
    for (nkt = 0; CPUtime < CPUTEST; nkt++) {
      reducep(a);		/* Timing run for about CPUTEST seconds */
      if (a == a1) {
        interlvf(a1, a0, r);		/* Forward interleave */
        a = a0;				/* Adjust pointer to data */
        }
      else {
        interlvr(a0, a1, r);		/* Reverse interleave */
        a = a1;
        }
      if ((nkt & 0x7F) == 0)		/* Reduce overhead of clockd calls */
        CPUtime += clockd(&cstart, false);
      }
    CPUtime = CPUtime/(double)nkt;	/* Normalise CPU time */
    if ((nkt == 0) || (CPUtime < CPUbest)) {
      CPUbest = CPUtime;		/* Save best result yet */
      bestkt = kt;
      }
    if (CPUtime > CPUworst)
      CPUworst = CPUtime;  		/* and (for comparison) worst too */
    }
  for (kt = nits; kt > 0; kt--) {	/* Free all but the "best" space */
    if (kt != bestkt) 	
      free(savea0[kt]);
    }
  CPUtime = ((double)r)*CPUtime;	/* Estimate of overall time */
#if VERBOSE
  if (nits > 1)
    printf("Worst/best ratio in fastmem %1.2f\n", CPUworst/CPUbest);
  printf("Estimated CPU time (not sieving) %2.2f sec = %2.2f r^2 nsec\n\n",
  	CPUtime, 1.0e9*CPUtime/r/r);
  fflush(stdout); 	
#endif
  *CPUest = CPUtime;			/* Return CPU time estimate */
  return savea0[bestkt];		/* and the "best" pointer */
  }

void hex8(n, str)
unsigned int n;
char *str;

/* Return str as 8-character hex string representing n,
   with leading zeros rather than blanks */

  {
  snprintf(str, 8+1, "%08x", n);
  }

bool prime(n)
int n;

/* Returns true if n is prime. Simple and not intended to be efficient. */

  {
  long j;
  if (n <= 1) return(false);		/* n < 2 */
  if (n <= 3) return(true);		/* n = 2 or 3 */
  if ((n & 1) == 0) return(false);	/* n even, not 2 */
  for (j = 3; j <= (n/j); j += 2) {	/* n > 3 */
    if ((n%j) == 0) return(false);
    }
  return(true);
  }

int main(argc, argv)
int argc;
char *argv[];
  {
  FILE *fp;
  double CPUtime1, CPUlast, CPUest;
  double CPUtotal = 0;
  double CPUtotal1 = 0;
  int minutes = 0;		/* minutes before stopping */
  int j, k, rv;
  int s = -1;
  int s1 = 0, s2 = 0;
  int n, sizeah, sizep;
  int skt = 0;
  unsigned int temp;
  uint64_t *a;			/* For polynomial of degree (r-1) */
  uint64_t *a0, *a1;		/* a = a0 or a1 */
  uint64_t *p, *q;			/* For sieving */
  uint64_t new;

  char line[MC];		/* Line buffer */
  char log[] = " LOG";		/* Identifier for log 
  				   (if no log file specified) */
  char str8[9];			/* 8 characters plus null */
  struct skip *skiplist;	/* The head of the skip list */
  struct skip *skiprec;
  int slow, shigh;
  bool done, found, g, swan;
  
  printf("\nThis is irred version 3.15\n");	  /* Date 20030328 */
  
#if GNU
  printf("\nCopyright (C) 2003 R. P. Brent.\n");
  printf("irred comes with absolutely no warranty. This is free software,\n");
  printf("and you are welcome to redistribute it under the conditions of\n");
  printf("the GNU General Public License, version 2, June 1991.\n");
  printf("See http://www.gnu.org/copyleft/gpl.html for further details.\n\n");
#endif

  printf("Options ");				  /* Print relevant options */
  if (IBMPC) printf("IBMPC, ");
  if (ALPHA) printf("ALPHA, ");
  if (SPARC) printf("SPARC, ");
  if (SGI)   printf("SGI, ");
  if ((IBMPC + ALPHA + SPARC + SGI) != 1) {
    printf("\nPlease set one machine type to true, ");
    printf("others to false, and recompile\n");
    return(ERROR);
    }
  if (UNROLL) 
    printf("UNROLL, ");
  else if (ULTRA) 
    printf("ULTRA, ");
  if (FASTTRY > 0) printf("FASTTRY = %d, ", FASTTRY);
  printf("NEXTRA = %d\n", NEXTRA);

  printf("\n");
      
  CPUtime = clockd(&cstart, true); 		/* Initialise timer */
  a0 = NULL;
  skiplist = NULL;
  found = false;
  r = 0;
  line[0] = '\0';

  if (argc > 4) {				/* Process skip file */	
    if ((fp = fopen(argv[4], "r")) == NULL)
      printf("Warning - can not open skip file\n"); 
    else {					/* Ignore if can't open */
      while (fgets(line, MC, fp) != NULL) {
        skiprec = (struct skip *)mymalloc((int)sizeof(struct skip));
        slow = shigh = 0;  
        if (sscanf(line, "%d %d", &slow, &shigh) < 2) break;
        if ((slow <= 0) || (slow > shigh)) break; /* >= changed to > 20010307 */
#if false				  	   /* Print skip list */
        printf("Will skip s in [%d, %d]\n", slow, shigh);
#endif      
        skiprec->low = slow;
        skiprec->high = shigh;
        skiprec->next = skiplist;
        skiplist = skiprec;
        }
      fclose(fp);	
      }
    }  

  if (argc > 3) {   /* Run-time bound if nonzero, ignored if zero */
    if (sscanf(argv[3], "%d", &minutes) <= 0) minutes = 0;
    }
    	  
  if (argc > 1) {
    fp = myfopen(argv[1], "r");
    fgets(line, MC, fp);   
    fclose(fp);	
    }
  else {
    printf ("r, s1, s2 ?\n"); 	  	  /* Interactive input */
    fgets(line, MC, stdin);
    }
  rv = sscanf(line, "%d %d %d", &r, &s1, &s2);
  if (rv < 3) s2 = 0;
  if (rv < 2) s1 = 0;
  if (rv < 1) r = 0;
  if (! prime(r)) {			  /* See comments above */
    printf("\nWarning: r = %d is not prime, so irreducibility test", r);
    printf(" is incomplete\n\n");
    }	
    					  /* Searches up from initial s1
    					     to s2-1 (or down to s2+1) */ 
  if (s1 < 0) return(OK);
  if (s1 >= r) s1 = r/2;		  /* Reasonable defaults */
  if (s2 < 0) s2 = 0;
  if (s2 > r) s2 = (s1 > r/2) ? r - 64: r/2 + 1;  

  /* Swan's theorem says trinomial is reducible if r = +- 3 (mod 8)
     and s is not 2 or r-2 */
     
  swan = ((r%8) == 3) || ((r%8) == 5);
  if (swan) {
    printf("Swan's theorem applies, only need to test s = 2\n");
    s1 = 2;
    s2 = 1;
    }   

  if (s1 < s2)
    printf("Searching from %d to %d\n", s1, s2-1);
  else if (s1 > s2)
    printf("Searching from %d down to %d\n", s1, s2+1);
  else
    printf("Empty search interval at %d\n", s1);
  if (minutes != 0) printf("Time limit %d minutes\n", minutes);
  fflush(stdout);  
  while ((s1 != s2) && (! found) && 
    ((minutes == 0) || (CPUtotal < 60*minutes))) {

    done = false;
    for (;;) {		/* Look for next s in range but not in skiplist */
      s = s1; 
      if (s1 == s2) {
        done = true;
        break;
        }
      if (s1 < s2)
        s1++;
      else
        s1--;
      if (! skips(skiplist, s)) break;		  
      }   
    if (done) break;  

    sodd = ODD(s) ? s: r-s;		/* whichever of s, r-s is odd */

#if VERBOSE       
    printf ("\nr %d, s %d, r-s %d\n", r, s, r-s);
#endif

    if ((sodd <= 0) || (r <= sodd) || 
      (! ODD(r)) || (! ODD(sodd)) || ((r-sodd) < LIM)) {
        printf ("Illegal parameters\n");
        if (! ODD(s))
          printf("For small even s try version 1.35 or earlier\n");
        return(ERROR);
        }            

    CPUtime += clockd(&cstart, false);
    CPUtotal += CPUtime;		/* CPUtotal is for all (r, s) */
    CPUtime = 0;			/* CPUtime is for current (r, s) */
    CPUtime1 = 0;			/* CPUtime1 is for current sieving */
    CPUlast = 0;
    skt++;

    /* Set up variables depending only on r */

    q1  = r >> WD;			/* q1 > 0 */

    sizeah = (q1 >> 1) + 12;		/* r/(2*WLEN) + 4 = half size of
    					   a/b array plus safety margin */
    sizep  = q1 + 4;			/* sizep = r/WLEN + 4
    					   = size of p/q array */
    	
    sizeah += sizeah & 1;		/* Round up to even values */
    sizep  += sizep  & 1;		/* to preserve 8-byte boundaries */					   		

#if IBMPC
    sizeah = (((sizeah+7)>>3)<<3);	/* Make sizeah = 0 (mod 8)
    					   because 8 4-byte words per cache
    					   line on IBM PC (added v. 2.71)
    					   (does this help ?) */
#endif

    if (a0 == NULL) {	/* Allocate space for a, b, p and q arrays */

      CPUtime += clockd(&cstart, false);
      CPUtotal += CPUtime;
      CPUtime = 0;			/* Don't count in non-sieving time */

      a0 = fastmem(r, sodd, sizeah, &CPUest); /* Find a "good" a0 on the heap 
      					   	 and estimate CPU time */
      					   	 
      CPUtime += clockd(&cstart, false);/* Count call to fastmem in total */
      CPUtotal += CPUtime;		/* but not in sieving time */
      CPUtime = 0;		

      /* The "+ 4"s were added in version 3.15 to avoid valgrind diagnostics
         caused by reading p[-1] or q[-1] */
         
      p  = (uint64_t *)mymalloc((sizep + 4)*(int)sizeof(uint64_t)) + 4;
      q  = (uint64_t *)mymalloc((sizep + 4)*(int)sizeof(uint64_t)) + 4;
      	
      a0 += 2;		/* For optimisation in interlvr, may access a[-1] */	
      			/* See comments in fastmem re version 3.15 change */
      			
      /* Following is to align a0, a1 on 8-byte boundary */

      a1 = a0 + sizeah;
      a = a1;		/* Will cycle through a1, a0 */
      }

    g = true;

    CPUtime += clockd(&cstart, false);

    if (g) { /* Phase 1 sieve no longer done */

      setupx (a);
      for (k = 0; g && (k < r); k++) {

      reducep(a);			/* Reduce (square of) a */

      /* Interleave in forward/reverse directions alternately */

      if (a == a1) {
        interlvf(a1, a0, r);		/* Forward interleave */
        a = a0;				/* Adjust pointer to data */
        }
      else {
        interlvr(a0, a1, r);		/* Reverse interleave */
        a = a1;
        }

      if ((k & 0x7FFFL) == 0) { 	 	/* Avoid timer overflow by */
        CPUtime += clockd(&cstart, false);	/* calling clockd sometimes */
        }
        
	if (! g) break;  
        }

      CPUtime += clockd(&cstart, false);

      if (g) {       
       if (comparex(a)) {
        found = (! CONTINUE);
        printf("%d %d irreducible (primitive if 2^%d-1 is prime)\n", r, s, r);
        fflush(stdout);
        if (argc > 2) {
          fp = myfopen(argv[2], "a");
	  fprintf(fp, "%d %d irreducible (may be primitive)\n", r, s);
 	  fclose(fp);
          }
        else {
	  printf("%d %d irreducible (may be primitive)%s\n", r, s, log);
	  fflush(stdout);
	  }
        }
      else {

	temp = (unsigned int)((a[0]<<32)>>32);

        hex8(temp, str8);
#if VERBOSE
        printf("Not irreducible/primitive, low word %s (hex)\n", str8); 
#endif
        if (argc > 2) {
          fp = myfopen(argv[2], "a");
          if (s == sodd)
   	    fprintf(fp, "%d %d x%s\n", r, s, str8);
 	  else
   	    fprintf(fp, "%d %d y%s\n", r, s, str8);
 	  fclose(fp);
          }
        else {
          if (s == sodd)
   	    printf("%d %d x%s%s\n", r, s, str8, log);
	  else
   	    printf("%d %d y%s%s\n", r, s, str8, log);
 	  fflush(stdout);
 	  }
        }
      }
      
#if VERBOSE        
      if (CPUtime > CPUTOL) {  
        printf("CPU time (not sieving) %.2f sec = ", CPUtime);
        printf("%2.2f r^2 nsec\n", 1.0e9*CPUtime/r/r);
        }
      if (CPUtime1 > CPUTOL) {
	printf("CPU time (sieving)     %.2f sec", CPUtime1);
	if (CPUtime > CPUtime1)
          printf(", ratio %2.0f", CPUtime/CPUtime1);
        printf("\n");   
        }
#endif
      CPUtotal1 += CPUtime1;  
      CPUtotal += CPUtime;
      CPUtime = 0;
#if VERBOSE
      if (CPUtotal > CPUTOL) {
        printf("\nOverall time %.2f sec = ", CPUtotal);
        printf("%.4f sec per trinomial\n", CPUtotal/skt);
        if ((CPUtotal1 > CPUTOL) && (CPUtotal1 < CPUtotal))
          printf("Overall %1.2f percent of time spent sieving\n", 
        	100.0*CPUtotal1/CPUtotal); 
	printf("\n");
        }
      fflush(stdout);
#endif        
      }
  
    /* Try again until irreducible/primitive trinomial found
       (or until all s checked if CONTINUE is true) */

    if (argc > 1) {
      fp = myfopen(argv[1], "w");
      fprintf(fp, "%d %d %d\n", r, s1, s2);
      fclose(fp);
      }  

    if (swan) break;
    }

  if (argc > 1) {
    fp = myfopen(argv[1], "w");
    fprintf(fp, "%d %d %d\n", r, s1, s2);
    fclose(fp);
    }

  if ((! found) && (p != NULL) && (s != -1))
    printf("Searched to %d\n", s);
#if VERBOSE
  printf("Working set after sieving about %d bytes\n", (3*r)>>4);
  printf("Dynamic storage allocation %d bytes\n", space);
  printf("%d clock calls\n", clockcalls);
  if (syscalls > 0) printf("%d system/sleep calls\n", syscalls);
#endif
  return(OK);				/* Normal exit */
  }
