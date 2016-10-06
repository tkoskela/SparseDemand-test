C *******************************************************
C * Copyright (c) 2006 by Ziena Optimization, Inc.      *
C * All Rights Reserved                                 *
C *******************************************************

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  This file contains an example program that calls KNITRO to solve
C  a simple test problem.  The example is compiled with knitro_fortran.c
C  to gain access to the KNITRO application programming interface (API).
C  The problem is defined in problemQCQP.f.
C
C  This example uses the reverse communications mode of calling KNITRO.
C  Both exact second derivatives and Hessian-vector products are
C  implemented, although KNITRO needs just one Hessian option (in fact,
C  neither implementation is required if KNITRO is told to use a
C  quasi-Newton Hessian approximation).
C
C  The example illustrates only part of the KNITRO API.  Consult the
C  KNITRO Users Manual for full details.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      program          example

      implicit none

C---- DEFINE PARAMETERS THAT BOUND THE PROBLEM SIZE.
C---- THE SIZE OF problemQCQP.f IS QUITE SMALL, BUT A RANGE OF SIZES
C---- ARE DEFINED HERE FOR GENERALITY.
      integer  max_n, max_m, max_nnzj, max_nnzh

C---- GUIDELINES FOR MAXIMUM ARRAY SIZES, BASED ON MEMORY CAPACITY
C---- OF TYPICAL NONLINEAR PROGRAMMING PROBLEMS SOLVED ON A STANDARD
C---- WORKSTATION.  ACTUAL MEMORY REQUIREMENTS VARY, DEPENDING ON
C---- SPARSITY OF THE DERIVATIVES.
C----    VERY SMALL PROBLEMS   n AND m <=    200
C----    SMALL PROBLEMS        n AND m <=   1000
C----    MEDIUM PROBLEMS       n AND m <=  10000
C----    LARGE PROBLEMS        n AND m <=  50000
C----    VERY LARGE PROBLEMS   n AND m <= 250000
C---- WHERE
C----    n IS THE NUMBER OF VARIABLES
C----    m IS THE NUMBER OF CONSTRAINTS

C---- UNCOMMENT ONE OF THESE SECTIONS.

C---- VERY SMALL
C     parameter (max_n       = 200)
C     parameter (max_m       = 200)
C     parameter (max_nnzj    = 40000)
C     parameter (max_nnzh    = 40000)
C
C---- SMALL
C     parameter (max_n       = 1000)
C     parameter (max_m       = 1000)
C     parameter (max_nnzj    = 100000)
C     parameter (max_nnzh    = 100000)
C
C---- MEDIUM
      parameter (max_n       = 10000)
      parameter (max_m       = 10000)
      parameter (max_nnzj    = 200000)
      parameter (max_nnzh    = 200000)
C
C---- LARGE
C     parameter (max_n       = 50000)
C     parameter (max_m       = 50000)
C     parameter (max_nnzj    = 1000000)
C     parameter (max_nnzh    = 1000000)
C
C---- VERY LARGE
C     parameter (max_n       = 250000)
C     parameter (max_m       = 250000)
C     parameter (max_nnzj    = 5000000)
C     parameter (max_nnzh    = 5000000)


C---- DEFINE PARAMETERS TAKEN FROM "knitro.h".
      integer           KTR_RC_EVALFC, KTR_RC_EVALGA,
     $                  KTR_RC_EVALH, KTR_RC_EVALHV
      parameter        (KTR_RC_EVALFC = 1,
     $                  KTR_RC_EVALGA = 2,
     $                  KTR_RC_EVALH  = 3,
     $                  KTR_RC_EVALHV = 7)

C---- DECLARE VARIABLES THAT ARE PASSED TO KNITRO.
      integer           n, m, nnzh, nnzj,
     $                  objgoal, objtype, ctype(max_m),
     $                  jacindvar(max_nnzj), jacindcon(max_nnzj),
     $                  hcol(max_nnzh), hrow(max_nnzh)
      double precision  obj, x(max_n), xinitial(max_n),
     $                  xlobnds(max_n), xupbnds(max_n),
     $                  c(max_m), clobnds(max_m), cupbnds(max_m),
     $                  objgrad(max_n), jac(max_nnzj),
     $                  lambda(max_m + max_n), hess(max_nnzh),
     $                  hvector(max_n)

C---- OTHER LOCAL VARIABLES.
      integer           status, evalstatus


C---- FETCH THE SIZES OF THE PROBLEM TO BE SOLVED.
      call get_problem_sizes (n, m, nnzj, nnzh)

C---- CHECK THE SIZES.
      if (n .gt. max_n) then
         print*,'ERROR: max_n too small.'
         stop
      elseif (m .gt. max_m) then
         print*,'ERROR: max_m too small.'
         stop
      elseif (nnzj .gt. max_nnzj) then
         print*,'ERROR: max_nnzj too small.'
         stop
      elseif (nnzh .gt. max_nnzh) then
         print*,'ERROR: max_nnzh too small.'
         stop
      endif

C---- FETCH THE DEFINITION OF THE PROBLEM TO BE SOLVED.
      call load_problem_def (n, m, nnzj, nnzh, objtype, objgoal,
     $                       xlobnds, xupbnds, xinitial,
     $                       ctype, clobnds, cupbnds,
     $                       jacindvar, jacindcon, hrow, hcol)

C---- CREATE A NEW KNITRO SOLVER INSTANCE.
C---- ANY OLD INSTANCE MUST BE CLOSED BEFORE CREATING A NEW ONE.
      call ktrF_open_instance

C---- READ THE USER OPTIONS FILE, APPLYING ANY OPTIONS SPECIFIED
C---- BY THE USER.
      call ktrF_load_param_file

C---- INITIALIZE KNITRO WITH THE PROBLEM DEFINITION.
      call ktrF_init_problem (n, objgoal, objtype,
     $                        xlobnds, xupbnds,
     $                        m, ctype, clobnds, cupbnds,
     $                        nnzj, jacindvar, jacindcon,
     $                        nnzh, hrow, hcol, xinitial)


C---- SOLVE THE PROBLEM.  IN REVERSE COMMUNICATIONS MODE, KNITRO
C---- RETURNS WHENEVER IT NEEDS MORE PROBLEM INFORMATION.  THE CALLING
C---- PROGRAM MUST INTERPRET KNITRO'S RETURN STATUS AND CONTINUE
C---- SUPPLYING PROBLEM INFORMATION UNTIL KNITRO IS COMPLETE.
C----
C---- RETURN STATUS CODES ARE DEFINED IN "knitro.h" AND DESCRIBED
C---- IN THE KNITRO MANUAL.

 100  continue

      evalstatus = 0
      call ktrF_solve (x, lambda,
     $                 evalstatus, obj, c,
     $                 objgrad, jac, hess, hvector, status)

      if (status .EQ. KTR_RC_EVALFC) then
          call compute_fc (x, obj, c)
      endif

      if (status .EQ. KTR_RC_EVALGA) then
          call compute_ga (x, objgrad, jac)
      endif

      if (status .EQ. KTR_RC_EVALH) then
          call compute_h (x, lambda, hess)
      endif

      if (status .EQ. KTR_RC_EVALHV) then
          call compute_hv (n, x, lambda, hvector)
      endif

      if (status .GT. 0) goto 100

      write (6,600) status
 600  format ('Final KNITRO status is ', I4)


C---- DELETE THE KNITRO SOLVER INSTANCE.
      call ktrF_close_instance


      stop
      end


C----- End of source code -------------------------------------------
