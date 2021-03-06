! Routines for computing rankings of an array
!
! Uses the routines `I_unirnk`, `I_unirnk`, `D_nearless` and `I_nearless`
! from ORDERPACK 2.0 by by Michel Olagnon (http://www.fortran-2000.com/rank).
!
! The only updates are to wrap this subset of routines in to a module and to use the
! local type definitions in ConstantsModule.

module Ranking

  use ConstantsModule

  implicit none

  public unirank
  private I_unirank, D_unirank, I_nearless, D_nearless, nearless

  interface unirank
    module procedure D_unirank, I_unirank
  end interface unirank

  interface nearless
    module procedure D_nearless, I_nearless
  end interface nearless

  contains

  subroutine D_unirank (XVALT, IRNGT, NUNI)
  ! __________________________________________________________
  !   UNIRNK = Merge-sort ranking of an array, with removal of
  !   duplicate entries.
  !   The routine is similar to pure merge-sort ranking, but on
  !   the last pass, it discards indices that correspond to
  !   duplicate entries.
  !   For performance reasons, the first 2 passes are taken
  !   out of the standard loop, and use dedicated coding.
  ! __________________________________________________________
  ! __________________________________________________________
        Real(dp), Dimension (:), Intent (In) :: XVALT
        Integer(i4b), Dimension (:), Intent (Out) :: IRNGT
        Integer(i4b), Intent (Out) :: NUNI
  ! __________________________________________________________
        Integer(i4b), Dimension (SIZE(IRNGT)) :: JWRKT
        Integer(i4b) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
        Integer(i4b) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
        Real(dp) :: XTST, XVALA, XVALB
  !
  !
        NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
        NUNI = NVAL
  !
        Select Case (NVAL)
        Case (:0)
           Return
        Case (1)
           IRNGT (1) = 1
           Return
        Case Default
           Continue
        End Select
  !
  !  Fill-in the index array, creating ordered couples
  !
        Do IIND = 2, NVAL, 2
           If (XVALT(IIND-1) < XVALT(IIND)) Then
              IRNGT (IIND-1) = IIND - 1
              IRNGT (IIND) = IIND
           Else
              IRNGT (IIND-1) = IIND
              IRNGT (IIND) = IIND - 1
           End If
        End Do
        If (Modulo(NVAL, 2) /= 0) Then
           IRNGT (NVAL) = NVAL
        End If
  !
  !  We will now have ordered subsets A - B - A - B - ...
  !  and merge A and B couples into     C   -   C   - ...
  !
        LMTNA = 2
        LMTNC = 4
  !
  !  First iteration. The length of the ordered subsets goes from 2 to 4
  !
        Do
           If (NVAL <= 4) Exit
  !
  !   Loop on merges of A and B into C
  !
           Do IWRKD = 0, NVAL - 1, 4
              If ((IWRKD+4) > NVAL) Then
                 If ((IWRKD+2) >= NVAL) Exit
  !
  !   1 2 3
  !
                 If (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) Exit
  !
  !   1 3 2
  !
                 If (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) Then
                    IRNG2 = IRNGT (IWRKD+2)
                    IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                    IRNGT (IWRKD+3) = IRNG2
  !
  !   3 1 2
  !
                 Else
                    IRNG1 = IRNGT (IWRKD+1)
                    IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                    IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                    IRNGT (IWRKD+2) = IRNG1
                 End If
                 Exit
              End If
  !
  !   1 2 3 4
  !
              If (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) Cycle
  !
  !   1 3 x x
  !
              If (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) Then
                 IRNG2 = IRNGT (IWRKD+2)
                 IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                 If (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) Then
  !   1 3 2 4
                    IRNGT (IWRKD+3) = IRNG2
                 Else
  !   1 3 4 2
                    IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                    IRNGT (IWRKD+4) = IRNG2
                 End If
  !
  !   3 x x x
  !
              Else
                 IRNG1 = IRNGT (IWRKD+1)
                 IRNG2 = IRNGT (IWRKD+2)
                 IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                 If (XVALT(IRNG1) <= XVALT(IRNGT(IWRKD+4))) Then
                    IRNGT (IWRKD+2) = IRNG1
                    If (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) Then
  !   3 1 2 4
                       IRNGT (IWRKD+3) = IRNG2
                    Else
  !   3 1 4 2
                       IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                       IRNGT (IWRKD+4) = IRNG2
                    End If
                 Else
  !   3 4 1 2
                    IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                    IRNGT (IWRKD+3) = IRNG1
                    IRNGT (IWRKD+4) = IRNG2
                 End If
              End If
           End Do
  !
  !  The Cs become As and Bs
  !
           LMTNA = 4
           Exit
        End Do
  !
  !  Iteration loop. Each time, the length of the ordered subsets
  !  is doubled.
  !
        Do
           If (2*LMTNA >= NVAL) Exit
           IWRKF = 0
           LMTNC = 2 * LMTNC
  !
  !   Loop on merges of A and B into C
  !
           Do
              IWRK = IWRKF
              IWRKD = IWRKF + 1
              JINDA = IWRKF + LMTNA
              IWRKF = IWRKF + LMTNC
              If (IWRKF >= NVAL) Then
                 If (JINDA >= NVAL) Exit
                 IWRKF = NVAL
              End If
              IINDA = 1
              IINDB = JINDA + 1
  !
  !  One steps in the C subset, that we create in the final rank array
  !
  !  Make a copy of the rank array for the iteration
  !
              JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
              XVALA = XVALT (JWRKT(IINDA))
              XVALB = XVALT (IRNGT(IINDB))
  !
              Do
                 IWRK = IWRK + 1
  !
  !  We still have unprocessed values in both A and B
  !
                 If (XVALA > XVALB) Then
                    IRNGT (IWRK) = IRNGT (IINDB)
                    IINDB = IINDB + 1
                    If (IINDB > IWRKF) Then
  !  Only A still with unprocessed values
                       IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                       Exit
                    End If
                    XVALB = XVALT (IRNGT(IINDB))
                 Else
                    IRNGT (IWRK) = JWRKT (IINDA)
                    IINDA = IINDA + 1
                    If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                    XVALA = XVALT (JWRKT(IINDA))
                 End If
  !
              End Do
           End Do
  !
  !  The Cs become As and Bs
  !
           LMTNA = 2 * LMTNA
        End Do
  !
  !   Last merge of A and B into C, with removal of duplicates.
  !
        IINDA = 1
        IINDB = LMTNA + 1
        NUNI = 0
  !
  !  One steps in the C subset, that we create in the final rank array
  !
        JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
        If (IINDB <= NVAL) Then
          XTST = NEARLESS (Min(XVALT(JWRKT(1)), XVALT(IRNGT(IINDB))))
        Else
          XTST = NEARLESS (XVALT(JWRKT(1)))
        Endif
        Do IWRK = 1, NVAL
  !
  !  We still have unprocessed values in both A and B
  !
           If (IINDA <= LMTNA) Then
              If (IINDB <= NVAL) Then
                 If (XVALT(JWRKT(IINDA)) > XVALT(IRNGT(IINDB))) Then
                    IRNG = IRNGT (IINDB)
                    IINDB = IINDB + 1
                 Else
                    IRNG = JWRKT (IINDA)
                    IINDA = IINDA + 1
                 End If
              Else
  !
  !  Only A still with unprocessed values
  !
                 IRNG = JWRKT (IINDA)
                 IINDA = IINDA + 1
              End If
           Else
  !
  !  Only B still with unprocessed values
  !
              IRNG = IRNGT (IWRK)
           End If
           If (XVALT(IRNG) > XTST) Then
              XTST = XVALT (IRNG)
              NUNI = NUNI + 1
              IRNGT (NUNI) = IRNG
           End If
  !
        End Do
  !
        Return
  !
  end subroutine D_unirank

  subroutine I_unirank (XVALT, IRNGT, NUNI)
  ! __________________________________________________________
  !   UNIRNK = Merge-sort ranking of an array, with removal of
  !   duplicate entries.
  !   The routine is similar to pure merge-sort ranking, but on
  !   the last pass, it discards indices that correspond to
  !   duplicate entries.
  !   For performance reasons, the first 2 passes are taken
  !   out of the standard loop, and use dedicated coding.
  ! __________________________________________________________
  ! __________________________________________________________
        Integer(i4b), Dimension (:), Intent (In) :: XVALT
        Integer(i4b), Dimension (:), Intent (Out) :: IRNGT
        Integer(i4b), Intent (Out) :: NUNI
  ! __________________________________________________________
        Integer(i4b), Dimension (SIZE(IRNGT)) :: JWRKT
        Integer(i4b) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
        Integer(i4b) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
        Integer(i4b) :: XTST, XVALA, XVALB
  !
  !
        NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
        NUNI = NVAL
  !
        Select Case (NVAL)
        Case (:0)
           Return
        Case (1)
           IRNGT (1) = 1
           Return
        Case Default
           Continue
        End Select
  !
  !  Fill-in the index array, creating ordered couples
  !
        Do IIND = 2, NVAL, 2
           If (XVALT(IIND-1) < XVALT(IIND)) Then
              IRNGT (IIND-1) = IIND - 1
              IRNGT (IIND) = IIND
           Else
              IRNGT (IIND-1) = IIND
              IRNGT (IIND) = IIND - 1
           End If
        End Do
        If (Modulo(NVAL, 2) /= 0) Then
           IRNGT (NVAL) = NVAL
        End If
  !
  !  We will now have ordered subsets A - B - A - B - ...
  !  and merge A and B couples into     C   -   C   - ...
  !
        LMTNA = 2
        LMTNC = 4
  !
  !  First iteration. The length of the ordered subsets goes from 2 to 4
  !
        Do
           If (NVAL <= 4) Exit
  !
  !   Loop on merges of A and B into C
  !
           Do IWRKD = 0, NVAL - 1, 4
              If ((IWRKD+4) > NVAL) Then
                 If ((IWRKD+2) >= NVAL) Exit
  !
  !   1 2 3
  !
                 If (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) Exit
  !
  !   1 3 2
  !
                 If (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) Then
                    IRNG2 = IRNGT (IWRKD+2)
                    IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                    IRNGT (IWRKD+3) = IRNG2
  !
  !   3 1 2
  !
                 Else
                    IRNG1 = IRNGT (IWRKD+1)
                    IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                    IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                    IRNGT (IWRKD+2) = IRNG1
                 End If
                 Exit
              End If
  !
  !   1 2 3 4
  !
              If (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) Cycle
  !
  !   1 3 x x
  !
              If (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) Then
                 IRNG2 = IRNGT (IWRKD+2)
                 IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                 If (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) Then
  !   1 3 2 4
                    IRNGT (IWRKD+3) = IRNG2
                 Else
  !   1 3 4 2
                    IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                    IRNGT (IWRKD+4) = IRNG2
                 End If
  !
  !   3 x x x
  !
              Else
                 IRNG1 = IRNGT (IWRKD+1)
                 IRNG2 = IRNGT (IWRKD+2)
                 IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                 If (XVALT(IRNG1) <= XVALT(IRNGT(IWRKD+4))) Then
                    IRNGT (IWRKD+2) = IRNG1
                    If (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) Then
  !   3 1 2 4
                       IRNGT (IWRKD+3) = IRNG2
                    Else
  !   3 1 4 2
                       IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                       IRNGT (IWRKD+4) = IRNG2
                    End If
                 Else
  !   3 4 1 2
                    IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                    IRNGT (IWRKD+3) = IRNG1
                    IRNGT (IWRKD+4) = IRNG2
                 End If
              End If
           End Do
  !
  !  The Cs become As and Bs
  !
           LMTNA = 4
           Exit
        End Do
  !
  !  Iteration loop. Each time, the length of the ordered subsets
  !  is doubled.
  !
        Do
           If (2*LMTNA >= NVAL) Exit
           IWRKF = 0
           LMTNC = 2 * LMTNC
  !
  !   Loop on merges of A and B into C
  !
           Do
              IWRK = IWRKF
              IWRKD = IWRKF + 1
              JINDA = IWRKF + LMTNA
              IWRKF = IWRKF + LMTNC
              If (IWRKF >= NVAL) Then
                 If (JINDA >= NVAL) Exit
                 IWRKF = NVAL
              End If
              IINDA = 1
              IINDB = JINDA + 1
  !
  !  One steps in the C subset, that we create in the final rank array
  !
  !  Make a copy of the rank array for the iteration
  !
              JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
              XVALA = XVALT (JWRKT(IINDA))
              XVALB = XVALT (IRNGT(IINDB))
  !
              Do
                 IWRK = IWRK + 1
  !
  !  We still have unprocessed values in both A and B
  !
                 If (XVALA > XVALB) Then
                    IRNGT (IWRK) = IRNGT (IINDB)
                    IINDB = IINDB + 1
                    If (IINDB > IWRKF) Then
  !  Only A still with unprocessed values
                       IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                       Exit
                    End If
                    XVALB = XVALT (IRNGT(IINDB))
                 Else
                    IRNGT (IWRK) = JWRKT (IINDA)
                    IINDA = IINDA + 1
                    If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                    XVALA = XVALT (JWRKT(IINDA))
                 End If
  !
              End Do
           End Do
  !
  !  The Cs become As and Bs
  !
           LMTNA = 2 * LMTNA
        End Do
  !
  !   Last merge of A and B into C, with removal of duplicates.
  !
        IINDA = 1
        IINDB = LMTNA + 1
        NUNI = 0
  !
  !  One steps in the C subset, that we create in the final rank array
  !
        JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
        If (IINDB <= NVAL) Then
          XTST = NEARLESS (Min(XVALT(JWRKT(1)), XVALT(IRNGT(IINDB))))
        Else
          XTST = NEARLESS (XVALT(JWRKT(1)))
        Endif
        Do IWRK = 1, NVAL
  !
  !  We still have unprocessed values in both A and B
  !
           If (IINDA <= LMTNA) Then
              If (IINDB <= NVAL) Then
                 If (XVALT(JWRKT(IINDA)) > XVALT(IRNGT(IINDB))) Then
                    IRNG = IRNGT (IINDB)
                    IINDB = IINDB + 1
                 Else
                    IRNG = JWRKT (IINDA)
                    IINDA = IINDA + 1
                 End If
              Else
  !
  !  Only A still with unprocessed values
  !
                 IRNG = JWRKT (IINDA)
                 IINDA = IINDA + 1
              End If
           Else
  !
  !  Only B still with unprocessed values
  !
              IRNG = IRNGT (IWRK)
           End If
           If (XVALT(IRNG) > XTST) Then
              XTST = XVALT (IRNG)
              NUNI = NUNI + 1
              IRNGT (NUNI) = IRNG
           End If
  !
        End Do
  !
        Return
  !
  end subroutine I_unirank

  function D_nearless (XVAL) result (D_nl)
  !  Nearest value less than given value
  ! __________________________________________________________
        Real(dp), Intent (In) :: XVAL
        Real(dp) :: D_nl
  ! __________________________________________________________
        D_nl = nearest (XVAL, -1.0_dp)
        return
  !
  end function D_nearless

  function I_nearless (XVAL) result (I_nl)
  !  Nearest value less than given value
  ! __________________________________________________________
        Integer(i4b), Intent (In) :: XVAL
        Integer(i4b) :: I_nl
  ! __________________________________________________________
        I_nl = XVAL - 1
        return
  !
  end function I_nearless

end module