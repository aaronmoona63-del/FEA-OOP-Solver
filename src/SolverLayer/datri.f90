!-----------------------------------------------------------------------
SUBROUTINE datri (al, au, ad, jp, neq, flg, jfile)
!-----------------------------------------------------------------------
!
!.... triangular decomposition of a matrix stored in profile form
!
!.... input parameters
!         al(jp(neq)) - lower triangular part of matrix
!         au(jp(neq)) - upper part of triangular matrix
!         ad(neq)     - diagonals of triangular matrix
!         jp(neq)     - pointers to bottom of colums of al and au arrays
!         neq         - number of equations to be solved
!         flg         - if true equations are unsymmetric
!                       if false equations are symmetric and calling
!                       address of al may be same as that for au
!                       (i.e., au and al share same memory)
!         jfile       - unit number for printed output of warning
!                       messages.
!.... output parameters
!         al(jp(neq)) - lower triangular factor of matrix
!         au(jp(neq)) - upper triangular factor of matrix
!         ad(neq)     - inverse of diagonal matrix in triangular factor
!
!-----------------------------------------------------------------------
  IMPLICIT DOUBLE PRECISION (a - h, o - z)
  LOGICAL flg
  DIMENSION al ( * ), au ( * ), ad ( * ), jp ( * )
  DATA tol / 1.d-7 /
!
!.... n.b.  tol should be set to approximate half-word precision.
!.... loop through the columns to perform the triangular decomposition
!
      jd = 1
      DO j = 1, neq
        jr = jd+1
        jd = jp (j)
        jh = jd-jr
        IF ( jh > 0 ) then
          is = j - jh
          ie = j - 1
!
!.... if diagonal is 0.0 compute a norm for singularity test
!
          IF ( ad (j) == 0.d0 ) CALL datest (au (jr), jh, daval)
          DO i = is, ie
            jr = jr + 1
            id = jp (i)
            ih = min0 (id-jp (i - 1), i - is + 1)
            IF ( ih > 0 ) then
              jrh = jr - ih
              idh = id-ih + 1
              au (jr) = au (jr) - dot (au (jrh), al (idh), ih)
              IF ( flg ) al (jr) = al (jr) &
                                 - dot (al (jrh), au (idh), ih)
            END IF
          END DO
        END IF
!
!.... reduce the diagonal
!
        IF ( jh >= 0 ) then
          dd = ad (j)
          jr = jd-jh
          jrh = j - jh - 1
          CALL dredu (al (jr), au (jr), ad (jrh), jh + 1, flg, ad (j) )
!
!.... check for possible errors and print warnings
!
          IF ( abs (ad (j)) < tol * abs (dd) ) write (jfile, 2000) j
            2000 FORMAT('datri: Solver warning: loss of at least 7 ', &
             & ' digits in reducing diagonal: ',i5)
          IF ( dd < 0.d0 .and. ad (j) > 0.d0 ) write (jfile, 2001) j
          IF ( dd > 0.d0 .and. ad (j) < 0.d0 ) write (jfile, 2001) j
            2001 FORMAT('datri: Solver warning: sign of diagonal ', &
            & 'changed when reducing equation: ',i5)
          IF ( ad (j) == 0.d0 ) write (jfile, 2002) j
            2002 FORMAT('datri: Solver warning: reduced diagonal ', &
            &  'is zero for equation:',i5)
!
!.... complete rank test for a 0.0 diagonal case
!
          IF ( dd == 0.d0 .and. jh > 0 ) then
            IF ( abs (ad (j) ) < tol * daval) write (jfile, 2003) j
              2003 FORMAT('datri: Solver warning: rank failure for ', &
              & 'zero unreduced diagonal in equation:',i5)
          END IF
        END IF
!
!.... store reciprocal of diagonal
!
        IF ( ad (j) /= 0.d0 ) ad (j) = 1.d0 / ad (j)
      END DO
!
END SUBROUTINE datri
 
!-----------------------------------------------------------------------
SUBROUTINE dasol (al, au, ad, b, u, jp, neq, jfile, energy)
!-----------------------------------------------------------------------
!
!.... solution of equations stored in profile form: a * u = b
!.... coefficient matrix must be decomposed into its triangular
!.... factors using datri before using dasol.
!
!.... input parameters
!         al(jp(neq)) - lower triangular factor of matrix
!         au(jp(neq)) - upper triangular factor of matrix
!                       (au and al have same calling address for
!                        symmetric matrices)
!         ad(neq)     - reciprocal of diagonal of triangular factor
!         b(neq)      - right hand side vector in equations
!         jp(neq)     - pointer array to bottom of columns of al and au
!         neq         - number of equations to be solved
!         jfile       - unit number for printed warning message.
!         energy      - energy norm for equations (rhs * soln)
!
!.... output parameter
!         u(neq)      - solution of equations
!
!-----------------------------------------------------------------------
  IMPLICIT DOUBLE PRECISION (a - h, o - z)
  DIMENSION al ( * ), au ( * ), ad ( * ), b ( * ), jp ( * ),       &
      u ( * )
!
!.... find the first non zero entry in the right hand side
!
      DO is = 1, neq
        u (is) = b (is)
        IF ( b (is) /= 0.d0 ) goto 200
      END DO
      WRITE (jfile, 2000)
      2000 FORMAT('dasol: Solver warning: zero right-hand-side vector')
      RETURN
!
!.... reduce the right hand side
!
  200 IF ( is < neq ) then
        DO j = is + 1, neq
          jr = jp (j - 1)
          jh = jp (j) - jr
          u (j) = b (j)
          IF ( jh > 0 ) then
            u (j) = u (j) - dot (al (jr + 1), u (j - jh), jh)
          END IF
        END DO
      END IF
!
!.... multiply by inverse of diagonal elements
!
      energy = 0.d0
      DO j = is, neq
        bd = u (j)
        u (j) = u (j) * ad (j)
        energy = energy + bd * u (j)
      END DO
!
!.... backsubstitution
!
      IF ( neq > 1 ) then
        j = neq
        DO WHILE (j > 1)
          jr = jp (j - 1)
          jh = jp (j) - jr
          IF ( jh > 0 ) then
            CALL colred (au (jr + 1), u (j), jh, u (j - jh) )
          END IF
          j = j - 1
        END DO
      END IF
!
END SUBROUTINE dasol