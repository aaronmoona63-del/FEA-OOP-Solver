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