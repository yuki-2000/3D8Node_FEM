!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     3次元弾性FEMプログラム      !!!
!!!  8節点アイソパラメトリック要素  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM MAIN

IMPLICIT NONE


INTEGER i, j, k

INTEGER NUM_NODE, NUM_ELEME, MATERIAL
INTEGER NUM_FIX, NUM_FORCE
DOUBLE PRECISION AMP

DOUBLE PRECISION, ALLOCATABLE :: node(:,:) !節点座標
INTEGER, ALLOCATABLE :: eleme(:,:) !各要素のコネクティビティ
INTEGER, ALLOCATABLE :: mat_num(:) !各要素の材料種類
INTEGER, ALLOCATABLE :: fix_pnt(:,:) !変位境界条件
INTEGER, ALLOCATABLE :: force_pnt(:,:) !力学的境界条件
DOUBLE PRECISION, ALLOCATABLE :: fix(:) !変位境界条件の値
DOUBLE PRECISION, ALLOCATABLE :: force(:) !力学的境界条件の値
DOUBLE PRECISION, ALLOCATABLE :: Bmat(:,:,:,:) !Bマトリックス（形状関数の微分）
DOUBLE PRECISION, ALLOCATABLE :: Dmat(:,:,:) !弾性剛性マトリックス
DOUBLE PRECISION, ALLOCATABLE :: Kmat(:,:) !全体剛性マトリックス
DOUBLE PRECISION, ALLOCATABLE :: Fmat(:) !節点荷重ベクトル
DOUBLE PRECISION, ALLOCATABLE :: Umat(:) !節点変位ベクトル
DOUBLE PRECISION, ALLOCATABLE :: det(:,:) !ヤコビアン
DOUBLE PRECISION gauss(8,3), polar(8,3) !ガウスの積分点，要素座標(xi,eta,zeta)における節点座標
INTEGER, ALLOCATABLE :: known_DOF(:) !既知節点変位ベクトルの自由度
INTEGER, ALLOCATABLE :: unknown_DOF(:) !未知節点変位ベクトルの自由度
DOUBLE PRECISION, ALLOCATABLE :: K11(:,:) !変位境界条件付加後の小行列
DOUBLE PRECISION, ALLOCATABLE :: K12(:,:) !変位境界条件付加後の小行列
DOUBLE PRECISION, ALLOCATABLE :: K22(:,:) !変位境界条件付加後の小行列
DOUBLE PRECISION, ALLOCATABLE :: F1(:) !変位境界条件付加後の小行列
DOUBLE PRECISION, ALLOCATABLE :: F2(:) !変位境界条件付加後の小行列
DOUBLE PRECISION, ALLOCATABLE :: U1(:) !変位境界条件付加後の小行列
DOUBLE PRECISION, ALLOCATABLE :: U2(:) !変位境界条件付加後の小行列

OPEN(10,file='input_AnalysisConditions.txt',form='formatted')
	READ(10,*) NUM_NODE !モデル節点数
	READ(10,*) NUM_ELEME !モデル要素数
	READ(10,*) MATERIAL !材料種類数
	READ(10,*) NUM_FIX !拘束点数
	READ(10,*) NUM_FORCE !荷重点数
	READ(10,*) AMP !変形図倍率
CLOSE(10)

ALLOCATE(node(NUM_NODE,3))
ALLOCATE(eleme(NUM_ELEME,8))
ALLOCATE(mat_num(NUM_ELEME))
ALLOCATE(fix_pnt(NUM_FIX,2))
ALLOCATE(force_pnt(NUM_FORCE,2))
ALLOCATE(fix(NUM_FIX))
ALLOCATE(force(NUM_FORCE))
ALLOCATE(Bmat(6,24,NUM_ELEME,8))
ALLOCATE(Dmat(6,6,MATERIAL))
ALLOCATE(Kmat(3*NUM_NODE,3*NUM_NODE))
ALLOCATE(Fmat(3*NUM_NODE))
ALLOCATE(Umat(3*NUM_NODE))
ALLOCATE(det(NUM_ELEME,8))
ALLOCATE(known_DOF(NUM_FIX))
ALLOCATE(unknown_DOF(3*NUM_NODE-NUM_FIX))
ALLOCATE(K11(3*NUM_NODE-NUM_FIX,3*NUM_NODE-NUM_FIX))
ALLOCATE(K12(3*NUM_NODE-NUM_FIX,NUM_FIX))
ALLOCATE(K22(NUM_FIX,NUM_FIX))
ALLOCATE(F1(3*NUM_NODE-NUM_FIX))
ALLOCATE(F2(NUM_FIX))
ALLOCATE(U1(3*NUM_NODE-NUM_FIX))
ALLOCATE(U2(NUM_FIX))

CALL input (NUM_NODE, NUM_ELEME, NUM_FIX, NUM_FORCE, node, eleme, mat_num, fix_pnt, force_pnt, fix, force)

CALL makeDmat (MATERIAL, Dmat)

CALL makeBmat (NUM_NODE, NUM_ELEME, node, eleme, Bmat, det, gauss, polar)

CALL makeKmat (NUM_NODE, NUM_ELEME, MATERIAL, eleme, mat_num, Bmat, Dmat, Kmat, det)

CALL makeFmat (NUM_NODE, NUM_FORCE, Fmat, force_pnt, force)

CALL makeUmat (NUM_NODE, NUM_FIX, Umat, fix_pnt, fix)

CALL makeSUBmatrix (NUM_NODE, NUM_FIX, Kmat, Fmat, Umat, fix_pnt, known_DOF, unknown_DOF, K11, K12, K22, &
&                   F1, F2, U1, U2)

CALL solveUmat (NUM_NODE, NUM_FIX, Umat, K11, K12, F1, U1, U2, unknown_DOF)

CALL solveFmat (NUM_NODE, NUM_FIX, Fmat, K12, K22, F2, U1, U2, known_DOF)

CALL displacement (NUM_NODE, AMP, node, Umat)

CALL distribution (NUM_NODE, NUM_ELEME, MATERIAL, eleme, mat_num, Bmat, Dmat, Umat, gauss, polar)

DEALLOCATE(node)
DEALLOCATE(eleme)
DEALLOCATE(mat_num)
DEALLOCATE(fix_pnt)
DEALLOCATE(force_pnt)
DEALLOCATE(fix)
DEALLOCATE(force)
DEALLOCATE(Bmat)
DEALLOCATE(Dmat)
DEALLOCATE(Kmat)
DEALLOCATE(Fmat)
DEALLOCATE(Umat)
DEALLOCATE(det)
DEALLOCATE(known_DOF)
DEALLOCATE(unknown_DOF)
DEALLOCATE(K11)
DEALLOCATE(K12)
DEALLOCATE(K22)
DEALLOCATE(F1)
DEALLOCATE(F2)
DEALLOCATE(U1)
DEALLOCATE(U2)


STOP
END PROGRAM

!!!!!=======================================!!!!!

SUBROUTINE input (NUM_NODE, NUM_ELEME, NUM_FIX, NUM_FORCE, node, eleme, mat_num, fix_pnt, force_pnt, fix, force)

IMPLICIT NONE


INTEGER i, j, k
INTEGER NUM_NODE, NUM_ELEME, NUM_FIX, NUM_FORCE
DOUBLE PRECISION node(NUM_NODE,3)
INTEGER eleme(NUM_ELEME,8)
INTEGER mat_num(NUM_ELEME)
INTEGER fix_pnt(NUM_FIX,2)
INTEGER force_pnt(NUM_FORCE,2)
DOUBLE PRECISION fix(NUM_FIX)
DOUBLE PRECISION force(NUM_FORCE)

OPEN(110,file='input_point.txt',form='formatted')
	DO i=1, NUM_NODE
		READ(110,*) k, (node(i,j), j=1, 3)
	END DO
CLOSE(110)

OPEN(120,file='input_eleme.txt',form='formatted')
	DO i=1, NUM_ELEME
		READ(120,*) k, (eleme(i,j), j=1, 8)
	END DO
CLOSE(120)

OPEN(130,file='input_material.txt',form='formatted')
	DO i=1, NUM_ELEME
		READ(130,*) k, mat_num(i)
	END DO
CLOSE(130)

OPEN(140,file='input_fixednodes.txt',form='formatted')
	DO i=1, NUM_FIX
		READ(140,*) k, (fix_pnt(i,j), j=1, 2), fix(i)
	END DO
CLOSE(140)

OPEN(150,file='input_forcednodes.txt',form='formatted')
	DO i=1, NUM_FORCE
		READ(150,*) k, (force_pnt(i,j), j=1, 2), force(i)
	END DO
CLOSE(150)


END SUBROUTINE

!!!!!=======================================!!!!!

SUBROUTINE makeDmat (MATERIAL, Dmat)

IMPLICIT NONE


INTEGER i, j, k, num, mat_type
INTEGER MATERIAL
DOUBLE PRECISION Dmat(6,6,MATERIAL)
DOUBLE PRECISION Young1, Young2
DOUBLE PRECISION Poisson1, Poisson2, Poisson3
DOUBLE PRECISION GL

DO i=1, 6
	DO j=1, 6
		DO k=1, MATERIAL
			Dmat(i,j,k) = 0.0d0
		END DO
	END DO
END DO

OPEN(210,file='input_matinfo.txt',form='formatted')

	DO k=1, MATERIAL
		READ(210,*) num
		WRITE(*,*) 'MATERIAL ID :', num
		IF (k /= num) THEN
			WRITE(*,*) 'MATERIAL ID IS NOT APPROPRIATE.'
			STOP
		END IF
		
		READ(210,*) mat_type
		
		IF (mat_type == 1) THEN !Isotropic Material
			WRITE(*,*) 'MATERIAL TYPE IS ISOTROPIC MATERIAL.'
			
			READ(210,*) Young1
			READ(210,*) Poisson1
			
			WRITE(*,*) 'YOUNG''S MODULUS [MPa] :', Young1
			WRITE(*,*) 'POISSON''S RATIO :', Poisson1
			
			Young1 = Young1 * (10.0d0 ** 6)
			
				Dmat(1,1,k) = Young1 * (1.0d0 - Poisson1) / ( 1.0d0 + Poisson1) / (1.0d0 - 2.0d0 * Poisson1)
				Dmat(2,2,k) = Young1 * (1.0d0 - Poisson1) / ( 1.0d0 + Poisson1) / (1.0d0 - 2.0d0 * Poisson1)
				Dmat(3,3,k) = Young1 * (1.0d0 - Poisson1) / ( 1.0d0 + Poisson1) / (1.0d0 - 2.0d0 * Poisson1)
				Dmat(1,2,k) = Young1 * Poisson1 / ( 1.0d0 + Poisson1) / (1.0d0 - 2.0d0 * Poisson1)
				Dmat(2,1,k) = Young1 * Poisson1 / ( 1.0d0 + Poisson1) / (1.0d0 - 2.0d0 * Poisson1)
				Dmat(2,3,k) = Young1 * Poisson1 / ( 1.0d0 + Poisson1) / (1.0d0 - 2.0d0 * Poisson1)
				Dmat(3,2,k) = Young1 * Poisson1 / ( 1.0d0 + Poisson1) / (1.0d0 - 2.0d0 * Poisson1)
				Dmat(3,1,k) = Young1 * Poisson1 / ( 1.0d0 + Poisson1) / (1.0d0 - 2.0d0 * Poisson1)
				Dmat(1,3,k) = Young1 * Poisson1 / ( 1.0d0 + Poisson1) / (1.0d0 - 2.0d0 * Poisson1)
				Dmat(4,4,k) = Young1 / 2.0d0 / (1.0d0 + Poisson1)
				Dmat(5,5,k) = Young1 / 2.0d0 / (1.0d0 + Poisson1)
				Dmat(6,6,k) = Young1 / 2.0d0 / (1.0d0 + Poisson1)
			
		ELSE IF (mat_type == 2) THEN !Transeversely Isotropic Material
			WRITE(*,*) 'MATERIAL TYPE IS TRANSEVERSELY ISOTROPIC MATERIAL.'
			WRITE(*,*) 'TRANSEVERSELY ISOTROPIC MATERIAL IS NOT APPLIED YET.'
			STOP
			
		ELSE
			WRITE(*,*) 'MATERIAL TYPE IS NOT APPROPRIATE.'
			STOP
		END IF
		
	END DO
	
CLOSE(210)

OPEN(220,file='output_Dmat.dat',form='formatted')
	DO k=1, MATERIAL
		WRITE(220,'(A15,I5)') 'MATERIAL ID :', k
		DO i=1, 6
			WRITE(220,'(6E15.5)') (Dmat(i,j,k), j=1, 6)
		END DO
	END DO
	WRITE(*,*) 'OUTPUT FILE: output_Dmat.dat'
CLOSE(220)

WRITE(*,*) 'MAKE D-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!

SUBROUTINE makeBmat (NUM_NODE, NUM_ELEME, node, eleme, Bmat, det, gauss, polar)

IMPLICIT NONE


INTEGER i, j, k, l, m
INTEGER NUM_NODE, NUM_ELEME
DOUBLE PRECISION node(NUM_NODE,3)
INTEGER eleme(NUM_ELEME,8)
DOUBLE PRECISION Bmat(6,24,NUM_ELEME,8)
DOUBLE PRECISION Hmat(3,8) !dN/d(xi),dN/d(eta),dN/d(zeta)を成分に持つ行列（xi,eta,zetaにはガウスの積分点を代入）
DOUBLE PRECISION det(NUM_ELEME,8) !ガウスの積分点におけるヤコビアン
DOUBLE PRECISION gauss(8,3), polar(8,3)
DOUBLE PRECISION Jacobi(3,3) !ヤコビ行列
DOUBLE PRECISION Jacobiinv(3,3) !ヤコビ行列の逆行列
DOUBLE PRECISION dNdxy(3,8) !dN/dx,dN/dyを成分に持つ行列
DOUBLE PRECISION e_point(8,3) !要素内節点番号での要素座標
DOUBLE PRECISION, PARAMETER :: ga = 0.577350269189626d0 !1/sqrt(3)

gauss(1,1:3) = (/ -ga, -ga, -ga /)
gauss(2,1:3) = (/  ga, -ga, -ga /)
gauss(3,1:3) = (/  ga,  ga, -ga /)
gauss(4,1:3) = (/ -ga,  ga, -ga /)
gauss(5,1:3) = (/ -ga, -ga,  ga /)
gauss(6,1:3) = (/  ga, -ga,  ga /)
gauss(7,1:3) = (/  ga,  ga,  ga /)
gauss(8,1:3) = (/ -ga,  ga,  ga /)

polar(1,1:3) = (/ -1.0d0, -1.0d0, -1.0d0 /)
polar(2,1:3) = (/  1.0d0, -1.0d0, -1.0d0 /)
polar(3,1:3) = (/  1.0d0,  1.0d0, -1.0d0 /)
polar(4,1:3) = (/ -1.0d0,  1.0d0, -1.0d0 /)
polar(5,1:3) = (/ -1.0d0, -1.0d0,  1.0d0 /)
polar(6,1:3) = (/  1.0d0, -1.0d0,  1.0d0 /)
polar(7,1:3) = (/  1.0d0,  1.0d0,  1.0d0 /)
polar(8,1:3) = (/ -1.0d0,  1.0d0,  1.0d0 /)

DO i=1, 6
	DO j=1, 24
		DO k=1, NUM_ELEME
			DO l=1, 8
				Bmat(i,j,k,l) = 0.0d0
			END DO
		END DO
	END DO
END DO

DO i=1, NUM_ELEME

	DO j=1, 8
		DO k=1, 3
			e_point(j,k) = node(eleme(i,j),k)
		END DO
	END DO
	
	DO j=1, 8
		DO k=1, 8
			Hmat(1,k) = polar(k,1) * (1.0d0 + polar(k,2) * gauss(j,2)) * (1.0d0 + polar(k,3) * gauss(j,3)) / 8.0d0
			Hmat(2,k) = polar(k,2) * (1.0d0 + polar(k,1) * gauss(j,1)) * (1.0d0 + polar(k,3) * gauss(j,3)) / 8.0d0
			Hmat(3,k) = polar(k,3) * (1.0d0 + polar(k,1) * gauss(j,1)) * (1.0d0 + polar(k,2) * gauss(j,2)) / 8.0d0
		END DO
		
		DO k=1, 3
			DO l=1, 3
				Jacobi(k,l) = 0.0d0
			END DO
		END DO
		
		DO k=1, 3
			DO l=1, 3
				DO m=1, 8
					Jacobi(k,l) = Jacobi(k,l) + Hmat(k,m) * e_point(m,l)
				END DO
			END DO
		END DO
		
		det(i,j) = Jacobi(1,1) * Jacobi(2,2) * Jacobi(3,3) + Jacobi(1,2) * Jacobi(2,3) * Jacobi(3,1) &
		       & + Jacobi(1,3) * Jacobi(2,1) * Jacobi(3,2) - Jacobi(1,3) * Jacobi(2,2) * Jacobi(3,1) &
		       & - Jacobi(1,2) * Jacobi(2,1) * Jacobi(3,3) - Jacobi(1,1) * Jacobi(2,3) * Jacobi(3,2)
		
		Jacobiinv(1,1) = (Jacobi(2,2) * Jacobi(3,3) - Jacobi(2,3) * Jacobi(3,2)) / det(i,j)
		Jacobiinv(1,2) = (Jacobi(1,3) * Jacobi(3,2) - Jacobi(1,2) * Jacobi(3,3)) / det(i,j)
		Jacobiinv(1,3) = (Jacobi(1,2) * Jacobi(2,3) - Jacobi(1,3) * Jacobi(2,2)) / det(i,j)
		Jacobiinv(2,1) = (Jacobi(2,3) * Jacobi(3,1) - Jacobi(2,1) * Jacobi(3,3)) / det(i,j)
		Jacobiinv(2,2) = (Jacobi(1,1) * Jacobi(3,3) - Jacobi(1,3) * Jacobi(3,1)) / det(i,j)
		Jacobiinv(2,3) = (Jacobi(1,3) * Jacobi(2,1) - Jacobi(1,1) * Jacobi(2,3)) / det(i,j)
		Jacobiinv(3,1) = (Jacobi(2,1) * Jacobi(3,2) - Jacobi(2,2) * Jacobi(3,1)) / det(i,j)
		Jacobiinv(3,2) = (Jacobi(1,2) * Jacobi(3,1) - Jacobi(1,1) * Jacobi(3,2)) / det(i,j)
		Jacobiinv(3,3) = (Jacobi(1,1) * Jacobi(2,2) - Jacobi(1,2) * Jacobi(2,1)) / det(i,j)
		
		DO k=1, 3
			DO l=1, 8
				dNdxy(k,l) = 0.0d0
			END DO 
		END DO
		
		DO k=1, 3
			DO l=1, 8
				DO m=1, 3
					dNdxy(k,l) = dNdxy(k,l) + Jacobiinv(k,m) * Hmat(m,l)
				END DO
			END DO
		END DO
		
		DO k=1, 8
			Bmat(1,3*(k-1)+1,i,j) = dNdxy(1,k)
			Bmat(2,3*(k-1)+2,i,j) = dNdxy(2,k)
			Bmat(3,3*(k-1)+3,i,j) = dNdxy(3,k)
			Bmat(4,3*(k-1)+2,i,j) = dNdxy(3,k)
			Bmat(4,3*(k-1)+3,i,j) = dNdxy(2,k)
			Bmat(5,3*(k-1)+1,i,j) = dNdxy(3,k)
			Bmat(5,3*(k-1)+3,i,j) = dNdxy(1,k)
			Bmat(6,3*(k-1)+1,i,j) = dNdxy(2,k)
			Bmat(6,3*(k-1)+2,i,j) = dNdxy(1,k)
		END DO
	END DO
	
END DO

WRITE(*,*) 'MAKE B-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!

SUBROUTINE makeKmat (NUM_NODE, NUM_ELEME, MATERIAL, eleme, mat_num, Bmat, Dmat, Kmat, det)

IMPLICIT NONE


INTEGER i, j, k, l, m, pt1, pt2
INTEGER NUM_NODE, NUM_ELEME, MATERIAL
INTEGER eleme(NUM_ELEME,8)
INTEGER mat_num(NUM_ELEME)
DOUBLE PRECISION Bmat(6,24,NUM_ELEME,8)
DOUBLE PRECISION Dmat(6,6,MATERIAL)
DOUBLE PRECISION Kmat(3*NUM_NODE,3*NUM_NODE)
DOUBLE PRECISION det(NUM_ELEME,8)
DOUBLE PRECISION e_Kmat(24,24) !要素剛性マトリックス
DOUBLE PRECISION BTD(24,6)
DOUBLE PRECISION BTDB(24,24)

DO i=1, 3*NUM_NODE
	DO j=1, 3*NUM_NODE
		Kmat(i,j) = 0.0d0
	END DO
END DO

DO i=1, NUM_ELEME

	DO j=1, 24
		DO k=1, 24
			e_Kmat(j,k) = 0.0d0
		END DO
	END DO
	
	DO j=1, 8
		
		DO k=1, 24
			DO l=1, 6
				BTD(k,l) = 0.0d0
			END DO
			DO l=1, 24
				BTDB(k,l) = 0.0d0
			END DO
		END DO
		
		DO k=1, 24
			DO l=1, 6
				DO m=1, 6
					BTD(k,l) = BTD(k,l) + Bmat(m,k,i,j) * Dmat(m,l,mat_num(i))
				END DO
			END DO
		END DO
		
		DO k=1, 24
			DO l=1, 24
				DO m=1, 6
					BTDB(k,l) = BTDB(k,l) + BTD(k,m) * Bmat(m,l,i,j)
				END DO
			END DO
		END DO
		
		DO k=1, 24
			DO l=1, 24
				e_Kmat(k,l) = e_Kmat(k,l) + BTDB(k,l) * det(i,j)
			END DO
		END DO
		
	END DO
	
	DO j=1, 8
		DO k=1, 8
			pt1 = eleme(i,j)
			pt2 = eleme(i,k)
			DO l=1, 3
				DO m=1, 3
					Kmat(3*(pt1-1)+l,3*(pt2-1)+m) = Kmat(3*(pt1-1)+l,3*(pt2-1)+m) + e_Kmat(3*(j-1)+l,3*(k-1)+m)
				END DO
			END DO
		END DO
	END DO
	
END DO

OPEN(310,file='output_Kmat.dat',form='formatted')
	DO i=1, 3*NUM_NODE
		WRITE(310,*) (Kmat(i,j), j=1, 3*NUM_NODE)
	END DO
	WRITE(*,*) 'OUTPUT FILE: output_Kmat.dat'
CLOSE(310)

WRITE(*,*) 'MAKE K-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!

SUBROUTINE makeFmat (NUM_NODE, NUM_FORCE, Fmat, force_pnt, force)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_FORCE

INTEGER i, j, k
INTEGER force_pnt(NUM_FORCE,2)
DOUBLE PRECISION force(NUM_FORCE)
DOUBLE PRECISION Fmat(3*NUM_NODE)

DO i=1, 3*NUM_NODE
	Fmat(i) = 0.0d0
END DO

DO i=1, NUM_FORCE
	Fmat(3*(force_pnt(i,1)-1)+force_pnt(i,2)) = force(i)
	
	IF (force_pnt(i,2) /= 1 .AND. force_pnt(i,2) /= 2 .AND. force_pnt(i,2) /=3) THEN
		WRITE(*,*) 'INPUT DATA "input_forcednodes.txt" IS NOT APPROPRIATE.'
		STOP
	END IF
END DO

WRITE(*,*) 'MAKE F-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!

SUBROUTINE makeUmat (NUM_NODE, NUM_FIX, Umat, fix_pnt, fix)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_FIX
INTEGER i, j, k
INTEGER fix_pnt(NUM_FIX,2)
DOUBLE PRECISION fix(NUM_FIX)
DOUBLE PRECISION Umat(3*NUM_NODE)

DO i=1, 3*NUM_NODE
	Umat(i) = 0.0d0
END DO

DO i=1, NUM_FIX
	Umat(3*(fix_pnt(i,1)-1)+fix_pnt(i,2)) = fix(i)
	
	IF(fix_pnt(i,2) /= 1 .AND. fix_pnt(i,2) /= 2 .AND. fix_pnt(i,2) /= 3) THEN
		WRITE(*,*) 'INPUT DATA "input_fixednodes.txt" IS NOT APPROPRIATE.'
		STOP
	END IF
END DO

WRITE(*,*) 'MAKE U-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!

SUBROUTINE makeSUBmatrix (NUM_NODE, NUM_FIX, Kmat, Fmat, Umat, fix_pnt, known_DOF, unknown_DOF, K11, K12, K22, &
&                         F1, F2, U1, U2)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_FIX
INTEGER i, j, k
INTEGER TEMP
DOUBLE PRECISION Kmat(3*NUM_NODE,3*NUM_NODE)
DOUBLE PRECISION Fmat(3*NUM_NODE)
DOUBLE PRECISION Umat(3*NUM_NODE)
INTEGER fix_pnt(NUM_FIX,2)
INTEGER known_DOF(NUM_FIX)
INTEGER unknown_DOF(3*NUM_NODE-NUM_FIX)
DOUBLE PRECISION K11(3*NUM_NODE-NUM_FIX,3*NUM_NODE-NUM_FIX)
DOUBLE PRECISION K12(3*NUM_NODE-NUM_FIX,NUM_FIX)
DOUBLE PRECISION K22(NUM_FIX,NUM_FIX)
DOUBLE PRECISION F1(3*NUM_NODE-NUM_FIX)
DOUBLE PRECISION F2(NUM_FIX)
DOUBLE PRECISION U1(3*NUM_NODE-NUM_FIX)
DOUBLE PRECISION U2(NUM_FIX)

DO i=1, NUM_FIX
	F2(i) = 0.0d0
END DO

DO i=1, 3*NUM_NODE-NUM_FIX
	U1(i) = 0.0d0
END DO

DO i=1, NUM_FIX
	known_DOF(i) = 3*(fix_pnt(i,1)-1) + fix_pnt(i,2)
END DO

!バブルソート
DO i=NUM_FIX-1, 1, -1
	DO j=1, i
		IF (known_DOF(j) > known_DOF(j+1)) THEN
			TEMP = known_DOF(j)
			known_DOF(j) = known_DOF(j+1)
			known_DOF(j+1) = TEMP
		END IF
	END DO
END DO

DO j=1, known_DOF(1)-1
	unknown_DOF(j) = j
END DO

DO i=2, NUM_FIX
	DO j=known_DOF(i-1)+1, known_DOF(i)-1
		unknown_DOF(j-(i-1)) = j
	END DO
END DO

DO j=known_DOF(NUM_FIX)+1, 3*NUM_NODE
	unknown_DOF(j-NUM_FIX) = j
END DO

DO i=1, 3*NUM_NODE-NUM_FIX
	DO j=1, 3*NUM_NODE-NUM_FIX
		K11(i,j) = Kmat(unknown_DOF(i),unknown_DOF(j))
	END DO
	
	DO j=1, NUM_FIX
		K12(i,j) = Kmat(unknown_DOF(i),known_DOF(j))
	END DO
	
	F1(i) = Fmat(unknown_DOF(i))
	U1(i) = Umat(unknown_DOF(i))
END DO

DO i=1, NUM_FIX
	DO j=1, NUM_FIX
		K22(i,j) = Kmat(known_DOF(i),known_DOF(j))
	END DO
	
	F2(i) = Fmat(known_DOF(i))
	U2(i) = Umat(known_DOF(i))
END DO

WRITE(*,*) 'MAKE SUB-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!

SUBROUTINE solveUmat (NUM_NODE, NUM_FIX, Umat, K11, K12, F1, U1, U2, unknown_DOF)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_FIX
INTEGER i, j, k
DOUBLE PRECISION Umat(3*NUM_NODE)
DOUBLE PRECISION K11(3*NUM_NODE-NUM_FIX,3*NUM_NODE-NUM_FIX)
DOUBLE PRECISION K12(3*NUM_NODE-NUM_FIX,NUM_FIX)
DOUBLE PRECISION F1(3*NUM_NODE-NUM_FIX)
DOUBLE PRECISION U1(3*NUM_NODE-NUM_FIX)
DOUBLE PRECISION U2(NUM_FIX)
DOUBLE PRECISION FKU(3*NUM_NODE-NUM_FIX)
INTEGER unknown_DOF(3*NUM_NODE-NUM_FIX)
DOUBLE PRECISION P, Q, S

DO i=1, 3*NUM_NODE-NUM_FIX
	FKU(i) = 0.0d0
END DO

!Gaussの消去法を用いて，U1を求める．
DO i=1, 3*NUM_NODE-NUM_FIX
	DO j=1, NUM_FIX
		FKU(i) = FKU(i) + K12(i,j) * U2(j)
	END DO
	FKU(i) = F1(i) - FKU(i)
END DO

DO i=1, 3*NUM_NODE-NUM_FIX-1
	P = K11(i,i)
	DO j=i+1, 3*NUM_NODE-NUM_FIX
		K11(i,j) = K11(i,j) / P
	END DO
	FKU(i) = FKU(i) / P
	DO j=i+1, 3*NUM_NODE-NUM_FIX
		Q = K11(j,i)
		DO k=i+1, 3*NUM_NODE-NUM_FIX
			K11(j,k) = K11(j,k) - Q * K11(i,k)
		END DO
		FKU(j) = FKU(j) - Q * FKU(i)
	END DO
END DO

U1(3*NUM_NODE-NUM_FIX) = FKU(3*NUM_NODE-NUM_FIX) / K11(3*NUM_NODE-NUM_FIX,3*NUM_NODE-NUM_FIX)

DO i=3*NUM_NODE-NUM_FIX-1, 1, -1
	S = FKU(i)
	DO j=i+1, 3*NUM_NODE-NUM_FIX
		S = S - K11(i,j) * U1(j)
	END DO
	U1(i) = S
END DO

DO i=1, 3*NUM_NODE-NUM_FIX
	Umat(unknown_DOF(i)) = U1(i)
END DO

OPEN(410,file='output_Umat.dat',form='formatted')
	DO i=1, 3*NUM_NODE
		WRITE(410,'(I5,E15.5)') i, Umat(i)
	END DO
	WRITE(*,*) 'OUTPUT FILE: output_Umat.dat'
CLOSE(410)

WRITE(*,*) 'SOLVE U-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!

SUBROUTINE solveFmat (NUM_NODE, NUM_FIX, Fmat, K12, K22, F2, U1, U2, known_DOF)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_FIX
INTEGER i, j, k
DOUBLE PRECISION Fmat(3*NUM_NODE)
DOUBLE PRECISION K12(3*NUM_NODE-NUM_FIX,NUM_FIX)
DOUBLE PRECISION K22(NUM_FIX,NUM_FIX)
DOUBLE PRECISION F2(NUM_FIX)
DOUBLE PRECISION U1(3*NUM_NODE-NUM_FIX)
DOUBLE PRECISION U2(NUM_FIX)
INTEGER known_DOF(NUM_FIX)

DO i=1, NUM_FIX
	DO j=1, 3*NUM_NODE-NUM_FIX
		F2(i) = F2(i) + K12(j,i) * U1(j)
	END DO
	
	DO j=1, NUM_FIX
		F2(i) = F2(i) + K22(i,j) * U2(j)
	END DO
END DO

DO i=1, NUM_FIX
	Fmat(known_DOF(i)) = F2(i)
END DO

OPEN(420,file='output_Fmat.dat',form='formatted')
	DO i=1, 3*NUM_NODE
		WRITE(420,'(I5,F15.5)') i, Fmat(i)
	END DO
	WRITE(*,*) 'OUTPUT FILE: output_Fmat.dat'
CLOSE(420)

WRITE(*,*) 'SOLVE F-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!

SUBROUTINE displacement (NUM_NODE, AMP, node, Umat)

IMPLICIT NONE


INTEGER NUM_NODE
DOUBLE PRECISION AMP
INTEGER i, j, k
DOUBLE PRECISION node(NUM_NODE,3)
DOUBLE PRECISION disp(NUM_NODE,3)
DOUBLE PRECISION Umat(3*NUM_NODE)

DO i=1, NUM_NODE
	DO j=1, 3
		disp(i,j) = 0.0d0
	END DO
END DO

DO i=1, NUM_NODE
	DO j=1, 3
		disp(i,j) = node(i,j) + Umat(3*(i-1)+j) * AMP
	END DO
END DO

OPEN(510,file='output_disp.dat',form='formatted')
	DO i=1, NUM_NODE
		WRITE(510,'(I5,3F15.10)') i, (disp(i,j), j=1, 3)
	END DO
	WRITE(*,*) 'OUTPUT FILE: output_disp.dat'
CLOSE(510)

WRITE(*,*) 'CALCULATE DISPLACEMENT'


END SUBROUTINE

!!!!!=======================================!!!!

SUBROUTINE distribution (NUM_NODE, NUM_ELEME, MATERIAL, eleme, mat_num, Bmat, Dmat, Umat, gauss, polar)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_ELEME, MATERIAL
INTEGER i, j, k, l, m
INTEGER eleme(NUM_ELEME,8)
INTEGER mat_num(NUM_ELEME)
DOUBLE PRECISION Bmat(6,24,NUM_ELEME,8)
DOUBLE PRECISION Dmat(6,6,MATERIAL)
DOUBLE PRECISION Umat(3*NUM_NODE)
DOUBLE PRECISION e_Umat(24)
DOUBLE PRECISION gauss(8,3), polar(8,3)
DOUBLE PRECISION Nmat(8,8)
DOUBLE PRECISION GAUSSstrain(NUM_ELEME,8,6), GAUSSstress(NUM_ELEME,8,6)
DOUBLE PRECISION NODALstrain(NUM_ELEME,8,6), NODALstress(NUM_ELEME,8,6)
DOUBLE PRECISION P, Q, S

DO i=1, NUM_ELEME
	DO j=1, 8
		DO k=1, 6
			GAUSSstrain(i,j,k) = 0.0d0
			GAUSSstress(i,j,k) = 0.0d0
			NODALstrain(i,j,k) = 0.0d0
			NODALstress(i,j,k) = 0.0d0
		END DO
	END DO
END DO

DO i=1, NUM_ELEME
	DO j=1, 8
		DO k=1, 3
			e_Umat(3*(j-1)+k) = Umat(3*(eleme(i,j)-1)+k)
		END DO
	END DO
	
	DO j=1, 8
		DO k=1, 6
			DO l=1, 24
				GAUSSstrain(i,j,k) = GAUSSstrain(i,j,k) + Bmat(k,l,i,j) * e_Umat(l)
			END DO
		END DO
	END DO
	
	DO j=1, 8
		DO k=1, 6
			DO l=1, 6
				GAUSSstress(i,j,k) = GAUSSstress(i,j,k) + Dmat(k,l,mat_num(i)) * GAUSSstrain(i,j,l)
			END DO
		END DO
	END DO
END DO

OPEN(610,file='output_GAUSSstrain.dat',form='formatted')
	DO i=1, NUM_ELEME
		DO j=1, 8
			WRITE(610,'(2I6,6F20.15)') i, eleme(i,j), (GAUSSstrain(i,j,k), k=1, 6)
		END DO
	END DO
	WRITE(*,*) 'OUTPUT FILE: output_GAUSSstrain.dat'
CLOSE(610)

OPEN(620,file='output_GAUSSstress.dat',form='formatted')
	DO i=1, NUM_ELEME
		DO j=1, 8
			WRITE(620,'(2I6,6F20.5)') i, eleme(i,j), (GAUSSstress(i,j,k), k=1, 6)
		END DO
	END DO
	WRITE(*,*) 'OUTPUT FILE: output_GAUSSstress.dat'
CLOSE(620)

WRITE(*,*) 'CALCULATE DISTRIBUTIONS'

DO i=1, NUM_ELEME

	DO j=1, 6
	
		DO k=1, 8
			DO l=1, 8
				Nmat(k,l) = (1.0d0 + polar(l,1) * gauss(k,1)) * (1.0d0 + polar(l,2) * gauss(k,2))&
				        & * (1.0d0 + polar(l,3) * gauss(k,3)) / 8.0d0
			END DO
		END DO
		
		DO k=1, 8-1
			P = Nmat(k,k)
			DO l=k+1, 8
				Nmat(k,l) = Nmat(k,l) / P
			END DO
			GAUSSstrain(i,k,j) = GAUSSstrain(i,k,j) / P
			DO l=k+1, 8
				Q = Nmat(l,k)
				DO m=k+1, 8
					Nmat(l,m) = Nmat(l,m) - Q * Nmat(k,m)
				END DO
				GAUSSstrain(i,l,j) = GAUSSstrain(i,l,j) - Q * GAUSSstrain(i,k,j)
			END DO
		END DO
		
		NODALstrain(i,8,j) = GAUSSstrain(i,8,j) / Nmat(8,8)
		
		DO k=8-1, 1, -1
			S = GAUSSstrain(i,k,j)
			DO l=k+1, 8
				S = S - Nmat(k,l) * NODALstrain(i,l,j)
			END DO
			NODALstrain(i,k,j) = S
		END DO
		
	END DO
	
END DO

DO i=1, NUM_ELEME

	DO j=1, 6
	
		DO k=1, 8
			DO l=1, 8
				Nmat(k,l) = (1.0d0 + polar(l,1) * gauss(k,1)) * (1.0d0 + polar(l,2) * gauss(k,2))&
				        & * (1.0d0 + polar(l,3) * gauss(k,3)) / 8.0d0
			END DO
		END DO
		
		DO k=1, 8-1
			P = Nmat(k,k)
			DO l=k+1, 8
				Nmat(k,l) = Nmat(k,l) / P
			END DO
			GAUSSstress(i,k,j) = GAUSSstress(i,k,j) / P
			DO l=k+1, 8
				Q = Nmat(l,k)
				DO m=k+1, 8
					Nmat(l,m) = Nmat(l,m) - Q * Nmat(k,m)
				END DO
				GAUSSstress(i,l,j) = GAUSSstress(i,l,j) - Q * GAUSSstress(i,k,j)
			END DO
		END DO
		
		NODALstress(i,8,j) = GAUSSstress(i,8,j) / Nmat(8,8)
		
		DO k=8-1, 1, -1
			S = GAUSSstress(i,k,j)
			DO l=k+1, 8
				S = S - Nmat(k,l) * NODALstress(i,l,j)
			END DO
			NODALstress(i,k,j) = S
		END DO
		
	END DO
	
END DO

OPEN(630,file='output_NODALstrain.dat',form='formatted')
	DO i=1, NUM_ELEME
		DO j=1, 8
			WRITE(630,'(2I6,6F20.15)') i, eleme(i,j), (NODALstrain(i,j,k), k=1, 6)
		END DO
	END DO
	WRITE(*,*) 'OUTPUT FILE: output_NODALstrain.dat'
CLOSE(630)

OPEN(640,file='output_NODALstress.dat',form='formatted')
	DO i=1, NUM_ELEME
		DO j=1, 8
			WRITE(640,'(2I6,6F20.5)') i, eleme(i,j), (NODALstress(i,j,k), k=1, 6)
		END DO
	END DO
	WRITE(*,*) 'OUTPUT FILE: output_NODALstress.dat'
CLOSE(640)

WRITE(*,*) 'CALCULATE DISTRIBUTIONS'


END SUBROUTINE