!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                          Ver.1.1 (2019.12.1) !!!
!!!    2次元弾性FEMプログラム（3節点CST要素）    !!!
!!! 計算力学 -有限要素法の基礎-（森下出版）5.2章 !!!
!!!                                              !!!
!!! 2019.12.1 Ver.1.1 モデル厚さを変数に追加     !!!
!!! 2017.7.11 Ver.1.0                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM main

IMPLICIT NONE


INTEGER i, j, k
INTEGER NUM_NODE, NUM_ELEME
INTEGER NUM_FIX, NUM_FORCE
DOUBLE PRECISION THICKNESS
DOUBLE PRECISION AMP

DOUBLE PRECISION, ALLOCATABLE :: node(:,:) !節点座標
INTEGER, ALLOCATABLE :: eleme(:,:) !各要素のコネクティビティ
INTEGER, ALLOCATABLE :: fix_pnt(:,:) !変位境界条件
INTEGER, ALLOCATABLE :: force_pnt(:,:) !力学的境界条件
DOUBLE PRECISION, ALLOCATABLE :: fix(:) !変位境界条件の値
DOUBLE PRECISION, ALLOCATABLE :: force(:) !力学的境界条件の値
DOUBLE PRECISION, ALLOCATABLE :: Bmat(:,:,:) !Bマトリックス（形状関数の偏微分）
DOUBLE PRECISION, ALLOCATABLE :: Dmat(:,:) !弾性剛性マトリックス
DOUBLE PRECISION, ALLOCATABLE :: Kmat(:,:) !全体剛性マトリックス
DOUBLE PRECISION, ALLOCATABLE :: Fmat(:) !節点荷重ベクトル
DOUBLE PRECISION, ALLOCATABLE :: Umat(:) !節点変位ベクトル
DOUBLE PRECISION, ALLOCATABLE :: Ae(:) !要素面積
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
  READ(10,*) THICKNESS !モデル厚さ
  READ(10,*) NUM_FIX !拘束点数
  READ(10,*) NUM_FORCE !荷重点数
  READ(10,*) AMP !変形図倍率
CLOSE(10)

ALLOCATE(node(NUM_NODE,2))
ALLOCATE(eleme(NUM_ELEME,3))
ALLOCATE(fix_pnt(NUM_FIX,2))
ALLOCATE(force_pnt(NUM_FORCE,2))
ALLOCATE(fix(NUM_FIX))
ALLOCATE(force(NUM_FORCE))
ALLOCATE(Bmat(3,6,NUM_ELEME))
ALLOCATE(Dmat(3,3))
ALLOCATE(Kmat(2*NUM_NODE,2*NUM_NODE))
ALLOCATE(Fmat(2*NUM_NODE))
ALLOCATE(Umat(2*NUM_NODE))
ALLOCATE(Ae(NUM_ELEME))
ALLOCATE(known_DOF(NUM_FIX))
ALLOCATE(unknown_DOF(2*NUM_NODE-NUM_FIX))
ALLOCATE(K11(2*NUM_NODE-NUM_FIX,2*NUM_NODE-NUM_FIX))
ALLOCATE(K12(2*NUM_NODE-NUM_FIX,NUM_FIX))
ALLOCATE(K22(NUM_FIX,NUM_FIX))
ALLOCATE(F1(2*NUM_NODE-NUM_FIX))
ALLOCATE(F2(NUM_FIX))
ALLOCATE(U1(2*NUM_NODE-NUM_FIX))
ALLOCATE(U2(NUM_FIX))

CALL input (NUM_NODE, NUM_ELEME, NUM_FIX, NUM_FORCE, node, eleme, fix_pnt, force_pnt, fix, force)

CALL makeDmat (Dmat)

CALL makeBmat (NUM_NODE, NUM_ELEME, node, eleme, Bmat, Ae)

CALL makeKmat (NUM_NODE, NUM_ELEME, eleme, Bmat, Dmat, Kmat, Ae, THICKNESS)

CALL makeFmat (NUM_NODE, NUM_FORCE, Fmat, force_pnt, force)

CALL makeUmat (NUM_NODE, NUM_FIX, Umat, fix_pnt, fix)

CALL makeSUBmatrix (NUM_NODE, NUM_FIX, Kmat, Fmat, Umat, fix_pnt, known_DOF, unknown_DOF, K11, K12, K22, &
&                   F1, F2, U1, U2)

CALL makeInverse(NUM_NODE, NUM_FIX, K11)

CALL solveUmat(NUM_NODE, NUM_FIX, Umat, K11, K12, F1, U1, U2, unknown_DOF)

CALL solveFmat(NUM_NODE, NUM_FIX, Fmat, K12, K22, F2, U1, U2, known_DOF)

CALL displacement (NUM_NODE, AMP, node, Umat)

CALL distribution (NUM_NODE, NUM_ELEME, eleme, Bmat, Dmat, Umat)

DEALLOCATE(node)
DEALLOCATE(eleme)
DEALLOCATE(fix_pnt)
DEALLOCATE(force_pnt)
DEALLOCATE(fix)
DEALLOCATE(force)
DEALLOCATE(Bmat)
DEALLOCATE(Dmat)
DEALLOCATE(Kmat)
DEALLOCATE(Fmat)
DEALLOCATE(Umat)
DEALLOCATE(Ae)
DEALLOCATE(known_DOF)
DEALLOCATE(unknown_DOF)
DEALLOCATE(K11)
DEALLOCATE(K12)
DEALLOCATE(K22)
DEALLOCATE(F1)
DEALLOCATE(F2)
DEALLOCATE(U1)
DEALLOCATE(U2)


END PROGRAM

!!!!!=======================================!!!!!
SUBROUTINE input (NUM_NODE, NUM_ELEME, NUM_FIX, NUM_FORCE, node, eleme, fix_pnt, force_pnt, fix, force)

IMPLICIT NONE


INTEGER i, j, k
INTEGER NUM_NODE, NUM_ELEME, NUM_FIX, NUM_FORCE
DOUBLE PRECISION node(NUM_NODE,2)
INTEGER eleme(NUM_ELEME,3)
INTEGER fix_pnt(NUM_FIX,2)
INTEGER force_pnt(NUM_FORCE,2)
DOUBLE PRECISION fix(NUM_FIX)
DOUBLE PRECISION force(NUM_FORCE)

OPEN(10,file='input_point.txt',form='formatted')
  DO i=1, NUM_NODE
    READ(10,*) k, (node(i,j),j=1,2)
  END DO
CLOSE(10)

OPEN(20,file='input_eleme.txt',form='formatted')
  DO i=1, NUM_ELEME
    READ(20,*) k, (eleme(i,j),j=1,3)
  END DO
CLOSE(20)

OPEN(30,file='input_fixednodes.txt',form='formatted')
  DO i=1, NUM_FIX
    READ(30,*) k, (fix_pnt(i,j),j=1,2), fix(i)
  END DO
CLOSE(30)

OPEN(40,file='input_forcednodes.txt',form='formatted')
  DO i=1, NUM_FORCE
    READ(40,*) k, (force_pnt(i,j),j=1,2), force(i)
  END DO
CLOSE(40)


END SUBROUTINE

!!!!!=======================================!!!!!
SUBROUTINE makeDmat (Dmat)

IMPLICIT NONE


INTEGER i, j, k
DOUBLE PRECISION Dmat(3,3)
DOUBLE PRECISION Young
DOUBLE PRECISION Poisson

!配列の初期化
DO i=1, 3
  DO j=1, 3
    Dmat(i,j) = 0.0d0
  END DO
END DO

OPEN(210,file='input_matinfo.txt',form='formatted')
  READ(210,*) Young
  READ(210,*) Poisson
CLOSE(210)

WRITE(*,*) 'YOUNG''S MODULUS [Pa] :', Young
WRITE(*,*) 'POISSON''S RATIO :', Poisson

!平面応力状態 P.121 式(5.53)
Dmat(1,1) = Young / (1.0d0 - (Poisson ** 2))
Dmat(1,2) = Young / (1.0d0 - (Poisson ** 2)) * Poisson
Dmat(2,1) = Young / (1.0d0 - (Poisson ** 2)) * Poisson
Dmat(2,2) = Young / (1.0d0 - (Poisson ** 2))
Dmat(3,3) = Young / (1.0d0 - (Poisson ** 2)) * (1.0d0 - Poisson) / 2.0d0

!平面ひずみ状態 P.122 式(5.54)
!Dmat(1,1) = Young * (1.0d0 - Poisson) / (1.0d0 - 2.0d0 * Poisson) / (1.0d0 + Poisson)
!Dmat(1,2) = Young / (1.0d0 - 2.0d0 * Poisson) / (1.0d0 + Poisson) * Poisson
!Dmat(2,1) = Young / (1.0d0 - 2.0d0 * Poisson) / (1.0d0 + Poisson) * Poisson
!Dmat(2,2) = Young * (1.0d0 - Poisson) / (1.0d0 - 2.0d0 * Poisson) / (1.0d0 + Poisson)
!Dmat(3,3) = Young / (1.0d0 + Poisson) / 2.0d0

OPEN(220,file='output_Dmat.dat',form='formatted')
  DO i=1, 3
    WRITE(220,'(3F20.5)') (Dmat(i,j),j=1,3)
  END DO
  WRITE(*,*) 'OUTPUT FILE: output_Dmat.dat'
CLOSE(220)

WRITE(*,*) 'MAKE D-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!

SUBROUTINE makeBmat (NUM_NODE, NUM_ELEME, node, eleme, Bmat, Ae)

IMPLICIT NONE


INTEGER i, j, k
INTEGER NUM_NODE, NUM_ELEME
DOUBLE PRECISION node(NUM_NODE,2)
INTEGER eleme(NUM_ELEME,3)
DOUBLE PRECISION Bmat(3,6,NUM_ELEME)
DOUBLE PRECISION Ae(NUM_ELEME)
DOUBLE PRECISION e_node(3,2)

!配列の初期化
DO i=1, 3
  DO j=1, 6
    DO k=1, NUM_ELEME
      Bmat(i,j,k) = 0.0d0
    END DO
  END DO
END DO

DO i=1, NUM_ELEME
  Ae(i) = 0.0d0
END DO

!各要素の面積を求める
DO i=1, NUM_ELEME
  DO j=1, 3
    e_node(j,1) = node(eleme(i,j),1)
    e_node(j,2) = node(eleme(i,j),2)
  END DO
  
  !P.102 式(5.19)
  Ae(i) = 0.50d0 * ((e_node(1,1) * (e_node(2,2) - e_node(3,2))) + (e_node(2,1) * (e_node(3,2) - e_node(1,2))) &
  &     + (e_node(3,1) * (e_node(1,2) - e_node(2,2))))
  
END DO

!各要素のB-matrixを求める
DO i=1, NUM_ELEME
  DO j=1, 3
    e_node(j,1) = node(eleme(i,j),1)
    e_node(j,2) = node(eleme(i,j),2)
  END DO
  
  !P.129 式(5.77)
  Bmat(1,1,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(2,2) - e_node(3,2))
  Bmat(1,3,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(3,2) - e_node(1,2))
  Bmat(1,5,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(1,2) - e_node(2,2))
  
  Bmat(2,2,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(3,1) - e_node(2,1))
  Bmat(2,4,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(1,1) - e_node(3,1))
  Bmat(2,6,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(2,1) - e_node(1,1))
  
  Bmat(3,1,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(3,1) - e_node(2,1))
  Bmat(3,2,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(2,2) - e_node(3,2))
  Bmat(3,3,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(1,1) - e_node(3,1))
  Bmat(3,4,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(3,2) - e_node(1,2))
  Bmat(3,5,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(2,1) - e_node(1,1))
  Bmat(3,6,i) = (1.0d0 / (2.0d0 * Ae(i))) * (e_node(1,2) - e_node(2,2))
END DO

OPEN(310,file='output_Bmat.dat',form='formatted')
  DO i=1, NUM_ELEME
    WRITE(310,'(A9,I5)') 'ELEMENT :', i
    DO j=1, 3
      WRITE(310,'(6F10.5)') (Bmat(j,k,i),k=1,6)
    END DO
  END DO
  WRITE(*,*) 'OUTPUT FILE: output_Bmat.dat'
CLOSE(310)

OPEN(320,file='output_Ae.dat',form='formatted')
  DO i=1, NUM_ELEME
    WRITE(320,'(A9,F10.5)') 'ELEMENT :', Ae(i)
  END DO
  WRITE(*,*) 'OUTPUT FILE: output_Ae.dat'
CLOSE(320)

WRITE(*,*) 'MAKE B-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!
SUBROUTINE makeKmat (NUM_NODE, NUM_ELEME, eleme, Bmat, Dmat, Kmat, Ae, THICKNESS)

IMPLICIT NONE


INTEGER i, j, k, l, pt1, pt2
INTEGER NUM_NODE, NUM_ELEME
DOUBLE PRECISION THICKNESS
INTEGER eleme(NUM_ELEME,3)
DOUBLE PRECISION Bmat(3,6,NUM_ELEME)
DOUBLE PRECISION Dmat(3,3)
DOUBLE PRECISION Kmat(2*NUM_NODE,2*NUM_NODE) !全体剛性マトリックス
DOUBLE PRECISION Ae(NUM_ELEME)
DOUBLE PRECISION e_Kmat(6,6) !要素剛性マトリックス
DOUBLE PRECISION BTD(6,3)
DOUBLE PRECISION BTDB(6,6)

!配列の初期化
DO i=1, 2*NUM_NODE
  DO j=1, 2*NUM_NODE
    Kmat(i,j) = 0.0d0
  END DO
END DO

DO i=1, NUM_ELEME
  
  DO j=1, 6
    DO k=1, 6
      e_Kmat(j,k) = 0.0d0
      BTDB(j,k) = 0.0d0
    END DO
    DO k=1, 3
      BTD(j,k) = 0.0d0
    END DO
  END DO
  
  !要素剛性マトリックスの構築 P.135 式(5.94)
  DO j=1, 6
    DO k=1, 3
      DO l=1, 3
        BTD(j,k) = BTD(j,k) + Bmat(l,j,i) * Dmat(k,l)
      END DO
    END DO
  END DO
  
  DO j=1, 6
    DO k=1, 6
      DO l=1, 3
        BTDB(j,k) = BTDB(j,k) + BTD(j,l) * Bmat(l,k,i)
      END DO
    END DO
  END DO
  
  DO j=1, 6
    DO k=1, 6
      e_Kmat(j,k) = Ae(i) * THICKNESS * BTDB(j,k)
    END DO
  END DO
  
  !全体剛性マトリックスへの組込み P.137 式(5.97)
  DO j=1, 3
    DO k=1, 3
      pt1 = eleme(i,j)
      pt2 = eleme(i,k)
      
      ![2x2]の成分ごとに組込み
      Kmat(2*(pt1-1)+1,2*(pt2-1)+1) = Kmat(2*(pt1-1)+1,2*(pt2-1)+1) + e_Kmat(2*(j-1)+1,2*(k-1)+1)
      Kmat(2*(pt1-1)+1,2*(pt2-1)+2) = Kmat(2*(pt1-1)+1,2*(pt2-1)+2) + e_Kmat(2*(j-1)+1,2*(k-1)+2)
      Kmat(2*(pt1-1)+2,2*(pt2-1)+1) = Kmat(2*(pt1-1)+2,2*(pt2-1)+1) + e_Kmat(2*(j-1)+2,2*(k-1)+1)
      Kmat(2*(pt1-1)+2,2*(pt2-1)+2) = Kmat(2*(pt1-1)+2,2*(pt2-1)+2) + e_Kmat(2*(j-1)+2,2*(k-1)+2)
    END DO
  END DO
  
END DO

WRITE(*,*) 'MAKE K-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!
SUBROUTINE makeFmat (NUM_NODE, NUM_FORCE, Fmat, force_pnt, force)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_FORCE
INTEGER i, j, k
INTEGER force_pnt(NUM_FORCE,2)
DOUBLE PRECISION force(NUM_FORCE)
DOUBLE PRECISION Fmat(2*NUM_NODE)

DO i=1, 2*NUM_NODE
  Fmat(i) = 0.0d0
END DO

DO i=1, NUM_FORCE
  Fmat(2*(force_pnt(i,1)-1)+force_pnt(i,2)) = force(i)
  
  IF(force_pnt(i,2) /= 1 .AND. force_pnt(i,2) /= 2)THEN
    WRITE(*,*) 'INPUT DATA "input_forcednodes.txt" IS NOT APPROPREATE.'
    STOP
  END IF
END DO

!OPEN(410,file='output_Fmat.dat',form='formatted')
!  DO i=1, 2*NUM_NODE
!    WRITE(410,'(I5,F15.5)') i, Fmat(i)
!  END DO
!  WRITE(*,*) 'OUTPUT FILE: output_Fmat.dat'
!CLOSE(410)

WRITE(*,*) 'MAKE F-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!
SUBROUTINE makeUmat (NUM_NODE, NUM_FIX, Umat, fix_pnt, fix)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_FIX
INTEGER i, j, k
INTEGER fix_pnt(NUM_FIX,2)
DOUBLE PRECISION fix(NUM_FIX)
DOUBLE PRECISION Umat(2*NUM_NODE)

DO i=1, 2*NUM_NODE
  Umat(i) = 0.0d0
END DO

DO i=1, NUM_FIX
  Umat(2*(fix_pnt(i,1)-1)+fix_pnt(i,2)) = fix(i)
  
  IF(fix_pnt(i,2) /= 1 .AND. fix_pnt(i,2) /= 2)THEN
    WRITE(*,*) 'INPUT DATA "input_fixednodes.txt" IS NOT APPROPREATE.'
    STOP
  END IF
END DO

!OPEN(420,file='output_Umat.dat',form='formatted')
!  DO i=1, 2*NUM_NODE
!    WRITE(420,'(I5,F15.5)') i, Umat(i)
!  END DO
!  WRITE(*,*) 'OUTPUT FILE: output_Umat.dat'
!CLOSE(420)

WRITE(*,*) 'MAKE U-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!
SUBROUTINE makeSUBmatrix (NUM_NODE, NUM_FIX, Kmat, Fmat, Umat, fix_pnt, known_DOF, unknown_DOF, K11, K12, K22, &
&                          F1, F2, U1, U2)

IMPLICIT NONE

INTEGER NUM_NODE, NUM_FIX
INTEGER i, j, k
DOUBLE PRECISION Kmat(2*NUM_NODE,2*NUM_NODE)
DOUBLE PRECISION Fmat(2*NUM_NODE)
DOUBLE PRECISION Umat(2*NUM_NODE)
INTEGER fix_pnt(NUM_FIX,2)
INTEGER known_DOF(NUM_FIX)
INTEGER unknown_DOF(2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION K11(2*NUM_NODE-NUM_FIX,2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION K12(2*NUM_NODE-NUM_FIX,NUM_FIX)
DOUBLE PRECISION K22(NUM_FIX,NUM_FIX)
DOUBLE PRECISION F1(2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION F2(NUM_FIX)
DOUBLE PRECISION U1(2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION U2(NUM_FIX)

!境界条件適用後の小行列を作成
DO i=1, NUM_FIX
  F2(i) = 0.0d0
END DO

DO i=1, 2*NUM_NODE-NUM_FIX
  U1(i) = 0.0d0
END DO

DO i=1, NUM_FIX
  known_DOF(i) = 2*(fix_pnt(i,1)-1) + fix_pnt(i,2)
END DO

DO j=1, known_DOF(1)-1
  unknown_DOF(j) = j
END DO
DO i=2, NUM_FIX
  DO j=known_DOF(i-1)+1, known_DOF(i)-1
    unknown_DOF(j-(i-1)) = j
  END DO
END DO
DO j=known_DOF(NUM_FIX)+1, 2*NUM_NODE
  unknown_DOF(j-NUM_FIX) = j
END DO

DO i=1, 2*NUM_NODE-NUM_FIX
  DO j=1, 2*NUM_NODE-NUM_FIX
    K11(i,j) = Kmat(unknown_DOF(i),unknown_DOF(j))
  END DO
  
  DO j=1, NUM_FIX
    K12(i,j) = Kmat(unknown_DOF(i),known_DOF(j))
  END DO
  F1(i) = Fmat(unknown_DOF(i))
  U1(i) = Umat(unknown_DOF(i)) !未知成分
END DO

DO i=1, NUM_FIX
  DO j=1, NUM_FIX
    K22(i,j) = Kmat(known_DOF(i),known_DOF(j))
  END DO
  F2(i) = Fmat(known_DOF(i)) !未知成分
  U2(i) = Umat(known_DOF(i))
END DO

WRITE(*,*) 'MAKE SUB-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!
SUBROUTINE makeInverse (NUM_NODE, NUM_FIX, K11)

IMPLICIT NONE

INTEGER NUM_NODE, NUM_FIX
INTEGER i, j, k
DOUBLE PRECISION K11(2*NUM_NODE-NUM_FIX,2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION inv(2*NUM_NODE-NUM_FIX,2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION pivot(2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION num

!（部分）ピボット交換
DO i=1, 2*NUM_NODE-NUM_FIX
  pivot(i) = i
END DO

DO i=1, 2*NUM_NODE-NUM_FIX
  num = K11(i,i)
  k = 0
  DO j=1, 2*NUM_NODE-NUM_FIX
    IF(num < K11(j,i))THEN
      num = K11(j,i)
      k = j
    END IF
  END DO
  
  IF(k /= 0)THEN
    pivot(i) = k
    pivot(k) = i
    
    DO j=1, 2*NUM_NODE-NUM_FIX
      num = K11(i,j)
      K11(i,j) = K11(k,j)
      K11(k,j) = num
    END DO
  END IF
END DO

!K11を単位行列化することで逆行列を計算
!数値計算的には，Gaussの消去法やLU分解の方が有用…
DO i=1, 2*NUM_NODE-NUM_FIX
  DO j=1, 2*NUM_NODE-NUM_FIX
    inv(i,j) = 0.0d0
  END DO
  inv(i,i) = 1.0d0
END DO

DO i=1, 2*NUM_NODE-NUM_FIX
  num = K11(i,i)
  DO j=1, 2*NUM_NODE-NUM_FIX
    K11(i,j) = K11(i,j) / num
    inv(i,j) = inv(i,j) / num
  END DO
  
  DO j=1, i-1
    num = K11(j,i)
    DO k=1, 2*NUM_NODE-NUM_FIX
      K11(j,k) = K11(j,k) - K11(i,k) * num
      inv(j,k) = inv(j,k) - inv(i,k) * num
    END DO
  END DO
  DO j=i+1, 2*NUM_NODE-NUM_FIX
    num = K11(j,i)
    DO k=1, 2*NUM_NODE-NUM_FIX
      K11(j,k) = K11(j,k) - K11(i,k) * num
      inv(j,k) = inv(j,k) - inv(i,k) * num
    END DO
  END DO
  
END DO

!K11の逆行列を格納
DO i=1, 2*NUM_NODE-NUM_FIX
  DO j=1, 2*NUM_NODE-NUM_FIX
    K11(i,j) = inv(i,j)
  END DO
END DO

!ピボット交換を元に戻す
DO i=1, 2*NUM_NODE-NUM_FIX
  IF(pivot(i) /= i)THEN
    DO j=1, 2*NUM_NODE-NUM_FIX
      num = K11(i,j)
      K11(i,j) = K11(k,j)
      K11(k,j) = num
    END DO
  END IF
END DO

WRITE(*,*) 'MAKE K11-INV-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!
SUBROUTINE solveUmat (NUM_NODE, NUM_FIX, Umat, K11, K12, F1, U1, U2, unknown_DOF)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_FIX
INTEGER i, j, k
DOUBLE PRECISION Umat(2*NUM_NODE)
DOUBLE PRECISION K11(2*NUM_NODE-NUM_FIX,2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION K12(2*NUM_NODE-NUM_FIX,NUM_FIX)
DOUBLE PRECISION F1(2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION U1(2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION U2(NUM_FIX)
DOUBLE PRECISION FKU(2*NUM_NODE-NUM_FIX)
INTEGER unknown_DOF(2*NUM_NODE-NUM_FIX)

DO i=1, 2*NUM_NODE-NUM_FIX
  U1(i) = 0.0d0
  FKU(i) = 0.0d0
END DO

!P.139 式(5.104)
DO i=1, 2*NUM_NODE-NUM_FIX
  DO j=1, NUM_FIX
    FKU(i) = FKU(i) + K12(i,j) * U2(j)
  END DO
  FKU(i) = F1(i) - FKU(i)
END DO

DO i=1, 2*NUM_NODE-NUM_FIX
  DO j=1, 2*NUM_NODE-NUM_FIX
    U1(i) = U1(i) + K11(i,j) * FKU(j)
  END DO
END DO

DO i=1, 2*NUM_NODE-NUM_FIX
  Umat(unknown_DOF(i)) = U1(i)
END DO

OPEN(510,file='output_Umat.dat',form='formatted')
  DO i=1, 2*NUM_NODE
    WRITE(510,'(I5,F15.5)') i, Umat(i)
  END DO
  WRITE(*,*) 'OUTPUT FILE: output_Umat.dat'
CLOSE(510)

WRITE(*,*) 'SOLVE U-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!
SUBROUTINE solveFmat(NUM_NODE, NUM_FIX, Fmat, K12, K22, F2, U1, U2, known_DOF)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_FIX
INTEGER i, j, k
DOUBLE PRECISION Fmat(2*NUM_NODE)
DOUBLE PRECISION K12(2*NUM_NODE-NUM_FIX,NUM_FIX)
DOUBLE PRECISION K22(NUM_FIX,NUM_FIX)
DOUBLE PRECISION F2(NUM_FIX)
DOUBLE PRECISION U1(2*NUM_NODE-NUM_FIX)
DOUBLE PRECISION U2(NUM_FIX)
INTEGER known_DOF(NUM_FIX)

!P.140
DO i=1, NUM_FIX
  DO j=1, 2*NUM_NODE-NUM_FIX
    F2(i) = F2(i) + K12(j,i) * U1(j)
  END DO
  
  DO j=1, NUM_FIX
    F2(i) = F2(i) + K22(i,j) * U2(j)
  END DO
END DO

DO i=1, NUM_FIX
  Fmat(known_DOF(i)) = F2(i)
END DO

OPEN(520,file='output_Fmat.dat',form='formatted')
  DO i=1, 2*NUM_NODE
    WRITE(520,'(I5,F15.5)') i, Fmat(i)
  END DO
  WRITE(*,*) 'OUTPUT FILE: output_Fmat.dat'
CLOSE(520)

WRITE(*,*) 'SOLVE F-MATRIX'


END SUBROUTINE

!!!!!=======================================!!!!!
SUBROUTINE displacement (NUM_NODE, AMP, node, Umat)

IMPLICIT NONE


INTEGER NUM_NODE
DOUBLE PRECISION AMP
INTEGER i, j, k
DOUBLE PRECISION node(NUM_NODE,2)
DOUBLE PRECISION disp(NUM_NODE,2)
DOUBLE PRECISION Umat(2*NUM_NODE)

DO i=1, NUM_NODE
  DO j=1, 2
    disp(i,j) = 0.0d0
  END DO
END DO

DO i=1, NUM_NODE
  DO j=1, 2
    disp(i,j) = node(i,j) + Umat(2*(i-1)+j) * AMP
  END DO
END DO

OPEN(610,file='output_disp.dat',form='formatted')
  DO i=1, NUM_NODE
    WRITE(610,'(I5,2F10.5)') i, (disp(i,j),j=1,2)
  END DO
  WRITE(*,*) 'OUTPUT FILE: output_disp.dat'
CLOSE(610)

WRITE(*,*) 'CALCULATE DISPLACEMENT'


END SUBROUTINE

!!!!!=======================================!!!!!
SUBROUTINE distribution (NUM_NODE, NUM_ELEME, eleme, Bmat, Dmat, Umat)

IMPLICIT NONE


INTEGER NUM_NODE, NUM_ELEME
INTEGER i, j, k, l
INTEGER eleme(NUM_ELEME,3)
DOUBLE PRECISION Bmat(3,6,NUM_ELEME)
DOUBLE PRECISION Dmat(3,3)
DOUBLE PRECISION Umat(2*NUM_NODE)
DOUBLE PRECISION e_Umat(6)
DOUBLE PRECISION strain(3,NUM_ELEME), stress(3,NUM_ELEME)

DO i=1, 3
  DO j=1, NUM_ELEME
    strain(i,j) = 0.0d0
    stress(i,j) = 0.0d0
  END DO
END DO

DO i=1, NUM_ELEME
  DO j=1, 3
    DO k=1, 2
      e_Umat(2*(j-1)+k) = Umat(2*(eleme(i,j)-1)+k)
    END DO
  END DO
  
  !P.140 式(5.105)
  DO j=1, 3
    DO k=1, 6
      strain(j,i) = strain(j,i) + Bmat(j,k,i) * e_Umat(k)
    END DO
  END DO
  
  !P.141
  DO j=1, 3
    DO k=1, 3
      stress(j,i) = stress(j,i) + Dmat(j,k) * strain(k,i)
    END DO
  END DO
END DO

OPEN(710,file='output_strain.dat',form='formatted')
  DO i=1, NUM_ELEME
    WRITE(710,'(I5,3F20.5)') i, (strain(j,i),j=1,3)
  END DO
  WRITE(*,*) 'OUTPUT FILE: output_strain.dat'
CLOSE(710)

OPEN(720,file='output_stress.dat',form='formatted')
  DO i=1, NUM_ELEME
    WRITE(720,'(I5,3F20.5)') i, (stress(j,i),j=1,3)
  END DO
  WRITE(*,*) 'OUTPUT FILE: output_stress.dat'
CLOSE(720)

WRITE(*,*) 'CALCULATE DISTRIBUTIONS'


END SUBROUTINE

