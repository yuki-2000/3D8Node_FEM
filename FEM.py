# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 09:47:20 2022

@author: yuki
"""

#FEM




"""

このプログラムは、後藤先生のfortranのプログラムをpythonで書き直したものです。

目標
numpyで行列計算を高速に
0での初期化を一瞬で
行列の内積をforで回さずに一瞬で
変数名は元のプとグラムを踏襲
変数型はこちらで指定



#今後
コメント多め
GPU
マルチスレッド
numba
疎行列 
class
わざわざ入れ替えずに直接行列式を解く
メッシュの可視化


参考
https://qiita.com/damyarou/items/8ca3432f9bba20a1f5ea
http://nkl.cc.u-tokyo.ac.jp/12w/CW-intro/CW-intro01.pdf
https://ipsj.ixsq.nii.ac.jp/ej/?action=repository_action_common_download&item_id=75552&item_no=1&attribute_id=1&file_no=1




#その他
!でのコメントアウトを#に変更、データ読み込み時も。
最初に配列を確保しなくていい
配列でデフォで行優先なのでfortranとメモリ確保の仕方が違う

fortranは配列は1始まりだが、pythonは0始まり。さらに[1:3]はfortranは[1,2,3]、pythonは[1,2]

floatで定義した配列にはstrを代入しても数値に代わるが、ただの変数はstrのままなので変換が必要

node番号!=配列番号

疑問点：境界以外の内部については、すべて荷重0かもしくはあたえられているという認識でいいのか？

"""


import numpy as np
from matplotlib import pyplot as plt
import time
import sys
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import inv, spsolve
from scipy.linalg import solve, det
#処理時間計測
start_time = time.time()
lap_time = time.time()



#デフォでmode=`r`
#pythonではすべて大文字の変数は定数扱いなので小文字に変換
#https://qiita.com/naomi7325/items/4eb1d2a40277361e898b    
#fortranは大文字小文字の区別なし

#withを使うと自動でcloseされる
#入力データに!を使ったコメントアウトがあるので分割して取得後、空白含むstrをintに変換
#pythonのfloatは64bit倍精度

#fortranでは単精度では1.23e4、倍精度では1.23d4とかくが、pythonはeのみ対応。よって置換
#https://docs.python.org/ja/3/library/functions.html#float

with open('input_AnalysisConditions.txt') as f:
#with open('benchmark_input_AnalysisConditions.txt') as f:
    l = f.readlines()
    num_node  = int(l[0].split('!')[0]) #モデル節点数
    num_eleme = int(l[1].split('!')[0]) #モデル要素数
    num_material = int(l[2].split('!')[0])#材料種類数
    num_fix   = int(l[3].split('!')[0]) #拘束点数
    num_force = int(l[4].split('!')[0]) #荷重点数
    amp       = float(l[5].split('!')[0].replace('d','e')) #変形図倍率














#input (NUM_NODE, NUM_ELEME, NUM_FIX, NUM_FORCE, node, eleme, fix_pnt, force_pnt, fix, force)




#numpy.zeros(shape, dtype=float, order='C')
#order : {‘C’, ‘F’},多次元配列を行優先 (C-style) か列優先 (Fortran-style) でメモリに格納するかを指定します。  デフォでC  
#dtypedata-type, optionalDesired output data-type for the array, e.g, numpy.int8. Default is numpy.float64.
#https://numpy.org/doc/stable/reference/generated/numpy.zeros.html


#node = np.zeros((num_node,2),dtype=np.float64)
#emptyのほうがより高速だが、初期化されない
node      = np.empty((num_node,3), dtype=np.float64) #節点座標
eleme     = np.empty((num_eleme,8),dtype=np.int32) #各要素のコネクティビティ #つまりある四角形elementを構成する接点node番号(1スタートに注意)
material  = np.empty((num_eleme),dtype=np.int32) #各要素の素材番号 (元と違うので注意)
fix_pnt   = np.empty((num_fix,2),  dtype=np.int32) #変位境界条件
force_pnt = np.empty((num_force,2),dtype=np.int32) #力学的境界条件 #接点番号と向きの配列
fix       = np.empty((num_fix),    dtype=np.float64) #変位境界条件の値
force     = np.empty((num_force),  dtype=np.float64) #力学的境界条件の値


#dを使った指数表現に対応、スペース区切りに変更
with open('input_point.txt') as f:
#with open('benchmark_input_point.txt') as f:
    l = f.readlines()
    for i, input_point in enumerate(l):
        node[i,0] = input_point.split()[1].replace('d','e')
        node[i,1] = input_point.split()[2].replace('d','e')
        node[i,2] = input_point.split()[3].replace('d','e')
        
#スペース区切りに変更
with open('input_eleme.txt') as f:
#with open('benchmark_input_eleme.txt') as f:
    l = f.readlines()
    for i, input_eleme in enumerate(l):
        eleme[i] = input_eleme.split()[1:9]

#追加   
with open('input_material.txt') as f:
    l = f.readlines()
    for i, input_material in enumerate(l):
        material[i] = input_material.split()[1]

      
#行の最後に文章があるので行番号を厳密に指定        
with open('input_fixednodes.txt') as f:
    l = f.readlines()
    for i in range(num_fix):
        fix_pnt[i] = l[i].split()[1:3]
        fix[i] = l[i].split()[3].replace('d','e')
        

#行の最後に文章があるので行番号を厳密に指定        
with open('input_forcednodes.txt') as f:
    l = f.readlines()
    for i in range(num_force):
        force_pnt[i] = l[i].split()[1:3]
        force[i] = l[i].split()[3].replace('d','e')



print("Finish reading input text")

print("経過時間:", time.time() - start_time)
print("処理時間:", time.time() - lap_time)
lap_time = time.time()














#makeDmat (Dmat)

#配列の初期化
Dmat = np.zeros((6,6,num_material),dtype=np.float64) #弾性剛性マトリックス

#最悪　めちゃくちゃわかりにくい
#ちゃんとしたデータセット形式のinput作れ

#readで読んでいるから
line_num = 0
  
with open('input_matinfo.txt') as f:
    l = f.readlines()
    
      
        
        
    for k in range(num_material):
        num = int(l[line_num].split('!')[0]) #MATERIAL ID  変数名変えるべき
        line_num += 1
        print('MATERIAL ID :', num)
        #MATERIAL ID は1始まりだが、kは0始まりなので合わせている。
        if num != k+1:
            print('MATERIAL ID IS NOT APPROPRIATE.')
            
        
        mat_type = int(l[line_num].split('!')[0]) #Isotropic Material: 1, Transversely Isotropic Material: 2
        line_num += 1
        
        if mat_type == 1: #Isotropic Material 等方材料
            print('MATERIAL TYPE IS ISOTROPIC MATERIAL.')
            
            Young1 = float(l[line_num].split('!')[0].replace('d','e'))
            line_num += 1
            Poisson1 = float(l[line_num].split('!')[0].replace('d','e'))
            line_num += 1

            print('YOUNG''S MODULUS [MPa] :', Young1)
            print('POISSON''S RATIO :', Poisson1)  
            
            Young1 *= 10**6 #MPaからPaに変更
            
            
            #配列がpythonでは0始まりなので[0,1]は1行2列、fortranでは[1,2]とかく。
            Dmat[0,0,k] = Young1 * (1 - Poisson1) / ( 1 + Poisson1) / (1 - 2 * Poisson1)
            Dmat[1,1,k] = Young1 * (1 - Poisson1) / ( 1 + Poisson1) / (1 - 2 * Poisson1)
            Dmat[2,2,k] = Young1 * (1 - Poisson1) / ( 1 + Poisson1) / (1 - 2 * Poisson1)
            Dmat[0,1,k] = Young1 * Poisson1 / ( 1 + Poisson1) / (1 - 2 * Poisson1)
            Dmat[1,0,k] = Young1 * Poisson1 / ( 1 + Poisson1) / (1 - 2 * Poisson1)
            Dmat[1,2,k] = Young1 * Poisson1 / ( 1 + Poisson1) / (1 - 2 * Poisson1)
            Dmat[2,1,k] = Young1 * Poisson1 / ( 1 + Poisson1) / (1 - 2 * Poisson1)
            Dmat[2,0,k] = Young1 * Poisson1 / ( 1 + Poisson1) / (1 - 2 * Poisson1)
            Dmat[0,2,k] = Young1 * Poisson1 / ( 1 + Poisson1) / (1 - 2 * Poisson1)
            Dmat[3,3,k] = Young1 / 2 / (1 + Poisson1)
            Dmat[4,4,k] = Young1 / 2 / (1 + Poisson1)
            Dmat[5,5,k] = Young1 / 2 / (1 + Poisson1)
        
        
        elif mat_type == 2: #Transeversely Isotropic Material　直交異方性材料、横等方性材料
            print('MATERIAL TYPE IS TRANSEVERSELY ISOTROPIC MATERIAL.')
            print('TRANSEVERSELY ISOTROPIC MATERIAL IS NOT APPLIED YET.')
   
                





#----------------------------------
#ouptputを今は書いていない


print('MAKE D-MATRIX')

print("経過時間:", time.time() - start_time)
print("処理時間:", time.time() - lap_time)
lap_time = time.time()

























# makeBmat (NUM_NODE, NUM_ELEME, node, eleme, Bmat, Ae)


#配列の初期化
#ξ:xi η:eta ζ:zeta
Bmat      = np.zeros((6,24,num_eleme,8), dtype=np.float64) #Bマトリックス（形状関数の偏微分） #ただしすでにガウス積分点代入済み 8はガウスの積分点の個数
Hmat      = np.zeros((3,8), dtype=np.float64) #dN/d(xi),dN/d(eta),dN/d(zeta)を成分に持つ行列（xi,eta,zetaにはガウスの積分点を代入）　あるガウス積分点における値
det_Jacobi       = np.zeros((num_eleme,8), dtype=np.float64) #ガウスの積分点におけるヤコビアン（ヤコビ行列の行列式)
gauss_nodes     = np.zeros((8,3), dtype=np.float64) #ガウスの積分点(polarと同じように左下から反時計)
polar     = np.zeros((8,3), dtype=np.float64) #要素座標(xi,eta,zeta)における節点座標　左下から反時計
Jacobi    = np.zeros((3,3), dtype=np.float64) #ヤコビ行列 (ガウスの積分点代入済み)
Jacobiinv = np.zeros((3,3), dtype=np.float64) #ヤコビ行列の逆行列
dNdxy     = np.zeros((3,8), dtype=np.float64) #dN/dx,dN/dy,dN/dzを成分に持つ行列 行がxyz列が形状関数N Hmat、Jacobiinvともにあるガウスの積分点代入済みのため、これもガウスの積分点代入済み
e_node    = np.zeros((8,3), dtype=np.float64) #ある四角形elementを構成する4接点のxy座標　#e_pointから戻した。


#ガウスの積分点 左下から
gauss_nodes = np.array([[-1/np.sqrt(3), -1/np.sqrt(3), -1/np.sqrt(3)],
                        [ 1/np.sqrt(3), -1/np.sqrt(3), -1/np.sqrt(3)],
                        [ 1/np.sqrt(3),  1/np.sqrt(3), -1/np.sqrt(3)],
                        [-1/np.sqrt(3),  1/np.sqrt(3), -1/np.sqrt(3)],
                        [-1/np.sqrt(3), -1/np.sqrt(3),  1/np.sqrt(3)],
                        [ 1/np.sqrt(3), -1/np.sqrt(3),  1/np.sqrt(3)],
                        [ 1/np.sqrt(3),  1/np.sqrt(3),  1/np.sqrt(3)],
                        [-1/np.sqrt(3),  1/np.sqrt(3),  1/np.sqrt(3)]])

#自然座標での４点、左下から
polar = np.array([[-1, -1, -1],
                  [ 1, -1, -1],
                  [ 1,  1, -1],
                  [-1,  1, -1],
                  [-1, -1,  1],
                  [ 1, -1,  1],
                  [ 1,  1,  1],
                  [-1,  1,  1]])





#各要素の面積を求める

#配列0始まりに変更
#eleme[i,j]は接点番号であり、pythonにおける配列位置にするためには-1する必要あり
#enodeは要素を構成する接点の座標

#各要素のB-matrixを求める
#配列0始まりに変更
#eleme[i,j]は接点番号であり、pythonにおける配列位置にするためには-1する必要あり
for i in range(num_eleme):
    for j in range(8): #節点の8
        e_node[j,0] = node[eleme[i,j]-1,0]
        e_node[j,1] = node[eleme[i,j]-1,1]
        e_node[j,2] = node[eleme[i,j]-1,2]
    
   
    for j in range(len(gauss_nodes)): #各ガウスの積分点を代入した時
        for k in range(8): #各接点の8
            #pythonは0スタート
            #Nは(ξ,η)て定義、4節点に対応するNの偏微分に、ある積分点を代入
            Hmat[0,k] = polar[k,0] * (1 + polar[k,1] * gauss_nodes[j,1]) * (1 + polar[k,2] * gauss_nodes[j,2]) /8
            Hmat[1,k] = polar[k,1] * (1 + polar[k,0] * gauss_nodes[j,0]) * (1 + polar[k,2] * gauss_nodes[j,2]) /8
            Hmat[2,k] = polar[k,2] * (1 + polar[k,0] * gauss_nodes[j,0]) * (1 + polar[k,1] * gauss_nodes[j,1]) /8
        
        

        
        #p220-221 ただしガウスの積分点代入済み
        Jacobi = Hmat @ e_node
        
        
        #p220(B.25)
        det_Jacobi[i,j] = np.linalg.det(Jacobi)
        if det_Jacobi[i,j] <= 0:
            print("error det_Jacobi<=0")
            print("element:", i)
            print("gauss:", j)
            
        
       
        Jacobiinv = np.linalg.inv(Jacobi)
        

        #p222らへん            
        dNdxy = Jacobiinv @ Hmat
        
        #fortranは1始まり、pythonは0始まり
        #p222
        #for k in range(4):
            #Bmat[0,2*k,i,j] = dNdxy[0,k]
            #Bmat[1,2*k+1,i,j] = dNdxy[1,k]
            #Bmat[2,2*k,i,j] = dNdxy[1,k]
            #Bmat[2,2*k+1,i,j] = dNdxy[0,k]
            
        Bmat[0,::3, i,j] = dNdxy[0,:]
        Bmat[1,1::3,i,j] = dNdxy[1,:]
        Bmat[2,2::3,i,j] = dNdxy[2,:]
        Bmat[3,1::3,i,j] = dNdxy[2,:]
        Bmat[3,2::3,i,j] = dNdxy[1,:]
        Bmat[4,::3, i,j] = dNdxy[2,:]
        Bmat[4,2::3,i,j] = dNdxy[0,:]
        Bmat[5,::3, i,j] = dNdxy[1,:]
        Bmat[5,1::3,i,j] = dNdxy[0,:]
        
        



#----------------------------------
#ouptputを今は書いていない

print('MAKE B-MATRIX')

print("経過時間:", time.time() - start_time)
print("処理時間:", time.time() - lap_time)
lap_time = time.time()










# makeKmat (NUM_NODE, NUM_ELEME, eleme, Bmat, Dmat, Kmat, Ae, THICKNESS)




#配列の初期化
Kmat   = np.zeros((3*num_node,3*num_node), dtype=np.float64) #全体剛性マトリックス


for i in range(num_eleme):
    
    e_Kmat = np.zeros((24,24), dtype=np.float64)  #要素剛性マトリックス

    for j in range(len(gauss_nodes)): #各ガウスの積分点を代入した時
    
        #要素剛性マトリックスの構築 P.135 式(5.94)
        #一発で、メモリのほんのちょっとの節約
        #material[i]は、i要素の素材番号1始まりだが、Dmatの格納場所は0なので注意
        #ガウスの積分点の個数だけ足されていることに注意
        e_Kmat += det_Jacobi[i,j] * Bmat[:,:,i,j].T @ Dmat[:,:,material[i]-1] @ Bmat[:,:,i,j]
    
    
    
    
    
    #全体剛性マトリックスへの組込み P.137 式(5.97)

    #ここもっと行列計算したい
    for j in range(8):
        for k in range(8):

            #eleme[i,j]は接点番号であり、pythonにおける配列位置にするためには-1する必要があると思ったが、
            #Kmatの式の作り方からやめておく
            pt1 = eleme[i,j] #-1 #行
            pt2 = eleme[i,k] #-1 #列
            

            #[3x3]の成分ごとに組込み
            #1行でできる
            Kmat[3*(pt1-1):3*(pt1-1)+3, 3*(pt2-1):3*(pt2-1)+3] += e_Kmat[3*j:3*j+3, 3*k:3*k+3]

#疎行列に変換
#Kmat = coo_matrix(Kmat).tolil()
Kmat = coo_matrix(Kmat).tocsr()
#Kmat = coo_matrix(Kmat).tocsc()


print( 'MAKE K-MATRIX')

print("経過時間:", time.time() - start_time)
print("処理時間:", time.time() - lap_time)
lap_time = time.time()







# makeFmat (NUM_NODE, NUM_FORCE, Fmat, force_pnt, force)

Fmat = np.zeros((3*num_node), dtype=np.float64) #節点荷重ベクトル

#unknown_DOFをつかってファンシーインデックスにしたほうが早い
for i in range(num_force):
    #force_pnt[i,1]は接点番号であり、pythonにおける配列位置にするために変更、
    #各接点のx,yの順に配列が並んでいるので、xは+1、yは+2が割り振られうまく位置を計算している。
    #pythonの配列番号0始まりに変更
    Fmat[3*(force_pnt[i,0]-1) + force_pnt[i,1] -1] = force[i]
    
    if (force_pnt[i,1] != 1 and force_pnt[i,1] != 2 and force_pnt[i,1] != 3):
        print('INPUT DATA "input_forcednodes.txt" IS NOT APPROPREATE.')
        print('load direction is now',force_pnt[i,1], ', not 1(x) or 2(y) or 3(z)' )
        break
 


#output省略













# makeUmat (NUM_NODE, NUM_FIX, Umat, fix_pnt, fix)



Umat = np.zeros((3*num_node), dtype=np.float64)


#known_DOFをつかってファンシーインデックスにしたほうが早い
for i in range(num_fix):
    #fix_pnt[i,1]は接点番号であり、pythonにおける配列位置にするために変更、
    #各接点のx,yの順に配列が並んでいるので、xは+1、yは+2が割り振られうまく位置を計算している。
    #pythonの配列番号0始まりに変更
    Umat[3*(fix_pnt[i,0]-1) + fix_pnt[i,1] -1] = fix[i]
    
    if (fix_pnt[i,1] != 1 and fix_pnt[i,1] != 2 and force_pnt[i,1] != 3):
        print('IINPUT DATA "input_fixednodes.txt" IS NOT APPROPREATE.')
        print('Fixed direction is now', fix_pnt[i,1], ', not 1(x) or 2(y) or 3(z)' )
        break


#output省略、というか解けていないから不必要











# makeSUBmatrix (NUM_NODE, NUM_FIX, Kmat, Fmat, Umat, fix_pnt, known_DOF, unknown_DOF, K11, K12, K22, F1, F2, U1, U2)



#境界条件適用後の小行列を作成
known_DOF   = np.empty(num_fix, dtype=np.int32)              #既知節点変位ベクトルの自由度  #既知接点変位の行番号であり、未知荷重行に対応
unknown_DOF = np.empty(3*num_node - num_fix, dtype=np.int32) #未知節点変位ベクトルの自由度

K11 = np.zeros((3*num_node-num_fix, 3*num_node-num_fix), dtype=np.float64) #変位境界条件付加後の小行列
K12 = np.zeros((3*num_node-num_fix, num_fix), dtype=np.float64)            #変位境界条件付加後の小行列 #K21の転置
K22 = np.zeros((num_fix, num_fix), dtype=np.float64)                       #変位境界条件付加後の小行列


F1  = np.zeros((3*num_node-num_fix), dtype=np.float64)                     #変位境界条件付加後の小行列 #与えられる
F2  = np.zeros(num_fix, dtype=np.float64)                                  #変位境界条件付加後の小行列
U1  = np.zeros((3*num_node-num_fix), dtype=np.float64)                     #変位境界条件付加後の小行列
U2  = np.zeros(num_fix, dtype=np.float64)                                  #変位境界条件付加後の小行列　#与えられる



##既知接点変位の行番号配列作成
#pythonの配列番号0始まりに変更
#各接点のx,yの順に配列が並んでいるので、xは+1、yは+2が割り振られうまく位置を計算している。
known_DOF = 3*(fix_pnt[:,0]-1) + fix_pnt[:,1] -1





#すべての行番号の中から、known_DOFの行番号のインデックスを削除
#unknown_DOFのインデックスと値が一致しているためこう書くが、本質はknown_DOFの行番号の値を削除。
unknown_DOF = np.array(range(3*num_node))
unknown_DOF = np.delete(unknown_DOF, known_DOF)
        

#zerosで作ったものを上書きしている？
#ファンシーインデックスはビュー（参照）ではなくコピーを返す。
K11 = Kmat[unknown_DOF, :]
K11 = K11[:, unknown_DOF]
K12 = Kmat[unknown_DOF, :]
K12 = K12[:, known_DOF]

#ファインシーインデックスはviewでなくcopyを返す        
F1 = Fmat[unknown_DOF]
U1 = Umat[unknown_DOF] #未知成分


K22 = Kmat[known_DOF, :]
K22 = K22[:, known_DOF]

#ファインシーインデックスはviewでなくcopyを返す              
F2 = Fmat[known_DOF] #未知成分
U2 = Umat[known_DOF] 



print('MAKE SUB-MATRIX')


print("経過時間:", time.time() - start_time)
print("処理時間:", time.time() - lap_time)
lap_time = time.time()

















# makeInverse (NUM_NODE, NUM_FIX, K11)

#solveUmatで連立方程式を直接解くので必要なしになった
# 逆行列のアルゴリズムがわからない。


#test ちゃんと単位行列になるか
#K11inv = np.linalg.inv(K11)
#a = K11inv @ K11
#originalK11 = K11.copy()

#K11を上書きして逆行列
#K11 = np.linalg.inv(K11)

#疎行列
#普通より遅い
#invはcsc_matrixを使わないと非効率
#K11 = coo_matrix(K11).tocsc()
#K11 = inv(K11)

print('MAKE K11-INV-MATRIX')

print("経過時間:", time.time() - start_time)
print("処理時間:", time.time() - lap_time)
lap_time = time.time()


















# solveUmat (NUM_NODE, NUM_FIX, Umat, K11, K12, F1, U1, U2, unknown_DOF)


#U1  = np.zeros((2*num_node-num_fix), dtype=np.float64)   #変位境界条件付加後の小行列 #前に作成済み
#fku = np.zeros((2*num_node-num_fix), dtype=np.float64)   #わからない   (F-Kd) #大文字から小文字に変更 いらない

#P.139 式(5.104)
#一気に計算する
#fku = F1 - K12 @ U2

#K11は逆行列をすでにとっている。
#U1は未知成分だったが、ここで判明
#U1 = K11 @ fku

#一気に、メモリの節約
#U1 = K11 @ (F1 - K12 @ U2)

#普通の連立方程式を解く。K11invより速いが、疎行列連立方程式のほうが早い。
#U1 = solve(K11, F1 - K12 @ U2)

#疎行列連立方程式を解く。　早い
#http://www.turbare.net/transl/scipy-lecture-notes/advanced/scipy_sparse/solvers.html#sparse-direct-solvers
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve.html#scipy.sparse.linalg.spsolve
#sparce@ndarry=ndarray, coo_matrix(1D)は[1,n]の2次元
K11 =  coo_matrix(K11).tocsr()
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve.html#scipy.sparse.linalg.spsolve
#use_umfpack：倍精度
U1 = spsolve(K11, F1 - K12 @ U2, use_umfpack=True)



#元の並びのUmatに、判明部分を代入
#view? copy?
Umat[unknown_DOF] = U1

#output省略

print('SOLVE U-MATRIX')

print("経過時間:", time.time() - start_time)
print("処理時間:", time.time() - lap_time)
lap_time = time.time()




# solveFmat(NUM_NODE, NUM_FIX, Fmat, K12, K22, F2, U1, U2, known_DOF)

#P.140

#K21=K12.T 対称性より
#F2は未知成分だったが、ここで判明
F2 = K12.T @ U1 + K22 @ U2


#元の並びのUmatに、判明部分を代入  
#view? copy?
Fmat[known_DOF] = F2

#output省略

print('SOLVE F-MATRIX')

print("経過時間:", time.time() - start_time)
print("処理時間:", time.time() - lap_time)
lap_time = time.time()





# displacement (NUM_NODE, AMP, node, Umat)


disp = np.zeros((num_node, 3), dtype=np.float64)  #amp倍した変位後の座標


disp[:,0] = node[:,0] + Umat[0::3] * amp
disp[:,1] = node[:,1] + Umat[1::3] * amp
disp[:,2] = node[:,2] + Umat[2::3] * amp

#output省略



print('CALCULATE DISPLACEMENT')

print("経過時間:", time.time() - start_time)
print("処理時間:", time.time() - lap_time)
lap_time = time.time()








# distribution (NUM_NODE, NUM_ELEME, eleme, Bmat, Dmat, Umat)

#念のため残しておく。
AVEstrain = np.zeros((num_eleme,3), dtype=np.float64)  #各四角形の平均ひずみ(εx,εy,γxy)
AVEstress = np.zeros((num_eleme,3), dtype=np.float64)  #各四角形の平均応力(σx,σy,τxy)

GAUSSstrain = np.zeros((num_eleme,4,3), dtype=np.float64) #各ガウスの積分点におけるひずみ。四角形では要素内で一定でない。
GAUSSstress = np.zeros((num_eleme,4,3), dtype=np.float64)
NODALstrain = np.zeros((num_eleme,4,3), dtype=np.float64) #各節点におけるひずみ
NODALstress = np.zeros((num_eleme,4,3), dtype=np.float64)
Nmat        = np.zeros((4,4), dtype=np.float64)   #行:各積分点　列:各形状関数N　に対応した行列。


e_Umat = np.empty(8, dtype=np.float64)               #ある四角形要素の変位


for i in range(num_eleme):
    for j in range(4): #四角形の4
        e_Umat[2*j]   = Umat[2*(eleme[i,j]-1)]     #四角形要素のx変位
        e_Umat[2*j+1] = Umat[2*(eleme[i,j]-1)+1]   #四角形要素のy変位
        
    for j in range(len(gauss_nodes)): #各ガウスの積分点を代入した時
        GAUSSstrain[i,j,:] = Bmat[:,:,i,j] @ e_Umat
        GAUSSstress[i,j,:] = Dmat[:,:,material[i]-1] @ GAUSSstrain[i,j,:]
        
        
    #ガウス求積より、積分点の平均をとれば要素全体の平均となる。
    for j in range(len(gauss_nodes)):
        AVEstrain[i,:] += GAUSSstrain[i,j,:]
        AVEstress[i,:] += GAUSSstress[i,j,:]
    AVEstrain[i,:] /= 4
    AVEstress[i,:] /= 4



#ガウスの積分点の形状関数と、ひずみ、応力から節点を求めている。
#ひずみ、応力の要素内の分布は同じ形状関数なの？
for i in range(num_eleme):
    for j in range(3): #εx,εy、τxyの3
        
        for k in range(len(gauss_nodes)):
            for l in range(4): #形状関数4こ
                Nmat[k,l] = 0.25 * (1 + polar[l,0] * gauss_nodes[k,0]) * (1 + polar[l,1] * gauss_nodes[k,1])
                
        NODALstrain[i,:,j] = solve(Nmat, GAUSSstrain[i,:,j])
        NODALstress[i,:,j] = solve(Nmat, GAUSSstress[i,:,j])

print('CALCULATE DISTRIBUTIONS')

print("経過時間:", time.time() - start_time)
print("処理時間:", time.time() - lap_time)
lap_time = time.time()


#output省略







#https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib


import matplotlib.pyplot as plt
import matplotlib.collections
import numpy as np


def showMeshPlot(nodes, elements, values, title):

    y = nodes[:,0]
    z = nodes[:,1]

    def quatplot(y,z, quatrangles, values, ax=None, **kwargs):

        if not ax: ax=plt.gca()
        yz = np.c_[y,z]
        verts= yz[quatrangles]
        pc = matplotlib.collections.PolyCollection(verts, **kwargs)
        pc.set_array(values)
        ax.add_collection(pc)
        ax.autoscale()
        return pc

    fig, ax = plt.subplots(dpi=500)
    ax.set_aspect('equal')

    pc = quatplot(y,z, np.asarray(elements), values, ax=ax, 
             edgecolor="black", cmap="rainbow",linewidths=(0.1,))
    fig.colorbar(pc, ax=ax)        
    #ax.plot(y,z, marker="o", ls="", color="black")

    ax.set(title=title, xlabel='X Axis', ylabel='Y Axis')

    plt.show()
    #fig.savefig(f'result_{title}.png')







#可視化
#https://qiita.com/itotomball/items/e63039d186fa1f564513


showMeshPlot(nodes=node, elements=eleme-1, values=np.zeros(num_eleme), title = 'mesh')


result_list = (('strain_x', AVEstrain[:,0]),('strain_y', AVEstrain[:,1]),('strain_xy', AVEstrain[:,2]),('stress_x', AVEstress[:,0]),('stress_y', AVEstress[:,1]),('stress_xy', AVEstress[:,2]))
for title, C in result_list:

    #接点番号は1から、pythonの行番号は0から始まるので修正
    showMeshPlot(nodes=disp, elements=eleme-1, values=C, title = title)
    

for matrix_name in["Kmat", "K11", "K12", "K22"] :
    #疎行列の可視化
    fig = plt.figure()
    ax = fig.add_subplot()
    fig.suptitle(matrix_name)
    ax.spy(eval(matrix_name))
    # アスペクト比を1対1に, レイアウトを調整
    #ax.set_aspect('equal')
    fig.tight_layout()
    plt.show()
    #fig.savefig('Kmat.png')

    
    
#メモリ確認
#http://harmonizedai.com/article/%E5%A4%89%E6%95%B0%E3%81%AE%E3%83%A1%E3%83%A2%E3%83%AA%E5%86%85%E5%AE%B9%E3%82%92%E4%B8%80%E8%A6%A7%E8%A1%A8%E7%A4%BA%E3%81%97%E3%81%A6/


print("{}{: >15}{}{: >15}{}".format('|','Variable Name','|','Memory[Byte]','|'))
print("|---------------|---------------|")
for var_name in dir():
    if not var_name.startswith("_"):
        print("{}{: >15}{}{: >15}{}".format('|',var_name,'|',sys.getsizeof(eval(var_name)),'|'))

