# 3D8Node_FEM
これは3次元8節点6面体アイソパラメトリック要素の有限要素法プログラムです。




# Verion 1.0.0 概要

少々バグがあるので非推奨

このプログラムは、
[2D4Node_FEM Version1.1.0](https://github.com/yuki-2000/2D4Node_FEM/releases/tag/v1.1.0)
をもとに、
[3D-8NODE-ISOPARAMETRIC2.F90](https://github.com/yuki-2000/3D8Node_FEM/blob/main/3D-8NODE-ISOPARAMETRIC2.F90)
に従ってpythonに書き換えたプログラムです。

## What's Changed
* 3D-8NODE-ISOPARAMETRIC-FEM2に差し替え by @yuki-2000 in https://github.com/yuki-2000/3D8Node_FEM/pull/1
* Conform to 3 d8node by @yuki-2000 in https://github.com/yuki-2000/3D8Node_FEM/pull/4

## New Contributors
* @yuki-2000 made their first contribution in https://github.com/yuki-2000/3D8Node_FEM/pull/1

**Full Changelog**: https://github.com/yuki-2000/3D8Node_FEM/commits/v1.0.0



# プログラム概要

これは8節点6面体アイソパラメトリック要素を使用した有限要素法プログラムです。
材料は、等方材料にのみ対応しています。

可視化については、smartgraphを使用することを前提としています。
モデルの作成は、femapを使用して、アバカス形式に変換後、中身をうまくコピペします。



[2D4Node_FEM Version1.1.0](https://github.com/yuki-2000/2D4Node_FEM/releases/tag/v1.1.0)
との変更点は、2Dから3Dに対応するために行列のサイズなどをいじったり、しただけで、ほとんだ変わっていません。
ただし、結果の可視化は不可能となっています。



# Fortranからの変数名の変更


|  意味 |Fortran | このプログラム |
|---|-----|---|
|    素材数    |  MATERIAL | num_material   |
|    素材識別番号    |  mat_num | material   |
| 要素構成節点座標   | e_point  | e_node  |
| ヤコビ行列式 | det  |  det_Jacobi |
| 積分点 | gauss  | gauss_nodes|


# 数式

$$\begin{equation}
\boldsymbol{J}=
\begin{bmatrix}
\frac{\partial x}{\partial \xi}  &\frac{\partial y}{\partial \xi} &\frac{\partial z}{\partial \xi}  \\
\frac{\partial x}{\partial \eta}  &\frac{\partial y}{\partial \eta}  &\frac{\partial z}{\partial \eta} \\
\frac{\partial x}{\partial \zeta}  &\frac{\partial y}{\partial \zeta}  &\frac{\partial z}{\partial \zeta} \\
\end{bmatrix}
\end{equation}$$






# Version 1.0.1変更点


 入力ファイルを1000要素の梁の曲げに変更

平均ひずみ、平均応力を出力するように変更

メモリ表示を文字に拡大

#6 バグフィックス

アバカス形式からコピペで作るためにカンマ区切りに変更



## What's Changed
* Chnage input file by @yuki-2000 in https://github.com/yuki-2000/3D8Node_FEM/pull/8


**Full Changelog**: https://github.com/yuki-2000/3D8Node_FEM/compare/v1.0.0...v1.0.1








