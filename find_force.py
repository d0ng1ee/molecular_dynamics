#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author： donglee 
# environment：python2.7
import numpy as np


def f_force(N, D, NN, NL, L, pbc, r):
    """

    :param N: number of atoms in the system 系统里的总原子数
    :param D: space-dimension, which should be 3 维度，为3
    :param NN: NN(i) is the number of neighbors of atom i 表示第n个原子近邻列表中原子的数目
    :param NL: NL(i,k) is the index of the k-th neighbor of atom i 表示与第n个原子相邻的原子序号所构成的矩阵，因为去掉自身，因此与每一个原子相邻的原子数目不超过N-1
	
			第0个原子近邻原子序号	[序号0 序号1 序号2 ... 序号NN[0] (数目为NN[0])
			第1个原子近邻原子序号	 序号0 序号1 序号2 ... 序号NN[1] (数目为NN[1])
			.....................	 ...   ...   ...   ... ...			
			第n个原子近邻原子序号    序号0 序号1 序号2 ... 序号NN[n]](数目为NN[n])
    :param L: L[d] is the box length in the d-th direction L[d]为研究系统的某一方向的长度
    :param pbc: pbc[d]=1[0] means periodic [free] in the d-th direction pbc[d]为研究系统的某一方向的自由度,这里为[1,1,1]
    :param r: r[i,d] is the position of atom i in the d-th direction r[i,3]为第i个原子的三维坐标所构成的矩阵
																			    [r0(x) r0(y) r0(z)
																				 r1(x) r1(y) r1(z)
																				 ...   ...   ...
																				 rn(x) rn(y) rn(z)]
    :return: f: f[i,d] is the total force on atom i in the d-th direction f[i,3]为第i个原子在各个方向所受到的力所构成的矩阵
																			    [f0(x) f0(y) f0(z)
																				 f1(x) f1(y) f1(z)
																				 ...   ...   ...
																				 fn(x) fn(y) fn(z)]
			 U: total potential energy of the system 系统的总势能
    """
    EPSILON = 1.032e-2      											# in units of eV(only for Argon)
    SIGMA = 3.405         												# in units of Angstrom(only for Argon)
    sigma_6 = SIGMA**6
    sigma_12 = SIGMA**12
    L_times_pbc = L*pbc
    U = 0																# 初始化系统的总势能
    f = np.zeros((N, D))												# 初始化每个原子上的总力
    for n1 in range(N):													# 最外层循环遍历系统中的所有原子n1
        m = np.arange(NN[n1], dtype=np.int32)
        n2 = NL[n1, m]
        r12 = np.zeros((NN[n1], D))
        r12 = r[n2] - r[n1]
        r12 = r12 - np.round(r12 / L) * L_times_pbc
        d12_square = np.sum(r12**2, axis=1)
        d12_square.shape = (NN[n1], 1)
        d_6 = d12_square ** 3
        d_8 = d12_square ** 4
        d_14 = d12_square ** 7
        d_12 = d_6 ** 2
        f12 = (sigma_6 / d_8 - 2.0 * sigma_12 / d_14) * 24.0 * EPSILON
        f[n1] = np.sum(np.dot(f12.T, r12),  axis=0)
        U += 2.0 * EPSILON * np.sum(sigma_12 / d_12 - sigma_6 / d_6)


    return [f, U]				    										# 返回一个list
# just for a test
#print f_force(4, 3, np.array([[3], [3], [3], [3]]), np.array([[1,2,3],[0,2,3],[0,1,3],[0,1,2]]), np.array([5.45,5.45,5.45]), np.array([1,1,1]),
#                np.array([[0,0,0],[0,2.725,2.725],[2.725,0,2.725],[2.725,2.725,0]]))


