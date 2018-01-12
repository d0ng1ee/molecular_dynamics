#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author： donglee 
# environment：python2.7
import numpy as np


def f_neighbor(N, L, pbc, rc, r):					# 寻找近邻列表
    """

    :param N:number of atoms in the system N为系统里总的原子数目
    :param L:L[d] is the box length in the d-th direction L[d]为研究系统的某一方向的长度
    :param pbc:pbc[d]=1[0] means periodic [free] in the d-th direction pbc[d]为研究系统的某一方向的自由度,这里为[1,1,1]
    :param rc:cutoff distance rc为设定的截断距离用于寻找近邻列表
    :param r: r[i,d] is the position of atom i in the d-th direction r[i,3]为第i个原子的三维坐标所构成的矩阵
																			    [r0(x) r0(y) r0(z)
																				 r1(x) r1(y) r1(z)
																				 ...   ...   ...
																				 rn(x) rn(y) rn(z)]
    :return:NN[N,1]表示第n个原子近邻列表中原子的数目
			NL[N,N-1]表示与第n个原子相邻的原子序号所构成的矩阵，因为去掉自身，因此与每一个原子相邻的原子数目不超过N-1
					
			第0个原子近邻原子序号	[序号0 序号1 序号2 ... 序号NN[0] (数目为NN[0])
			第1个原子近邻原子序号	 序号0 序号1 序号2 ... 序号NN[1] (数目为NN[1])
			.....................	 ...   ...   ...   ... ...			
			第n个原子近邻原子序号    序号0 序号1 序号2 ... 序号NN[n]](数目为NN[n])
    """
    NN = np.zeros((N, 1), dtype=np.int32)
    NL = np.zeros((N, N-1), dtype=np.int32)
    L_times_pbc = L * pbc
    rc_square = rc * rc								
    for n1 in range(N-1):							# 外层循环对象为序号为n1的原子
        for n2 in range(n1+1, N):     				# 固定n1，遍历其他原子与n1所固定原子的距离（只遍历序号大于n1的原子避免重复遍历）
            r12 = r[n2] - r[n1]						# 两个原子相对位置的矢量
            r12 = r12 - np.round(r12/L)*L_times_pbc	# 解决周期性边界条件：如果x_ij<-Lx/2,则将x_ij换为x_ij+Lx；如果x_ij>+Lx/2,则将x_ij换为x_ij-Lx。
            d12_square = np.sum(r12 * r12)
            if d12_square < rc_square:              # 若被判断为近邻原子
                NL[n1, int(NN[n1])] = n2            # 第n1个原子近邻原子序号增加一个n2
                NN[n1] += 1                         # 第n1个原子近邻列表中原子的数目加1
                NL[n2, int(NN[n2])] = n1            # 第n2个原子近邻原子序号增加一个n1
                NN[n2] += 1                         # 第n2个原子近邻列表中原子的数目加1
    return [NN, NL]                                # 返回一个list



