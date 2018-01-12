#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author： donglee 
# environment：python2.7
import numpy as np


def ini_position(N, D, n0, nxyz, a):    	# 位置初始化函数，这里按照面心立方结构排列
    """
    :param N:number of atoms in the system, which is n0*nxyz(0)*nxyz(1)*nxyz(2) N系统里总的原子数目=一个晶胞中的原子数*x方向晶胞数目*y方向晶胞数目*z方向晶胞数目
    :param D:dimension, which should be 3 D维度，这里为3维
    :param n0:number of atoms in a cubic unit cell, which is 4 in this case n0一个晶胞中的原子数，这里为4个
    :param nxyz:nxyz[d] is the number of unit cells in the d-th direction nxyz[d]由各个方向的晶胞数目构成的矩阵
    :param a:a[d] is the lattice constant in the d-th direction a[d]为d方向的晶格常数
    :return:r[N,3]: r[i,d] is the position of atom i along the d-th direction r[i,3]为第i个原子的三维坐标所构成的矩阵
																			    [r0(x) r0(y) r0(z)
																				 r1(x) r1(y) r1(z)
																				 ...   ...   ...
																				 rn(x) rn(y) rn(z)]
    """
    r0 = np.array([[0, 0, 0],
                  [0, 0.5, 0.5],
                  [0.5, 0, 0.5],
                  [0.5, 0.5, 0]])          	# 元胞内的相对坐标，参考面心立方结构的排列
    r = np.zeros((N, D))                   	# 创建一个N*3的矩阵以保存N个原子的初始位置坐标
    n = 0
# 循环迭代出所模拟的所有原子的位置初始坐标
    for nx in range(nxyz[0]):             	# 外层循环——x方向
        for ny in range(nxyz[1]):			# 外层循环——y方向
            for nz in range(nxyz[2]):		# 外层循环——z方向
                for m in range(n0):			# 最内层循环最基本的单元（元胞），一个元胞里共四个原子
                    r[n] = a * ((np.array([nx, ny, nz])) + r0[m])
                    # just for a test
                    # print r[n]
                    n += 1
    return r
# just for a test
# print ini_position(32, 3, 4, np.array([2, 2, 2]), np.array([1, 1, 1]))

