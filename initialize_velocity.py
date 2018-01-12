#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author： donglee 
# environment：python2.7
import numpy as np
import math


def ini_velocity(K_B, N, D, T, m):										# 初始化速度函数；速度分布要求符合玻尔兹曼分布条件且满足能量均分定理
    """
    :param K_B: Boltzmann's constant K_B波尔兹曼常数
    :param N: number of atoms in the system N系统里的总原子数
    :param D: dimension, which is 3 D维度，为3
    :param T: temperature prescribed T所设定的温度
    :param m: m[i] is the mass of atom i m[i]第i个原子的质量
    :return: v[i,d] is the velocity of atom i in the d-th direction v[i,d]第i个原子在各个维度上的速度
																	[v0(x) v0(y) v0(z)
																	 v1(x) v1(y) v1(z)
																	 ...   ...   ...
																	 vn(x) vn(y) vn(z)]
    """
    v = np.random.random((N, 3)) - 0.5									# 随机生成速度矩阵
    momentum_average_t = (np.sum(v*m, axis=0))/N						# 计算出x，y，z三个方向的动量平均值，为下面设置总动量为0做准备[vx均 vy均 vz均]
    momentum_average = momentum_average_t.T								# 矩阵转置
    for n in range(N):
        v[n] = v[n] - momentum_average/m[n]								# 总动量为0
    v *= math.sqrt(T * D * K_B * N / np.sum(m * np.sum(v**2, axis=1)))	# 速度分布满足能量均分定理
    return v
# just for a test
# print ini_velocity(8.625e-5, 5, 3, 270, np.array([1, 1, 1, 1, 1]))
