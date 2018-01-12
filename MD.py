#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author： donglee 
# environment：python2.7
# 模拟的为单原子分子所构成的系统
# My unit system: energy-eV; length-Angstrom; mass-atomic mass unit
import numpy as np
import math

import initialize_position as i_p
import initialize_velocity as i_i
import find_neighbor as f_n
import find_force as f_f


def MD(D, n0, nxyz, a, pbc, Ne, Np, Ns, rc, dt, T):
    """

    :param D:dimension, which should be 3 D维度，为3
    :param n0:number of atoms in a cubic unit cell, which is 4 in this case 一个晶胞中的原子数，这里为4个
    :param nxyz:nxyz(d) is the number of unit cells in the d-th direction nxyz[d]由各个方向的晶胞数目构成的矩阵
    :param a:a(d) is the lattice constant in the d-th direction a[d]为d方向的晶格常数
    :param pbc:pbc(d)=1(0) means periodic (free) in the d-th direction pbc[d]为研究系统的某一方向的自由度,这里为[1,1,1]
    :param Ne: number of time steps in the equilibration stage 在平衡阶段中的时间步数
    :param Np: number of time steps in the production stage 生产阶段的时间步数 
    :param Ns: sampling interval in the production stage 在生产阶段的采样间隔
    :param rc: cutoff distance 截断距离
    :param dt: time step of integration in units of fs 以fs为单位的时间间隔
    :param T: temperature prescribed in units of K 温度规定单位为K
    :return:E[:, 2]: total eneryg per atom at regular time points 在规则时间点的每原子的总能量
			E[:, 1]: kinetic eneryg per atom at regular time points 在规则时间点的每原子的总动能
			E[:, 0]: potenital eneryg per atom at regular time points 在规则时间点的每原子的总势能
    """
    K_B = 8.625e-5																# Boltzmann's constant in my unit system 
    TIME_UNIT_CONVERSION = 10.18												# from fs to my unit system
    N = n0*nxyz[0]*nxyz[1]*nxyz[2]												# number of atoms
    L = a*nxyz																	# box size (Angstrom)
    dt /= TIME_UNIT_CONVERSION												    # time step in my unit system
    m = np.ones((N,1))*40														# mass for Argon atom in my unit system
    r = i_p.ini_position(N, D, n0, nxyz, a)										# 初始化位置坐标
    v = i_i.ini_velocity(K_B, N, D, T, m)										# 初始化速度坐标
    [NN, NL] = f_n.f_neighbor(N, L, pbc, rc, r)									# 寻找近邻粒子
    [f, U] = f_f.f_force(N, D, NN, NL, L, pbc, r)								# 初始化各个粒子的受力
    E = np.zeros([(Np+Ne)/Ns, 3])													# 能量数据
    print 0.5 * np.sum(m * np.sum(v ** 2, axis=1))/N                            # just for a test
    i = 0
# 这里采用Velocity-Verlet算法
    for step in range(Ne + Np):
        if step % 100 == 0:
            i += 1
            print i,															# 进度条orz
        for d in range(D):
            v[:, d] += (f[:, d]/m[0])*(dt*0.5)
            r[:, d] += v[:, d]*dt
        [f, U] = f_f.f_force(N, D, NN, NL, L, pbc, r)
        for d in range(D):
            v[:, d] += (f[:, d]/m[0])*(dt*0.5)									# Velocity-Verlet算法进行迭代
        if step <= Ne:
            v *= math.sqrt(T * D * K_B * N / np.sum(m * np.sum(v ** 2, axis=1)))# 平衡（控制温度）过程用1000步
        if step % Ns == 0:													# 产出（不控制温度）过程也用1000步，并每隔20步抽样
            E[(step/Ns), 0] = U/N
            E[(step/Ns), 1] = 0.5*np.sum(m * np.sum(v ** 2, axis=1))/N
    E[:, 2] = E[:, 1] + E[:, 0]
    return E[:, 2], E[:, 1], E[:, 0]
print MD(3, 4, np.array([3,3,3]), np.array([5.45,5.45,5.45]), [1,1,1], 5000, 5000, 20, 10, 5, 80)



