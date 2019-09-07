#!/usr/bin/python
# -*- coding: utf-8 -*

# python3

##################################################################################
#
# Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
# All rights researved.
#
##################################################################################

import os
import sys
import os.path
import getopt
from operator import itemgetter
import time
import copy

from multiprocessing import Process, Value



######################################
def usage():
	print('$ npt2scat_mp.py [num_proc]')



######################################
def parse(elements, f1, f2, f3, lst, nPart):
	q, w = getElem2(elements, 'step_time')

	# f1は1ファイルにつき1回だけ1になる == 'step_time'は1回だけ
	if q >= 0:
		f1 += 1
		f2 = 0 # 初期化

	if f1 == 1:
		f2, f3, nPart = read_chunk(elements, f2, f3, lst, nPart)

	return w, f1, f2, f3, nPart


######################################
def read_chunk(elements, f2, f3, lst, nPart):
	
	q = getElem(elements, 'particles')
	if q >= 0:
		nPart = q
		f2 += 1
		f3 = 0 # 初期化

	if f2 == 1:
		f2, f3 = read_body(elements, f2, f3, lst, nPart)

	return f2, f3, nPart


######################################
def read_body(elements, f2, f3, lst, nPart):

	if elements[0] == "1" or elements[0] == "0":
		actv =   int(elements[0])
		x    = float(elements[1])
		y    = float(elements[2])
		z    = float(elements[3])
		life =   int(elements[4])
		u    = float(elements[5])
		v    = float(elements[6])
		w    = float(elements[7])
		uid  =   int(elements[8])

		a = [x, y, z, life, u, v, w, uid]
		lst.append(a)
		f3 += 1

	if f3 == nPart:
		f2 = 0
		f3 = 0

	return f2, f3


######################################
# リストの確認
def show():
	n = len(lst)
	for i in range(n):
		print(lst[i])


######################################
# keyの値を抜き出す
def getElem(elements, key):
  if key in elements:
    return int(elements[elements.index(key)+1])
  else:
    return -1


######################################
# keyの値を抜き出す 2要素
def getElem2(elements, key):
	if key in elements:
		return int(elements[elements.index(key)+1]), float(elements[elements.index(key)+2])
	else:
		return -1, -1.0



######################################
# scatter formatで書き出し
def write_scatter(fn, step, time, lst):
	with open(fn, mode='w') as f:

		f.write('# Scatter format\n')
		f.write('#\n')
		f.write('#TS %d %f\n' % (step, time) )
		f.write('#\n')

		# 書き出しデータ数は8-3=5
		n = len(lst)
		f.write('%d 5\n' % n)

		for i in range(n):
			elements = lst[i]
			x    = float(elements[0])
			y    = float(elements[1])
			z    = float(elements[2])
			life =   int(elements[3])
			u    = float(elements[4])
			v    = float(elements[5])
			w    = float(elements[6])
			uid  =   int(elements[7])

			f.write('%.07f %.07f %.07f %d %.07f %.07f %.07f %d\n'
							%( x, y, z, life, u, v, w, uid) )


######################################
def makeParticleList(stp, rank, lst):
	
	for w in rank:
		fn = 'pt_' + str(stp).zfill(8) + '_' + str(w).zfill(6) + '.npt'
		f1 = 0
		f2 = 0
		f3 = 0
		nPart = 0
		
		with open(fn, 'r') as f:

			for line in f:
				elements = line.split(' ')
				tm, f1, f2, f3, nPart = parse(elements, f1, f2, f3, lst, nPart)

	return tm



######################################
# 粒子数を積算する
def accumlateNumPart(fn):

	with open(fn, 'r') as f:

		ens = 0

		for line in f:
			elements = line.split(' ')

			y = getElem(elements, 'total_particles')
			if y > 0:
				ens = 1
				return y, ens


######################################
# 粒子数リストを作成
def scan_nPart(stepList, rankList, nparList):
	
	for stp in stepList:
		q = stepList.index(stp) # stp を含むのは何番目か？
		r = rankList[q]
		c = 0
		
		for w in r:
			fn = 'pt_' + str(stp).zfill(8) + '_' + str(w).zfill(6) + '.npt'
			q, s = accumlateNumPart(fn)
			
			if s == 0:
				print('File %s does not have \"total_particles\"' % fn)
				sys.exit(0)
			else:
				c += q
		
		#if len(nparList) == 0:
		#	nparList.insert(0, c)
		#else:
		nparList.append(c)



######################################
def rw_particle(np, s, r, n, tp, rp):

	for q in range(np):
		stp = s[q]
		rnk = r[q]
		
		rp.acquire()
		rp.value -= n[q]
		rp.release()
		
		fp = float(rp.value) / float(tp) * 100.0
		print('processing %.01f' % fp )
		lst = []
		tm = makeParticleList(stp, rnk, lst)
			
		# 全ランクのデータを読み込み後、複合ソート  uidでソートし、次いでlifeでソートする
		lst.sort(key=itemgetter(7,3))
			
		# 書き出し
		fn = 'streak_' + str(stp).zfill(8) + '.scat'
		write_scatter(fn, stp, tm, lst)
			
		# lstを削除
		del lst[:]



######################################
def main(num_proc):
	
	start = time.time()
	
	stepList = []
	rankList = []
	nparList = []


	# ディレクトリ直下にあるnptファイルのリストを作る
	cdir = './'
	flist = []
	npt_files = []
	for x in os.listdir(cdir):
		if os.path.isfile(cdir + x):
			flist.append(x)

	for y in flist:
		if (y[-4:] == '.npt'):
			npt_files.append(y)


	npt_files.sort()
	#print('\n'.join(npt_files))
	print('Number of files = %d\n' % len(npt_files))


	# stepListにはステップ数
	# ranklistにはあるステップで現れるランク番号を保存

	a = []

	for x in npt_files:
		
		# step
		s = int(x[3:11])

		# rank
		r = int(x[12:18])

		if len(stepList) == 0:
			stepList.insert(0, s)
			a.insert(0,r)
			rankList.insert(0, a)
		else:
			if not s in stepList:
				stepList.append(s)
				del a[:]
				a.insert(0, r)
				rankList.append(a)
			else:
				t = stepList.index(s)
				b = copy.deepcopy(rankList[t])
				b.append(r)
				rankList[t] = b

	num_step =  len(stepList)
	print('Number of steps = %d\n' % num_step)

	#print(stepList)
	
	
	# 各ステップに表れる粒子数をnparListに保存
	scan_nPart(stepList, rankList, nparList)
	total_pt = sum(nparList)
	print('Total particles = %d\n' % total_pt)


	# 並列処理時の担当範囲を決める
	sz = []
	st = []
	div = int(num_step/num_proc)
	if num_step % div != 0:
		div += 1

	if num_step % div == 0:
		for i in range(num_proc):
			sz.append(div)
	else:
		for i in range(num_proc-1):
			sz.append(div)
		sz.append(num_step % div)

	st.append(0)
	for i in range(1, num_proc):
		st.append(st[i-1]+sz[i-1])


	elapsed_time = time.time() - start
	print ("pre-process : {0}".format(elapsed_time) + " [sec]")

	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	start = time.time()


	# 各プロセス毎のリストを作成
	sList = []
	rList = []
	nList = []

	for i in range(num_proc):
		s = st[i]
		e = st[i]+sz[i]
		sList.append( stepList[s:e] )
		rList.append( rankList[s:e] )
		nList.append( nparList[s:e] )

	# 並列処理
	jobs = []
	ratio = Value('I', total_pt)

	for i in range(num_proc):
		s = sList[i]
		r = rList[i]
		n = nList[i]
		np = sz[i]
		tp = total_pt
		job = Process(target=rw_particle, args=(np, s, r, n, tp, ratio, ))
		jobs.append(job)
		job.start()

	[job.join() for job in jobs]


	elapsed_time = time.time() - start
	print ("main process : {0}".format(elapsed_time) + " [sec]")



######################################
if __name__ == '__main__':
	argv = sys.argv
	
	ag = argv[1]
	if ag.isalnum(): # 半角英数字のときtrue
		num_proc = int(ag)
		main(num_proc)
	else:
		usage()

