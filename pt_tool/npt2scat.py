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


# Global
lst = [[]]
stepList = []
rankList = []
nparList = []

idx =  0
nPart = 0
f1 = 0
f2 = 0
f3 = 0




def usage():
	print('Usage:')
	print('       -h help')



def parse_1(elements):
	global f1, f2

	q, w = getElem2(elements, 'step_time')

	# f1は1ファイルにつき1回だけ1にする
	if q >= 0:
		f1 += 1
		f2 = 0 # 初期化

	if f1 == 1:
		read_chunk(elements)

	return w



def read_chunk(elements):
	global f2, f3, nPart

	q = getElem(elements, 'particles')
	if q >= 0:
		nPart = q
		f2 += 1
		f3 = 0 # 初期化

	if f2 == 1:
		read_body(elements)




def read_body(elements):
	global f2, f3, idx, lst
	#print(elements)

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

		#print(idx)
		lst[idx] = a

		idx += 1
		f3 += 1


	if f3 == nPart:
		f2 = 0
		f3 = 0




# リストの確認
def show():
	n = len(lst)
	for i in range(n):
		print(lst[i])



# keyの値を抜き出す
def getElem(elements, key):
  if key in elements:
    return int(elements[elements.index(key)+1])
  else:
    return -1



# keyの値を抜き出す 2要素
def getElem2(elements, key):
	if key in elements:
		return int(elements[elements.index(key)+1]), float(elements[elements.index(key)+2])
	else:
		return -1, -1.0




# scatter formatで書き出し
def write_scatter(fn, step, time):
	with open(fn, mode='w') as f:

		f.write('# Scatter format\n')
		f.write('#\n')
		f.write('#TS %d %f\n' % (step, time) )
		f.write('#\n')

		# 書き出しデータ数は8
		n = len(lst)
		f.write('%d 8\n' % n)

		for i in range(n):

			l_str = [str(j) for j in lst[i]]

			f.write(' '.join(l_str))
			f.write('\n')



def makeParticleList(stp, rank):
	global f1
	
	for w in rank:
		fn = 'pt_' + str(stp).zfill(8) + '_' + str(w).zfill(6) + '.npt'

		with open(fn, 'r') as f:

			f1 = 0
			for line in f:
				elements = line.split(' ')
				tm = parse_1(elements)

	return tm




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



# 粒子数リストを作成
def scan_nPart():
	
	for stp in stepList:
		q = stepList.index(stp)
		r = rankList[q]

		c = 0
		for w in r:

			fn = 'pt_' + str(stp).zfill(8) + '_' + str(w).zfill(6) + '.npt'
			#print(fn)
			q, w = accumlateNumPart(fn)
			if w == 0:
				print('File %s does not have \"total_particles\"' % fn)
				sys.exit(0)
			else:
				c += q

		
		if len(nparList) == 0:
			nparList.insert(0, c)
		else:
			nparList.append(c)

	# check
	#for x in nparList:
	#	print(x)




def main():

	global npt_files, stepList, rankList, lst, idx


	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	start = time.time()

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
			#print(a)
			rankList.insert(0, a)
			#print(rankList)

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
				#print(s, rankList[t])


	# check
	#for x in stepList:
	#	r = stepList.index(x)
	#	print(stepList[r], rankList[r])	
	#	if s > 100:
	#		sys.exit(0)


	elapsed_time = time.time() - start
	print ("make file list : {0}".format(elapsed_time) + " [sec]")
	

	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	start = time.time()

	# 各ステップに表れる粒子数をnparListに保存
	scan_nPart()

	# 各ステップの粒子データをリストに登録しソート
	for stp in stepList:
		q = stepList.index(stp)
		r = rankList[q]
		print('processing %d' % stp)

		# 利用するリストの型で必要数だけ初期化しておく
		# a = [x, y, z, life, u, v, w, uid]
		n = nparList[q]
		lst = [[0.0, 0.0, 0.0, -1, 0.0, 0.0, 0.0, -1] for j in range(n) ]
		#print('lst size = %d' % len(lst))

		# lst[idx] のインデクス初期化
		idx = 0

		tm = makeParticleList(stp, r)

		# 全ランクのデータを読み込み後、複合ソート  uidでソートし、次いでlifeでソートする
		lst.sort(key=itemgetter(7,3))


		# 書き出し
		fn = 'streak_' + str(stp).zfill(8) + '.scat'
		write_scatter(fn, stp, tm)

		# lstを削除
		del lst[:]


	elapsed_time = time.time() - start
	print ("scan files : {0}".format(elapsed_time) + " [sec]")



if __name__ == '__main__':
	main()
