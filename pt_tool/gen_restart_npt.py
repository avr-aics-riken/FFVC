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
lst = []
stepList = []
rankList = []

idx =  0
nPart = 0
f1 = 0
f2 = 0
f3 = 0
tm = 0.0



def usage():
	print('Usage: gen_reatart_npt.py [-s restart step]')
	print('       -s restat step')
	print('       -h help')



def parse_1(elements):
	global f1, f2, tm

	q, w = getElem2(elements, 'step_time')

	# f1は1ファイルにつき1回だけ1にする
	if q >= 0:
		f1 += 1
		tm = w
		f2 = 0 # 初期化

	if f1 == 1:
		read_chunk(elements)




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
	global f2, f3, lst
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

		if len(lst) == 0:
			lst.insert(0, a)
		else:
			lst.append(a)

		f3 += 1


	if f3 == nPart:
		f2 = 0
		f3 = 0




# リストの確認
def show():
	n = len(lst)
	for i in range(n):
		print(lst[i])

 
def getOnly(q):
	n = len(lst)
	c = 0
	for i in range(n):
		w = lst[i]
		e = w[3]
		if e == q:
		#	print(lst[i])
			c += 1

	return c



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




# npt formatで書き出し
def writeRestartNpt(fn, step, time, npart, q):

	with open(fn, mode='w') as f:

		f.write('%d step\n' % step)
		f.write('%f time\n' % time)
		f.write('%d number_of_particles\n' % npart )

		for i in range(len(lst)):
			w = lst[i]
			e = w[3]

			if e == q:
				l_str = [str(j) for j in lst[i]]
				f.write(' '.join(l_str))
				f.write('\n')



def makeParticleList(stp, w):
	global f1
	
	fn = 'pt_' + str(stp).zfill(8) + '_' + str(w).zfill(6) + '.npt'
	print(fn)

	with open(fn, 'r') as f:

		f1 = 0
		for line in f:
			elements = line.split(' ')
			tm = parse_1(elements)

	return tm




def writeRankListFile(rank):
	with open('rpt_rank.txt', mode='w') as f:

		f.write('%d number_of_ranks\n' % len(rank))

		for j in rank:
			f.write('%d\n' % j)
	



def make_nptFileList():

	global npt_files, stepList, rankList

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


def getMaxLife():
	lmax = 0

	for j in range(len(lst)):
		w = lst[j]
		e = w[3] # life
		if e > lmax:
			lmax = e
	
	return lmax




def main():

	global lst

	nStep = 0

  # スクリプト名を除いたコマンドライン引数
	argv = sys.argv[1:]
	opts, args = getopt.getopt(argv, "hs:")

	for name, value in opts:

		if name in ('-s'):
			nStep = int(value)
			print('Restat step = %s' % value)
    
		elif name in ('-h'):
			usage()

        
	for arg in args:
		print('Argument="%s"' % arg)




	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	start = time.time()


	make_nptFileList()


	if nStep == 0:
		m = len(stepList)
		nStep = stepList[m-1]
		print(m, nStep)
	else:
		print(stepList.index(nStep), nStep)


	

	stp = nStep
	q = stepList.index(stp)
	r = rankList[q]

	writeRankListFile(r)

	# a = [x, y, z, life, u, v, w, uid]

	for d in r:
		makeParticleList(stp, d)

		# データを読み込み後、複合ソート  uidでソートし、次いでlifeでソートする
		lst.sort(key=itemgetter(7,3))

		# 粒子のライフカウントの最大値
		q = getMaxLife()
		#print(q)

		# ライフカウントがqである個数
		z = getOnly(q)

	
		fn = 'restart_' + str(stp).zfill(8) + '_' + str(d).zfill(6) + '.nrpt'	
		writeRestartNpt(fn, stp, tm, z, q)

		del lst[:]

	elapsed_time = time.time() - start
	print ("exec time : {0}".format(elapsed_time) + " [sec]")



if __name__ == '__main__':
	main()
