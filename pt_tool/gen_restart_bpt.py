#!/usr/bin/python
# -*- coding: utf-8 -*

# python3

##################################################################################
#
# Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
# All rights researved.
#
##################################################################################

import struct
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



def usage():
	print('Usage: gen_reatart_npt.py [-s restart step]')
	print('       -s restat step')
	print('       -h help')



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




# アスキーで書き出し
def writeRestartAscii(fn, step, time, npart, q):

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


# バイナリで書き出し
def writeRestartBinary(fn, step, time, npart, q):
	
	with open(fn, mode='wb') as f:
		
		f.write(struct.pack("I", step))
		f.write(struct.pack("f", time))
		f.write(struct.pack("I", npart))
		
		for i in range(len(lst)):
			w = lst[i]
			e = w[3]
			
			if e == q:
				px = w[0]
				py = w[1]
				pz = w[2]
				lc = w[3]
				vx = w[4]
				vy = w[5]
				vz = w[6]
				uid= w[7]
				f.write(struct.pack("f", px))
				f.write(struct.pack("f", py))
				f.write(struct.pack("f", pz))
				f.write(struct.pack("I", lc))
				f.write(struct.pack("f", vx))
				f.write(struct.pack("f", vy))
				f.write(struct.pack("f", vz))
				f.write(struct.pack("I", uid))




def readChunkHeader(stp, d):
	
	global lst
	
	fn = 'pt_' + str(stp).zfill(8) + '_' + str(d).zfill(6) + '.bpt'
		
	fid = open(fn, 'rb')
	stp = struct.unpack('I', fid.read(4))[0]  # unsigned, 4bytes
	tm  = struct.unpack('d', fid.read(8))[0]  # double, 8bytes
	np  = struct.unpack('I', fid.read(4))[0]  # unsigned, 4bytes
	csz = struct.unpack('I', fid.read(4))[0]  # unsigned, 4bytes
		
	for j in range (csz):
		psz = struct.unpack('I', fid.read(4))[0]  # unsigned, 4bytes
		uid = struct.unpack('I', fid.read(4))[0]  # unsigned, 4bytes
		eox = struct.unpack('f', fid.read(4))[0]  # float32, 4bytes
		eoy = struct.unpack('f', fid.read(4))[0]  # float32, 4bytes
		eoz = struct.unpack('f', fid.read(4))[0]  # float32, 4bytes
		est = struct.unpack('I', fid.read(4))[0]  # unsigned, 4bytes
			
		for k in range (psz):
			act = struct.unpack('i', fid.read(4))[0]  # int, 4bytes
			px  = struct.unpack('f', fid.read(4))[0]  # float32, 4bytes
			py  = struct.unpack('f', fid.read(4))[0]  # float32, 4bytes
			pz  = struct.unpack('f', fid.read(4))[0]  # float32, 4bytes
			lc  = struct.unpack('I', fid.read(4))[0]  # unsigned, 4bytes
			vx  = struct.unpack('f', fid.read(4))[0]  # float32, 4bytes
			vy  = struct.unpack('f', fid.read(4))[0]  # float32, 4bytes
			vz  = struct.unpack('f', fid.read(4))[0]  # float32, 4bytes
			uid = struct.unpack('i', fid.read(4))[0]  # int, 4bytes
				
			a = [px, py, pz, lc, vx, vy, vz, uid]
				
			if len(lst) == 0:
				lst.insert(0, a)
			else:
				lst.append(a)

	return tm




def writeRankListFile(rank):
	with open('rpt_rank.txt', mode='w') as f:

		f.write('%d number_of_ranks\n' % len(rank))

		for j in rank:
			f.write('%d\n' % j)
	



def make_bptFileList():

	global npt_files, stepList, rankList

	# ディレクトリ直下にあるnptファイルのリストを作る
	cdir = './'
	flist = []
	npt_files = []
	for x in os.listdir(cdir):
		if os.path.isfile(cdir + x):
			flist.append(x)

	for y in flist:
		if (y[-4:] == '.bpt'):
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


	make_bptFileList()


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

		tm = readChunkHeader(stp, d)

		# データを読み込み後、複合ソート  uidでソートし、次いでlifeでソートする
		lst.sort(key=itemgetter(7,3))

		# 粒子のライフカウントの最大値
		q = getMaxLife()

		# ライフカウントがqである個数
		z = getOnly(q)

	
		fn = 'restart_' + str(stp).zfill(8) + '_' + str(d).zfill(6) + '.brpt'
		#writeRestartAscii(fn, stp, tm, z, q)
		writeRestartBinary(fn, stp, tm, z, q)

		del lst[:]

	elapsed_time = time.time() - start
	print ("exec time : {0}".format(elapsed_time) + " [sec]")



if __name__ == '__main__':
	main()
