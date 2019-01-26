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
from operator import itemgetter
import time
import copy


# Global
lst = []
stepList = []
rankList = []
nparList = []





# リストの確認
def show():
	n = len(lst)
	for i in range(n):
		print(lst[i])




# scatter formatで書き出し
def write_scatter(fn, step, time):
	with open(fn, mode='w') as f:

		f.write('# Scatter format\n')
		f.write('#\n')
		f.write('#TS %d %f\n' % (step, time) )
		f.write('#\n')

		# 書き出しデータ数は8-3=5
		n = len(lst)
		f.write('%d 5\n' % n)

		for i in range(n):

			l_str = [str(j) for j in lst[i]]

			f.write(' '.join(l_str))
			f.write('\n')




# Chunk::write_binary()
#	ofs.write((char*)&pchkSz, sizeof(unsigned));
#	ofs.write((char*)&uid, sizeof(unsigned));
#	ofs.write((char*)&EmitOrigin.x, sizeof(REAL_TYPE));
#	ofs.write((char*)&EmitOrigin.y, sizeof(REAL_TYPE));
#	ofs.write((char*)&EmitOrigin.z, sizeof(REAL_TYPE));
#	ofs.write((char*)&est, sizeof(unsigned));

#		ofs.write((char*)&actv, sizeof(int));
#		ofs.write((char*)&p.x, sizeof(REAL_TYPE));
#		ofs.write((char*)&p.y, sizeof(REAL_TYPE));
#		ofs.write((char*)&p.z, sizeof(REAL_TYPE));
#		ofs.write((char*)&lc, sizeof(unsigned));
#		ofs.write((char*)&v.x, sizeof(REAL_TYPE));
#		ofs.write((char*)&v.y, sizeof(REAL_TYPE));
#		ofs.write((char*)&v.z, sizeof(REAL_TYPE));
#		ofs.write((char*)&foo, sizeof(int));

def readChunkHeader(stp, rank):
	
	global lst

	for w in rank:
		fn = 'pt_' + str(stp).zfill(8) + '_' + str(w).zfill(6) + '.bpt'

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




# ヘッダを読む
# Cloud::write_binary()
# 	ofs.write((char*)&stp, sizeof(unsigned));
#   ofs.write((char*)&tm,	sizeof(double));
#  	ofs.write((char*)&nParticle, sizeof(unsigned));

def readFileHeader(fn):

	fid = open(fn, 'rb')

	stp = struct.unpack('I', fid.read(4))  # unsigned, 4bytes
	tm  = struct.unpack('d', fid.read(8))  # double, 8bytes
	np  = struct.unpack('I', fid.read(4))  # unsigned, 4bytes

	#print(stp, tm, np)
	#print(np[0])
	return np[0]




# 粒子数リストを作成
def scan_nPart():
	
	for stp in stepList:
		q = stepList.index(stp)
		r = rankList[q]

		c = 0
		for w in r:

			fn = 'pt_' + str(stp).zfill(8) + '_' + str(w).zfill(6) + '.bpt'

			c += readFileHeader(fn)

		
		if len(nparList) == 0:
			nparList.insert(0, c)
		else:
			nparList.append(c)

	# check
	#for x in nparList:
	#	print(x)




def main():

	global stepList, rankList, lst


	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	start = time.time()

	# ディレクトリ直下にあるbptファイルのリストを作る
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


	elapsed_time = time.time() - start
	print ("make file list : {0}".format(elapsed_time) + " [sec]")
	

	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	start = time.time()

	# 各ステップに表れる粒子数をnparListに保存
	#scan_nPart()


	# 各ステップの粒子データをリストに登録しソート
	for stp in stepList:
		q = stepList.index(stp)
		r = rankList[q]
		print('processing %d' % stp)

		# 利用するリストの型で必要数だけ初期化しておく
		# a = [x, y, z, life, u, v, w, uid]
		#n = nparList[q]
		#lst = [[0.0, 0.0, 0.0, -1, 0.0, 0.0, 0.0, -1] for j in range(n) ]
		#print('lst size = %d' % len(lst))


		tm = readChunkHeader(stp, r)


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
