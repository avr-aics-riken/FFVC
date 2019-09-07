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

from multiprocessing import Process, Value


######################################
def usage():
	print('$ bpt2scat_mp.py [num_proc]')



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
def readChunkHeader(stp, rank, lst):

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
				px  = float(struct.unpack('f', fid.read(4))[0])  # float32, 4bytes
				py  = float(struct.unpack('f', fid.read(4))[0])  # float32, 4bytes
				pz  = float(struct.unpack('f', fid.read(4))[0])  # float32, 4bytes
				lc  = struct.unpack('I', fid.read(4))[0]  # unsigned, 4bytes
				vx  = float(struct.unpack('f', fid.read(4))[0])  # float32, 4bytes
				vy  = float(struct.unpack('f', fid.read(4))[0])  # float32, 4bytes
				vz  = float(struct.unpack('f', fid.read(4))[0])  # float32, 4bytes
				uid = struct.unpack('i', fid.read(4))[0]  # int, 4bytes

				a = [px, py, pz, lc, vx, vy, vz, uid]

				lst.append(a)

	return tm



######################################
# ヘッダを読む
# Cloud::write_binary()
# 	ofs.write((char*)&stp, sizeof(unsigned));
#   ofs.write((char*)&tm,	sizeof(double));
#  	ofs.write((char*)&nParticle, sizeof(unsigned));
def readFileHeader(step, rank):
	
	tp = 0
	
	for w in rank:
		fn = 'pt_' + str(step).zfill(8) + '_' + str(w).zfill(6) + '.bpt'
	
		fid = open(fn, 'rb')
	
		stp = struct.unpack('I', fid.read(4))  # unsigned, 4bytes
		tm  = struct.unpack('d', fid.read(8))  # double, 8bytes
		np  = struct.unpack('I', fid.read(4))  # unsigned, 4bytes
		tp += np[0]

	return tp


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
		tm = readChunkHeader(stp, rnk, lst)
		
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


	# 各ステップに表れる粒子数をnparListに保存
	for j in range(len(stepList)):
		stp = stepList[j]
		rank= rankList[j]
		num_p = readFileHeader(stp, rank)
		nparList.append(num_p)
	
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
	ratio = Value('I', total_pt) # unsigned 桁落ちしないため

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
