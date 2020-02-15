#!/usr/bin/python
# -*- coding: utf-8 -*

# python3

##################################################################################
#
# Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
# All rights researved.
#
# 2020−02−14 Ver. 1.4 add pathtrace,
#                     change the order of the output list from [x, y, z, life, u, v, w, uid] to [x, y, z, u, v, w, uid, life]
#
##################################################################################

import os
import sys
import os.path
import getopt
from operator import itemgetter
import time
import copy
import struct

from multiprocessing import Process, Value



######################################
def usage():
	print('$ xpt2scat.py')
	print('\t$ xpt2scat.py -i <ascii | binary> // input  mode, if omit, binary input mode is selected.')
	print('\t              -o <ascii | binary> // output mode, if omit, binary input mode is selected.')
	print('\t              -n <number of processes> // default is 1 process')
	#print('\t$             -p <number of emitting points to be processed>')
	print('\t              -s <number of steps to be processed from first>')



######################################
# リストの確認
def show(lst):
	n = len(lst)
	for i in range(n):
		print(lst[i])


######################################
# リストelementsに含まれる要素keyのn要素後のの値を抜き出す
# @param [in]  elements  検査対象リスト
# @param [in]  key       検索文字列
# @param [in]  n         返却要素のオフセットインデクス
# @usage getElem(elements, `total_particles`, 1)
#        ex) total_particles 20000 >> 20000を返す、リストになければ-1
def getElemInt(elements, key, n):

	flag = 0
	
	if key in elements:
		return int(elements[elements.index(key)+n]), 1
	else:
		return -1, 0


######################################
def getElemFloat(elements, key, n):

	flag = 0
	
	if key in elements:
		return float(elements[elements.index(key)+n]), 1
	else:
		return -1.0, 0



######################################
# scatter formatで書き出し
def write_scat(fn, step, time, lst):
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
			u    = float(elements[3])
			v    = float(elements[4])
			w    = float(elements[5])
			uid  =   int(elements[6])
			life =   int(elements[7])

			f.write('%.07f %.07f %.07f %.07f %.07f %.07f %d %d\n'
							%( x, y, z, u, v, w, uid, life ) )



######################################
# scab formatで書き出し
def write_scab(fn, step, time, lst):

	np = len(lst)
	nvar = 5

	# open scab file
	try:
		ofp = open(fn, "wb")
	except:
		print("open failed: %s" % fn)
		return False

	# write step, time, np, nvar
	try:
		ofp.write(struct.pack('ifii', step, time, np, nvar))
	except:
		print("write header failed: %s" % fn)
		return False

	for i in range(np):
		elements = lst[i]
		x    = float(elements[0])
		y    = float(elements[1])
		z    = float(elements[2])
		u    = float(elements[3])
		v    = float(elements[4])
		w    = float(elements[5])
		uid  =   int(elements[6])
		life =   int(elements[7])
		
		# body
		try:
			ofp.write(struct.pack('3f', x, y, z))
			ofp.write(struct.pack('3f', u, v, w))
			ofp.write(struct.pack('i', uid))
			ofp.write(struct.pack('I', life))
		except:
			print("write data failed: %s" % fn)
			return False

	ofp.close()

	return True


######################################
# アスキーファイルの粒子数を取得
def accumlateNumPart(fn):
	
	with open(fn, 'r') as f:

		for line in f:
			elements = line.split(' ')

			di, e = getElemInt(elements, 'total_particles', 1)
			if e==1:
				return di
	
	# ここに到達するとtoral_particlesを含まない	
	print('File %s does not have \"total_particles\"' % fn)
	sys.exit(1)


######################################
# アスキーファイルの指定ステップの粒子数を積算する 
def scan_nParticle(step, rank):

	c = 0
		
	for w in rank:
		fn = 'pt_' + str(step).zfill(8) + '_' + str(w).zfill(6) + '.npt'
		q = accumlateNumPart(fn)
		c += q
	
	return c


######################################
# バイナリファイルのヘッダを読み、粒子数を取得
# Cloud::write_binary()
# 	ofs.write((char*)&stp, sizeof(unsigned));
#   ofs.write((char*)&tm,	sizeof(double));
#  	ofs.write((char*)&nParticle, sizeof(unsigned));
def readFileHeaderBinary(step, rank):
	
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
# リストelementsに含まれる要素を並べ替えlstに追加
# @param [in]     elements  検査行のリスト
# @param [in,out] lst       粒子リスト
# @retval cnt1  粒子数
def readChunkAscii(elements, lst, cnt1):
	
	if elements[0] == "1" or elements[0] == "0":
		actv =   int(elements[0])
		x    = float(elements[1])
		y    = float(elements[2])
		z    = float(elements[3])
		u    = float(elements[4])
		v    = float(elements[5])
		w    = float(elements[6])
		uid  =   int(elements[7])
		life =   int(elements[8])

		a = [x, y, z, u, v, w, uid, life]

		if len(lst) == 0:
			lst.insert(0, a)
		else:
			lst.append(a)
		
		cnt1 += 1
		
	return cnt1
		


######################################
# 指定ステップ数のときの、粒子リスト（全ランク）を作成
# @param [in]     stp  ステップ数
# @param [in]     rank ランクリスト
# @param [in,out] lst  粒子リスト
# @param [out]    tm   ファイルが含む粒子の時刻
def readChunkHeaderAscii(stp, rank, lst):
	
	tm = 0.0

	for w in rank:
		fn = 'pt_' + str(stp).zfill(8) + '_' + str(w).zfill(6) + '.npt'
		
		csz = 0
		cnt1 = 0 # ファイル単位でリセット
		
		with open(fn, 'r') as f:
			for line in f:
				elements = line.split(' ')
				
				df, e = getElemFloat(elements, 'step_time', 2) # time
				if e==1:
					tm = df
					
				di, e = getElemInt(elements, 'total_particles', 1)
				if e==1:
					csz = di

				cnt1 = readChunkAscii(elements, lst, cnt1)
		
		if not csz == cnt1:
			print("total_particles is not match %d %d" % (csz, cnt1))
			sys.exit(0)

	return tm


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
#		ofs.write((char*)&v.x, sizeof(REAL_TYPE));
#		ofs.write((char*)&v.y, sizeof(REAL_TYPE));
#		ofs.write((char*)&v.z, sizeof(REAL_TYPE));
#		ofs.write((char*)&foo, sizeof(int));
#    ofs.write((char*)&lc, sizeof(unsigned));
#
def readChunkHeaderBinary(stp, rank, lst):

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
				vx  = float(struct.unpack('f', fid.read(4))[0])  # float32, 4bytes
				vy  = float(struct.unpack('f', fid.read(4))[0])  # float32, 4bytes
				vz  = float(struct.unpack('f', fid.read(4))[0])  # float32, 4bytes
				uid = struct.unpack('i', fid.read(4))[0]  # int, 4bytes
				lc  = struct.unpack('I', fid.read(4))[0]  # unsigned, 4bytes

				a = [px, py, pz, vx, vy, vz, uid, lc]

				lst.append(a)
		fid.close()

	return tm



######################################
def rw_particle(np, s, r, n, tp, rp, imode, omode):

	lst = []
	
	for q in range(np):
		stp = s[q]
		rnk = r[q]
		
		rp.acquire()
		rp.value -= n[q]
		rp.release()
		
		fp = float(rp.value) / float(tp) * 100.0
		print('processing %.01f' % fp )
		
		if imode == 'ascii' :
			tm = readChunkHeaderAscii(stp, rnk, lst)
		else:
			tm = readChunkHeaderBinary(stp, rnk, lst)
			
		# 全ランクのデータを読み込み後、複合ソート  uidでソートし、次いでlifeでソートする
		# [x, y, z, life, u, v, w, uid]のときは key=itemgetter(7,3)
		# [x, y, z, u, v, w, uid, life] => key=itemgetter(6,7)
		lst.sort(key=itemgetter(6,7))
			
		# 書き出し
		if omode == 'ascii' :
			fn = 'streak_' + str(stp).zfill(8) + '.scat'
			write_scat(fn, stp, tm, lst)
		else:
			fn = 'streak_' + str(stp).zfill(8) + '.scab'
			if not write_scab(fn, stp, tm, lst):
				return False
			
		# lstを削除
		del lst[:]
	
	return True


######################################
def main():
	
	stepList = []
	rankList = []
	nparList = []
	num_proc = 1      # 並列数
	Nstep = 0         # number of steps to be processed from first
	stp_flag = False  # flag
	ext_i = ''        # 拡張子
	ext_o = ''        # 拡張子
	i_mode = 'binary' # 読み込みモード
	o_mode = 'binary' # 出力モード

	# スクリプト名を除いたコマンドライン引数
	argv = sys.argv[1:]
	opts, args = getopt.getopt(argv, "e:i:n:o:s:h")

	for name, value in opts:

		if name in ('-i'):
			if value == 'ascii':
				i_mode = 'ascii'

		elif name in ('-o'):
			if value == 'ascii':
				o_mode = 'ascii'
				
		elif name in ('-n'):
			num_proc = int(value)
			if not value.isalnum() or num_proc<1 : # 半角英数字のときtrue
				print('Invalid number of processors')
				return False

		#elif name in ('-p'):
		#	nPts = int(value)
		#	if not value.isalnum() or nPts<1 : # 半角英数字のときtrue
		#		print('Invalid number of emitting points to be processed')
		#		sys.exit()

		elif name in ('-s'):
			Nstep = int(value)
			stp_flag = True
			if not value.isalnum() or Nstep<1 : # 半角英数字のときtrue
				print('Invalid number of steps to be processed')
				return False

		elif name in ('-h'):
			usage()

	for arg in args:
		print('Argument="%s"' % arg)

	if i_mode == 'ascii' :
		print('\tInput  : ASCII (*.npt)')
		ext_i = '.npt'
	else:
		print('\tInput  : BINARY (*.bpt)')
		ext_i = '.bpt'
		
	if o_mode == 'ascii' :
		print('\tOutput : ASCII (*.scat)')
		ext_o = '.scat'
	else:
		print('\tOutput : BINARY (*.scab)')
		ext_o = '.scab'

	print('\tNumber of processes = %d' % num_proc)
	
	
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	start = time.time()
	

	# ディレクトリ直下にあるファイルのリストを作る
	cdir = './'
	flist = []
	xpt_files = []
	for x in os.listdir(cdir):
		if os.path.isfile(cdir + x):
			flist.append(x)

	for y in flist:
		if (y[-4:] == ext_i):
			xpt_files.append(y)


	xpt_files.sort()
	#print('\n'.join(npt_files))
	num_f = len(xpt_files)
	print('Number of files = %d\n' % num_f)

	if num_f < 1:
		print('No *%s files' % ext_i)
		return False


	# stepListにはステップ数
	# ranklistにはあるステップで現れるランク番号を保存
	for x in xpt_files:
		
		# step
		s = int(x[3:11])

		# rank
		r = int(x[12:18])

		# 一番最初の、stepList, rankListに最初の値を入れる
		if len(stepList) == 0:
			stepList.insert(0, s)
			a = []
			a.insert(0,r)
			rankList.insert(0, a)
		else:
			# sがstepListにない場合にsを追加、ランク番号を新たに挿入
			if not s in stepList:
				stepList.append(s)
				b = []
				b.insert(0,r)
				rankList.append(b)

        # sがstepListに既に存在する場合
			else:
				i = len(rankList)-1 # リストの最後の要素のインデクス
				b = rankList.pop(i) # 最後の要素を取得し、リストから削除してから、
				b.append(r)         # 追加
				rankList.append(b)

        # ここまでで、stepListにはステップ数、rankListには対応するステップでのランク番号リストが作成される

	# 処理するステップ数を決定（先頭からnum_step数）
	if stp_flag == False:
		num_step = len(stepList)
		print('Number of steps = %d\n' % num_step)
	else:
		num_step = Nstep
		print('Number of steps read = %d, but specified = %d' % (len(stepList), num_step) )


	
	# 各ステップに現れる粒子数をnparListに保存
	for j in range(num_step):
		stp = stepList[j]
		rank= rankList[j]
		
		if i_mode == 'ascii' :
			num_p = scan_nParticle(stp, rank)
		else:
			num_p = readFileHeaderBinary(stp, rank)
		nparList.append(num_p)
	
	total_pt = sum(nparList)
	print('Total particles = %d\n' % total_pt)


	# 並列処理時の担当範囲を決める
	sz = []
	st = []
	
	# 基本分割数
	div = int(num_step/num_proc)
	
	# 余りを付与するプロセス数
	mnp = num_step - div*num_proc 

	# 割り切れるとき、等分
	if num_step == div*num_proc:
		for i in range(num_proc):
			sz.append(div)
	
	# 割り切れないとき、余りを最初のプロセスから追加
	else:
		for i in range(mnp):
			sz.append(div+1)
			
		for i in range(mnp, num_proc):
			sz.append(div)
			
	# 開始点
	st.append(0)
	for i in range(1, num_proc):
		st.append(st[i-1]+sz[i-1])
	

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
		
	elapsed_time = time.time() - start
	print ("pre-process : {0}".format(elapsed_time) + " [sec]")



	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	start = time.time()

	# 並列処理
	jobs = []
	ratio = Value('I', total_pt) # unsigned 桁落ちしないため

	for i in range(num_proc):
		s = sList[i]
		r = rList[i]
		n = nList[i]
		np = sz[i]
		tp = total_pt
		job = Process(target=rw_particle, args=(np, s, r, n, tp, ratio, i_mode, o_mode, ))
		jobs.append(job)
		job.start()

	[job.join() for job in jobs]


	elapsed_time = time.time() - start
	print ("main process : {0}".format(elapsed_time) + " [sec]")

	return True


######################################
if __name__ == '__main__':
	if not main():
		sys.exit(1)

	sys.exit(0)
