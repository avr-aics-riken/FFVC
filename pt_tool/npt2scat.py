#!/usr/bin/python
# -*- coding: utf-8 -*

# python3

##################################################################################
#
# Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
# All rights researved.
#
##################################################################################


import sys
import os.path
import getopt
from operator import itemgetter


# Global
lst = [[]]
idx =  0
nPart = 0
gid = 0
pid = 0
f1 = 0
f2 = 0
f3 = 0
np_all = 0
part_sum = 0
file_step = 0
file_time = 0.0


def usage():
	print('Usage:')
	print('       -h help')
	print('       -n number of procs at calculation')
	print('       -p number of start points')
	print('       -s start step')
	print('       -e end step')
	print('options marked * are mandatory')




# nptファイルの中をパース
def read_npt(fn):

	with open(fn, 'r') as f:
		for line in f:
			elements = line.split(' ')

			parse_1(elements)



# parse 1st level
def parse_1(elements):
	global file_step, file_time, f1, f2

	q, w = getElem2(elements, 'step_time')

	if q >= 0:
		file_step = q
		file_time = w
		f1 += 1
		f2 = 0 # 初期化

	if f1 == 1:
		parse_2(elements)

	


# parse 2nd level
def parse_2(elements):
	global f2, f3, nPart, gid, pid, part_sum

	if f2 < 4:

		q = getElem(elements, 'particles')
		if q >= 0:
			nPart = q
			part_sum += nPart
			f2 += 1

		q = getElem(elements, 'group_id')
		if q >= 0:
			gid = q
			f2 += 1

		q = getElem(elements, 'emit_pnt_id')
		if q >= 0:
			pid = q
			f2 += 1
			f3 = 0 # 初期化

		if f2 == 3:
			parse_3(elements)




# paese 3rd level
def parse_3(elements):
	global f1, f2, f3, idx, lst

	# nPartだけ行を読む, +4はemit_pnt_id, emit_origin, stat_origin, start_emitの行を含む意味
	if elements[0] == "1" or elements[0] == "0":
		actv =   int(elements[0])
		x    = float(elements[1])
		y    = float(elements[2])
		z    = float(elements[3])
		life =   int(elements[4])
		u    = float(elements[5])
		v    = float(elements[6])
		w    = float(elements[7])
		foo  =   int(elements[8])
		a = [x, y, z, life, u, v, w, pid]

		# リスト要素があれば置換、なければ追加
		if idx < np_all:
			lst[idx] = a
		else:
			lst.append(a)

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
def write_scatter(fn):
	with open(fn, mode='w') as f:

		f.write('# Scatter format\n')
		f.write('#\n')
		f.write('#TS %d %f\n' % (file_step, file_time) )
		f.write('#\n')
		# 書き出しデータ数は8
		n = len(lst)
		f.write('%d 8\n' % n)

		for i in range(n):

			l_str = [str(j) for j in lst[i]]

			f.write(' '.join(l_str))
			f.write('\n')





def main():

	global lst, np_all, f1, idx
	nProc = 0

	# スクリプト名を除いたコマンドライン引数
	argv = sys.argv[1:]
	opts, args = getopt.getopt(argv, "e:n:hs:p:")

	for name, value in opts:

    
		if name in ('-n'):
			nProc = int(value)
			print('Number of processes at calulation = %s' % value)
    
		elif name in ('-p'):
			nEmit = int(value)
			print('Number of start points = %s' % value)

		elif name in ('-s'):
			nStart = int(value)
			print('start = %s' % value)
 
		elif name in ('-e'):
			nEnd = int(value)
			print('end = %s' % value)
 
		elif name in ('-h'):
			usage()

        
	for arg in args:
		print('Argument="%s"' % arg)


	nStep = nEnd - nStart + 1


	if nProc == 0 or nEmit == 0 or nStep == 0 :
		usage()
		sys.exit(0)




	# 各ステップ毎にファイルを作成
	for step in range(nStart, nEnd+1):

		np_all = nEmit * (step + 1)

  	# 利用するリストの型で初期化しておく サイズは適当
		# a = [x, y, z, life, u, v, w, pid]
		lst = [[0.0, 0.0, 0.0, -1, 0.0, 0.0, 0.0,  -1] for j in range(np_all) ]

		s = str(step)
		idx = 0

		# プロセス数ループ
		for rank in range(0, nProc):
    
			# ファイルリスト名を生成
			r = str(rank)
			fn = 'pt_' + s.zfill(8) + '_' + r.zfill(6) + '.npt'
			print(fn)

			f1 = 0

			# ファイルが存在するならデータ読み込み
			if os.path.exists(fn):
				read_npt(fn)
			else:
				print('file %s does not exist' % fn)


	
		# 全ランクのデータを読み込み後、複合ソート  pidでソートし、次いでlifeでソートする
		lst.sort(key=itemgetter(7,3))

		#show()

		# 書き出し
		fn = 'streak_' + s.zfill(8) + '.scat'
		write_scatter(fn)

		print('length of list = %s' % len(lst))
		print('sum of particles = %s' % part_sum)

		# lstをクリア
		lst.clear()




if __name__ == '__main__':
	main()
