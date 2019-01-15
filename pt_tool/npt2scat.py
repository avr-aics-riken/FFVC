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
import getopt
from operator import itemgetter


# Global
lst = [[]]
idx =  0
startEmit = 0
nStart = 0
nStep = 0
nProc = 0
Interval = 0


def usage():
  print('Usage:')
  print('       -h help')
  print('       -i intreval *')
  print('       -n number of procs at calculation *')
  print('       -p number of start points *')
  print('       -r generate restart files')
  print('       -s steps *')
  print('options marked * are mandatory')


# pt_file_xxxxx.lstを1行ずつ読む
def readtxt(fn):
  with open(fn, 'r') as f:
    for line in f:
      line = line.strip()
      read_npt(line)


# リストの示すnptファイルの中をパース
def read_npt(fn):
  global lst, idx, startEmit

  f1 = 0
  f2 = 0
  f3 = 0
  step_time = 0
  no_of_chunk = 0
  nChunk = 0
  nPart = 0
  gid = 0
  pid = 0
  startEmit = 0

  with open(fn, 'r') as f:
    for line in f:
      elements = line.split(' ')
      #print(len(elements))

      q = getElem(elements, 'step_time')
      if q >= 0 :
        step_time = q
        f1 += 1
	
        q = getElem(elements, 'no_of_chunks')
        if q >= 0 :
          no_of_chunk = q
          f1 += 1
          f2 = 0 # 初期化

          # f1は1ファイルにつき2
          if f1 == 2 and f2 < 6:
            q = getElem(elements, 'Chunk')
            if q >= 0 :
              nChunk = q
              f2 += 1

              q = getElem(elements, 'particles')
              if q >= 0 :
                nPart = q
                f2 += 1

              q = getElem(elements, 'group_id')
              if q >= 0 :
                gid = q
                f2 += 1

              q = getElem(elements, 'particle_id')
              if q >= 0 :
                pid = q
                f2 += 1

              q = getElem(elements, 'start_emit')
              if q >= 0 :
                startEmit = q
                f2 += 1
                f3 = 0 # 初期化
                #print(elements)

              #print(step_time, no_of_chunk, nChunk, nPart, gid, pid, startEmit)

              # nPart数だけ行を読む
              if f2 == 5 and f3 < nPart :
                #print(elements, f3, nPart)
                if elements[0] == "1" or elements[0] == "0":
                  actv =   int(elements[0])
                  x    = float(elements[1])
                  y    = float(elements[2])
                  z    = float(elements[3])
                  life =   int(elements[4])
                  u    = float(elements[5])
                  v    = float(elements[6])
                  w    = float(elements[7])
                  a = [pid, life, x, y, z, actv, u, v, w, gid, startEmit]

                if idx < np_all:
                  lst[idx] = a
                else:
                  lst.append(a)

                idx += 1
                f3 += 1
                #print(a)
					
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
    #print(elements)
    return int(elements[elements.index(key)+1])
  else:
    return -1




def main():

  global lst, nProc, nStart, nStep, Interval

  # variables
  nProc = 0
  nStart = 0
  nStep = 0
  Interval = 0

  # スクリプト名を除いたコマンドライン引数
  argv = sys.argv[1:]
  opts, args = getopt.getopt(argv, "i:r:n:hs:p:")

  for name, value in opts:

    if name in ('-i'):
      Interval = int(value)
      print('Interval = %s' % value)
    
    elif name in ('-r'):
      print('generate restart file=%s' % value)
    
    elif name in ('-n'):
      nProc = int(value)
      print('Number of processes at calulation = %s' % value)
    
    elif name in ('-p'):
      nStart = int(value)
      print('Number of start points = %s' % value)

    elif name in ('-s'):
      nStep = int(value)
      print('Number of steps = %s' % value)
    
    elif name in ('-h'):
      usage()

        
	for arg in args:
      print('Argument="%s"' % arg)


	if nProc == 0 or nStart == 0 or nStep == 0 or Interval == 0:
      usage()
      sys.exit(0)


  # 利用するリストの型で初期化しておく
  #a = [pid, life, x, y, z, actv, u, v, w, gid, startEmit]
  lst = [[-1, -1, 0.0, 0.0, 0.0, -1, 0.0, 0.0, 0.0, -1, -1] for j in range(np_all) ]


  for i in range(0, nProc):
    
    # ファイルリスト名を生成
    s = str(i)
    fn = 'pt_files_' + s.zfill(6) + '.lst'
    print(fn)

    #　データ読み込み
    readtxt(fn)

	
	# 全データを読み込み、複合ソート  pidでソートし、次いでlifeでソートする
	lst.sort(key=itemgetter(0,1))

	show()
    #print(len(lst))


if __name__ == '__main__':
	main()
