#! /usr/bin/env python
#
# scat2scab
#

import sys
import struct

def convert(srcf, dstf):
    # open srcf
    try:
        ifp = open(srcf, "r")
    except:
        print("open failed: %s" % srcf)
        return False

    # read head comment lines
    st = 0
    tm = 0.0
    line = ifp.readline()
    while line:
        if len(line) < 1:
            line = ifp.readline()
            continue
        if line[0] != '#': break
        if line.startswith('#TS'):
            toks = line[3:].split()
            try:
                st = int(toks[0])
                tm = float(toks[1])
            except:
                pass
        line = ifp.readline()
        continue # end of while(line)

    # parse np, nvar
    np = 0
    nvar = 0
    toks = line.split()
    try:
        np = int(toks[0])
        nvar = int(toks[1])
    except:
        print("can not get np nor nvar: %s" % srcf)
        return False
    if np < 1 or nvar < 1:
        print("invalid np(%d) or nvar(%d)" % (np, nvar))
        return False

    # open dstf
    try:
        ofp = open(dstf, "wb")
    except:
        print("open failed: %s" % dstf)
        return False

    # write step, time, np, nvar
    try:
        ofp.write(struct.pack('ifii', st, tm, np, nvar))
    except:
        print("write header failed: %s" % dstf)
        return False

    # read/write datas
    n = 0
    line = ifp.readline()
    vals = [0.0] * nvar
    vfmt = '%df' % nvar
    while line:
        toks = line.split()
        try:
            x = float(toks[0])
            y = float(toks[1])
            z = float(toks[2])
            for i in range(nvar):
                vals[i] = float(toks[i+3])
        except:
            line = ifp.readline()
            continue
        try:
            ofp.write(struct.pack('3f', x, y, z))
            ofp.write(struct.pack(vfmt, *vals))
        except:
            print("write data failed: %s" % dstf)
            return False
        n = n + 1
        if n >= np:
            break
        line = ifp.readline()
        continue # end of while(line)

    # epilogue
    ifp.close()
    ofp.close()
    return True


def usage():
    print("usage: scat2scab.py infile.scat outfile.scab")
    

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] in ('-h', '--help'):
        usage()
        sys.exit(0)
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)

    if not convert(sys.argv[1], sys.argv[2]):
        sys.exit(1)

    sys.exit(0)
