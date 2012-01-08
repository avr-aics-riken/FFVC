# Makefile

DIST_DIR=libtt-gfx-1.5
DATE=`LC_TIME=C date '+%Y%m%d'`

SRC_FILES=\
	$(DIST_DIR) \

all:
	cd src && make
	cd test && make

clean:
	cd src && make clean
	cd test && make clean
	cd work && rm -f out*

dist-src:
	cd .. && tar czvf $(DIST_DIR)-src-$(DATE).tar.gz $(SRC_FILES)

