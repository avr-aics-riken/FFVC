###
###########################################
# CPM - Cartesian Partition Manager       #
#                                         #
# Copyright (C) 2012 University of Tokyo. #
###########################################
###

SUBDIRS = @MAKE_SUB_DIRS@

bin_SCRIPTS=cpm-config
EXTRA_DIST = cpm-config.in cpm-uname

clean-local:
	find ./ -name *.obj | xargs rm -f
	find ./ -name *.exe | xargs rm -f
	find ./ -name *.lib | xargs rm -f

distclean-local: clean-local
	rm -f ./bin/*
	rm -f ./lib/*

