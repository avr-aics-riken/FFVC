###################################################################################
#
# FFV-C : Frontflow / violet Cartesian
#
# Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
# All rights reserved.
#
# Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2020 Research Institute for Information Technology(RIIT),
# Kyushu University.
# All rights reserved.
#
###################################################################################

# Install files

install(FILES
        2Dcavity/cavity.tp
        DESTINATION example/2Dcavity
)

install(FILES
        2Dcyl/cyl.tp
        DESTINATION example/2Dcyl
)

install(FILES
        3Dcavity/cavity.tp
        3Dcavity/comb.tp
        3Dcavity/layout.tp
        DESTINATION example/3Dcavity
)
