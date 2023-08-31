# Copyright (C) 2001-2020 Artifex Software, Inc.
# All Rights Reserved.
#
# This software is provided AS-IS with no warranty, either express or
# implied.
#
# This software is distributed under license and may not be copied,
# modified or distributed except as expressly authorized under the terms
# of the license contained in the file LICENSE in this distribution.
#
# Refer to licensing information at http://www.artifex.com or contact
# Artifex Software, Inc.,  1305 Grant Avenue - Suite 200, Novato,
# CA 94945, U.S.A., +1(415)492-9861, for further information.
#

# makefile for the lcms2mt library code.
# Users of this makefile must define the following:
#	SHARE_LCMS - whether to compile in or link to the library
#	LCMS2MTSRCDIR - the library source directory
#
# gs.mak and friends define the following:
#	LCMS2OBJDIR - the output obj directory
#	LCMS2GENDIR - generated (.dev) file directory
#	LCMS2I_ LCMS2_CFLAGS - include and cflags for compiling the lib

# We define the lcms2.dev target and its dependencies
#
# This partial makefile compiles the lcms2mt library for use in
# Ghostscript.


# Define the name of this makefile.
LCMS2_MAK=$(GLSRC)lcms2mt.mak $(TOP_MAKEFILES)

LCMS2SRC=$(LCMS2MTSRCDIR)$(D)src$(D)
LCMS2GEN=$(LCMS2GENDIR)$(D)
LCMS2OBJ=$(LCMS2OBJDIR)$(D)

# This makefile was stolen from the one for lcms2.

lcms2_OBJS=\
	$(LCMS2OBJ)cmscam02.$(OBJ) \
	$(LCMS2OBJ)cmscgats.$(OBJ) \
	$(LCMS2OBJ)cmscnvrt.$(OBJ) \
	$(LCMS2OBJ)cmserr.$(OBJ) \
	$(LCMS2OBJ)cmsgamma.$(OBJ) \
	$(LCMS2OBJ)cmsgmt.$(OBJ) \
	$(LCMS2OBJ)cmshalf.$(OBJ) \
	$(LCMS2OBJ)cmsintrp.$(OBJ) \
	$(LCMS2OBJ)cmsio0.$(OBJ) \
	$(LCMS2OBJ)cmsio1.$(OBJ) \
	$(LCMS2OBJ)cmslut.$(OBJ) \
	$(LCMS2OBJ)cmsmd5.$(OBJ) \
	$(LCMS2OBJ)cmsmtrx.$(OBJ) \
	$(LCMS2OBJ)cmsnamed.$(OBJ) \
	$(LCMS2OBJ)cmsopt.$(OBJ) \
	$(LCMS2OBJ)cmspack.$(OBJ) \
	$(LCMS2OBJ)cmspcs.$(OBJ) \
	$(LCMS2OBJ)cmsplugin.$(OBJ) \
	$(LCMS2OBJ)cmsps2.$(OBJ) \
	$(LCMS2OBJ)cmssamp.$(OBJ) \
	$(LCMS2OBJ)cmstypes.$(OBJ) \
	$(LCMS2OBJ)cmsvirt.$(OBJ) \
	$(LCMS2OBJ)cmswtpnt.$(OBJ) \
	$(LCMS2OBJ)cmsxform.$(OBJ) \
	$(LCMS2OBJ)cmsalpha.$(OBJ)

LCMS2_DEPS=\
        $(LCMS2MTSRCDIR)$(D)include$(D)lcms2mt.h \
	$(GLSRC)icc34.h $(LCMS2_MAK) $(MAKEDIRS)

lcms2.clean : lcms2.config-clean lcms2.clean-not-config-clean

lcms2.clean-not-config-clean :
	$(EXP)$(ECHOGS_XE) $(LCMS2MTSRCDIR) $(LCMS2OBJDIR)
	$(RM_) $(lcms2_OBJS)

lcms2.config-clean :
	$(RMN_) $(LCMS2GEN)$(D)lcms2mt*.dev

# NB: we can't use the normal $(CC_) here because msvccmd.mak
# adds /Za which conflicts with the lcms source.
LCMS2_CC=$(CC) $(CFLAGS_VISIBILITY) $(D_)SHARE_LCMS=$(SHARE_LCMS)$(_D) $(GENOPT) $(CAPOPT) $(CFLAGS) $(LCMS2_CFLAGS) $(I_)$(LCMS2MTSRCDIR)$(D)include $(LCMS2CF_)
LCMS2O_=$(O_)$(LCMS2OBJ)

# switch in the version of lcms2mt.dev we're actually using
$(LCMS2GEN)lcms2mt.dev : $(LCMS2GEN)lcms2mt_$(SHARE_LCMS).dev $(MAKEDIRS)
	$(CP_) $(LCMS2GEN)lcms2mt_$(SHARE_LCMS).dev $(LCMS2GEN)lcms2mt.dev

# dev file for shared (separately built) lcms library
$(LCMS2GEN)lcms2mt_1.dev : $(LCMS2_MAK) $(ECHOGS_XE) $(MAKEDIRS)
	$(SETMOD) $(LCMS2GEN)lcms2mt_1 -lib lcms2

# dev file for compiling our own from source
$(LCMS2GEN)lcms2mt_0.dev : $(LCMS2_MAK) $(ECHOGS_XE) $(lcms2_OBJS) $(LCMS2_DEPS)
	$(SETMOD) $(LCMS2GEN)lcms2mt_0 $(lcms2_OBJS)

# explicit rules for building the source files.

$(LCMS2OBJ)cmscam02.$(OBJ) : $(LCMS2SRC)cmscam02.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmscam02.$(OBJ) $(C_) $(LCMS2SRC)cmscam02.c

$(LCMS2OBJ)cmscgats.$(OBJ) : $(LCMS2SRC)cmscgats.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmscgats.$(OBJ) $(C_) $(LCMS2SRC)cmscgats.c

$(LCMS2OBJ)cmscnvrt.$(OBJ) : $(LCMS2SRC)cmscnvrt.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmscnvrt.$(OBJ) $(C_) $(LCMS2SRC)cmscnvrt.c

$(LCMS2OBJ)cmserr.$(OBJ) : $(LCMS2SRC)cmserr.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmserr.$(OBJ) $(C_) $(LCMS2SRC)cmserr.c

$(LCMS2OBJ)cmsgamma.$(OBJ) : $(LCMS2SRC)cmsgamma.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsgamma.$(OBJ) $(C_) $(LCMS2SRC)cmsgamma.c

$(LCMS2OBJ)cmsgmt.$(OBJ) : $(LCMS2SRC)cmsgmt.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsgmt.$(OBJ) $(C_) $(LCMS2SRC)cmsgmt.c

$(LCMS2OBJ)cmshalf.$(OBJ) : $(LCMS2SRC)cmshalf.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmshalf.$(OBJ) $(C_) $(LCMS2SRC)cmshalf.c

$(LCMS2OBJ)cmsintrp.$(OBJ) : $(LCMS2SRC)cmsintrp.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsintrp.$(OBJ) $(C_) $(LCMS2SRC)cmsintrp.c

$(LCMS2OBJ)cmsio0.$(OBJ) : $(LCMS2SRC)cmsio0.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsio0.$(OBJ) $(C_) $(LCMS2SRC)cmsio0.c

$(LCMS2OBJ)cmsio1.$(OBJ) : $(LCMS2SRC)cmsio1.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsio1.$(OBJ) $(C_) $(LCMS2SRC)cmsio1.c

$(LCMS2OBJ)cmslut.$(OBJ) : $(LCMS2SRC)cmslut.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmslut.$(OBJ) $(C_) $(LCMS2SRC)cmslut.c

$(LCMS2OBJ)cmsmd5.$(OBJ) : $(LCMS2SRC)cmsmd5.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsmd5.$(OBJ) $(C_) $(LCMS2SRC)cmsmd5.c

$(LCMS2OBJ)cmsmtrx.$(OBJ) : $(LCMS2SRC)cmsmtrx.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsmtrx.$(OBJ) $(C_) $(LCMS2SRC)cmsmtrx.c

$(LCMS2OBJ)cmsnamed.$(OBJ) : $(LCMS2SRC)cmsnamed.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsnamed.$(OBJ) $(C_) $(LCMS2SRC)cmsnamed.c

$(LCMS2OBJ)cmsopt.$(OBJ) : $(LCMS2SRC)cmsopt.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsopt.$(OBJ) $(C_) $(LCMS2SRC)cmsopt.c

$(LCMS2OBJ)cmspack.$(OBJ) : $(LCMS2SRC)cmspack.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmspack.$(OBJ) $(C_) $(LCMS2SRC)cmspack.c

$(LCMS2OBJ)cmspcs.$(OBJ) : $(LCMS2SRC)cmspcs.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmspcs.$(OBJ) $(C_) $(LCMS2SRC)cmspcs.c

$(LCMS2OBJ)cmsplugin.$(OBJ) : $(LCMS2SRC)cmsplugin.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsplugin.$(OBJ) $(C_) $(LCMS2SRC)cmsplugin.c

$(LCMS2OBJ)cmsps2.$(OBJ) : $(LCMS2SRC)cmsps2.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsps2.$(OBJ) $(C_) $(LCMS2SRC)cmsps2.c

$(LCMS2OBJ)cmssamp.$(OBJ) : $(LCMS2SRC)cmssamp.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmssamp.$(OBJ) $(C_) $(LCMS2SRC)cmssamp.c

$(LCMS2OBJ)cmstypes.$(OBJ) : $(LCMS2SRC)cmstypes.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmstypes.$(OBJ) $(C_) $(LCMS2SRC)cmstypes.c

$(LCMS2OBJ)cmswtpnt.$(OBJ) : $(LCMS2SRC)cmswtpnt.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmswtpnt.$(OBJ) $(C_) $(LCMS2SRC)cmswtpnt.c

$(LCMS2OBJ)cmsvirt.$(OBJ) : $(LCMS2SRC)cmsvirt.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsvirt.$(OBJ) $(C_) $(LCMS2SRC)cmsvirt.c

$(LCMS2OBJ)cmsxform.$(OBJ) : $(LCMS2SRC)cmsxform.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsxform.$(OBJ) $(C_) $(LCMS2SRC)cmsxform.c

$(LCMS2OBJ)cmsalpha.$(OBJ) : $(LCMS2SRC)cmsalpha.c $(LCMS2_DEPS)
	$(LCMS2_CC) $(LCMS2O_)cmsalpha.$(OBJ) $(C_) $(LCMS2SRC)cmsalpha.c
