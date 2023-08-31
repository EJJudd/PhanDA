#
#
# CUPS driver makefile for Ghostscript.
#
# Copyright 2001-2005 by Easy Software Products.
# Copyright 2007 Artifex Software, Inc.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#

# define the name of this makefile
CUPS_MAK=$(LCUPSSRCDIR)$(D)cups.mak $(TOP_MAKEFILES)

### ----------------- CUPS Ghostscript Driver ---------------------- ###

cups_=	$(GLOBJ)gdevcups.$(OBJ)

# These are set in the toplevel Makefile via autoconf(1)
# CUPSCFLAGS=`cups-config --cflags`
# CUPSSERVERBIN=`cups-config --serverbin`
# CUPSSERVERROOT=`cups-config --serverroot`
# CUPSDATA=`cups-config --datadir`
# CUPSPDFTORASTER= 1 if CUPS is new enough (cups-config --version)
# CUPSDIR

install:
