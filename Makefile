#/*
#* Copyright (C) 2012-2017 Jianxing Feng
#*
#* This program is free software; you can redistribute it and/or modify
#* it under the terms of the GNU General Public License as published by
#* the Free Software Foundation; either version 3 of the License, or (at
#* your option) any later version.
#*
#* This program is distributed in the hope that it will be useful, but
#* WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#* General Public License for more details.
#*
#* You should have received a copy of the GNU General Public License
#* along with this program; if not, write to the Free Software
#* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#*/

VERSION = V1.1.4
package=gfold.$(VERSION)
packagename=gfold.$(VERSION).tar.gz 

all: program

debug: DataProcessor.hpp GFOLD.hpp Utility.hpp GeneInfo.hpp main.cc
	g++ -Wall -g main.cc -o gfold -lgsl -lgslcblas 

program: DataProcessor.hpp GFOLD.hpp Utility.hpp GeneInfo.hpp main.cc
	g++ -O3 -Wall -g main.cc -o gfold -lgsl -lgslcblas 

docu: doc/gfold.pod
	pod2man doc/gfold.pod > doc/gfold.man
	pod2html -css gfold.css --noindex --header --title="GFOLD $(VERSION)" doc/gfold.pod > doc/gfold_doc.html
	head -n 12 doc/gfold_doc.html > doc/gfold.html
	echo '<tr><td width="10%" class="block" align="left"> <img align=center src="gfold.png"></img> </td>'>> doc/gfold.html
	echo '<td class="block" valign="middle">' >> doc/gfold.html
	awk '{if (NR > 13) print}' doc/gfold_doc.html >> doc/gfold.html
	rm doc/gfold_doc.html


dist: DataProcessor.hpp GFOLD.hpp Utility.hpp GeneInfo.hpp main.cc doc/gfold.pod Makefile README
	mkdir ${package}
	cp -r DataProcessor.hpp GFOLD.hpp Utility.hpp GeneInfo.hpp main.cc doc Makefile README ${package}
	tar cvzf ${packagename} ${package}
	rm -rf ${package}
	cp ${packagename} doc/gfold.html doc/gfold.css doc/gfold.png web/
	tar cvzf gfoldweb.tar.gz web
