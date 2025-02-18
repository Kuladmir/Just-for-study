# Gromacs Test Summary

From: http://www.mdtutorials.com/gmx/complex/01_pdb2gmx.html
## 0.Gromacs says

### >  Efraim Racker

>Don't waste pure thoughts on dirty enzymes.不要把纯粹的想法浪费在肮脏的酶上。

#### 1913.6.28--1991.9.9  奥地利生物化学家
> ?-1938   维也纳大学  --> 药物方向 

> 1941-1942  明尼苏达大学  生物学副研

> 1944  纽约大学医学院  微生物学副教授  -->糖酵解

> 1952-1954  耶鲁医学院

> 1954-1966  纽约市公共卫生研究所营养与生理学系  主任

> 1966  康奈尔大学  -->化学生物学

#### 主要成就

>第一个氧化磷酸化酶的发现和纯化

>ATP合酶的第一部分F1的鉴定和纯化

### > Irene Joliot-Curie

>The farther the experiment is from theory, the closer it is to the Nobel Prize.实验离理论越远，则离诺贝尔奖越近。

#### 1897.9.12-1956.3.17  法国原子物理学家

> 1918  镭研究所  玛丽·居里（其母）的助手

> 1925.3  获学位  --> 关于钋的α粒子的研究

> 1934.1  人工制造出了放射能  --> 诺贝尔奖

## 1.准备拓扑文件

### 1.1  初始化蛋白的pdb文件

需要获取蛋白的pdb文件，并且尽可能消除其中的水分子、PO4(-3)、BME（B-巯基乙醇）等共结晶溶剂等。
但是不是所有情况下都需要如此处理。

### 1.2  准备蛋白的拓扑文件

**1）获得小分子信息**

先将蛋白中的小分子信息提取出来，存入jz4.pdb文件。之后删除3HTB_clean.pdb中的小分子部分
```
# 案例
grep JZ4 3HTB_clean.pdb > jz4.pdb
```
**2）获得立场**

From: https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs

下载最新的CHARMM36立场文件（charmm36-date.ff.tgz）；

下载对应python版本的python文件（cgenff_charmm2gmx.py）

**3）解压**

解压后会产生一个charmm36-jul2022.ff的文件夹
```
mkdir working_file $$ cd working_file

tar -zxvf charmm36-jul2022.ff.tgz
```

**4）拓扑转化**

```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx pdb2gmx -f 3HTB_clean.pdb -o 3HTB_processed.gro -ter
```
之后会出现一些功能选择，执行完该命令后，会出现三个功能选择（其实准确来说是四个）：Force field  -->  Water model  -->  Terminus type(Start terminus MET-1 and End teminus ASN-163)

>Force field --> No.1

1: CHARMM all-atom force field
From '/pubhome/soft/gromacs/2020.5/thread_mpi/share/gromacs/top':

>Water model --> No.1

1: TIP3P      CHARMM-modified TIP3P water model (recommended over original TIP3P)

>Start terminus MET-1 --> No.1

1: NH3+

>End teminus ASN-163 --> No.0

0: COO-

***<文中解释：由于3HTB这个蛋白以MET为开头，这种相互作用力是必要的，会导致pdb2gmx功能选择一种不相容的碳水化合物末端类型>***

### 1.3  准备配体的拓扑文件

**1）加氢**

由于charmm是全原子立场，但是一般配体中会隐藏H，因此需要加氢。同时将pdb文件转化为mol2文件。

使用Avogrado完成这个工作。From: https://sourceforge.net/projects/avogadro/

***<使用Avogrado打开jz4.pdb，选择Build --> Add Hydrogen，选择File --> Save as --> .mol2 >***

**2）预处理**

使用记事本/Notepad++/etc等文本工具打开.mol2。前两行会显示这样：
>@<TRIPOS>MOLECULE
>### *****  
修改成这样：
>@<TRIPOS>MOLECULE
>
>JZ4

然后在ATOM部分，会显示这样（部分），应该给它进行修改：
>10 OAB        23.4120  -23.5360   -4.3420 O.3   167  JZ4167      -0.5065

>11 H          25.3133  -24.3619    0.1509 H       1  UNL1         0.0230

修改成这样：
>10 OAB        23.4120  -23.5360   -4.3420 O.3     1  JZ4     -0.5065

>11 H          25.3133  -24.3619    0.1509 H       1  JZ4      0.0230 

### **注意列对齐！！！**

***<如果键没有按升序列出，则在构建具有匹配坐标的正确拓扑时会出现问题。>***

创建一个sort_mol2_bonds.pl文件：
```
#!/usr/bin/perl

use strict;

# sort_mol2_bonds.pl - a script to reorder the listing in a .mol2 @<TRIPOS>BOND
# section so that the following conventions are preserved:
#   1. Atoms on each line are in increasing order (e.g. 1 2 not 2 1)
#   2. The bonds appear in order of ascending atom number
#   3. For bonds involving the same atom in the first position, the bonds appear
#       in order of ascending second atom
#
# Written by: Justin Lemkul (jalemkul@vt.edu)
#
# Distributed under the GPL-3.0 license

unless (scalar(@ARGV)==2)
{
    die "Usage: perl sort_mol2_bonds.pl input.mol2 output.mol2\n";
}

my $input = $ARGV[0];
my $output = $ARGV[1];

open(IN, "<$input") || die "Cannot open $input: $!\n";
my @in = <IN>;
close(IN);

# test for header lines that some scripts produce
unless($in[0] =~ /TRIPOS/)
{
    die "Nonstandard header found: $in[0]. Please delete header lines until the TRIPOS molecule definition.\n"; 
}

open(OUT, ">$output") || die "Cannot open $output: $!\n";

# get number of atoms and number of bonds from mol2 file
my @tmp = split(" ", $in[2]);
my $natom = $tmp[0];
my $nbond = $tmp[1];

# check
print "Found $natom atoms in the molecule, with $nbond bonds.\n";

# print out everything up until the bond section
my $i=0;
while (!($in[$i] =~ /BOND/))
{
    print OUT $in[$i];
    $i++;
}

# print the bond section header line to output
print OUT $in[$i];
$i++;

# read in the bonds and sort them
my $bondfmt = "%6d%6d%6d%5s\n";
my @tmparray; 

# sort the bonds - e.g. the one that has the
# lowest atom number in the first position and then the
# lowest atom number in the second position (swap if necessary)
for (my $j=0; $j<$nbond; $j++)
{
    my @tmp = split(" ", $in[$i+$j]);
    # parse atom numbers
    my $ai = $tmp[1];
    my $aj = $tmp[2];
    # reorder if second atom number < first
    if ($aj < $ai)
    {
        $ai = $tmp[2];
        $aj = $tmp[1];
    }
    # store new lines in a temporary array
    $tmparray[$j] = sprintf($bondfmt, $tmp[0], $ai, $aj, $tmp[3]); 
}

# loop over tmparray to find each atom number
my $nbond = 0;
for (my $x=1; $x<=$natom; $x++)
{
    my @bondarray;
    my $ntmp = scalar(@tmparray);
    for (my $b=0; $b<$ntmp; $b++)
    {
        my @tmp = split(" ", $tmparray[$b]);
        if ($tmp[1] == $x)
        {
            push(@bondarray, $tmparray[$b]);
            splice(@tmparray, $b, 1);
            $ntmp--;
            $b--;
        }
    }

    if (scalar(@bondarray) > 0) # some atoms will only appear in $aj, not $ai
    {
        my $nbondarray = scalar(@bondarray);
        if ($nbondarray > 1)
        {
            # loop over all bonds, find the one with lowest $aj
            # and then print it
            for (my $y=0; $y<$nbondarray; $y++)
            {
                my @tmp2 = split(" ", $bondarray[$y]);
                my $tmpatom = $tmp[2];
                my $lowindex = 0;
                if ($tmp2[2] < $tmpatom)
                {
                    $lowindex = $y; 
                }
                my $keep = splice(@bondarray, $lowindex, 1);
                $y--;
                $nbondarray--;
                my @sorted = split(" ", $keep);
                $nbond++;
                printf OUT $bondfmt, $nbond, $sorted[1], $sorted[2], $sorted[3]; 
            }
        }
        else
        {
            $nbond++;
            my @tmp2 = split(" ", $bondarray[0]);
            printf OUT $bondfmt, $nbond, $tmp2[1], $tmp2[2], $tmp2[3];
        }
    }
}

close(OUT);

exit;

```
然后执行以下命令：
```
perl sort_mol2_bonds.pl jz4.mol2 jz4_fix.mol2
```
**3）拓扑转化**

使用CGenFF将mol2文件转化。
```
/pubhome/soft/silcsbio.2024.1/cgenff/cgenff -f jz4.mol2 > jz4.str
```

**4）文件再处理**

强烈建议使用conda创建一个环境，并且使用python=2.x的版本，和networkx=1.11的版本进行处理，避免报错。
```
conda create -n gromacs python=2.7

conda activate gormacs

pip install numpy

pip install networkx==1.11

python cgenff_charmm2gmx_py2.py JZ4 jz4_fix.mol2 jz4.str charmm36-jul2022.ff

conda deactivate
```
cgenff_charmm2gmx_py2.py文件如下：
```
#!/usr/bin/env python

# USAGE: python cgenff_charmm2gmx.py DRUG drug.mol2 drug.str charmm36.ff
# Tested with Python 2.7.3 and 2.7.12. Requires numpy and networkx
# The networkx version MUST be in the 1.x series. Tested version: 1.11

# Copyright (C) 2014 E. Prabhu Raman prabhu@outerbanks.umaryland.edu
#
# Modified 11/6/2018 by Justin Lemkul to add lone pair support
# needed for CGenFF >= 4.0 halogens
#
# For help/questions/bug reports, please contact Justin Lemkul jalemkul@vt.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU Affero General Public License for more details.
# <http://www.gnu.org/licenses/>

# EXAMPLE: You have a drug-like molecule in drug.mol2 file
# ParamChem returns a CHARMM stream file drug.str with topology and parameters
# INPUT
# The program needs four inputs:
#	(1) The first argument (resname) is found in the RESI entry in the CHARMM stream file; for example
#		RESI DRUG		  0.000  ! here DRUG is the resname
#	(2) drug.mol2 is the .mol2 which you supplied to ParamChem
#	(3) drug.str contains the topology entry and CGenFF parameters that ParamChem generates for you
#	(4) charmm36.ff should contain the CHARMM force-field in GROMACS format
#		Download it from: http://mackerell.umaryland.edu/CHARMM_ff_params.html

# OUTPUT
# The program will generate 4 output files ("DRUG" is converted to lowercase and the files are named accordingly):
#	(1) drug.itp - contains GROMACS itp
#	(2) drug.prm - contains parameters obtained from drug.str which are converted to GROMACS format and units
#	(3) drug.top - A Gromacs topology file which incorporates (1) and (2)
#	(4) drug_ini.pdb - Coordinates of the molecule obtained from drug.mol2

# The program has been tested only on CHARMM stream files containing topology and parameters of a single molecule.

import string
import re
import sys
import os
import math
import numpy as np
import networkx as nx

#=================================================================================================================
def check_versions(str_filename,ffdoc_filename):
	ffver = 0	# CGenFF version in force field directory
	strver = 0	# CGenFF version in stream file
	f = open(str_filename, 'r')
	for line in f.readlines():
		if line.startswith("* For use with CGenFF version"):
			entry = re.split('\s+', string.lstrip(line))
			strver = entry[6]
			print "--Version of CGenFF detected in ",str_filename,":",strver
	f.close()
	f = open(ffdoc_filename, 'r')
	for line in f.readlines():
		if line.startswith("Parameters taken from CHARMM36 and CGenFF"):
			entry = re.split('\s+', string.lstrip(line))
			ffver = entry[6]
			print "--Version of CGenFF detected in ",ffdoc_filename,":",ffver
	f.close()

	# warn the user about version mismatch 
	if strver != ffver:
		print "\nWARNING: CGenFF versions are not equivalent!\n"

	# in case something has gone horribly wrong
	if (strver == 0) or (ffver == 0):
		print "\nERROR: Could not detect CGenFF version. Exiting.\n"
		exit()

#-----------------------------------------------------------------------
## jal
def is_lp(s):
	if ((s[0]=='L') and (s[1]=='P')):
		return True
	return False

#-----------------------------------------------------------------------
## jal
def is_lp_host_atom(self,name):
	for ai in range (0,self.nvsites):
		if (name==self.G.node[ai]['at1']):
			return True
	return False

#-----------------------------------------------------------------------
## jal - only for COLINEAR lone pairs, since CGenFF only needs this now
def construct_lp(x1,y1,z1,x2,y2,z2,dist):
	dx = x1-x2
	dy = y1-y2
	dz = z1-z2
	dr = math.sqrt(dx*dx+dy*dy+dz*dz)
	dr = dist/dr
	# set LP coords
	xlp = x1+dr*dx
	ylp = y1+dr*dy
	zlp = z1+dr*dz

	return xlp,ylp,zlp

#-----------------------------------------------------------------------
## jal
def find_vsite(self, atnum):
	for i in range (0, self.nvsites): 
		# if we find the LP host, find the LP atom index
		if (self.G.node[i]['at1'] == self.G.node[atnum]['name']):
			for j in range (0, self.natoms):
				if (self.G.node[i]['vsite'] == self.G.node[j]['name']):
					return j 

#-----------------------------------------------------------------------
def read_gmx_atomtypes(filename):
	atomtypes = []
	f = open(filename, 'r')
	for line in f.readlines():
		if line.startswith(";"):
			continue
		if line == '\n':
			continue
		entry = re.split('\s+', string.lstrip(line))
		var = [entry[0],entry[1]]
		atomtypes.append(var)
	f.close()
	return atomtypes
#-----------------------------------------------------------------------
def get_filelist_from_gmx_forcefielditp(ffdir,ffparentfile):
	filelist=[]
	f = open(ffdir+"/"+ffparentfile, 'r')
	for line in f.readlines():
		if line.startswith("#include"):
			entry = re.split('\s+', string.lstrip(line))
			filename = ffdir + "/" + entry[1].replace("\"","")
			filelist.append(filename)
	return filelist
#-----------------------------------------------------------------------
def read_gmx_anglpars(filename):

	angllines = []
	f = open(filename, 'r')
	section="NONE"
	for line in f.readlines():
		if line.startswith(";"):
			continue
		if line.startswith("\n"):
			continue
		if line.startswith("["):
			section="NONE"
		if(section=="ANGL"):
			angllines.append(line)
		if line.startswith("[ angletypes ]"):
			section="ANGL"

	anglpars = []
	anglpar = {}
	for line in angllines:
		entry = re.split('\s+', string.lstrip(line))
		ai, aj, ak, eq = entry[0],entry[1],entry[2],float(entry[4])
		anglpars.append([ai,aj,ak,eq])

	return anglpars
#-----------------------------------------------------------------------
def get_charmm_rtp_lines(filename,molname):
	foundmol=0
	store=0
	rtplines=[]
	f = open(filename, 'r')
	section="NONE"
	for line in f.readlines():
		if(store==1) and line.startswith("RESI"):
			store=0

		if line.startswith("RESI"):
			entry = re.split('\s+', string.lstrip(line))
			rtfmolname=entry[1]
			if(rtfmolname == molname):
				store=1

		if line.startswith("END"):
			store=0

		if(store==1):
			rtplines.append(line)

	return rtplines
#-----------------------------------------------------------------------
def get_charmm_prm_lines(filename):
	foundmol=0
	store=0
	prmlines=[]
	f = open(filename, 'r')
	section="NONE"
	for line in f.readlines():

		if line.startswith("END"):
			section="NONE"
			store=0

		if(store):
			prmlines.append(line)

		if line.startswith("read para"):
			section="PRM"
			store=1


	return prmlines
#-----------------------------------------------------------------------
def parse_charmm_topology(rtplines):
	topology = {}
	noblanks = filter(lambda x: len(x.strip())>0, rtplines)
	nocomments = filter(lambda x: x.strip()[0] not in ['*','!'], noblanks)
	section = "BONDS"	# default
	state = "free"
	for line in nocomments:
		if state == "free":
			if line.find("MASS") == 0:
				if "ATOMS" not in topology.keys():
					topology["ATOMS"] = {}
				s = line.split()
				idx,name,mass,type = int(s[1]),s[2],float(s[3]),s[4]
				if line.find("!"):
					comment = line[line.find("!")+1:].strip()
				else:
					comment = ""
				topology["ATOMS"][name] = [idx,mass,type,comment]
			elif line.find("DECL") == 0:
				if "DECL" not in topology.keys():
					topology["DECL"] = []
				decl = line.split()[1]
				topology["DECL"].append(decl)
			elif line.find("DEFA") == 0:
				topology["DEFA"] = line[4:]
			elif line.find("AUTO") == 0:
				topology["AUTO"] = line[4:]
			elif line.find("RESI") == 0:
				if "RESI" not in topology.keys():
					topology["RESI"] = {}
				state = "resi"
				s = line.split()
				resname, charge = s[1],float(s[2])
				topology["RESI"][resname] = {}
				topology["RESI"][resname]["charge"] = charge
				topology["RESI"][resname]["cmaps"] = []
				topology["RESI"][resname]["vsites"] = []
				topology["RESI"][resname]["bonds"] = []
				topology["RESI"][resname]["impropers"] = []
				topology["RESI"][resname]["double_bonds"] = []
				group = -1 
			elif line.find("PRES") == 0:
				state = "pres"
				s = line.split()
				presname, charge = s[1],float(s[2])
			elif line.find("END") == 0:
				return topology
		elif state == "resi":
			if line.find("RESI")==0:
				state = "resi"
				s = line.split()
				resname, charge = s[1],float(s[2])
				topology["RESI"][resname] = {}
				topology["RESI"][resname]["charge"] = charge
				topology["RESI"][resname]["cmaps"] = []
				topology["RESI"][resname]["vsites"] = []
				topology["RESI"][resname]["bonds"] = []
				topology["RESI"][resname]["impropers"] = []
				topology["RESI"][resname]["double_bonds"] = []
				#topology["RESI"][resname]["groups"] = []
				group = -1 

			elif line.find("GROU")==0:
				group += 1
				topology["RESI"][resname][group] = []
			elif line.find("ATOM")==0: 
				if line.find('!'):
					line = line[:line.find('!')]
				s = line.split()
				name,type,charge = s[1],s[2],float(s[3])
				topology["RESI"][resname][group].append((name,type,charge))
			## jal - adding lone pair support
			elif line.find("LONE")==0:
				if line.find('!'):
					line = line[:line.find('!')]
				s = line.split()
				name,at1,at2,dist = s[2],s[3],s[4],(float(s[6])*0.1)
				topology["RESI"][resname]["vsites"].append((name,at1,at2,dist))
			elif line.find("BOND")==0: 
				if line.find('!'):
					line = line[:line.find('!')]
				s = line.split()
				nbond = (len(s)-1)/2
				for i in range(nbond):
					p,q = s[1+2*i],s[2+2*i]
					## jal - ignore "bonds" to lone pairs
					if ((is_lp(p)==False) and (is_lp(q)==False)):
						topology["RESI"][resname]["bonds"].append((p,q))
			elif line.find("DOUB")==0: 
				if line.find('!'):
					line = line[:line.find('!')]
				s = line.split()
				ndouble = (len(s)-1)/2
				for i in range(ndouble):
					p,q = s[1+2*i],s[2+2*i]
					topology["RESI"][resname]["double_bonds"].append((p,q))
			elif line.find("IMPR")==0: 
				if line.find('!'):
					line = line[:line.find('!')]
				s = line.split()
				nimproper = (len(s)-1)/4
				for i in range(nimproper):
					impr = s[1+4*i],s[2+4*i],s[3+4*i],s[4+4*i]
					topology["RESI"][resname]["impropers"].append(impr)
			elif line.find("CMAP")==0: 
				if line.find('!'):
					line = line[:line.find('!')]
				s = line.split()
				#nimproper = (len(s)-1)/4
				#for i in range(nimproper):
				cmap = s[1:9]
				topology["RESI"][resname]["cmaps"].append(cmap)
			elif line.find("DONOR")==0: 
				continue	# ignore for now
			elif line.find("ACCEPTOR")==0: 
				continue
			elif line.find("IC")==0:
				continue

	return topology
#-----------------------------------------------------------------------
def parse_charmm_parameters(prmlines):

	parameters = {}
	cmapkey = ()
	noblanks = filter(lambda x: len(x.strip())>0, prmlines)
	nocomments = filter(lambda x: x.strip()[0] not in ['*','!'], noblanks)
	section = "ATOM"	# default
	for line in nocomments:
				#print line
		sectionkeys = [ "BOND", "ANGL", "DIHE", \
				"IMPR", "CMAP", "NONB", "HBON", "NBFI" ]
		key = line.split()[0]

				#exit()

		if key[0:4] in sectionkeys:
			section = key[0:4]
			continue

		if section not in parameters.keys():
			parameters[section] = []

				#print line
		if section == "BOND":
			if line.find('!'):
				line = line[:line.find('!')]
			s = line.split()
			ai, aj, kij, rij = s[0],s[1],float(s[2]),float(s[3])
			parameters["BOND"].append((ai,aj,kij,rij))
		elif section == "ANGL":
			if line.find('!'):
				line = line[:line.find('!')]
			s = line.split()
			ai, aj, ak = s[0],s[1],s[2]
			other = map(float,s[3:])
			parameters["ANGL"].append([ai,aj,ak]+other)
		elif section == "DIHE":
			if line.find('!'):
				line = line[:line.find('!')]
			s = line.split()
			ai, aj, ak, al, k, n, d = s[0],s[1],s[2],s[3],float(s[4]),int(s[5]),float(s[6])
			parameters["DIHE"].append([ai,aj,ak,al,k,n,d])
		elif section == "IMPR":
			if line.find('!'):
				line = line[:line.find('!')]
			s = line.split()
			ai, aj, ak, al, k, d = s[0],s[1],s[2],s[3],float(s[4]),float(s[6])
			parameters["IMPR"].append([ai,aj,ak,al,k,d])
		elif section == "CMAP":
			if line.find('!'):
				line = line[:line.find('!')]
			if cmapkey == ():
				s = line.split()
				a,b,c,d,e,f,g,h = s[0:8]
				N = int(s[8])
				cmapkey = (a,b,c,d,e,f,g,h,N)
				cmaplist = []
			else:
				cmapdata = map(float,line.split())
				cmaplist += cmapdata
				if len(cmaplist) == N**2:
					parameters["CMAP"].append([cmapkey,cmaplist])
					cmapkey = ()
		elif section == "NONB":
			if line.find("cutnb")>=0 or line.find("wmin")>=0 or line.find("CUTNB")>=0 or line.find("WMIN")>=0 :
				continue
			bang = line.find('!')
			if bang>0:
				comment = line[bang+1:]
				prm = line[:bang].split()
			else:
				comment = ""
				prm = line.split()
			atname = prm[0]
			epsilon = -float(prm[2])
			half_rmin = float(prm[3])
			parameters["NONB"].append((atname,epsilon,half_rmin))

			if len(prm)>4:
				epsilon14 = -float(prm[5])
				half_rmin14 = float(prm[6])
				if "NONBONDED14" not in parameters.keys():
					parameters["NONBONDED14"] = []
				parameters["NONBONDED14"].append((atname,epsilon14,half_rmin14))


	return parameters
#-----------------------------------------------------------------------
def write_gmx_bon(parameters,header_comments,filename):
	kcal2kJ = 4.18400

	outp = open(filename,"w")
	outp.write("%s\n"%(header_comments))
	outp.write("[ bondtypes ]\n")
	kbond_conversion = 2.0*kcal2kJ/(0.1)**2		# [kcal/mol]/A**2 -> [kJ/mol]/nm**2
						# factor of 0.5 because charmm bonds are Eb(r)=Kb*(r-r0)**2
	rbond_conversion = .1			# A -> nm
	outp.write(";%7s %8s %5s %12s %12s\n"%("i","j","func","b0","kb"))
	if(parameters.has_key("BOND")):
		for p in parameters["BOND"]:
			ai,aj,kij,rij = p
			rij *= rbond_conversion
			kij *= kbond_conversion
			outp.write("%8s %8s %5i %12.8f %12.2f\n"%(ai,aj,1,rij,kij))

	kangle_conversion = 2.0*kcal2kJ		# [kcal/mol]/rad**2 -> [kJ/mol]/rad**2
										# factor of 0.5 because charmm angles are Ea(r)=Ka*(a-a0)**2
	kub_conversion = 2.0*kcal2kJ/(0.1)**2		# [kcal/mol]/A**2 -> [kJ/mol]/nm**2prm
	ub0_conversion = 0.1						# A -> nm

	outp.write("\n\n[ angletypes ]\n")
	outp.write(";%7s %8s %8s %5s %12s %12s %12s %12s\n"\
			%("i","j","k","func","theta0","ktheta","ub0","kub"))
	if(parameters.has_key("ANGL")):
		for p in parameters["ANGL"]:
			if len(p) == 5:
				ai,aj,ak,kijk,theta = p
				kub = 0.0
				ub0 = 0.0
			else:
				ai,aj,ak,kijk,theta,kub,ub0 = p
		
			kijk *= kangle_conversion
			kub *= kub_conversion
			ub0 *= ub0_conversion
			outp.write("%8s %8s %8s %5i %12.6f %12.6f %12.8f %12.2f\n"\
					%(ai,aj,ak,5,theta,kijk,ub0,kub))

	kdihe_conversion = kcal2kJ
	outp.write("\n\n[ dihedraltypes ]\n")
	outp.write(";%7s %8s %8s %8s %5s %12s %12s %5s\n"\
			%("i","j","k","l","func","phi0","kphi","mult"))
	#parameters["DIHEDRALS"].sort(demote_wildcards)
	if(parameters.has_key("DIHE")):
		for p in parameters["DIHE"]:
			ai,aj,ak,al,k,n,d = p
			k *= kdihe_conversion
			outp.write("%8s %8s %8s %8s %5i %12.6f %12.6f %5i\n"\
					%(ai,aj,ak,al,9,d,k,n))

	kimpr_conversion = kcal2kJ*2	# see above
	outp.write("\n\n[ dihedraltypes ]\n")
	outp.write("; 'improper' dihedrals \n")
	outp.write(";%7s %8s %8s %8s %5s %12s %12s\n"\
			%("i","j","k","l","func","phi0","kphi"))
	if(parameters.has_key("IMPR")):
	#parameters["IMPROPERS"].sort(demote_wildcards)
		for p in parameters["IMPR"]:
			ai,aj,ak,al,k,d = p
			k *= kimpr_conversion
			outp.write("%8s %8s %8s %8s %5i %12.6f %12.6f\n"\
					%(ai,aj,ak,al,2,d,k))

	outp.close()
	return
#-----------------------------------------------------------------------
def write_gmx_mol_top(filename,ffdir,prmfile,itpfile,molname):
	outp = open(filename,"w")
	outp.write("#include \"%s/forcefield.itp\"\n" % (ffdir))
	outp.write("\n")
	outp.write("; additional params for the molecule\n")
	outp.write("#include \"%s\"\n" % (prmfile))
	outp.write("\n")
	outp.write("#include \"%s\"\n" % (itpfile))
	outp.write("\n")
	outp.write("#include \"%s/tip3p.itp\"\n" % (ffdir))
	outp.write("#ifdef POSRES_WATER\n")
	outp.write("; Position restraint for each water oxygen\n")
	outp.write("[ position_restraints ]\n")
	outp.write(";  i funct		 fcx		fcy		   fcz\n")
	outp.write("   1	1		1000	   1000		  1000\n")
	outp.write("#endif\n")
	outp.write("\n")
	outp.write("; Include topology for ions\n")
	outp.write("#include \"%s/ions.itp\"\n" % (ffdir))
	outp.write("\n")
	outp.write("[ system ]\n")
	outp.write("; Name\n")
	outp.write("mol\n")
	outp.write("\n")
	outp.write("[ molecules ]\n")
	outp.write("; Compound		  #mols\n")
	outp.write("%s			1\n" % (molname))
	outp.write("\n")

	outp.close()
#=================================================================================================================
class atomgroup:
	"""
	A class that contains the data structures and functions to store and process
	data related to groups of atoms (read molecules)

	USAGE: m = atomgroup()
	"""

	def __init__(self):
		self.G = nx.Graph()
		self.name = ""
		self.natoms = 0
		self.nvsites = 0
		self.nbonds = 0
		self.angles = []
		self.nangles = 0
		self.dihedrals = []
		self.ndihedrals = 0
		self.impropers = []
		self.nimpropers = 0
		#self.coord=np.zeros((self.natoms,3),dtype=float)

	#-----------------------------------------------------------------------
	def read_charmm_rtp(self,rtplines,atomtypes):
		"""
		Reads CHARMM rtp
		Reads atoms, bonds, impropers
		Stores connectivity as a graph
		Autogenerates angles and dihedrals

		USAGE: m = atomgroup() ; m.read_charmm_rtp(rtplines,atomtypes)

		"""
		#initialize everything
		self.G = nx.Graph()
		self.name = ""
		self.natoms = 0
		self.nvsites = 0
		self.nbonds = 0
		self.angles = []
		self.nangles = 0
		self.dihedrals = []
		self.ndihedrals = 0
		self.impropers = []
		self.nimpropers = 0

		atm = {}

		for line in rtplines:
			if line.find('!'):
				line = line[:line.find('!')]

				if line.startswith("RESI"):
					entry = re.split('\s+', string.lstrip(line))
					self.name=entry[1]

				if line.startswith("ATOM"):
					entry = re.split('\s+', string.lstrip(line))
					atm[self.natoms] = {'type':entry[2], 'resname':self.name, 'name':entry[1],
						  'charge':float(entry[3]),'mass':float(0.00), 'beta':float(0.0),
							'x':float(9999.9999),'y':float(9999.9999),'z':float(9999.9999),'segid':self.name, 'resid':'1' }

					for typei in atomtypes:
						if(typei[0] == atm[self.natoms]['type']):
							atm[self.natoms]['mass'] = float(typei[1])
							break

					self.G.add_node(self.natoms, atm[self.natoms])
					self.natoms=self.natoms+1

				## jal - adding lone pair support
				if line.startswith("LONE"):
					entry = re.split('\s+', string.rstrip(string.lstrip(line)))
					atm[self.nvsites] = {'vsite':entry[2], 'at1':entry[3], 'at2':entry[4],
							'dist':(float(entry[6])*0.1) }
					# DEBUG
					# print "Found lone pair in RTF: %s %s %s %.3f\n" % (atm[self.nvsites]['vsite'], atm[self.nvsites]['at1'], atm[self.nvsites]['at2'], atm[self.nvsites]['dist'])

					self.G.add_node(self.nvsites, atm[self.nvsites])
					self.nvsites=self.nvsites+1

				if line.startswith("BOND") or line.startswith("DOUB"):
					entry = re.split('\s+', string.rstrip(string.lstrip(line)))
					numbonds = int((len(entry)-1)/2)
					for bondi in range(0,numbonds):
						found1 = False
						found2 = False
						for i in range(0,self.natoms):
							if(atm[i]['name'] == entry[(bondi*2)+1]):
								found1 = True
								break
						for j in range(0,self.natoms):
							if(atm[j]['name'] == entry[(bondi*2)+2]):
								found2 = True
								break
						if(not found1):
							print "Error:atomgroup:read_charmm_rtp> Atomname not found in top",entry[(bondi*2)+1]
						if(not found2):
							print "Error:atomgroup:read_charmm_rtp> Atomname not found in top",entry[(bondi*2)+2]
						## jal - ignore "bonds" to lone pairs
						if ((is_lp(atm[i]['name'])==False) and (is_lp(atm[j]['name'])==False)):
							self.G.add_edge(i,j)
							self.G[i][j]['order']='1' # treat all bonds as single for now
							self.nbonds=self.nbonds+1

				if line.startswith("IMP"):
					entry = re.split('\s+', string.lstrip(line))
					numimpr = int((len(entry)-2)/4)
					for impi in range(0,numimpr):
						for i in range(0,self.natoms):
							if(atm[i]['name'] == entry[(impi*4)+1]):
								break
						for j in range(0,self.natoms):
							if(atm[j]['name'] == entry[(impi*4)+2]):
								break
						for k in range(0,self.natoms):
							if(atm[k]['name'] == entry[(impi*4)+3]):
								break
						for l in range(0,self.natoms):
							if(atm[l]['name'] == entry[(impi*4)+4]):
								break
						var = [i,j,k,l]
						self.impropers.append(var)

		self.nimpropers = len(self.impropers)
		if(self.ndihedrals > 0 or self.nangles > 0):
			print "WARNING:atomgroup:read_charmm_rtp> Autogenerating angl-dihe even though they are preexisting",self.nangles,self.ndihedrals
		self.autogen_angl_dihe()
		self.coord = np.zeros((self.natoms,3),dtype=float)
#-----------------------------------------------------------------------
	def autogen_angl_dihe(self):
		self.angles = []
		for atomi in range(0,self.natoms):
			nblist = []
			for nb in self.G.neighbors(atomi):
				nblist.append(nb)
			for i in range(0,len(nblist)-1):
				for j in range(i+1,len(nblist)):
					var = [nblist[i],atomi,nblist[j]]
					self.angles.append(var)
		self.nangles = len(self.angles)
		self.dihedrals = []
		for i,j in self.G.edges_iter():
			nblist1 = []
			for nb in self.G.neighbors(i):
				if(nb != j):
					nblist1.append(nb)
			nblist2 = []
			for nb in self.G.neighbors(j):
				if(nb != i):
					nblist2.append(nb)
			if(len(nblist1) > 0 and len(nblist2) > 0 ):
				for ii in range(0,len(nblist1)):
					for jj in range(0,len(nblist2)):
						var = [nblist1[ii],i,j,nblist2[jj]]
						if(var[0] != var[3]):
							self.dihedrals.append(var)
		self.ndihedrals = len(self.dihedrals)
#-----------------------------------------------------------------------
	def get_nonplanar_dihedrals(self,angl_params):
		nonplanar_dihedrals=[]
		cutoff=179.9
		for var in self.dihedrals:
			d1=self.G.node[var[0]]['type']			  
			d2=self.G.node[var[1]]['type']
			d3=self.G.node[var[2]]['type']
			d4=self.G.node[var[3]]['type']
			keep=1
			for angl_param in angl_params:
				p1=angl_param[0]
				p2=angl_param[1]
				p3=angl_param[2]
				eq=angl_param[3]
				if( d2==p2 and ( ( d1==p1 and d3==p3) or (d1==p3 and d3==p1))):
					if(eq > cutoff):
						keep=-1
						break
				if( d3==p2 and ( (d2==p1 and d4==p3) or (d2==p3 and d4==p1))):
					if(eq > cutoff):
						keep=-1
						break

				 
			if(keep==1):
				nonplanar_dihedrals.append(var)

		return nonplanar_dihedrals
#-----------------------------------------------------------------------
	def write_gmx_itp(self,filename,angl_params):
		f = open(filename, 'w')
		f.write("; Created by cgenff_charmm2gmx.py\n")
		f.write("\n")
		f.write("[ moleculetype ]\n")
		f.write("; Name			   nrexcl\n")
		f.write("%s				 3\n" % self.name)
		f.write("\n")
		f.write("[ atoms ]\n")
		f.write(";	 nr		  type	resnr residue  atom   cgnr	   charge		mass  typeB    chargeB		massB\n")
		f.write("; residue	 1 %s rtp %s q	qsum\n" % (self.name,self.name))
		pairs14 = nx.Graph()
		for atomi in range(0,self.natoms):
			pairs14.add_node(atomi)
			f.write("%6d %10s %6s %6s %6s %6d %10.3f %10.3f   ;\n" % 
			   ( atomi+1,self.G.node[atomi]['type'],
			   self.G.node[atomi]['resid'],self.name,self.G.node[atomi]['name'],atomi+1,
			   self.G.node[atomi]['charge'],self.G.node[atomi]['mass'] ) )
		f.write("\n")
		f.write("[ bonds ]\n")
		f.write(";	ai	  aj funct			  c0			c1			  c2			c3\n")
		for i,j in self.G.edges_iter():
			f.write("%5d %5d	 1\n" % (i+1,j+1) )
		f.write("\n")
		f.write("[ pairs ]\n")
		f.write(";	ai	  aj funct			  c0			c1			  c2			c3\n")
		for var in self.dihedrals:
			if (len(nx.dijkstra_path(self.G,var[0],var[3])) == 4): #this is to remove 1-2 and 1-3 included in dihedrals of rings
				pairs14.add_edge(var[0],var[3])
		for i,j in pairs14.edges_iter():
			f.write("%5d %5d	 1\n" % (i+1,j+1) )
			## jal - add LP pairs, same as parent atom
			## Use is_lp_host_atom() to test each index, then find associated vsite
			if ((is_lp_host_atom(self,self.G.node[i]['name'])==True)):
				k = find_vsite(self, i)
				f.write("%5d %5d	 1\n" % (k+1,j+1) )
			if ((is_lp_host_atom(self,self.G.node[j]['name'])==True)):
				k = find_vsite(self, j)
				f.write("%5d %5d	 1\n" % (k+1,i+1) )
		f.write("\n")
		f.write("[ angles ]\n")
		f.write(";	ai	  aj	ak funct			c0			  c1			c2			  c3\n")
		for var in self.angles:
			f.write("%5d %5d %5d	5\n" % (var[0]+1,var[1]+1,var[2]+1) )
		f.write("\n")
		f.write("[ dihedrals ]\n")
		f.write(";	ai	  aj	ak	  al funct			  c0			c1			  c2			c3			  c4			c5\n")
		nonplanar_dihedrals=self.get_nonplanar_dihedrals(angl_params)
		for var in nonplanar_dihedrals:
			f.write("%5d %5d %5d %5d	 9\n" % (var[0]+1,var[1]+1,var[2]+1,var[3]+1) )
		f.write("\n")
		if(self.nimpropers > 0):
			f.write("[ dihedrals ]\n")
			f.write(";	ai	  aj	ak	  al funct			  c0			c1			  c2			c3\n")
			for var in self.impropers:
				f.write("%5d %5d %5d %5d	 2\n" % (var[0]+1,var[1]+1,var[2]+1,var[3]+1) )
			f.write("\n")
		## jal - add vsite directive
		## we use 3fd construction with a bit of a hack
		##	1. we manually set the value of a = 0
		##	2. constructing atoms j and k are intentionally the same
		##	3. and the distance becomes negative to mean "outside the bond"
		if (self.nvsites > 0):
			func=2
			a=0
			f.write("[ virtual_sites3 ]\n")
			f.write("; Site   from				funct a	   d\n")
			for atomi in range (0,self.nvsites):
				vsite = 0
				at1 = 0
				at2 = 0
				# find atom name matches
				for ai in range (0, self.natoms):
					if (self.G.node[ai]['name'] == self.G.node[atomi]['vsite']):
						vsite = ai
					if (self.G.node[ai]['name'] == self.G.node[atomi]['at1']):
						at1 = ai
					if (self.G.node[ai]['name'] == self.G.node[atomi]['at2']):
						at2 = ai
				dist=self.G.node[atomi]['dist']*-1
				f.write("%5d %5d %5d %5d %5d %5d %8.3f\n" % (vsite+1, at1+1, at2+1, at2+1, func, a, dist))
			f.write("\n")

		## jal - add exclusions for vsite
		if (self.nvsites > 0):
			f.write("[ exclusions ]\n")
			f.write(";	ai	  aj\n")
			## jal - explicitly add all 1-2, 1-3, and 1-4 exclusions
			## for the lone pair, which are the same as 1-2, 1-3, 1-4
			## exclusions for the host (bonds, angles, pairs) 
			# first, exclude any LP from its host
			for i in range (0, self.natoms):
				if ((is_lp_host_atom(self,self.G.node[i]['name'])==True)):
					# find the LP attached to this host, not necessarily consecutive
					# in the topology
					j = find_vsite(self, i)
					f.write("%5d %5d	 1\n" % (i+1,j+1) )
			# first neighbors: 1-2
			for i,j in self.G.edges_iter():
				if ((is_lp_host_atom(self,self.G.node[i]['name'])==True)):
					k = find_vsite(self, i)
					f.write("%5d %5d	 1\n" % (k+1,j+1) )
				if ((is_lp_host_atom(self,self.G.node[j]['name'])==True)):
					k = find_vsite(self, j)
					f.write("%5d %5d	 1\n" % (k+1,i+1) )
			# second neighbors: 1-3
			for var in self.angles:
				# only need to consider ends of the angle, not middle atom
				ai = var[0]
				ak = var[2]
				if ((is_lp_host_atom(self,self.G.node[ai]['name'])==True)):
					l = find_vsite(self, ai)
					f.write("%5d %5d	 1\n" % (l+1,ak+1) )
				if ((is_lp_host_atom(self,self.G.node[ak]['name'])==True)):
					l = find_vsite(self, ak)
					f.write("%5d %5d	 1\n" % (l+1,ai+1) )
			# third neighbors: 1-4
			for i,j in pairs14.edges_iter():
				if ((is_lp_host_atom(self,self.G.node[i]['name'])==True)):
					k = find_vsite(self, i)
					f.write("%5d %5d	 1\n" % (k+1,j+1) )
				if ((is_lp_host_atom(self,self.G.node[j]['name'])==True)):
					k = find_vsite(self, j)
					f.write("%5d %5d	 1\n" % (k+1,i+1) )
			f.write("\n")

		f.close()

#-----------------------------------------------------------------------
	def read_mol2_coor_only(self,filename):
		check_natoms = 0
		check_nbonds = 0
		f = open(filename, 'r')
		atm = {}
		section="NONE"
		for line in f.readlines():
			secflag=False
			if line.startswith("@"):
				secflag=True
				section="NONE"

			if((section=="NATO") and (not secflag)):
				entry = re.split('\s+', string.lstrip(line))
				check_natoms=int(entry[0])
				check_nbonds=int(entry[1])
				if(check_natoms != self.natoms):
					# jal - if there are lone pairs, these will not be in the mol2 file
					if (self.nvsites == 0):
						print "Error in atomgroup.py: read_mol2_coor_only: no. of atoms in mol2 (%d) and top (%d) are unequal" % (check_natoms, self.natoms)
						print "Usually this means the specified residue name does not match between str and mol2 files"
						#print check_natoms,self.natoms
						exit()
					else:
						print ""
						print "NOTE 4: %d lone pairs found in topology that are not in the mol2 file. This is not a problem, just FYI!\n" % (self.nvsites)
				# jal - if we have correctly ignored bonds to LP then there is no need
				# for any check here
				if(check_nbonds != self.nbonds):
					print "Error in atomgroup.py: read_mol2_coor_only: no. of bonds in mol2 (%d) and top (%d) are unequal" % (check_nbonds, self.nbonds)
					#print check_nbonds,self.nbonds
					exit()

				section="NONE"	

			if((section=="MOLE") and (not secflag)):
				self.name=line.strip()
				section="NATO" #next line after @<TRIPOS>MOLECULE contains atom, bond numbers

			if((section=="ATOM") and (not secflag)):
				entry = re.split('\s+', string.lstrip(line))
				## guard against blank lines
				if (len(entry) > 1):
					## jal - if there are lone pairs, these are not in mol2
					## and are not necessarily something we can just tack on at the
					## end of the coordinate section. Here, check the atom to see if it is
					## the first constructing atom, and if so, we put in a dummy LP entry.
					atomi = int(entry[0])-1
					self.G.node[atomi]['x'] = float(entry[2])
					self.G.node[atomi]['y'] = float(entry[3])
					self.G.node[atomi]['z'] = float(entry[4])
					self.coord[atomi][0] = float(entry[2])
					self.coord[atomi][1] = float(entry[3])
					self.coord[atomi][2] = float(entry[4])
					## jal - if we have an atom that is the host for a LP, insert
					## the LP into the list 
					if (is_lp_host_atom(self,self.G.node[atomi]['name'])):
						atomj = find_vsite(self, atomi)
						# insert dummy entry for LP
						self.G.node[atomj]['x'] = float(9999.99)
						self.G.node[atomj]['y'] = float(9999.99)
						self.G.node[atomj]['z'] = float(9999.99)
						self.coord[atomj][0] = float(9999.99)
						self.coord[atomj][1] = float(9999.99)
						self.coord[atomj][2] = float(9999.99)

			if line.startswith("@<TRIPOS>MOLECULE"):
				section="MOLE"
			if line.startswith("@<TRIPOS>ATOM"):
				section="ATOM"
			if line.startswith("@<TRIPOS>BOND"):
				section="BOND"
#-----------------------------------------------------------------------
	def write_pdb(self,f):
		for atomi in range(0,self.natoms):
			if(len(self.G.node[atomi]['name']) > 4):
				print "error in atomgroup.write_pdb(): atom name > 4 characters"
				exit()
			if (len(self.name) > 4):
				resn = self.name[:4]
			else:
				resn = self.name
			## jal - construct LP sites
			if (is_lp(self.G.node[atomi]['name'])):
				# DEBUG
				# print "Found LP in write_pdb: %s\n" % self.G.node[atomi]['name']
				# find constructing atoms, get their coordinates and construction distance*10
				atn1 = "dum" 
				atn2 = "dum"
				dist = 0
				# loop over vsites
				for ai in range (0,self.nvsites):
					if (self.G.node[ai]['vsite'] == self.G.node[atomi]['name']):
						atn1 = self.G.node[ai]['at1']	 # atom name
						atn2 = self.G.node[ai]['at2']	 # atom name
						dist = self.G.node[ai]['dist']*10 # Angstrom for PDB, was saved as *0.1 for GMX

				# get atom indices
				at1 = 0
				at2 = 0
				for ai in range (0, self.natoms):
					if (self.G.node[ai]['name'] == atn1):
						at1 = ai
					if (self.G.node[ai]['name'] == atn2):
						at2 = ai

				# in case of failure
				if ((at1==0) and (at2==0)):
					print "Failed to match LP-constructing atoms in write_pdb!\n"
					exit()

				# DEBUG
				# print "Found LP in write_pdb: %d %s %s with dist: %.3f\n" % ((atomi+1,at1+1,at2+1,dist))

				# at1, at2, and dist only exist in vsite structure!
				x1=self.coord[at1][0]
				y1=self.coord[at1][1]
				z1=self.coord[at1][2]
				x2=self.coord[at2][0]
				y2=self.coord[at2][1]
				z2=self.coord[at2][2]

				xlp,ylp,zlp = construct_lp(x1,y1,z1,x2,y2,z2,dist)
				self.coord[atomi][0] = xlp
				self.coord[atomi][1] = ylp
				self.coord[atomi][2] = zlp
			f.write("%-6s%5d %-4s %-4s%5s%12.3f%8.3f%8.3f%6.2f%6.2f\n" %
				("ATOM",atomi+1,self.G.node[atomi]['name'],resn,self.G.node[atomi]['resid'],self.coord[atomi][0],
				self.coord[atomi][1],self.coord[atomi][2],1.0,self.G.node[atomi]['beta']))
		f.write("END\n")

#=================================================================================================================


if(len(sys.argv) != 5):
	print "Usage: RESNAME drug.mol2 drug.str charmm36.ff"
	exit()

# check for compatible NetworkX version
if(float(nx.__version__) > 1.11):
	print "Your NetworkX version is: ",nx.__version__
	print "This script requires a version no higher than 1.11."
	print "Your NetworkX package is incompatible with this conversion script and cannot be used."
	exit()

if(sys.version_info > (3,0)):
	print("You are using a Python version in the 3.x series. This script requires Python 2.x.")
	print("Please visit http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs to get a script for Python 3.x")
	exit()

mol_name = sys.argv[1]
mol2_name = sys.argv[2]
rtp_name = sys.argv[3]
ffdir = sys.argv[4]
atomtypes_filename = ffdir + "/atomtypes.atp"

print "NOTE 1: Code tested with python 2.7.12. Your version:",sys.version
print ""
print "NOTE 2: Please be sure to use the same version of CGenFF in your simulations that was used during parameter generation:"
check_versions(rtp_name,ffdir + "/forcefield.doc")
print ""
print "NOTE 3: To avoid duplicated parameters, do NOT select the 'Include parameters that are already in CGenFF' option when uploading a molecule into CGenFF."


#for output
itpfile = mol_name.lower() + ".itp"
prmfile = mol_name.lower() + ".prm"
initpdbfile = mol_name.lower() + "_ini.pdb"
topfile = mol_name.lower() +".top"

atomtypes = read_gmx_atomtypes(atomtypes_filename)

angl_params = []  #needed for detecting triple bonds
filelist = get_filelist_from_gmx_forcefielditp(ffdir,"forcefield.itp")
for filename in filelist:
	anglpars = read_gmx_anglpars(filename)
	angl_params = angl_params + anglpars


m = atomgroup()
rtplines=get_charmm_rtp_lines(rtp_name,mol_name)
m.read_charmm_rtp(rtplines,atomtypes)


m.read_mol2_coor_only(mol2_name)
f = open(initpdbfile, 'w')
m.write_pdb(f)
f.close()


prmlines=get_charmm_prm_lines(rtp_name)
params = parse_charmm_parameters(prmlines)
write_gmx_bon(params,"",prmfile)
anglpars = read_gmx_anglpars(prmfile)
angl_params = angl_params + anglpars # append the new angl params


m.write_gmx_itp(itpfile,angl_params)
write_gmx_mol_top(topfile,ffdir,prmfile,itpfile,mol_name)

print "============ DONE ============"
print "Conversion complete."
print "The molecule topology has been written to %s" % (itpfile)
print "Additional parameters needed by the molecule are written to %s, which needs to be included in the system .top" % (prmfile)
print "\nPLEASE NOTE: lone pair construction requires duplicate host atom numbers, which will make grompp complain"
print "To produce .tpr files, the user MUST use -maxwarn 1 to circumvent this check"
print "============ DONE ============"

exit()

```
## 2.构建复合物

**1）准备gro文件**
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx  editconf -f jz4_ini.pdb -o jz4.gro
```
将3HTB_processed.gro文件复制到新的complex.gro文件中。
```
cp 3HTB_processed.gro complex.gro
```
然后将jz4.gro中，包含了jz4行的内容复制到complex.gro中（紧跟在蛋白行后），保持列对齐 。

 **注意！！！添加配体信息后，一定要修改complex.gro文件最上方的原子数量！！！**。

**2）构建拓扑**

打开topol.top文件，在最前部分出现`#include "./charmm36-jul2022.ff/forcefield.itp"`时，在这行后添加以下内容：

```
; Include ligand parameters
#include "jz4.prm"
```

翻到最下面到出现`; Include water topology`时，在这行前添加以下内容：
```
; Include ligand topology
#include "jz4.itp"
```

然后再[molecules]的最后，添加以下内容（对齐上文）：`JZ4         1`
## 3.溶剂化处理

```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx  editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0

/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx  solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
```
## 4.添加离子

查看电荷，可以从topol.top的[atoms]部分的最后一行看到一句`qtot 6`。

***<由于生命不是由网电荷组成的，因此要在此添加电荷>***

先创建一个ions.mpd文件
```
; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
coulombtype	    = cutoff	; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; long range electrostatic cut-off
rvdw		    = 1.0		; long range Van der Waals cut-off
pbc             = xyz 		; Periodic Boundary Conditions
```

***<使用.mdp文件运行能量最小化，有着最少参数、更易维护的特点>***

```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx  grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr

/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```
如果出现一个功能选择：【Select a continuous group of solvent molecules】，应该选择【Group    15 (            SOL) has 30882 elements】。

***<这一命令意思是，在这个体系里添加Na+和Cl-来中和电荷，在#4.里说提到了qtot 6，因此本例会添加6个Cl-。>***

***<通过将生成的水替换为Cl-，以此来实现电荷平衡，因此也应该选择group 15 (SOL)>***

## 5.能量最小化

创建一个文件em.mdp，写入以下内容
```
; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		        ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		    ; Method to determine neighbor list (simple, grid)
rlist		    = 1.2		    ; Cut-off for making neighbor list (short range forces)
coulombtype	    = PME		    ; Treatment of long range electrostatic interactions
rcoulomb	    = 1.2		    ; long range electrostatic cut-off
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw		    = 1.2		    ; long range Van der Waals cut-off
pbc             = xyz 		    ; Periodic Boundary Conditions
DispCorr        = no
```
然后运行：
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr

/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx mdrun -v -deffnm em
```
执行完上述命令后，会开始计算一个能量最低状态，就像这样：

>Steepest Descents converged to Fmax < 1000 in 153 steps

>Potential Energy  = -4.9159522e+05

>Maximum force     =  8.9174823e+02 on atom 27

>Norm of force     =  5.5987542e+01

## 6.平衡蛋白-受体
在此过程，有两个特殊考虑：1.限制配体；2.温度耦合组的处理

### 6.1  限制配体
通过给配体产生一个位置限制拓扑。
**2）创建一个不含氢的配体结构的索引组**
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx make_ndx -f jz4.gro -o index_jz4.ndx
```
之后会出现一些选择，输入(>后)
```
> 0 & ! a H*
…………………………………………………………
> q
```
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx genrestr -f jz4.gro -n index_jz4.ndx -o posre_jz4.itp -fc 1000 1000 1000
```
之后会产生一个选项，选择3：
> Select group to position restrain --> No.3

Group     3 (   System_&_!H*) has    10 elements

在拓扑中添加何种信息，完全取决于所需要的条件。

如果想在蛋白也受到抑制的时候也抑制配体，打开topol.top文件，在最后部分出现`; Include water topology`时，在这行前添加以下内容：
```
; Ligand position restraints
#ifdef POSRES
#include "posre_jz4.itp"
#endif
```

如果你想在平衡过程中有更多的控制，即独立抑制蛋白质和配体。打开topol.top文件，在最后部分出现`; Include water topology`时，在这行前添加以下内容：
```
; Ligand position restraints
#ifdef POSRES_LIG
#include "posre_jz4.itp"
#endif
```

如果要同时限制配体和蛋白，需要在.mdp文件里制定` define = -DPOSRES -DPOSRES_LIG`。

### 6.2  温度耦合

不要将系统中的每一个物种都单独配对。例如`tc-grps = Protein JZ4 SOL CL`这样会使得系统被破坏。

经典的方法是：设定`tc-grps = Protein JZ4 SOL CL`并继续。不幸的是，“非蛋白质”组也包括JZ4。由于JZ4和蛋白质在物理上紧密相连，最好将它们视为一个单一的实体。也就是说，JZ4与蛋白质结合用于温度耦合。同样地，之前插入的少数Cl-离子被认为是溶剂的一部分。因此需要一个特殊的索引组。
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx make_ndx -f em.gro -o index.ndx
```
之后会出现一些选择，输入(>后)
```
> 1 | 13
…………………………………………………………
> q
```
创建nvt.mdp文件，写入以下内容：
```
title                   = Protein-ligand complex NVT equilibration 
define                  = -DPOSRES  ; position restrain the protein and ligand
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 500   ; save energies every 1.0 ps
nstlog                  = 500   ; update log file every 1.0 ps
nstxout-compressed      = 500   ; save coordinates every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_JZ4 Water_and_ions    ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 300   300                     ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 300       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
```
注意`; Temperature coupling`部分，里面设置了`tc-grps                 = Protein_JZ4 Water_and_ions    `

然后执行以下命令：
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr

/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx mdrun -deffnm nvt
```
### 6.3  再平衡
当nvt模拟完成后，开始npt模拟。
创建npt.mdp文件，写入以下内容：
```
title                   = Protein-ligand complex NPT equilibration 
define                  = -DPOSRES  ; position restrain the protein and ligand
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
nstxout-compressed      = 500       ; save coordinates every 1.0 ps
; Bond parameters
continuation            = yes       ; continuing from NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_JZ4 Water_and_ions    ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 300   300                     ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = Berendsen                     ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; velocity generation off after NVT 
```
然后运行以下内容：
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr

/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx mdrun -deffnm npt
```
## 7.开始MD
经过两次平衡后，温度和压力已经稳定，通过释放位置限制，并开始MD。

创建一个md.mdp文件，并写入以下内容：
```
title                   = Protein-ligand complex MD simulation 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 5000000   ; 2 * 5000000 = 10000 ps (10 ns)
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save coordinates every 10.0 ps
; Bond parameters
continuation            = yes       ; continuing from NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_JZ4 Water_and_ions    ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 300   300                     ; reference temperature, one for each group, in K
; Pressure coupling 
pcoupl                  = Parrinello-Rahman             ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; continuing from NPT equilibration 
```
然后执行以下命令：
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx rompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr

/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx mdrun -deffnm md_0_10
```
## 8.分析
### 8.1  重定位和坐标重绘制
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact

```
计算涉及跨越周期边界的许多跳跃的长模拟的中心复合物，是比较困难的。因此可能需要创建一个自定义的索引组用于定中心，对应于复合物中一种蛋白质的活性位点或一种单体的界面残基。

提取轨迹的第一帧：
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx  trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o start.pdb -dump 0
```
为了更顺滑地可视化，对其做旋转和平移拟合：
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o md_0_10_fit.xtc -fit rot+trans
```
选择“骨架”对蛋白质骨架进行最小二乘拟合，选择“系统”进行输出。注意，同时进行PBC重写和坐标拟合在数学上是不兼容的。如果希望执行拟合，则需要单独进行。

### 8.2  分析蛋白-配体相互作用力和配体动力学
##案例

该命令可以计算轨迹中的距离。
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx angle -f md_0_10_center.xtc -n index.ndx -ov angle.xvg
```
量化配体结合姿势在模拟过程中发生了多大变化。要计算JZ4的重原子RMSD。
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx make_ndx -f em.gro -n index.ndx
```
之后会出现一些选择，输入(>后)
```
> 13 & ! H*
…………………………………………………………
> name 26 JZ4_heavy
…………………………………………………………
> q 
```
选择“骨架”进行最小二乘拟合，选择“JZ4_Heavy”进行RMSD计算。通过拟合去除了蛋白质的整体旋转和平移，RMSD报告了JZ4位置相对于蛋白质的变化程度，这很好地表明了在模拟过程中结合姿势的保留程度。
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx rms -s em.tpr -f md_0_10_center.xtc -n index.ndx -tu ns -o rmsd_jz4.xvg
```
### 8.3  蛋白-配体相互作用能

创建一个ie.mdp文件，写入以下内容：
```
title                   = Protein-ligand complex MD simulation 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 5000000   ; 2 * 5000000 = 10000 ps (10 ns)
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save coordinates every 10.0 ps
energygrps              = Protein JZ4
; Bond parameters
continuation            = yes       ; continuing from NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_JZ4 Water_and_ions    ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 300   300                     ; reference temperature, one for each group, in K
; Pressure coupling 
pcoupl                  = Parrinello-Rahman             ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; continuing from NPT equilibration 
```
然后执行：
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o ie.tpr
```
之后从已有的模拟轨迹里重新计算能量
```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx mdrun -deffnm ie -rerun md_0_10.xtc -nb cpu
```

```
/pubhome/soft/gromacs/2020.5/thread_mpi/bin/gmx energy -f ie.edr -o interaction_energy.xvg
```
