# DOCK Test Summary

## 1.安装DOCK37和环境
(需要先安装miniconda)
 ```
 mkdir dock_learn && cd dock_learn

 installpath="."

 rsync -a n101:/home/soft/ucsfdock/bela30c $installpath
 ```
### 1.1  设置环境
创建dock37.yml环境调谐文件
```
name: dock37
channels:
  - conda-forge
  - free
dependencies:
  - python=2
  # from $DOCKBASE/install/environ/python/requirements.txt
  - numpy<1.9
  # - scipy
  - biopython
  - psycopg2
  - MySQL-python
  - sh
  - matplotlib
  - lockfile
  - psutil
  # from practices with scripts in $DOCKBASE/analysis
  - bsddb
  - requests
  - pip
  - pip:
    - scipy # along with libgfortran, in case the later not in system
  # Other packages may cause conflicts, raising warning like
  # 
  # Found conflicts! Looking for incompatible packages.
  # This can take several minutes.  Press CTRL-C to abort.
  # 
  # You can try something anyway.
```
### 1.2  创建环境
`conda env create -f dock37.yml`
### 1.3  获得文件
可先创一个文件夹，用以存放必要的安装包等。
例如已创建一个文件夹test(.等同于在此地址)

将这三个文件放到./bela30c/proteins/blastermaster下

**1.fmt_lig.py**
```
# If you really want to reproduce, put this script under the path below
# $DOCKBASE/proteins/blastermaster/fmt_lig.py

import os
import pdbMoveColumns,pdb

if True:
  class options:
    workingDir = "working"
    ligand = "xtal-lig.pdb"
  pdbMoveColumns.pdbMoveColumns(
      options.ligand,
      os.path.join(options.workingDir, 'pre-rename.' + options.ligand))
  pdbLigand = pdb.pdbData(
      os.path.join(options.workingDir, 'pre-rename.' + options.ligand))
  pdbLigand.replaceHETATMwithATOM()
  pdbLigand.write(os.path.join(options.workingDir, options.ligand))

```
**2.fmt_rec.py**
```
# If you really want to reproduce, put this script under the path below
# $DOCKBASE/proteins/blastermaster/fmt_rec.py

import os
import pdbMoveColumns,pdbMostOccupied,pdb

if True:
  class options:
    workingDir = "working"
    receptor = "rec.pdb"
  if not os.path.exists(options.workingDir):
    os.mkdir(options.workingDir)  # make working directory if it isn't there
  pdbMoveColumns.pdbMoveColumns(
      options.receptor,
      os.path.join(options.workingDir, 'pre-most.occ.' + options.receptor))
  pdbMostOccupied.pdbMostOccupied(
      os.path.join(options.workingDir, 'pre-most.occ.' + options.receptor),
      os.path.join(options.workingDir, 'post-most.occ.' + options.receptor))
  pdbPreRenameResidues = pdb.pdbData(
      os.path.join(options.workingDir, 'post-most.occ.' + options.receptor))
  pdbPreRenameResidues.replaceAltChars(' ')
  pdbPreRenameResidues.deleteInsertionCodes()
  pdbPreRenameResidues.fixChainIds()
  pdbPreRenameResidues.write(os.path.join(options.workingDir, options.receptor))
```
**3.polar_pdb.py**
```
# If you really want to reproduce, put this script under the path below
# $DOCKBASE/proteins/blastermaster/polarpdb.py

import os
import pdb

if True:
    chargedReceptor = "rec.crg.pdb"
    chargedReceptorPolarH = chargedReceptor + ".polarH"
    chargedReceptorFullH = chargedReceptor + ".fullh"
    workingDir = "."
    
    pdbD = pdb.pdbData(
        os.path.join(workingDir, chargedReceptorFullH), ignoreWaters=False)
    # ignoreWaters=True can only ignore waters with name "HOH"
    pdbD.removeApolarHydrogen() # A dictionary pdb.keepPolarH is used
    pdbD.write(os.path.join(workingDir, chargedReceptorPolarH))
    pdbD = pdb.pdbData(
        os.path.join(workingDir, chargedReceptorPolarH), ignoreWaters=False)
    # if initial name != "HIS"/"CYS", then will not rename
    pdbD.renameHistidines() # NAME: HIS -> HID/HIE/HIP
    pdbD.renameCysteines() # NAME: CYS -> CYS/CYX(no HG on -SH)
    pdbD.write(os.path.join(workingDir, chargedReceptor))
```
### 1.4  安装db2_converter
```
installpath="."

rsync -av k143:/pubhome/qcxia02/git-repo/db2_converter $installpath
```
### 1.5  安装db2_converter必须环境
```
conda create -n test python=3 -y

conda activate test

conda install -c conda-forge ambertools openbabel rdkit -y

installpath="."

rsync -av k143:/pubhome/qcxia02/2hnlab/share/soft/ZBH2022 $installpath

cd $installpath/ZBH2022

tar -xf conformator_1.2.1_CentOS-6.10-64bit.tar.gz

tar -xf unicon_1.4.2_CentOS-6.10-64bit.tar.gz
# 激活 unicon 和 conformator
./unicon_1.4.2/unicon --license=AAAAAAAlix8AAAAUcyTuOIbQiuYSaO37g+OLKBPN60Y=
# 应该不会报错
./conformator_1.2.1/conformator --license=AAAAAAAljSUAAAAU9AVqNnEVm+Hufq7HJxI/5Vyp8mo=
```

## 2.运行blastermaster
准备rec.pdb文件。rec.pdb即为蛋白受体分子（receipt.pdb）
注意：蛋白的cofactor的标识应该从HETATM改为ATOM.

准备xtal-lig.pdb文件。xtal-lig.pdb即为共晶小分子
可以通过选择蛋白活性口袋附近的坐标，对该小分子进行位置标识。
注意：小分子的pdb应该满足以下要求：
HETATM标识，原子顺序为cofactor+1

先执行：
```
export DOCKBASE=$(readlink -f "$installpath/DOCK-3.7-bbe1a30c")

conda activate dock37 #激活dock37环境
# conda deactivate  #退出环境
```
***直接执行***
```
mkdir test && cd test

wget -c https://files.rcsb.org/download/6N2W.pdb
# obtain rec.pdb

grep "^ATOM" 6N2W.pdb | grep " B " | grep -v "HOH" > rec.pdb

grep "^HETATM" 6N2W.pdb | grep "FE"  | grep "B" | sed "s/HETATM/ATOM  /" >> rec.pdb
# obtain xtal-lig.pdb

grep "30Z B 702" 6N2W.pdb > xtal-lig.pdb
```
```
# 对于直接执行的情况来说，一个默认解释器是 python2 的环境是必需的，否则 close_sph.py 会报错
python2.7 $DOCKBASE/proteins/blastermaster/blastermaster.py --addhOptions=" -HIS -FLIPs " -v
```
如果出现错误：
>WARNING:root:Could not find numpy or numeric

>CRITICAL:root:Could not find any matrix library (Numeric, numpy, pMatrix)! Cannot proceed

>Traceback (most recent call last):

>File "/home/shianine/dock_learn/DOCK-3.7-bbe1a30c/proteins/blastermaster/blastermaster.py", line 56, in <module>

>import blasterAddHydrogens

>File "/home/shianine/dock_learn/DOCK-3.7-bbe1a30c/proteins/blastermaster/blasterAddHydrogens.py", line 11, in <module>

>import pdb

>File "/home/shianine/dock_learn/DOCK-3.7-bbe1a30c/proteins/blastermaster/pdb.py", line 8, in <module>

>import geometry  # distance function

>File "/home/shianine/dock_learn/DOCK-3.7-bbe1a30c/proteins/blastermaster/geometry.py", line 27, in <module>

>import pMatrix

>ImportError: No module named pMatrix

则应在环境（dock37.yml）中安装
```
pip2 install numpy

pip2 install scipy
```
***分步执行***
### 2.1  准备必要文件
准备rec.pdb文件

```
mkdir test_steps && cd test_steps 
#将受体蛋白的pdb文件放到test_steps路径下

grep "^ATOM" name.pdb | grep -v "HOH" > rec.pdb
#去除水

grep "^HETATM" 6N2W.pdb | grep "FE"  | grep "B" | sed "s/HETATM/ATOM  /" >> rec.pdb
#将cofactor和金属离子等的标识改为ATOM(如有)
xtal-lig.pdb文件可以通过之前介绍的方式获得，或者直接从现成文件中提取。
```
按位读取
```
python2 $DOCKBASE/proteins/blastermaster/fmt_rec.py

python2 $DOCKBASE/proteins/blastermaster/fmt_lig.py
```
### 2.2  加氢

上述指令完成后，会在test_steps下产生一个working文件夹

---
#### 目的：给蛋白加氢；
#### 输入文件：rec.pdb  /{home}/bela30c/proteins/default/reduce_wwPDB_het_dict.txt
#### 输出文件：rec.crg.pdb.fullh    addh.log   rec.crg.pdb  rec.crg.pdb.polarH
---
```
cd working

cp /{home}/bela30c/proteins/defaults/reduce_wwPDB_het_dict.txt .

/{home}/bela30c/proteins/Reduce/reduce \
    -db reduce_wwPDB_het_dict.txt \
    -HIS -FLIPS \
    rec.pdb > rec.crg.pdb.fullh 2> addh.log

#    -OH -HIS -ALLALT -NOROTNH3 -FLIPS \
#   -Keep -METALBump1.5 -NONMETALBump-5.0 \

# format
sed -i 's/\s*new\s*//g' rec.crg.pdb.fullh

sed -i '/^USER.*/d' rec.crg.pdb.fullh

python2 /{home}/bela30c/proteins/blastermaster/polar_pdb.py
```
### 2.3  基于共晶小分子做结合位点筛选（filter）
---
#### 目的：筛选出结合位点；
#### 输入文件：rec.pdb  xtal-lig.pdb  /{home}/bela30c/proteins/defaults/filt.params
#### 输出文件：filter.log   rec.site
---
```
/{home}/bela30c/proteins/filt/bin/filt < /{home}/bela30c/proteins/defaults/filt.params > filter.log
```
### 2.4  生成表面（dms）

根据溶剂可及表面（Solvent Accessible Surface），默认用半径1.4 A的水探针小球在该面下滚动，生成一个溶剂排除表面（Solvent Excluded Surface）。

---
#### 目的：生成表面（覆盖外层原子，用于后续球面的生成和评估）；
#### 输入文件：rec.crg.pdb.dms  rec.site.dms  /{home}/bela30c/proteins/defaults/radii
#### 输出文件：rec.ms  dms.log  rec.ts.ms  dms.ts.log
---
```
grep -a -v HOH rec.crg.pdb > rec.crg.pdb.dms
grep -a -v HOH rec.site > rec.site.dms

# -a 使得grep命令将带查找的文件视为纯文本，保证查找不会出现问题
# -v 反向匹配，将排除包含目标词的内容

/{home}/bela30c/proteins/dms/bin/dms rec.crg.pdb.dms -a -d 1.0 -i rec.site.dms -g dms.log -p -n -o rec.ms # density 1.0

/{home}/bela30c/proteins/dms/bin/dms rec.crg.pdb.dms -a -d 1.0 -i rec.site.dms -g dms.ts.log -p -n -o rec.ts.ms

# -a all
# -d density 1.0 means 5 surface points per angstrom^2 on average
# -i 设置目标文件
# -g 生成日志.log文件
# -o 输出文件，后跟文件名
# -n include the unit normal vectors in surface point file
```
### 2.5  球面生成（sphgen）

sphgen会产生描述分子和分子表面形状的重叠球面，对于受体：会产生负电的表面内陷；对于配体：会产生全分子的正像。

上步生成溶剂排除表面后，球面生成时，会接触两处溶剂排除表面（类似相切），后形成一个球体。这个球面会在整个溶剂排除表面生成，每个曲面点约产生一个球体，之后会过滤较密集的，只保留与每个受体表面原子相关的最大球体。

---
#### 目的：生成描述结合口袋的负电图像；
#### 输入文件：/{home}/working/INSPH  (rec.ms)
#### 输出文件：OUTSPH  all_spheres.sph
---
```
echo "rec.ms
R  #1
X  #2
0.  #3
5.0  #4
1.4  #5
all_spheres.sph" > INSPH
/{home}/bela30c/proteins/sphgen/bin/sphgen
sed -i '1d' all_spheres.sph

#1.在表面外生成（R），内生成（L）；
#2.曲面点的使用，全（X）；
#3.放置产生与表面有紧密接触的大球体（默认0.0）；
#4.最大球半径（默认5.0）
#5.最小球半径（默认1.4）
```
### 2.6  细球体生成（thin_spheres）
---
#### 目的：生成低介电常数和配体去溶剂化的边界层球；
#### 输入文件：rec.ts.ms
#### 输出文件：low_die_thinspheres.sph  lig_die_thinspheres.sph
---
```
python2 /{home}/bela30c/proteins/thinspheres/thin_spheres.py -i rec.ts.ms -o low_die_thinspheres.sph -d 1.0 -s 1.0

python2 /{home}/bela30c/proteins/thinspheres/thin_spheres.py -i rec.ts.ms -o lig_des_thinspheres.sph -d 1.0 -s 1.0

# -d: sphere distance from receptor surface
# -s: sphere radius
# 默认都是1.0 A
```
### 2.7  选择共晶小分子距离内的边界球
---
#### 目的：选择共晶小分子距离内的边界球；
#### 输入文件：low_die_thinspheres.sph  lig_des_thinspheres.sph
#### 输出文件：low_die_thinspheres.sph.close  lig_des_thinspheres.sph.close   low_die_thinspheres.sph.close.log   lig_des_thinspheres.sph.close.log
---
```
python2 /{home}/bela30c/proteins/thinspheres/close_sph.py low_die_thinspheres.sph xtal-lig.pdb  low_die_thinspheres.sph.close 2.0 1.0

python2 /{home}/bela30c/proteins/thinspheres/close_sph.py lig_des_thinspheres.sph xtal-lig.pdb  lig_des_thinspheres.sph.close 2.0 1.0
# 设定球半径为1.0 A
# 任意球体和配体的截止距离为2.0 A
```
### 2.8  生成正像
---
#### 目的：将共晶小分子转化为球面
#### 输入文件：xtal-lig.pdb
#### 输出文件：xtal-lig.match.sph
---
```
/{home}/bela30c/proteins/pdbtosph/bin/pdbtosph xtal-lig.pdb xtal-lig.match.sph
```
### 2.9  筛选低介电常数球面
makespheres1.cli.pl

---
#### 目的：根据共晶小分子球面和描述结合口袋负电的球面，筛选出用于后续计算的球面
#### 输入文件：xtal-lig.match.sph  all_spheres.sph  rec.crg.pdb
#### 输出文件：lowdielectric.sph  lowdielectric_spheres.log
---
```
/{home}/bela30c/proteins/makespheres1/makespheres1.cli.pl xtal-lig.match.sph all_spheres.sph rec.crg.pdb  lowdielectric.sph 25 >& lowdielectric.spheres.log

# 设定至少产生25个球面，最大产生120个
# CONSTANTS defined here
# $M = 12;              # 保证球面靠近配体中心的最值为12 A
# $R = 7;               # 保证球面靠近受体原子的最值为7 A
# $Rclose = 1.2;        # 保证球面远离受体原子的边缘值为1.2 A
# $polardist = 3.3;     # 球体和极性重受体原子的极性距离
# $Hpolardist = 2.5;    # 球体和极性H受体原子的极性距离
# $nonpolardist = 4.5;  # 球体和C受体原子的非极性距离
# $gridsize = 1.5;      # 球面点之间的典型距离，即每个网格单位区域中只保留一个点
# $tooclose = 0.8;      # 无偏网格法后第二遍消除过于接近的球体
# $continuity = 3.0;    # 在delphi中，点之间的距离多远时，可以并且仍然被认为是连续的
# $contincrement = 0.5; # 如果根据连续性发现的球体少于$minspheres，则增加连续性
# $maxcontinuity = 4.5; # 当需要增加时，允许最大值连续性  
# $numrepeats = 10;     # 检查点连续性要多少次以保持球体生长
# $minspheres = $ARGV[4];     # 输出球体的最小数量，niu
# $maxspheres = 120;     # delphi输入的最大球体数量，基于DOCK的ARBITRARY，niu
# $polarbenefit =-0.25; # 最终选择球体的平均重量收益百分比
```
### 2.10  转化文件
提取球体坐标，将其转化为C，并写入到pdb文件

---
#### 目的：将.sph文件转化为.pdb，并且添加到对应文件里
#### 输入文件：rec.crg.pdb  lowdielectric.sph  low_die_thinspheres.sph.close  lig_des_thinspheres.sph.close
#### 输出文件：lowdielectric.sph.pdb    lowdielectric.sph.pdb.log  low_die_thinspheres.sph.close.pdb low_die_thinspheres.sph.close.pdb.log  lig_des_thinspheres.sph.close.pdb  lig_des_thinspheres.sph.close.pdb.log  receptor.crg.lowdielectric.pdb
---
```
/{home}/bela30c/proteins/showsphere/doshowsph.csh lowdielectric.sph 1 lowdielectric.sph.pdb >&  lowdielectric.sph.pdb.log

/{home}/bela30c/proteins/showsphere/doshowsph.csh low_die_thinspheres.sph.close 1  low_die_thinspheres.sph.close.pdb >& low_die_thinspheres.sph.close.log

/{home}/bela30c/proteins/showsphere/doshowsph.csh lig_des_thinspheres.sph.close 1  lig_des_thinspheres.sph.close.pdb >& lig_des_thinspheres.sph.close.log

cat rec.crg.pdb lowdielectric.sph.pdb > receptor.crg.lowdielectric.pdb
```
### 2.11  产生匹配球
makespheres3.cli.pl

---
#### 目的：产生匹配球
#### 输入文件：xtal-lig.match.sph  all_spheres.sph  rec.crg.pdb
#### 输出文件：matching_spheres.sph  matching_spheres.log
---
```
$DOCKBASE/proteins/makespheres3/makespheres3.cli.pl 1.5 0.8 45 xtal-lig.match.sph all_spheres.sph  rec.crg.pdb matching_spheres.sph >& matching_spheres.log

# CONSTANTS defined here
# $M = 10;              # 保证球面靠近配体中心的最值为10 A
# $R = 7;               # 保证球面靠近受体原子的最值为7 A
# $Rclose = 1.2;        # 保证球面远离受体原子的边缘值为1.2 A
# $polardist = 3.3;     # 球体和极性重受体原子的极性距离
# $Hpolardist = 2.5;    # 球体和极性H受体原子的极性距离
# $nonpolardist = 4.5;  # 球体和C受体原子的非极性距离
# $gridsize = $ARGV[0];      # 球面点之间的典型距离，即每个网格单位区域中只保留一个点
# $tooclose = $ARGV[1];      # 无偏网格法后第二遍消除过于接近的球体
# $continuity = 3.0;    # 在docking中，点之间的距离多远时，可以并且仍然被认为是连续的
# $contincrement = 0.5; # 如果根据连续性发现的球体少于$minspheres，则增加连续性
# $maxcontinuity = 4.5; # 当需要增加时，允许最大值连续性
# $numrepeats = 10;     # 检查点连续性要多少次以保持球体生长
# $minspheres = 20;     # 输出球体的最小数量，niu
# $maxspheres = $ARGV[2];     # 最大球体数量，基于DOCK的ARBITRARY，niu
# $polarbenefit =-0.25; # 最终选择球体的平均重量收益百分比
# $LIG = $ARGV[3]; #input sph/match.sph or sph/match0.sph file
# $SPH = $ARGV[4]; #input sph/sph file
# $REC = $ARGV[5]; #input grids/rec.crg file
# $MATCHSPH = $ARGV[6]; #output .sph file
```
## 3.生成网格
### 3.1  生成dock打分需要的网格
---
#### 目的：生成打分需要的网格
#### 输入文件：xtal-lig.match.sph  rec.crg.pdb
#### 输出文件：box  makebox.log
---
```
/{home}/bela30c/proteins/makebox/makebox.smallokay.pl xtal-lig.match.sph rec.crg.pdb box 10.0 >&  makebox.log
```
### 3.2  用数值解法生成静电网格
---
#### 目的：用数值解法生成静电网格
#### 输入文件：/{home}/bela30c/proteins/defaults/qnifft.parm  (receptor.crg.lowdielectric.pdb,amb.crg.oxt,vdw.siz)
#### 输出文件：qnifft.atm  qnifft.electrostatics.phi  qnifft.eps  qnifft.log  qnifft_sas.usr  
---
```
cp /{home}/bela30c/proteins/defaults/delphi.def qnifft.parm

echo "grid=193" >> qnifft.parm

echo "charge=amb.crg.oxt" >> qnifft.parm

echo "radius=vdw.siz" >> qnifft.parm # default extended atom radii based loosely on mike connolly's MS program

echo "pdb_input=receptor.crg.lowdielectric.pdb" >> qnifft.parm

echo "pdb_output_file=qnifft.atm" >> qnifft.parm

echo "phi_output_file=qnifft.electrostatics.phi" >> qnifft.parm

echo "sizing=border border_solvent=10." >> qnifft.parm

cp /{home}/bela30c/proteins/defaults/amb.crg.oxt . # charge info for each atom in each residue type

cp /{home}/bela30c/proteins/defaults/vdw.siz . # atom radii

/{home}/bela30c/proteins/qnifft/bin/qnifft22_193_pgf_32 qnifft.parm >& qnifft.log
```
### 3.3  休整phi-map
---
#### 目的：休整phi-map
#### 输入文件：box  qnifft.electrostatics.phi
#### 输出文件：trim.electrostatics.phi 
---
```
python2 /{home}/bela30c/proteins/blastermaster/phiTrim.py qnifft.electrostatics.phi box trim.electrostatics.phi
```
### 3.4  产生范德华打分网格
---
#### 目的：产生范德华打分网格
#### 输入文件：/{home}/working/INCHEM  (rec.crg.pdb  /{home}/bela30c/proteins/defaults/prot.table.ambcrg.ambH  /{home}/bela30c/proteins/defaults/vdm.parms.amb.mindock, box)
#### 输出文件：vdm.bmp  vdm.esp  vdm.log  vdm.vdw
---
```
echo "rec.crg.pdb" >> INCHEM ## receptor pdb file

echo "prot.table.ambcrg.ambH" >> INCHEM ## receptor parameters

echo "vdw.parms.amb.mindock" >> INCHEM ## van der Waals parameter file

echo "box" >> INCHEM

echo "0.2" >> INCHEM ## grid spacing (angstrom)

echo "1" >> INCHEM ## option 1: determine a distance-dependent dielectric

echo "4" >> INCHEM ## the dielectric function will be multiplied by   4.00

echo "10" >> INCHEM ## cutoff distance for vdw energy calculations

echo "2.3 2.6" >> INCHEM ## receptor polar atoms bump distance (2.3), carbon atom bump distance (2.6)

echo "vdw" >> INCHEM ## prefix

cp $DOCKBASE/proteins/defaults/prot.table.ambcrg.ambH . ## dock atom typing (tt)

cp $DOCKBASE/proteins/defaults/vdw.parms.amb.mindock . ## atom type -> A,B parameters (26 vdw atom types)

$DOCKBASE/proteins/chemgrid/bin/chemgrid >& vdw.log ## most output will be in OUTCHEM instead of vdw.log
```
### 3.5  产生非溶剂化配体打分网格
solvemap。运行两次，第一次对于配体重原子，第二次对于配体氢原子

---
#### 目的：产生非溶剂化配体打分网格
#### 输入文件：/{home}/working/heavy/INSEV  /{home}/working/hydrogen/INSEV  box  rec.crg.lds.pdb
#### 输出文件：/{home}/working/OUTSEV  /{home}/working/solvemap.log  /{home}/working/dielec_sev.box  /{home}/working/solv_sev.box  /{home}/working/solv_sev.plt  /{home}/working/heavy/ligand.desolv.heavy  /{home}/working/hydrogen/ligand.desolv.hydrogen  *heavy; hydrogen  
---
```
cp rec.crg.pdb rec.crg.lds.pdb

mkdir heavy && cd heavy

cp ../rec.crg.lds.pdb .

cp ../box .

echo "rec.crg.lds.pdb" >> INSEV # name of input pdb file

echo "ligand.desolv.heavy" >> INSEV # name of output contact grid file

echo "1.60,1.65,1.90,1.90,1.90,1.00" >> INSEV # atomic radii of O,N,C,S,P,Other

echo "1.4" >> INSEV # atomic radii of water probe (exclusion)

echo "2" >> INSEV # number of grid points per angstrom    

echo "box" >> INSEV # pdb file of limiting box

echo "1.8" >> INSEV # assumed ligand atomic radius

cd ..

mkdir hydrogen && cd hydrogen

cp ../rec.crg.lds.pdb .

cp ../box .

echo "rec.crg.lds.pdb" >> INSEV

echo "ligand.desolv.hydrogen" >> INSEV

echo "1.60,1.65,1.90,1.90,1.90,1.00" >> INSEV

echo "1.4" >> INSEV

echo "2" >> INSEV

echo "box" >> INSEV

echo "1.0" >> INSEV

cd ..

cd heavy; /{home}/bela30c/proteins/solvmap/bin/solvmap >& solvmap.log

cd ..

cd hydrogen; /{home}/bela30c/proteins/solvmap/bin/solvmap >& solvmap.log

cd ..
```
### 3.6  产生INDOCK文件
```
cd ..
delphi_nsize=57 # obtained from phiTrim@step13
echo """DOCK 3.7 parameter
#####################################################
# NOTE: split_database_index is reserved to specify a list of files
# defults for large scale docking.
ligand_atom_file               split_database_index
#####################################################
#                             OUTPUT
output_file_prefix            test.
#####################################################
#                             MATCHING
match_method                  2
distance_tolerance            0.05
match_goal                    1000
distance_step                 0.05
distance_maximum              0.5
timeout                       10.0
nodes_maximum                 4
nodes_minimum                 4
bump_maximum                  10.0
bump_rigid                    10.0
mol2_score_maximum            -10.0
#####################################################
#                             COLORING
chemical_matching             no
case_sensitive                no
#####################################################
#                             SEARCH MODE
atom_minimum                  4
atom_maximum                  25
number_save                   1
number_write                  1
flush_int                     100
#molecules_maximum            100000
check_clashes                 yes
do_premax                     no
do_clusters                   no
#####################################################
#                             SCORING
ligand_desolvation            volume
#vdw_maximum                   1.0e10
ligand_desolv_scale           1.0
electrostatic_scale           1.0
vdw_scale                     1.0
internal_scale                0.0
per_atom_scores               no
##################################################### 
#                             DOCKovalent 
dockovalent                   no
bond_len                      1.8
bond_ang1                     109.5
bond_ang2                     109.5
len_range                     0.0
len_step                      0.1
ang1_range                    10.0
ang2_range                    10.0
ang1_step                     2.5
ang2_step                     2.5
#####################################################
#                    MINIMIZATION
minimize                      yes
sim_itmax                     500
sim_trnstep                   0.2
sim_rotstep                   5.0
sim_need_to_restart           1.0
sim_cnvrge                    0.1
min_cut                       1.0e15
iseed                         777
##################################################### 
# INPUT FILES / THINGS THAT CHANGE
receptor_sphere_file          ../dockfiles/matching_spheres.sph
vdw_parameter_file            ../dockfiles/vdw.parms.amb.mindock
delphi_nsize                  $delphi_nsize
flexible_receptor             no
total_receptors               1
############## grids/data for one receptor
rec_number                    1
rec_group                     1
rec_group_option              1
solvmap_file                  ../dockfiles/ligand.desolv.heavy
hydrogen_solvmap_file         ../dockfiles/ligand.desolv.hydrogen
delphi_file                   ../dockfiles/trim.electrostatics.phi
chemgrid_file                 ../dockfiles/vdw.vdw
bumpmap_file                  ../dockfiles/vdw.bmp
############## end of INDOCK""" > INDOCK

mkdir dockfiles

cp working/matching_spheres.sph dockfiles/

cp working/vdw.parms.amb.mindock dockfiles/

cp working/heavy/ligand.desolv.heavy dockfiles/

cp working/hydrogen/ligand.desolv.hydrogen dockfiles/

cp working/trim.electrostatics.phi dockfiles/

cp working/vdw.vdw dockfiles/

cp working/vdw.bmp dockfiles/
```
## 4.DB2_converter
### 4.1  测试
```
conda activate test

which antechamber

which obabel

conda list rdkit

export DOCKBASE=$(readlink -f "/{home}/test_steps/db2_converter")
# 注意和之前的路径一致

pip install -e $installpath/db2_converter/ligand/generate/db2_converter
```
修改：/{home}/test_steps/db2_converter/ligand/generate/db2_converter/db2_converter/config.py
处的路径为自己的文件路径

测试
`echo $CONDA_PREFIX`

***TLDR处理***

https://tldr.docking.org/start。
使用build3d37进行，将共晶小分子的smiles文本内容，保存在input.smi文件里，build3d37可以对其进行多构象分析，产生一个多构象的.mol2文件。

### 4.2  DOCKing
```
/unicon -i xtal-lig.pdb -o xtal-lig.smi -t single -p single
#应该在unicon文件路径下执行，xtal-lig.pdb和xtal-lig.smi可以使用绝对路径。

max_conf=100 # default

workingpath=$(readlink -f "xtal-lig")

outputpath=$(readlink -f "xtal-lig")

method="conformator" # default

insmi="xtal-lig.smi"

build_ligand -i $insmi -n $max_conf --keep_max_conf --workingpath $workingpath --outputpath $workingpath --method $method -c --reseth --rotateh
```
如果有错误：一个表格中，某个组件的value部分，不是一个路径，available状态不是True。

则应在：/{home}/db2_converter/db2_converter/db2_converter/config.py 里修改对应组件的路径。如果该路径下的config.py无问题，应参考实际报错位置的config.py，修改和补全里面的组件的路径。

---
创建一个db2utils文件夹`mkdir db2utils`
在这个文件夹里，写入以下三个文件：
**1.check_TLDR.py**
```
from rdkit import Chem
from pathlib import Path
from tqdm import tqdm
import sys

def next_mol2_lines(infile):
    """Method to return one mol2 block once."""
    lines = list()

    for line in open(infile):
        if "@<TRIPOS>MOLECULE" in line:
            if len(lines) == 0:
                lines.append(line)
            else: # in case there are multiple mol2blocks in infile
                yield lines
                lines = list()
                lines.append(line)
        else:
            lines.append(line)

    yield lines

def checksmi(insmi, all_mols):
    canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(insmi), isomericSmiles=True)
    for i,mol in enumerate(all_mols):
        Chem.AssignAtomChiralTagsFromStructure(mol)
        Chem.AssignStereochemistryFrom3D(mol)
        smi = Chem.MolToSmiles(mol,isomericSmiles=True)
        if smi == canonical_smiles:
            return True
    return False

if True:
    insmifile = sys.argv[1]
    mol2path = sys.argv[2]
    successsums = []
    failsums = []
    for line in tqdm(Path(insmifile).read_text().split("\n"), total=len(Path(insmifile).read_text().split("\n"))):
        if line:
            statuss = []
            smi, name = line.split()
            try:
                mol2files = sorted(Path(mol2path).glob(f"{name}.*.mol2"))
                if mol2files:
                    for mol2file in mol2files:
                        all_blocks = [x for x in next_mol2_lines(mol2file)]
                        all_mols = [
                            Chem.MolFromMol2Block("".join(x), removeHs=True) for x in all_blocks
                        ]
                        statuss.append(checksmi(smi,all_mols))
                    summary = f"{name}"
                    for i in range(len(statuss)):
                        summary += f",({mol2files[i].name},{statuss[i]})"
                    successsums.append(summary)
                else:
                    failsums.append(name)
            except:
                print(name)
        # break
    with open("success.csv","w") as f:
        f.write("\n".join(successsums))
    with open("fail.csv","w") as f:
        f.write("\n".join(failsums))
```
**2.db2_to_mol2_teb.py**
```
#!/pubhome/jqzang02/miniconda3/envs/dock37_test/bin/python

#/home/xiaqc/mambaforge/envs/dock37/bin/python2.7

#################################################################################################################
##
## This libary for reading in db2 file
## 
#################################################################################################################
## Writen by Trent Balius in the Shoichet Lab, UCSF in 2015
#################################################################################################################


import math, sys
import os.path
import cmath
from math import sqrt

import mol2
#import mol2_debug as mol2

#################################################################################################################
#################################################################################################################
# data structure to store information about each residue with the docked ligand.
class db2Mol:
    def __init__(self,header,atom_list,bond_list,coord_list,seg_list,conformer_list):
        self.header         = str(header)
        #self.name           = str(name)
        self.atom_list      = atom_list
        self.bond_list      = bond_list
        self.coord_list     = coord_list
        self.seg_list       = seg_list
        self.conformer_list = conformer_list

class db2atom: # A
    def __init__(self,Q,type,name,num):
        #self.X = float(X)
        #self.Y = float(Y)
        #self.Z = float(Z)
        self.Q = float(Q)
        self.heavy_atom = False
        self.type = type
        self.name = name
        self.num  = int(num)
	#self.resnum  = int(resnum)
	#self.resname = resname
class db2bond: # B
     def __init__(self,a1_num,a2_num,num,type):
        self.a1_num = int(a1_num)
        self.a2_num = int(a2_num)
        self.num = int(num)
        self.type = type
class db2coords: # X
    def __init__(self,num,atomnum,segnum,X,Y,Z):
        self.num = int(num)
        self.atomnum = int(atomnum)
        self.segnum = int(segnum)
        self.X = float(X)
        self.Y = float(Y)
        self.Z = float(Z)
class db2segment: # 
    def __init__(self,num,start,stop):
        self.num = int(num)
        self.start = int(start)
        self.stop = int(stop)
class db2conformer: # C
    def __init__(self,num,seglist):
        self.num = int(num)
        self.seglist = seglist
     


#################################################################################################################
#################################################################################################################
#def read_Mol2_filehandel(filehandel,startline):
#    lines  =  filehandel.readlines()
#def read_Moldb2_lines(lines,startline):
def read_Moldb2_file(file):
    # reads in data from multi-Mol2 file.

#T ## namexxxx (implicitly assumed to be the standard 7)
#M zincname protname #atoms #bonds #xyz #confs #sets #rigid #Mlines #clusters
#M charge polar_solv apolar_solv total_solv surface_area
#M smiles
#M longname
#[M arbitrary information preserved for writing out]
#A stuff about each atom, 1 per line 
#B stuff about each bond, 1 per line
#X coordnum atomnum confnum x y z 
#R rigidnum color x y z
#C confnum coordstart coordend
#S setnum #lines #confs_total broken hydrogens omega_energy
#S setnum linenum #confs confs [until full column]
#D clusternum setstart setend matchstart matchend #additionalmatching
#D matchnum color x y z
#E 


    # reads in data from multi-Mol2 file.

    file1 = open(file,'r')
    lines  =  file1.readlines()
    file1.close()

    mol_list = []

    header = ''
    for line in lines:
         linesplit = line.split() #split on white space
         if(line[0] == "M"): # ATOM 
             header = header + line[1:-1]
             atomlist  = []
             bondlist  = []
             coordlist = []
             seglist   = []
             conflist  = []

         elif(line[0] == "A"): # ATOM 
            #print line
            #print line[0]
            atomnum    = linesplit[1]
            atomname   = linesplit[2]
            atomtype   = linesplit[3]
            atomcharge = linesplit[6]
            tempatom = db2atom(atomcharge,atomtype,atomname,atomnum)
            atomlist.append(tempatom)
            #print atomnum, atomname, atomtype, "q = ", atomcharge

         elif(line[0] == "B"): # BOND 
            #print line
            #print line[0]
            bondnum  = linesplit[1]
            atom1 = linesplit[2]
            atom2 = linesplit[3]
            bondtype = linesplit[4]
            #print "atom1,atom2,bondnum,bondtype = [" , atom1,atom2,bondnum,bondtype, "]"
            tempbond = db2bond(atom1,atom2,bondnum,bondtype)
            bondlist.append(tempbond)
            #print bondnum, atom1,atom2, bondtype
         elif(line[0] == "X"): # COORDS
            #exit()
            #print line
            #print line[0]
            coordnum = linesplit[1]
            atomnum  = linesplit[2]
            segnum   = linesplit[3]
            X        = linesplit[4]
            Y        = linesplit[5]
            Z        = linesplit[6]
            temp_coord = db2coords(coordnum,atomnum,segnum,X,Y,Z)
            coordlist.append(temp_coord)
            #print coordnum,X,Y,Z
         #elif(line[0] == "R"): # Rigid
         #   print line
         elif(line[0] == "C"): # Segment 
            #print line
            #print line[0]
            confnum    = linesplit[1]
            coordstart = linesplit[2]
            coordstop  = linesplit[3]
            #print confnum, coordstart, coordstop 
            tempseg = db2segment(confnum, coordstart, coordstop)
            seglist.append(tempseg)
            numold = 1
            fristflag = True
         elif(line[0] == "S"): # set -- Conformer 
            #print line
            num = int(linesplit[1])
            num2 = int(linesplit[2])
            #print numold, num
            if (fristflag):
                fristflag = False
                segnum_list = []
            elif (numold != num): # we know when it is a new conf when this number changes. 
                #print "new conformation" 
                tempconf = db2conformer(num,segnum_list)
                conflist.append(tempconf)
                segnum_list = []
                # This fist line does not contain the segment information
                # The second, and higher lines have more information
            else: # there may be multiple lines for enumerating sagments for one conformer. 
                #print "continueing, size of segnum_list = " + str(len(segnum_list))
                numofseg = linesplit[3]
                #print numofseg, len(linesplit) 
                for i in range(4,len(linesplit)):
                    segnum_list.append(int(linesplit[i]))
            numold = num
         elif(line[0] == "E"): # ATOM 
             #if (len(segnum_list) > 0): # this is to put the last conformation in the the list
             tempconf = db2conformer(num,segnum_list)
             conflist.append(tempconf)

             print "atomnum =", len(atomlist)
             print "bondnum =", len(bondlist)
             print "coordnum =", len(coordlist)
             print "segnum =", len(seglist)
             print "confnum =", len(conflist)
             tempmol = db2Mol(header, atomlist, bondlist, coordlist, seglist, conflist)  # this is an ensomble of conformation 
             header = ''
             mol_list.append(tempmol)
             #exit()
         else:
             print "Warrning: " + line[0] + " is not found in the if statments. "
             #exit()

    return mol_list


#################################################################################################################
#################################################################################################################

def convert_db2_to_mol2(db2mols):

    allmol2s = []
    # loop over each molecule
    for mol in db2mols: 
         # loop over each conformer in the molcule
         mol2mols = []
         for conf in mol.conformer_list:
              #print conf.seglist
              # the conformer is defined as a set of segement of the molecule
              mol2atomlist = []
              residue_list = {}
              for segint in conf.seglist:
                  segment =  mol.seg_list[segint-1]
                  print segment.num, segment.start, segment.stop
                  # the segement point to a bunch of coordenates, we know what atom the coordenate coresponds to. 
                  #print segment.start,segment.stop, range(segment.start,segment.stop)
                  for coordint in range(segment.start,segment.stop+1):
                      coord = mol.coord_list[coordint-1]
                      print coord.num, coord.atomnum, coord.segnum,coord.X,coord.Y,coord.Z
                      tempatom = mol.atom_list[coord.atomnum-1]
                      #X,Y,Z,Q,type,name,num,resnum,resname):
                      res_num = 1
                      resname = "lig"
                      mol2atom = mol2.atom(coord.X,coord.Y,coord.Z,tempatom.Q,tempatom.type,tempatom.name,tempatom.num,res_num,resname)
                      if residue_list.has_key(res_num):
                         residue_list[res_num].append(mol2atom)
                      else:
                         residue_list[res_num] = [mol2atom]
                      #residue_list[res_num] = [tempatom]
                      mol2atomlist.append(mol2atom)
              mol2bondlist = []
              for bond in mol.bond_list: 
                  #,a1_num,a2_num,num,type
                  mol2bond = mol2.bond(bond.a1_num,bond.a2_num,bond.num,bond.type)
                  mol2bondlist.append(mol2bond)
              mol2mol = mol2.Mol("","",mol2atomlist,mol2bondlist,residue_list)
              mol2mols.append(mol2mol)
         allmol2s.append(mol2mols)
         #exit()
         #return allmol2s
    return allmol2s
#################################################################################################################
#################################################################################################################
def main():
    if len(sys.argv) != 3: # if no input
        print "Syntax: python db2_to_mol2.py db2file.db2 outprefix ";
        return

    filename  = sys.argv[1];
    outfileprefix  = sys.argv[2];

    print "filename: " + filename 

    db2mols = read_Moldb2_file(filename)
    allmol2s = convert_db2_to_mol2(db2mols)
    filecount = 0
    for mol2mols in allmol2s:
       outfilename = outfileprefix + "." + str(filecount) + ".mol2"
       print "writing "+ outfilename 
       count = 0
       for mol2mol in mol2mols:
          if count == 0: 
             mol2.write_mol2(mol2mol,outfilename)
          else:
             mol2.append_mol2(mol2mol,outfilename)
          count = count + 1
       filecount = filecount + 1
    return;
#################################################################################################################
#################################################################################################################
main()
```
**3.db2gz_to_mol2.sh**
```
# Usage: /pubhome/qcxia02/db22mol2/db2gz_to_mol2.sh $1

zcat $1 > $(basename $1 .gz)
/{home}/test_steps/db2utils/db2_to_mol2_teb.py $(basename $1 .gz) $(basename $1 .db2.gz)
cat $(basename $1 .db2.gz).*.mol2 > $(basename $1 .db2.gz).mol2
# cat $(basename $1 .db2.gz).*.mol2 > $(basename $1 .db2.gz)-tot.mol2

rm $(basename $1 .db2.gz).*.mol2
rm $(basename $1 .gz)
```
**4.dock38todock37.py**
```
#!/home/xiaqc/mambaforge/envs/basic/bin/python3.10
import sys
import gzip

infile = sys.argv[1]
outfile = sys.argv[2]
with gzip.open(infile, "rb") as f:
    file_content = f.read().decode("utf-8")

lines = []
for line in file_content.split("\n"):
    if line.startswith("S"):
        # if isfloat(line.split()[-1]):
        if len(line.split()[-1]) > 4 and "." in line.split()[-1]:
            line = line[:line.rindex(line.split()[-1])].strip()
    lines.append(line)

dock37_file_content = "\n".join(lines).encode("utf-8")
with gzip.open(outfile, "wb") as f:
    f.write(dock37_file_content)
```
**5.mol2.py**
```
#!/user/bin/python

#################################################################################################################
##
## This libary was writen by Trent Balius and Sudipto Mukherjee in 
## the Rizzo Research Group at Stony Brook University released in 2012
## 
#################################################################################################################
## Modified by Trent Balius in the Shoichet Lab, UCSF in 2013
#################################################################################################################


import math, sys
import os.path
import cmath
from math import sqrt

#################################################################################################################
#################################################################################################################
# data structure to store information about each residue with the docked ligand.
class Mol:
    def __init__(self,header,name,atom_list,bond_list,residue_list):
        self.header       = str(header)
        self.name         = str(name)
        self.atom_list    = atom_list
        self.bond_list    = bond_list
        self.residue_list = residue_list

class atom:
    def __init__(self,X,Y,Z,Q,type,name,num,resnum,resname):
        self.X = float(X)
        self.Y = float(Y)
        self.Z = float(Z)
        self.Q = float(Q)
        self.heavy_atom = False
        self.type = type
        self.name = name
        self.num  = int(num)
        self.resnum  = int(resnum)
        self.resname = resname
class bond:
     def __init__(self,a1_num,a2_num,num,type):
        self.a1_num = int(a1_num)
        self.a2_num = int(a2_num)
        self.num = int(num)
        self.type = type
class residue:
     def __init__(self,atom_list,resnum,resname):
        self.atom_list = atom_list
        self.resnum  = int(resnum)
        self.resname = resname


#################################################################################################################
#################################################################################################################
#def read_Mol2_filehandel(filehandel,startline):
#    lines  =  filehandel.readlines()
def read_Mol2_lines(lines,startline):
    # reads in data from multi-Mol2 file.

    #print lines[startline]
    Name = ''
    atom_list = []
    bond_list = []
    residue_list = {}

    flag_atom    = False
    flag_bond    = False
    flag_substr  = False
    flag_nextmol = False
    flag_mol_set = False
    flag_mol     = False
    flag_getName = False
    

    #data = Mol('',[],[],[])

    i = 0  # i is the num of molecules read so far
    lnum = 0 
    #print len(lines)

    for lnum in range(startline,len(lines)):
         line = lines[lnum]
         linesplit = line.split() #split on white space
         if line[0] == "#":
            #print line
            continue
         if (len(linesplit) == 1):
            if(linesplit[0] == "@<TRIPOS>MOLECULE"):
               flag_mol_set = True

               if (flag_bond or flag_substr):
                   #print "I AM HERE"
                   flag_nextmol = True
               i = i + 1
               #print "READING IN MOL #" + str(i)
               #print "read in molecule info:"
               line_num = 0
               flag_mol = True
               flag_atom = False
               flag_bond = False
               flag_substr = False

            if(linesplit[0] == "@<TRIPOS>ATOM"):
               #print "read in atom info:"
               flag_atom = True
               flag_bond = False
               flag_substr = False
               flag_mol = False

            if(linesplit[0] == "@<TRIPOS>BOND"):
               #print "read in bond info:"
               flag_bond = True
               flag_substr = False
               flag_mol = False
               flag_atom = False

            if(linesplit[0] == "@<TRIPOS>SUBSTRUCTURE"):
               #print "read in substructure info:"
               flag_substr = True
               flag_mol = False
               flag_atom = False
               flag_bond = False
         if (flag_mol and (not flag_getName) and len(linesplit)==1 ):
             if (line_num == 1):
                #line_num = 0
                Name = linesplit[0]
                flag_getName = True
             line_num = line_num + 1

         if ((len(linesplit) >= 9 )and (flag_atom)):
             atom_num  = linesplit[0]
             atom_name = linesplit[1]
             X         = linesplit[2]
             Y         = linesplit[3]
             Z         = linesplit[4]
             atom_type = linesplit[5]
             res_num   = int(linesplit[6])
             res_name  = linesplit[7]
             Q         = linesplit[8]
             temp_atom = atom(X,Y,Z,Q,atom_type,atom_name,atom_num,res_num,res_name)
             atom_list.append(temp_atom)
             if residue_list.has_key(res_num):
                      residue_list[res_num].append(temp_atom)
             else:
                 residue_list[res_num] = [temp_atom]

         elif (len(linesplit) == 4 and flag_bond):
             #print line
             bond_num  = linesplit[0]
             a1_num    = linesplit[1]
             a2_num    = linesplit[2]
             bond_type = linesplit[3]
             temp_bond = bond(a1_num,a2_num,bond_num,bond_type)
             bond_list.append(temp_bond)

         #elif (flag_substr or flag_nextmol ):
         elif ( flag_nextmol ): ## we will braek the loop when we hit the next molecule
                 #ID_heavy_atoms(atom_list)
                 #data = Mol(Name,atom_list,bond_list,residue_list)
                 flag_getName = False
                 flag_substr = False
                 flag_nextmol = False
                 #atom_list = [];bond_list = []
                 #if (lnum != startline):
                 print (lnum, startline )
                 break
 
    ## we are reading in one molecule at a time
    ID_heavy_atoms(atom_list)
    data = Mol('',Name,atom_list,bond_list,residue_list)
    atom_list = [];bond_list = []
    print ("flag_mol_set", flag_mol_set )

    return flag_mol_set, data, lnum

#################################################################################################################
#################################################################################################################
def read_Mol2_file(file):
    # reads in data from multi-Mol2 file.

    file1 = open(file,'r')
    lines  =  file1.readlines()
    file1.close()

    atom_list = []
    bond_list = []
    residue_list = {}
    mol_list = []

    flag_atom      = False
    flag_bond      = False
    flag_substr    = False
    flag_mol       = False
    #flag_getName  = False
    flag_frist_mol = True

    i = 0  # i is the num of molecules read so far
    for line in lines:
         linesplit = line.split() #split on white space
         if line[0] == "#":
            #print line
            continue
         if (len(linesplit) == 1):
            if(linesplit[0] == "@<TRIPOS>MOLECULE"):
               if flag_frist_mol:
                  flag_frist_mol = False
               else: # when we have come to a new molecule put the pervious on on the list and reset arrays
                  ID_heavy_atoms(atom_list)
                  data = Mol('',Name,atom_list,bond_list,residue_list)
                  mol_list.append(data)
                  atom_list = [];bond_list = []
               i = i + 1
               #print "READING IN MOL #" + str(i)
               #print "read in molecule info:"
               line_num = 0
               flag_mol = True
               flag_atom = False
               flag_bond = False
               flag_substr = False

            if(linesplit[0] == "@<TRIPOS>ATOM"):
               #print "read in atom info:"
               flag_atom = True
               flag_bond = False
               flag_substr = False
               flag_mol = False

            if(linesplit[0] == "@<TRIPOS>BOND"):
               #print "read in bond info:"
               flag_bond = True
               flag_substr = False
               flag_mol = False
               flag_atom = False

            if(linesplit[0] == "@<TRIPOS>SUBSTRUCTURE"):
               #print "read in substructure info:"
               flag_substr = True
               flag_mol = False
               flag_atom = False
               flag_bond = False
         if (flag_mol and  len(linesplit) >= 1 ):
             if (line_num == 1):
                #print line 
                #print linesplit
                #line_num = 0
                Name = linesplit[0]
                #flag_getName = True
             line_num = line_num + 1

         if ((len(linesplit) >= 9 )and (flag_atom)):
             atom_num  = linesplit[0]
             atom_name = linesplit[1]
             X         = linesplit[2]
             Y         = linesplit[3]
             Z         = linesplit[4]
             atom_type = linesplit[5]
             res_num   = int(linesplit[6])
             res_name  = linesplit[7]
             Q         = linesplit[8]
             #print X,Y,Z,Q,atom_type,atom_name,atom_num,res_num,res_name
             temp_atom = atom(X,Y,Z,Q,atom_type,atom_name,atom_num,res_num,res_name)
             atom_list.append(temp_atom)
             if residue_list.has_key(res_num):
                      residue_list[res_num].append(temp_atom)
             else:
                 residue_list[res_num] = [temp_atom]

         elif (len(linesplit) == 4 and flag_bond):
             bond_num  = linesplit[0]
             a1_num    = linesplit[1]
             a2_num    = linesplit[2]
             bond_type = linesplit[3]
             #print a1_num,a2_num,bond_num,bond_type
             temp_bond = bond(a1_num,a2_num,bond_num,bond_type)
             bond_list.append(temp_bond)

         #elif (flag_substr):
         #        ID_heavy_atoms(atom_list)
         #        data = Mol(Name,atom_list,bond_list,residue_list)
         #        mol_list.append(data)
         #        #flag_getName = False
         #        flag_substr = False
         #        atom_list = [];bond_list = []
    # for the last molecule.
    ID_heavy_atoms(atom_list)
    data = Mol('',Name,atom_list,bond_list,residue_list)
    mol_list.append(data)
    atom_list = [];bond_list = []
    return mol_list
#################################################################################################################
#################################################################################################################
def read_Mol2_file_head(file):
    # reads in data from multi-Mol2 file.

    file1 = open(file,'r')
    lines  =  file1.readlines()
    file1.close()

    atom_list = []
    bond_list = []
    residue_list = {}
    mol_list = []

    flag_atom      = False
    flag_bond      = False
    flag_substr    = False
    flag_mol       = False
    #flag_getName  = False
    flag_frist_mol = True

    header1 = ''
    header2 = '' # we need two so we can recored the new header and print the perceding one. 
                 # we read in header i, then we read in mol i, then we read in header i+1 and then store mol i and header i.  

    i = 0  # i is the num of molecules read so far
    for line in lines:
         if (line[0] == '#'): 
             #print line
             header1 = header1 + line
             continue

         linesplit = line.split() #split on white space
         if (len(linesplit) == 1):
            if(linesplit[0] == "@<TRIPOS>MOLECULE"):
               if flag_frist_mol:
                  flag_frist_mol = False
               else: # when we have come to a new molecule put the pervious on on the list and reset arrays
                  #print header2
                  ID_heavy_atoms(atom_list)
                  data = Mol(header2,Name,atom_list,bond_list,residue_list)
                  mol_list.append(data)
                  atom_list = [];bond_list = []
               header2 = header1
               header1 = ''
               i = i + 1
               #print "READING IN MOL #" + str(i)
               #print "read in molecule info:"
               line_num = 0
               flag_mol = True
               flag_atom = False
               flag_bond = False
               flag_substr = False

            if(linesplit[0] == "@<TRIPOS>ATOM"):
               #print "read in atom info:"
               flag_atom = True
               flag_bond = False
               flag_substr = False
               flag_mol = False

            if(linesplit[0] == "@<TRIPOS>BOND"):
               #print "read in bond info:"
               flag_bond = True
               flag_substr = False
               flag_mol = False
               flag_atom = False

            if(linesplit[0] == "@<TRIPOS>SUBSTRUCTURE"):
               #print "read in substructure info:"
               flag_substr = True
               flag_mol = False
               flag_atom = False
               flag_bond = False
         if (flag_mol and  len(linesplit) >= 1 ):
             if (line_num == 1):
                #print line 
                #print linesplit
                #line_num = 0
                Name = linesplit[0]
                #flag_getName = True
             line_num = line_num + 1

         if ((len(linesplit) >= 9 )and (flag_atom)):
             atom_num  = linesplit[0]
             atom_name = linesplit[1]
             X         = linesplit[2]
             Y         = linesplit[3]
             Z         = linesplit[4]
             atom_type = linesplit[5]
             res_num   = int(linesplit[6])
             res_name  = linesplit[7]
             Q         = linesplit[8]
             #print X,Y,Z,Q,atom_type,atom_name,atom_num,res_num,res_name
             temp_atom = atom(X,Y,Z,Q,atom_type,atom_name,atom_num,res_num,res_name)
             atom_list.append(temp_atom)
             if residue_list.has_key(res_num):
                      residue_list[res_num].append(temp_atom)
             else:
                 residue_list[res_num] = [temp_atom]

         elif (len(linesplit) == 4 and flag_bond):
             bond_num  = linesplit[0]
             a1_num    = linesplit[1]
             a2_num    = linesplit[2]
             bond_type = linesplit[3]
             #print a1_num,a2_num,bond_num,bond_type
             temp_bond = bond(a1_num,a2_num,bond_num,bond_type)
             bond_list.append(temp_bond)

         #elif (flag_substr):
         #        ID_heavy_atoms(atom_list)
         #        data = Mol(Name,atom_list,bond_list,residue_list)
         #        mol_list.append(data)
         #        #flag_getName = False
         #        flag_substr = False
         #        atom_list = [];bond_list = []
    # for the last molecule.
    ID_heavy_atoms(atom_list)
    #print header2
    data = Mol(header2,Name,atom_list,bond_list,residue_list)
    mol_list.append(data)
    atom_list = [];bond_list = []
    return mol_list
#################################################################################################################
#################################################################################################################
def write_mol2(molecule,filename):

        # define a dictionary for help renumbering
        atom_dic = {}
        resid_dic = {}
        count = 1
        for atom in molecule.atom_list:
            if not atom_dic.has_key(atom.num):
               atom_dic[atom.num] = count
               #print atom.num, ",", count,",", atom_dic[atom.num]
               count=count+1
        count = 1
        for resnum in molecule.residue_list.keys():
            resid_dic[resnum] = count
            count=count+1

        outmol2 = open(filename,'w')
        #print molecule.header
        outmol2.write(molecule.header)      #dock info after #s 
        outmol2.write("@<TRIPOS>MOLECULE\n")      #start the MOLECULE RTI (Record Type Indicator)
        outmol2.write(molecule.name+'\n')         #print MOL2FILE name of the molecule
        outmol2.write("%-5d %-5d %-5d 0     0\n" % (len(molecule.atom_list), 
                len(molecule.bond_list), len(molecule.residue_list.keys()))) 
        # For now, the number of residues is hard-coded to 1. To be fixed.
        outmol2.write("SMALL\n")                   #mol_type
        outmol2.write("USER_CHARGES\n")           #charge_type

        outmol2.write("@<TRIPOS>ATOM\n")      #start the ATOM RTI (Record Type Indicator)
        for j in range(0,len(molecule.atom_list)):
                #print atom_dic[molecule.atom_list[j].num], molecule.atom_list[j].num
                outmol2.write("%-6d %-4s %9.4f %9.4f %9.4f %-5s %4s %6s %9.4f\n" % 
                (atom_dic[molecule.atom_list[j].num], molecule.atom_list[j].name, molecule.atom_list[j].X, molecule.atom_list[j].Y, 
                molecule.atom_list[j].Z, molecule.atom_list[j].type, resid_dic[molecule.atom_list[j].resnum], 
                molecule.atom_list[j].resname, molecule.atom_list[j].Q))

        outmol2.write("@<TRIPOS>BOND\n")
        count = 1
        for m in range(0,len(molecule.bond_list)):
                outmol2.write("%-5d %-5d %-5d %s\n" % (count, 
                atom_dic[molecule.bond_list[m].a1_num], atom_dic[molecule.bond_list[m].a2_num], molecule.bond_list[m].type))
                count = count + 1

        outmol2.write("@<TRIPOS>SUBSTRUCTURE\n")
        count = 1
        for resnum in molecule.residue_list.keys():
                #outmol2.write("%-3d %-5s %-5d RESIDUE    1   A     %-5s 1\n" % (resnum, 
                outmol2.write("%-3d %-5s %-5d RESIDUE    1   A     %-5s 1\n" % (resid_dic[resnum], 
                molecule.residue_list[resnum][0].resname, # residue name 
                atom_dic[molecule.residue_list[resnum][0].num], molecule.residue_list[resnum][0].resname[0:3]))   # atom num of first atom in this residue
                count = count + 1
        outmol2.close()
        return
#################################################################################################################
def append_mol2(molecule,filename):

        # define a dictionary for help renumbering
        atom_dic = {}
        resid_dic = {}
        count = 1
        for atom in molecule.atom_list:
            if not atom_dic.has_key(atom.num):
               atom_dic[atom.num] = count
               #print atom.num, ",", count,",", atom_dic[atom.num]
               count=count+1
        count = 1
        for resnum in molecule.residue_list.keys():
            resid_dic[resnum] = count
            count=count+1

        outmol2 = open(filename,'a')
        #print molecule.header
        outmol2.write(molecule.header)      #dock info after #s 
        outmol2.write("@<TRIPOS>MOLECULE\n")      #start the MOLECULE RTI (Record Type Indicator)
        outmol2.write(molecule.name+'\n')         #print MOL2FILE name of the molecule
        outmol2.write("%-5d %-5d %-5d 0     0\n" % (len(molecule.atom_list),
                len(molecule.bond_list), len(molecule.residue_list.keys())))
        # For now, the number of residues is hard-coded to 1. To be fixed.
        outmol2.write("SMALL\n")                  #mol_type
        outmol2.write("USER_CHARGES\n")           #charge_type

        outmol2.write("@<TRIPOS>ATOM\n")      #start the ATOM RTI (Record Type Indicator)
        for j in range(0,len(molecule.atom_list)):
                #print atom_dic[molecule.atom_list[j].num], molecule.atom_list[j].num
                outmol2.write("%-6d %-4s %9.4f %9.4f %9.4f %-5s %4s %6s %9.4f\n" %
                (atom_dic[molecule.atom_list[j].num], molecule.atom_list[j].name, molecule.atom_list[j].X, molecule.atom_list[j].Y,
                molecule.atom_list[j].Z, molecule.atom_list[j].type, resid_dic[molecule.atom_list[j].resnum],
                molecule.atom_list[j].resname, molecule.atom_list[j].Q))

        outmol2.write("@<TRIPOS>BOND\n")
        count = 1
        for m in range(0,len(molecule.bond_list)):
                outmol2.write("%-5d %-5d %-5d %s\n" % (count,
                atom_dic[molecule.bond_list[m].a1_num], atom_dic[molecule.bond_list[m].a2_num], molecule.bond_list[m].type))
                count = count + 1

        outmol2.write("@<TRIPOS>SUBSTRUCTURE\n")
        count = 1
        for resnum in molecule.residue_list.keys():
                #outmol2.write("%-3d %-5s %-5d RESIDUE    1   A     %-5s 1\n" % (resnum, 
                outmol2.write("%-3d %-5s %-5d RESIDUE    1   A     %-5s 1\n" % (resid_dic[resnum],
                molecule.residue_list[resnum][0].resname, # residue name 
                atom_dic[molecule.residue_list[resnum][0].num], molecule.residue_list[resnum][0].resname[0:3]))   # atom num of first atom in this residue
                count = count + 1
        outmol2.close()
        return

#################################################################################################################
# this fuction will convert the sybyl atomtypes to the DOCK3.7 atomtype.
# dictionary was obtained from Ryan's mol2db2.
#################################################################################################################

def convert_sybyl_to_dock (molecule):

  convertTypesDefault = {'C.3': 5,
                         'C.2': 1,
                         'C.ar': 1,
                         'C.1': 1,
                         'N.3': 10,
                         'N.2': 8,
                         'N.1': 8,
                         'O.3': 12,
                         'O.2': 11,
                         'S.3': 14,
                         'N.ar': 8,
                         'P.3': 13,
                         'H': 6,
                         'H-C': 7,
                         'Br': 17,
                         'Cl': 16,
                         'F': 15,
                         'I': 18,
                         'S.2': 14,
                         'N.pl3': 8,
                         'LP': 25,
                         'Na': 19,
                         'K': 19,
                         'Ca': 21,
                         'Li': 20,
                         'Al': 20,
                         'Du': 25,
                         'Du.C': 25,
                         'Si': 24,
                         'N.am': 8,
                         'S.o': 14,
                         'S.o2': 14,
                         'N.4': 9,
                         'O.co2': 11,
                         'C.cat': 1,
                         'H.spc': 6,
                         'O.spc': 11,
                         'H.t3p': 6,
                         'O.t3p': 11,
                         'ANY': 25,
                         'HEV': 25,
                         'HET': 25,
                         'HAL': 25,
                         'Mg': 20,
                         'Cr.oh': 25,
                         'Cr.th': 25,
                         'Se': 25,
                         'Fe': 25,
                         'Cu': 25,
                         'Zn': 26,
                         'Sn': 25,
                         'Mo': 25,
                         'Mn': 25,
                         'Co.oh': 25}
  #i = 0
  i = 1 # atoms labels in bonds start at one.
  docktype = []
  # loop over all atoms in molecule
  for atom in molecule.atom_list:
      type = atom.type
      if atom.type == 'H': 
         # if atom type is a hydrogen, loop over the bonds, see if it is attached to a carbon
         hflag = False
         for bond in molecule.bond_list:
             #print "looking at bond ", bond.num, " (", bond.a1_num, bond.a2_num,")"
             if (i == bond.a1_num):
                 #print i, "->", bond.a2_num
                 #print "bond (", molecule.atom_list[bond.a1_num-1].type, molecule.atom_list[bond.a2_num-1].type,")"
                 j = bond.a2_num
                 hflag = True
                 break # hydrogens are only attached to on atom. 
             elif (i == bond.a2_num):
                 #print i, "<-", bond.a1_num
                 #print "bond (", molecule.atom_list[bond.a1_num-1].type, molecule.atom_list[bond.a2_num-1].type,")"
                 j = bond.a1_num
                 hflag = True
                 break # hydrogens are only attached to on atom. 
         if hflag:
            if molecule.atom_list[j-1].type in [ 'C.1', 'C.3', 'C.2', 'C.ar', 'C.cat']:
               print (i, atom.type, j, molecule.atom_list[j-1].type)
               type = 'H-C'
            #else: 
            #   type = atom.type
         else: 
            print ("ERROR.")
            exit()
      #else:
      #   type = atom.type
      docktype.append(convertTypesDefault[type])
         
      i = i + 1   
  return docktype


#################################################################################################################
def get_pdbcode_list(filename):
    systems_list = open(filename,'r')
    lines  =  systems_list.readlines()
    return lines          
#################################################################################################################
def ID_heavy_atoms(atom_list):
    for i in range(len(atom_list)):
        if (atom_list[i].type[0] != 'H'):
            atom_list[i].heavy_atom = True
    return atom_list
#################################################################################################################
#################################################################################################################
def distance2_vec(vector1,vector2):
        if (len(vector1)!=len(vector2)):
                print ('function distance(): vectors differ in length')
                sys.exit(1)
        distance2 = 0
        for i in range(len(vector1)):
                distance2 += (vector1[i]-vector2[i])**2
        return distance2
##################################################################################################################
##################################################################################################################
#def norm(vector1):
#        norm = 0
#        for i in range(len(vector1)):
#                norm += (vector1[i])*(vector1[i])
#        return sqrt(norm)
##################################################################################################################
##################################################################################################################
def distance2(atom1,atom2):
    return (atom1.X - atom2.X )**2 + (atom1.Y - atom2.Y )**2 + (atom1.Z - atom2.Z )**2
#################################################################################################################
#################################################################################################################
# Make sure the heavy atoms are being declared as heavy
# i.e call ID_heavy atoms function
def heavy_atom_RMSD(ref,pose):
    if (len(ref.atom_list) != len(pose.atom_list)):
       return -1 # when atom numbers do not agree
    sum = 0.0
    num_hvy_atoms = 0
    for i in range(len(ref.atom_list)):
        if (ref.atom_list[i].heavy_atom and pose.atom_list[i].heavy_atom):
           sum += distance2(ref.atom_list[i],pose.atom_list[i])
           num_hvy_atoms+=1
    return sqrt(sum/num_hvy_atoms)

#################################################################################################################
#################################################################################################################
def formal_charge(molecule):
        total = 0
        for i in range(len(molecule.atom_list)):
                total += molecule.atom_list[i].Q
        return total
#################################################################################################################
def centre_of_mass(molecule):
        # Dictionary of atomic weights of elements
        atom_mass = {'O':15.9994 ,'N':14.00674 ,'C':12.011 ,'F':18.9984032 ,'Cl':35.4527 ,'Br':79.904
        ,'I':126.90447 ,'H':1.00794 ,'B':10.811 ,'S':32.066 ,'P':30.973762 ,'Li':6.941 ,'Na':22.98968
        ,'Mg':24.3050 ,'Al':26.981539 ,'Si':28.0855 ,'K':39.0983 ,'Ca':40.078 ,'Cr':51.9961 ,'Mn':54.93805
        ,'Fe':55.847 ,'Co':58.93320 ,'Cu':63.546 ,'Zn':65.39 ,'Se':78.96 ,'Mo':95.94 ,'Sn':118.710 ,'LP':0.0 }

        cmass = [0,0,0]
        centroid = [0,0,0]
        molecular_weight = 0
        for k in range(0,len(molecule.atom_list)):
                element = molecule.atom_list[k].type.split('.')[0]
                cmass[0] += molecule.atom_list[k].X * atom_mass[element]
                cmass[1] += molecule.atom_list[k].Y * atom_mass[element]
                cmass[2] += molecule.atom_list[k].Z * atom_mass[element]
                centroid[0] += molecule.atom_list[k].X
                centroid[1] += molecule.atom_list[k].Y
                centroid[2] += molecule.atom_list[k].Z
                molecular_weight += atom_mass[element]
        #print "Molecular Weight =",molecular_weight
        cmass[0] /= molecular_weight
        cmass[1] /= molecular_weight
        cmass[2] /= molecular_weight
        centroid[0] /= len(molecule.atom_list)
        centroid[1] /= len(molecule.atom_list)
        centroid[2] /= len(molecule.atom_list)
        #print 'Centroid =',centroid
        return cmass
#################################################################################################################
def molecular_weight(molecule):
        # Dictionary of atomic weights of elements
        atom_mass = {'O':15.9994 ,'N':14.00674 ,'C':12.011 ,'F':18.9984032 ,'Cl':35.4527 ,'Br':79.904
        ,'I':126.90447 ,'H':1.00794 ,'B':10.811 ,'S':32.066 ,'P':30.973762 ,'Li':6.941 ,'Na':22.98968
        ,'Mg':24.3050 ,'Al':26.981539 ,'Si':28.0855 ,'K':39.0983 ,'Ca':40.078 ,'Cr':51.9961 ,'Mn':54.93805
        ,'Fe':55.847 ,'Co':58.93320 ,'Cu':63.546 ,'Zn':65.39 ,'Se':78.96 ,'Mo':95.94 ,'Sn':118.710 ,'LP':0.0 }

        molecular_weight = 0
        for k in range(0,len(molecule.atom_list)):
                element = molecule.atom_list[k].type.split('.')[0]
                molecular_weight += atom_mass[element]
        return molecular_weight
#################################################################################################################
def calc_dipole_moment(molecule):
    uIsum=0
    uJsum=0
    uKsum=0
    dipolemoment=0
    conversion = 4.796 # Convert partialcharge*angstroms --> Coulombs*meters (Debye)

    cmass = centre_of_mass(molecule)
    #print "Centre of mass = ",cmass

    #cmass = [molecule.atom_list[0].X, molecule.atom_list[0].Y, molecule.atom_list[0].Z]
    for k in range(0,len(molecule.atom_list)):
        uIsum += molecule.atom_list[k].Q * (molecule.atom_list[k].X - cmass[0])
        uJsum += molecule.atom_list[k].Q * (molecule.atom_list[k].Y - cmass[1])
        uKsum += molecule.atom_list[k].Q * (molecule.atom_list[k].Z - cmass[2])

    umag          = sqrt( (uIsum*uIsum) + (uJsum*uJsum) + (uKsum*uKsum) )
    dipolemoment  = umag*conversion;
    uvector = [uIsum,uJsum,uKsum]

    return uvector, dipolemoment
#################################################################################################################
# Takes a single Mol object and returns a Mol object without the hydrogens
# Have to remove H from atom_list, bond_list and residue_list
def remove_hydrogens(m):
    atom_list = []
    bond_list = []
    residue_list = {}

    # Retain only heavy atoms in atom_list
    num_hvy_atoms = 0
    for i in range(len(m.atom_list)):
        if (m.atom_list[i].heavy_atom):
           atom_list.append(m.atom_list[i])
           num_hvy_atoms+=1

    # Retain only bonds containing heavy atoms
    for bond_id in range(len(m.bond_list)):
        retain_bond = True
        for atom_id in range(len(m.atom_list)):
           if (m.atom_list[atom_id].heavy_atom):
              continue  
           # Atoms down here are always hydrogen 
           if (m.bond_list[bond_id].a1_num == m.atom_list[atom_id].num):
              retain_bond = False 
           if (m.bond_list[bond_id].a2_num == m.atom_list[atom_id].num):
              retain_bond = False
        if (retain_bond):
            bond_list.append(m.bond_list[bond_id])

    # Assuming that residue list does not change

    data = Mol(m.header,m.name,atom_list,bond_list,m.residue_list)
    return data
#################################################################################################################
```
---

```
cd xtal-lig
#可以在test_steps路径下创建

./db2utils/db2gz_to_mol2.sh 30Z_B_702.db2.gz

mkdir docking

readlink -f xtal-lig/30Z_B_702.db2.gz > split_database_index
#注意xtal-lig文件夹的路径，可以使用绝对路径

cp ../INDOCK .

/{home}/bela30c/docking/DOCK/src/x86_64/dock64 INDOCK
```
