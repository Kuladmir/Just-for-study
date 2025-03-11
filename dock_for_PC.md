## DOCK3.7 For PC
*Kuladmir Present*
### 0.  前言
UCSF DOCK 3.7 是一个用于分子对接的计算工具，由 UCSF（加州大学旧金山分校）Shoichet 实验室开发。它是 DOCK 3 系列的最新版本，提供了许多改进功能。DOCK 3.7 主要用于蛋白质-配体对接研究，支持从 ZINC 数据库中获取配体。
作为开源软件，学术免费，但需要引用。

From: Coleman, Carchia, Sterling, Irwin & Shoichet, PLOS ONE 2013。

同时，过程中会使用unicon/ZBH2022等软件，可以在 https://www.zbh.uni-hamburg.de/forschung/amd/software/conformator.html 获取。
### 1.  安装WSL

Windows 10及以上开放了一个功能：Windows子系统（Windows Subsystem for Linux, WSL）。

#### 1.1 
在**菜单**里输入**启用或关闭 Windows 功能**，在弹出的页面中勾选**适用于 Linux 的 Windows 子系统**和**虚拟机平台**。
#### 1.2
重启系统
#### 1.3
打开终端，输入
```
wsl.exe --update
```
#### 1.4
打开Microsoft shop，搜索Ubuntu。
安装完成后，在终端页面处，可以查看Ubuntu的Linux终端。
#### 1.5
在Linux界面，创建账号和密码，**强烈建议备份账号和密码**。
在当前时间下，下载的应该是24.04的Ubuntu，以下命令可以查看版本。
```
lsb_release -a
```
实测结果：
>Distributor ID: Ubuntu

>Description:    Ubuntu 24.04.1 LTS

>Release:        24.04

>Codename:       noble
### 2.安装Conda环境

可以在本页面详细查看安装方法：
https://www.anaconda.com/docs/getting-started/miniconda/install#macos-linux-installation

在本页面下，找到Basic install instructions部分，选择Linux installation，然后选择Linux terminal installer。然后按照流程完成安装。

最终安装完成后，命令行界面一般是这样的:
> (bash)usr_name@usr_name:

### 3.安装DOCK3.7
可以直接在：https://pan.baidu.com/s/1gTBVGuuS3vLwLDeGahYSGQ?pwd=a2is 下载。
将.zip文件压缩后，确保bela30c为根目录，然后打开该目录就能看到其他所有文件。

### 4.使用DOCK3.7
#### 4.1 创建虚拟环境
使用dock37.yml文件，创建一个虚拟环境。
dock37.yml文件内容：
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
  # 
  # Other packages may cause conflicts, raising warning like
  # 
  # Found conflicts! Looking for incompatible packages.
  # This can take several minutes.  Press CTRL-C to abort.
  # 
  # You can try something anyway.
```
之后在Linux终端里输入命令：
```
conda env create -f dock37.yml
```
创建一个dock37环境。

#### 4.2 启动环境，开始运行
首先激活环境：
```
conda activate dock37
```
准备初始文件：
rec.pdb  --受体分子的pdb文件
xtal-lig.pdb  --共晶小分子的pdb文件

案例测试：
```
mkdir test && cd test

wget -c https://files.rcsb.org/download/6N2W.pdb

# obtain rec.pdb

grep "^ATOM" 6N2W.pdb | grep " B " | grep -v "HOH" > rec.pdb

grep "^HETATM" 6N2W.pdb | grep "FE"  | grep "B" | sed "s/HETATM/ATOM  /" >> rec.pdb

# obtain xtal-lig.pdb
grep "30Z B 702" 6N2W.pdb > xtal-lig.pdb
```
#### 4.3 一步运行
```
python2.7 /home/usr/bela30c/proteins/blastermaster/blastermaster.py --addhOptions=" -HIS -FLIPs " -v
```
如果出现以下错误：
```
WARNING:root:Could not find numpy or numeric
CRITICAL:root:Could not find any matrix library (Numeric, numpy, pMatrix)! Cannot proceed
Traceback (most recent call last):
  File "/home/shianine/dock_learn/DOCK-3.7-bbe1a30c/proteins/blastermaster/blastermaster.py", line 56, in <module>
    import blasterAddHydrogens
  File "/home/shianine/dock_learn/DOCK-3.7-bbe1a30c/proteins/blastermaster/blasterAddHydrogens.py", line 11, in <module>
    import pdb
  File "/home/shianine/dock_learn/DOCK-3.7-bbe1a30c/proteins/blastermaster/pdb.py", line 8, in <module>
    import geometry  # distance function
  File "/home/shianine/dock_learn/DOCK-3.7-bbe1a30c/proteins/blastermaster/geometry.py", line 27, in <module>
    import pMatrix
ImportError: No module named pMatrix
```
则应该使用如下命令：
```
pip2 install numpy && pip2 install scipy
```
此时，需要注意，测试一下numpy是否能正常使用：
```
python2 -c "import numpy as np; print(np.__version__)"
```
如果出现以下问题：
```
Traceback (most recent call last):
  File "<string>", line 1, in <module>
  File "/home/kuladmir/miniconda3/envs/dock37/lib/python2.7/site-packages/numpy/__init__.py", line 153, in <module>
    from . import add_newdocs
  File "/home/kuladmir/miniconda3/envs/dock37/lib/python2.7/site-packages/numpy/add_newdocs.py", line 13, in <module>
    from numpy.lib import add_newdoc
  File "/home/kuladmir/miniconda3/envs/dock37/lib/python2.7/site-packages/numpy/lib/__init__.py", line 18, in <module>
    from .polynomial import *
  File "/home/kuladmir/miniconda3/envs/dock37/lib/python2.7/site-packages/numpy/lib/polynomial.py", line 19, in <module>
    from numpy.linalg import eigvals, lstsq, inv
  File "/home/kuladmir/miniconda3/envs/dock37/lib/python2.7/site-packages/numpy/linalg/__init__.py", line 50, in <module>
    from .linalg import *
  File "/home/kuladmir/miniconda3/envs/dock37/lib/python2.7/site-packages/numpy/linalg/linalg.py", line 29, in <module>
    from numpy.linalg import lapack_lite, _umath_linalg
ImportError: libgfortran.so.3: cannot open shared object file: No such file or directory
```
则需要执行：
```
sudo apt update

sudo apt install libgfortran5
```
然后测试numpy是否运行：
```
python2 -c "import numpy as np; print(np.__version__)"
```
如果输出版本，则说明正常。
#### 4.4 分布运行（首次建议）
1. 执行完4.2的命令后，执行如下命令：
```
python2 $DOCKBASE/proteins/blastermaster/fmt_rec.py

python2 $DOCKBASE/proteins/blastermaster/fmt_lig.py
```
此时，会产生一个working文件夹。
```
cd working
```
2. 加氢
---
#### 目的：给蛋白加氢；
#### 输入文件：rec.pdb  /{home}/bela30c/proteins/default/reduce_wwPDB_het_dict.txt
#### 输出文件：rec.crg.pdb.fullh    addh.log   rec.crg.pdb  rec.crg.pdb.polarH
---
```
cp /home/usr/bela30c/proteins/defaults/reduce_wwPDB_het_dict.txt .

/home/usr/bela30c/proteins/Reduce/reduce \
    -db reduce_wwPDB_het_dict.txt \
    -HIS -FLIPS \
    rec.pdb > rec.crg.pdb.fullh 2> addh.log

sed -i 's/\s*new\s*//g' rec.crg.pdb.fullh

sed -i '/^USER.*/d' rec.crg.pdb.fullh
```
※如果执行上述命令时，显示reduce无法访问，则可以执行如下命令：
```
cd /home/usr/bela30c/proteins/Reduce

chmod a+x reduce
```
这是由于文件没有执行权，因此需要添加执行权，如果后面有类似的问题，可以用相同的方法解决。

※如果还出现问题，可以先检查其动态链接库依赖：
```
ldd /home/kuladmir/bela30c/proteins/Reduce/reduce
```
如果显示：
> not a dynamic executable

则需要执行：
```
file /home/kuladmir/bela30c/proteins/Reduce/reduce
```
如果输出以下内容：
>ELF 32-bit LSB executable, Intel 80386, version 1 (SYSV), dynamically linked, interpreter /lib/ld-linux.so.2, for GNU/Linux 2.6.9, BuildID[sha1]=3eee6c155a6bcf70b90118e2da0af37bf2762121, not stripped

则需要安装32位支持库：
```
sudo apt update

sudo dpkg --add-architecture i386

sudo apt update

sudo apt install libc6:i386 libstdc++6:i386

sudo apt install lib32gcc-s1
```
然后检查是否成功：
```
ldd /home/kuladmir/bela30c/proteins/Reduce/reduce
```
如果输出：
>linux-gate.so.1 (0xf7f7e000)

>libstdc++.so.6 => /lib/i386-linux-gnu/libstdc++.so.6 (0xf7cf3000)

>libm.so.6 => /lib/i386-linux-gnu/libm.so.6 (0xf7be9000)

>libgcc_s.so.1 => /lib/i386-linux-gnu/libgcc_s.so.1 (0xf7bb1000)

>libc.so.6 => /lib/i386-linux-gnu/libc.so.6 (0xf7975000)

>/lib/ld-linux.so.2 (0xf7f80000)

则说明问题解决。不再出现问题后，继续运行加氢：
```
sed -i 's/\s*new\s*//g' rec.crg.pdb.fullh

sed -i '/^USER.*/d' rec.crg.pdb.fullh

python2 /home/usr/bela30c/proteins/blastermaster/polar_pdb.py
```
3. 基于共晶小分子做结合位点筛选（filter）
---
#### 目的：筛选出结合位点；
#### 输入文件：rec.pdb  xtal-lig.pdb  /{home}/bela30c/proteins/defaults/filt.params
#### 输出文件：filter.log   rec.site
---
```
$DOCKBASE/proteins/filt/bin/filt < $DOCKBASE/proteins/defaults/filt.params > filter.log
```
4. 生成表面（dms）

根据溶剂可及表面（Solvent Accessible Surface），默认用半径1.4 A的水探针小球在该面下滚动，生成一个溶剂排除表面（Solvent Excluded Surface）。

```
grep -a -v HOH rec.crg.pdb > rec.crg.pdb.dms

grep -a -v HOH rec.site > rec.site.dms

$DOCKBASE/proteins/dms/bin/dms rec.crg.pdb.dms -a -d 1.0 -i rec.site.dms -g dms.log -p -n -o rec.ms

$DOCKBASE/proteins/dms/bin/dms rec.crg.pdb.dms -a -d 1.0 -i rec.site.dms -g dms.ts.log -p -n -o rec.ts.ms
```
如果显示错误：
>Bus error (core dumped)

则需要将/bela30c/proteins/dms/bin里的dms dmsd文件设置为可执行。

5. 球面生成（sphgen）

sphgen会产生描述分子和分子表面形状的重叠球面，对于受体：会产生负电的表面内陷；对于配体：会产生全分子的正像。

上步生成溶剂排除表面后，球面生成时，会接触两处溶剂排除表面（类似相切），后形成一个球体。这个球面会在整个溶剂排除表面生成，每个曲面点约产生一个球体，之后会过滤较密集的，只保留与每个受体表面原子相关的最大球体。

---
#### 目的：生成描述结合口袋的负电图像；
#### 输入文件：/{home}/working/INSPH  (rec.ms)
#### 输出文件：OUTSPH  all_spheres.sph
---
```
echo "rec.ms
R
X
0.
5.0
1.4
all_spheres.sph" > INSPH

$DOCKBASE/proteins/sphgen/bin/sphgen

sed -i '1d' all_spheres.sph
```
5. 生成微球
---
#### 目的：生成低介电常数和配体去溶剂化的边界层球；
#### 输入文件：rec.ts.ms
#### 输出文件：low_die_thinspheres.sph  lig_die_thinspheres.sph
---
```
python2 $DOCKBASE/proteins/thinspheres/thin_spheres.py -i rec.ts.ms -o low_die_thinspheres.sph -d 1.0 -s 1.0

python2 $DOCKBASE/proteins/thinspheres/thin_spheres.py -i rec.ts.ms -o lig_des_thinspheres.sph -d 1.0 -s 1.0
```
6. 选择配体附近的微球

---
#### 目的：选择共晶小分子距离内的边界球；
#### 输入文件：low_die_thinspheres.sph  lig_des_thinspheres.sph
#### 输出文件：low_die_thinspheres.sph.close  lig_des_thinspheres.sph.close   low_die_thinspheres.sph.close.log   lig_des_thinspheres.sph.close.log
---
```
python2 $DOCKBASE/proteins/thinspheres/close_sph.py low_die_thinspheres.sph xtal-lig.pdb  low_die_thinspheres.sph.close 2.0 1.0

python2 $DOCKBASE/proteins/thinspheres/close_sph.py lig_des_thinspheres.sph xtal-lig.pdb  lig_des_thinspheres.sph.close 2.0 1.0
```
7. 产生正像
---
#### 目的：将共晶小分子转化为球面
#### 输入文件：xtal-lig.pdb
#### 输出文件：xtal-lig.match.sph
---
```
$DOCKBASE/proteins/pdbtosph/bin/pdbtosph xtal-lig.pdb xtal-lig.match.sph
```
※如果执行上述命令后，出现错误：
>Segmentation fault (core dumped)

这是因为原来的pdbtosph是基于CentOS的，所以会导致在Ubuntu上产生问题。解决方案：
```
sudo apt update

sudo apt install gfortran

which gfortran

cd ~/bela30c/proteins/pdbtosph/src

make
```
8. 筛选低介电常数球面

makespheres1.cli.pl

---
#### 目的：根据共晶小分子球面和描述结合口袋负电的球面，筛选出用于后续计算的球面
#### 输入文件：xtal-lig.match.sph  all_spheres.sph  rec.crg.pdb
#### 输出文件：lowdielectric.sph  lowdielectric_spheres.log
---
```
$DOCKBASE/proteins/makespheres1/makespheres1.cli.pl xtal-lig.match.sph all_spheres.sph rec.crg.pdb  lowdielectric.sph 25 >& lowdielectric.spheres.log
```
9. 提取坐标并转化

提取球体坐标，将其转化为C，并写入到pdb文件。
在运行该步骤之前，需要先安装csh。
```
sudo apt-get install csh
```
然后执行以下命令：
```
$DOCKBASE/proteins/showsphere/doshowsph.csh lowdielectric.sph 1 lowdielectric.sph.pdb >&  lowdielectric.sph.pdb.log

$DOCKBASE/proteins/showsphere/doshowsph.csh low_die_thinspheres.sph.close 1  low_die_thinspheres.sph.close.pdb >& low_die_thinspheres.sph.close.log

$DOCKBASE/proteins/showsphere/doshowsph.csh lig_des_thinspheres.sph.close 1  lig_des_thinspheres.sph.close.pdb >& lig_des_thinspheres.sph.close.log

cat rec.crg.pdb lowdielectric.sph.pdb > receptor.crg.lowdielectric.pdb
```
执行这些命令时，需要保证/home/uer/bela30c/proteins/showsphere/doshowsph.sh文件和/home/uer/bela30c/proteins/showsphere/bin/showsphere文件都为可执行文件。

10. 产生匹配球

makespheres3.cli.pl

---
#### 目的：产生匹配球
#### 输入文件：xtal-lig.match.sph  all_spheres.sph  rec.crg.pdb
#### 输出文件：matching_spheres.sph  matching_spheres.log
---
```
$DOCKBASE/proteins/makespheres3/makespheres3.cli.pl 1.5 0.8 45 xtal-lig.match.sph all_spheres.sph  rec.crg.pdb matching_spheres.sph >& matching_spheres.log
```
### 5.打分
1. 产生打分网格
---
#### 目的：生成打分需要的网格
#### 输入文件：xtal-lig.match.sph  rec.crg.pdb
#### 输出文件：box  makebox.log
---
```
$DOCKBASE/proteins/makebox/makebox.smallokay.pl xtal-lig.match.sph rec.crg.pdb box 10.0 >&  makebox.log
```
2. 用数值解法生成静电网格
---
#### 目的：用数值解法生成静电网格
#### 输入文件：/{home}/bela30c/proteins/defaults/qnifft.parm  (receptor.crg.lowdielectric.pdb,amb.crg.oxt,vdw.siz)
#### 输出文件：qnifft.atm  qnifft.electrostatics.phi  qnifft.eps  qnifft.log  qnifft_sas.usr  
---
```
cp $DOCKBASE/proteins/defaults/delphi.def qnifft.parm

echo "grid=193" >> qnifft.parm

echo "charge=amb.crg.oxt" >> qnifft.parm

echo "radius=vdw.siz" >> qnifft.parm # default extended atom radii based loosely on mike connolly's MS program

echo "pdb_input=receptor.crg.lowdielectric.pdb" >> qnifft.parm

echo "pdb_output_file=qnifft.atm" >> qnifft.parm

echo "phi_output_file=qnifft.electrostatics.phi" >> qnifft.parm

echo "sizing=border border_solvent=10." >> qnifft.parm

cp $DOCKBASE/proteins/defaults/amb.crg.oxt .

cp $DOCKBASE/proteins/defaults/vdw.siz .

$DOCKBASE/proteins/qnifft/bin/qnifft22_193_pgf_32 qnifft.parm >& qnifft.log
```
3. 休整phi-map
---
#### 目的：休整phi-map
#### 输入文件：box  qnifft.electrostatics.phi
#### 输出文件：trim.electrostatics.phi 
---

```
python2 $DOCKBASE/proteins/blastermaster/phiTrim.py qnifft.electrostatics.phi box trim.electrostatics.phi
```
4. 产生范德华打分网格
---
#### 目的：产生范德华打分网格
#### 输入文件：/{home}/working/INCHEM  (rec.crg.pdb  /{home}/bela30c/proteins/defaults/prot.table.ambcrg.ambH  /{home}/bela30c/proteins/defaults/vdm.parms.amb.mindock, box)
#### 输出文件：vdm.bmp  vdm.esp  vdm.log  vdm.vdw
---
```
echo "rec.crg.pdb" >> INCHEM
echo "prot.table.ambcrg.ambH" >> INCHEM
echo "vdw.parms.amb.mindock" >> INCHEM
echo "box" >> INCHEM
echo "0.2" >> INCHEM
echo "1" >> INCHEM
echo "4" >> INCHEM
echo "10" >> INCHEM
echo "2.3 2.6" >> INCHEM
echo "vdw" >> INCHEM
cp $DOCKBASE/proteins/defaults/prot.table.ambcrg.ambH .
cp $DOCKBASE/proteins/defaults/vdw.parms.amb.mindock .

$DOCKBASE/proteins/chemgrid/bin/chemgrid >& vdw.log
```
5. 产生非溶剂化配体打分网格

solvemap需要运行两次，第一次对于配体重原子，第二次对于配体氢原子

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

echo "rec.crg.lds.pdb" >> INSEV

echo "ligand.desolv.heavy" >> INSEV

echo "1.60,1.65,1.90,1.90,1.90,1.00" >> INSEV

echo "1.4" >> INSEV

echo "2" >> INSEV

echo "box" >> INSEV

echo "1.8" >> INSEV

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

cd heavy

$DOCKBASE/proteins/solvmap/bin/solvmap >& solvmap.log

cd ..

cd hydrogen

$DOCKBASE/proteins/solvmap/bin/solvmap >& solvmap.log

cd ..
```
6. 产生INDOCK文件
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
收集对接所需文件：
```
mkdir dockfiles

cp working/matching_spheres.sph dockfiles/

cp working/vdw.parms.amb.mindock dockfiles/

cp working/heavy/ligand.desolv.heavy dockfiles/

cp working/hydrogen/ligand.desolv.hydrogen dockfiles/

cp working/trim.electrostatics.phi dockfiles/

cp working/vdw.vdw dockfiles/

cp working/vdw.bmp dockfiles/
```
### 6.安装DB2_converter
#### 1.准备环境
```
conda create -n test python=3.10

conda activate test

conda install -c conda-forge ambertools openbabel rdkit -y

conda deactivate
```
下载db2_converter包：

将其解压，放入主目录中：/home/usr，执行以下命令：
```
cd db2_converter

conda env create -f db2_converter.yml

conda activate db2_converter

pip install -e .

conda deactivate
```
对安装情况进行测试：
```
conda activate test

which antechamber

which obabel

conda list rdkit
# 在db2_converter里有一个setup.py文件，在与其相同的路径下，执行pip命令：
pip install -e ~/db2_converter

build-ligand -h
```
※如果产生错误：
>ModuleNotFoundError: No module named 'tabulate'

使用pip安装即可：
```
pip install tabulate
```

>ModuleNotFoundError: No module named 'posebusters'

使用pip安装即可
```
pip install posebusters
```
之后重新运行`build-ligand -h`
#### 2.安装unicon或ZBH2022
From: https://www.zbh.uni-hamburg.de/forschung/amd/software/conformator.html

安装并激活后，执行如下命令：
```
cd ~/test/working

~/unicon_version/unicon -i xtal-lig.pdb -o xtal-lig.smi -t single -p single

max_conf=100

workingpath=$(readlink -f "xtal-lig")

outputpath=$(readlink -f "xtal-lig")

method="conformator"

insmi="xtal-lig.smi"

build_ligand -i $insmi -n $max_conf --keep_max_conf --workingpath $workingpath --outputpath $workingpath --method $method -c --reseth --rotateh
```
※如果运行的时候出现了这样的错误：
> conformator | CONF_EXE    | PATH           | False
则说明这个组件的路径位置错误，需要进行修改。
```
cd ~/db2_converter/db2_converter/
```
然后打开config.py文件，找到对应组建的路径位置，将其修改为正确路径。
#### 3.产生db2文件
可以在 https://tldr.docking.org/results/all 上，使用smi文件获取不同构象状态的组合文件（db2.gz）。
下载db2_utils组件。
```
#当前位置：~/test/working
cd xtal-lig

~/db2utils/db2gz_to_mol2.sh file_name.db2.gz

mkdir docking

readlink -f xtal-lig/30Z_B_702.db2.gz > split_database_index

cp ../INDOCK .

$DOCKBASE/docking/DOCK/src/x86_64/dock64 INDOCK
```
---
Finish
