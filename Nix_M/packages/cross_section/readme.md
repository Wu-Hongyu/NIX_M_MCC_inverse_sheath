# 截面选取
结论参见 nix项目工作文档.docx，过程参见.\others\choose_Xsec_e_n_of_hydrogen_plasma（@pengchen）
# 文件结构
## 功能文件
### 获取截面
get_Xsec_e_H.m
get_Xsec_e_H2.m
test_get_Xsec.m
### 导入LXCat数据
ImportLXCat.m：从METHES移植到Nix_M的ImportLXCat类
test_importLXCat.m：对Nix_M-ImportLXCat的单元测试
### METHES-ImportLXCat文件
Ref：2016Rabie - METHES: A Monte Carlo collision code for the simulation of electron transport in low temperature plasmas  
  
License of METHES-ImportLXCat.txt：METHES-ImportLXCat的License  
e_Ar_METHES_test\：METHES-ImportLXCat中数据，test_importLXCat.m使用  
e_N2_METHES_test\：METHES-ImportLXCat中数据，test_importLXCat.m使用  
## 数据文件
文件夹命名方式：反应物A_反应物B_数据源_创建日期  
各子文件夹内，（唯一的）txt文件为截面数据文件，readme.md为详细信息，图片为截面绘图  
### e_H_CCC_20201214\
主要使用的e-H截面
### e_H2_Phelps_20201214\
主要使用的e-H2截面