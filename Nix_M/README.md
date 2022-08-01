# Nix_M
 HUST-Nix-Matlab  
 Matlab实现的 **灵活可靠的** 1/2D-PIC/MCC 等离子体模拟代码
# 使用
请以Nix_M/为MATLAB当前文件夹
## 运行已有算例  
运行main.m或example_examplename.m
## 开发新算例  
请在main或其他已有算例基础上，新建文件  
## 修改packages，见./packages/README.md  
# 文件结构
三级文件夹。请不要改动现有文件夹名与文件名。
## packages
功能模块所在文件夹，各模块含函数与单元测试
## examples
使用案例
## others
其他文件所在文件夹，数据文件、日志文件、临时文件等

2022.07.29吴鸿宇保存
根据张哲师姐课题复现逆鞘层仿真
路径：……\Nix_M_MCC_inverse_sheath\Nix_M\examples\example_sheath_potential_and_sheath_length_Benchmark
基本模型：H+离子碰撞到边界上，根据入射速度确定H-离子产生概率；H原子与碰撞到边界上，根据密度、温度确定流量、碰撞数目和产生H-概率；H-超过边界后消除；引入MCC碰撞，产生冷H+离子，冷H+离子速度服从麦克斯韦分布
目前情况：增大碰撞概率后可以产生逆鞘层，否则不行。一定条件下有产生逆鞘层的趋势，但是无法形成稳定的逆鞘层
