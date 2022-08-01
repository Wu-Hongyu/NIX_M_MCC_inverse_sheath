> 功能模块
# 简介
各packages包括功能实现代码和测试代码。  
测试代码不仅确保功能实现符合预期，并提供了packages使用的示例。
# 使用
请以Nix_M为MATLAB当前文件夹  
## main.m中使用
直接使用即可。详见main.m。
## 其他地方使用  
1. 设置文件路径  
在控制台执行以下指令：  
if exist('get_path_init','file')~=2
    addpath('./packages') % 添加./packages到路径
end
path_list=get_path_init('test000'); %文件路径
2. 运行test
运行test_all.m，或者在控制台执行以下指令（替换package_name）：  
results = runtests('test_package_name');disp(table(results))
3. 调用函数  
4. 如需修改package代码，请在修改后再次测试。可能需要补充测试代码。
# 文件结构  
每个文件夹是一个package，一般不支持再下级文件夹
## cross_section
截面数据文件，与导入代码。