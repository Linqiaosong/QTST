# QD 量子动力学计算软件
## 最新版本 v1.0
## 版本更新记录
* 支持基元反应量子隧穿效应透过系数的计算              
  * 基于Wigner方法的计算                          
  * 基于近似的Skodje-Truhlar方法的计算            
  * 基于完整的Skodje-Truhlar方法的计算            
* 支持单分子反应的量子过渡态理论动力学计算            
  * 基于自由能的计算                              
  * 基于配分函数的计算                            
* 支持双分子反应的量子过渡态理论动力学计算            
  * 基于自由能的计算                              
  * 基于配分函数的计算                            
## 一、简介
QD是基于C++开发的量子过渡态理论动力学计算软件，使用C++11标准，支持量子隧穿效应透过系数计算，单分子反应动力学计算，双分子反应动力学计算。
## 二、安装方法
iOS尚不支持
### 1. Windows
#### 方法1：下载二进制文件
下载Windows的可执行文件```QD_1_0_win64.exe```，直接运行。
#### 方法2：使用GCC编译器编译源码
下载mingw64-gcc: [下载地址](https://sourceforge.net/projects/mingw-w64/files/latest/download)

安装gcc 4.6以上版本

将mingw64/bin目录加入PATH：[方法教程](https://blog.csdn.net/Flood_Dragon/article/details/12363705)

下载```source```文件夹中的源码，将```makefile```文件夹中```makefile.win64.mk```复制到source文件夹，并修改文件名为```makefile```。

在当前目录下启动命令提示符或PowerShell，运行```mingw32-make```，编译完成将在当前目录下生成```QD.exe```二进制文件。

运行```QD.exe```即可。
#### 方法3：使用Microsoft Visual Studio编译源码
下载Microsoft Visual Studio 2019: [下载地址](https://visualstudio.microsoft.com/zh-hans/downloads/)

下载```source```文件夹中的源码，通过文本编辑器（如：Visual Studio Code, Vim, Notepad++等）将源码的文字编码由UTF-8改为GBK。

使用Microsoft Visual Studio新建空项目，将源码的.cpp文件和.h文件分别导入空项目，然后编译整个项目，编译完成将生成```main.exe```二进制文件。

运行```main.exe```即可。
### 2. Linux
#### 方法1：下载二进制文件
下载Linux的二进制文件```QD_1_0_linux```

使用```chmod +x QD_1_0_linux```赋予可执行权限

使用```./QD_1_0_linux```运行软件
#### 方法2：使用GCC编译器编译源码
### 3. MacOS
#### 方法：使用Clang编译器编译源码
### 4. Android
#### 方法：借助C4droid使用源码
