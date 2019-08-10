# 1. QD 量子动力学计算软件
<img src="http://www.forkosh.com/mathtex.cgi? \Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}">

最新版本 v1.3
--

<!-- TOC -->

- [1. QD 量子动力学计算软件](#1-qd-量子动力学计算软件)
    - [1.1. 版本更新记录](#11-版本更新记录)
    - [1.2. 简介](#12-简介)
    - [1.3. 安装使用方法](#13-安装使用方法)
        - [1.3.1. Windows](#131-windows)
            - [1.3.1.1. 下载二进制文件](#1311-下载二进制文件)
            - [1.3.1.2. 使用GCC编译器编译源码](#1312-使用gcc编译器编译源码)
            - [1.3.1.3. 使用Microsoft Visual Studio编译源码](#1313-使用microsoft-visual-studio编译源码)
        - [1.3.2. Linux](#132-linux)
            - [1.3.2.1. 下载二进制文件](#1321-下载二进制文件)
            - [1.3.2.2. 使用GCC编译器编译源码](#1322-使用gcc编译器编译源码)
        - [1.3.3. MacOS](#133-macos)
            - [1.3.3.1. 使用Clang编译器编译源码](#1331-使用clang编译器编译源码)
        - [1.3.4. Android](#134-android)
            - [1.3.4.1. 借助C4droid使用源码](#1341-借助c4droid使用源码)
- [2. QDynamic 量子动力学计算库](#2-qdynamic-量子动力学计算库)
    - [2.1. 获取和使用QDynamic库](#21-获取和使用qdynamic库)

<!-- /TOC -->

## 1.1. 版本更新记录
v1.3
* 更新了UI，更正了一些错误显示
* 更新了文件目录

v1.2
* 更新了UI，更正了透射系数的称呼
* 取消了输入反应路径简并数，改为输入转动对称数

v1.1
* 解决了makefile与文件夹重名的问题
* 修正了对转动对称数的描述


v1.0
* 支持基元反应量子隧穿效应透射系数的计算              
  * 基于Wigner方法的计算                          
  * 基于近似的Skodje-Truhlar方法的计算            
  * 基于完整的Skodje-Truhlar方法的计算            
* 支持单分子反应的量子过渡态理论动力学计算            
  * 基于自由能的计算                              
  * 基于配分函数的计算                            
* 支持双分子反应的量子过渡态理论动力学计算            
  * 基于自由能的计算                              
  * 基于配分函数的计算                            
## 1.2. 简介
QD是基于C++开发的量子过渡态理论动力学计算软件，使用C++11标准，支持量子隧穿效应透射系数计算，单分子反应动力学计算，双分子反应动力学计算。
## 1.3. 安装使用方法
iOS尚不支持
### 1.3.1. Windows
#### 1.3.1.1. 下载二进制文件
下载Windows的可执行文件```QD_1_3_win_binary.zip```，解压后直接运行。
#### 1.3.1.2. 使用GCC编译器编译源码
下载[mingw64-gcc](https://sourceforge.net/projects/mingw-w64/files/latest/download)。

安装gcc 4.6以上版本。

将```mingw64/bin```目录加入PATH：[方法教程](https://blog.csdn.net/Flood_Dragon/article/details/12363705)

下载```QD_1_3_source```文件夹中的源码，将```make_file```文件夹中```makefile.win64.mk```复制到```source```文件夹，并修改文件名为```makefile```。

在当前目录下启动命令提示符或PowerShell，运行```mingw32-make```，编译完成将在当前目录下生成```QD.exe```二进制文件。

运行```QD.exe```即可。
#### 1.3.1.3. 使用Microsoft Visual Studio编译源码
下载[Microsoft Visual Studio 2019](https://visualstudio.microsoft.com/zh-hans/downloads/)。

下载```QD_1_3_source```文件夹中的源码，通过文本编辑器（如：[Visual Studio Code](https://code.visualstudio.com/), [Vim](https://www.vim.org/), [Notepad++](https://notepad-plus-plus.org/)等）将源码的文字编码由UTF-8改为GBK。

使用Microsoft Visual Studio新建空项目，将源码的.cpp文件和.h文件分别导入空项目，然后编译整个项目，编译完成将生成二进制文件。
### 1.3.2. Linux
#### 1.3.2.1. 下载二进制文件
下载Linux的二进制文件```QD_1_3_linux_binary.tar.gz```

使用```tar -zxvf QD_1_3_linux_binary.tar.gz```解压。

进入文件夹赋予二进制文件可执行权限
```bash
cd QD_1_3_linux_binary
chmod +x QD
```

使用```./QD```启动软件

也可以将```.../QD_1_3_linux_binary```目录加入PATH，直接在终端使用```QD```启动软件。
#### 1.3.2.2. 使用GCC编译器编译源码
要求安装有GCC 4.6以上版本，可以在终端中通过```gcc -v```来查看系统自带GCC编译器版本。

下载```QD_1_3_source```文件夹中的源码，将```make_file```文件夹中```makefile.linux.mk```复制到```source```文件夹，并修改文件名为```makefile```。

在当前目录下启动终端，运行```make```，编译完成将在当前目录下生成```QD```二进制文件。

运行```./QD```即可。
### 1.3.3. MacOS
#### 1.3.3.1. 使用Clang编译器编译源码
要求安装有Clang 3.1以上版本，可以通过```clang --version```来查看系统自带Clang编译器版本。

下载```QD_1_3_source```文件夹中的源码，将```make_file```文件夹中```makefile.macos.mk```复制到```source```文件夹，并修改文件名为```makefile```。

在当前目录下启动终端，运行```make```，编译完成将在当前目录下生成```QD```二进制文件。

运行```./QD```即可。
### 1.3.4. Android
#### 1.3.4.1. 借助C4droid使用源码
下载安装```C4droid 5.96```及以上版本，```GCC for C4droid 6.1.0```及以上版本，```SDL plugin for C4droid 2.0.4```及以上版本。

修改C4droid默认编译器为g++编译器，修改g++编译器参数```-std=c++11```或更高C++语言标准。

下载```QD_1_3_android.cpp```文件，用C4droid打开，编译运行。


# 2. QDynamic 量子动力学计算库
QD使用了量子动力学计算库，计算库内提供了量子隧穿效应、过渡态理论的多个函数接口。

2.1. 最新版本v1.2
--

## 2.1. 获取和使用QDynamic库
下载```QD_1_3_source```源码中的```QDynamic.h```和```QDynamic.cpp```文件，使用时将其加入头文件中：
```c++
#include"QDynamic.h"
```

具体函数请参看```QDynamic.h```文件。
