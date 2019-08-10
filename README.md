# QD 量子动力学计算软件
## 最新版本 v1.2
## 版本更新记录
v1.2
* 更新了UI，更正了透射系数的称呼
* 取消了输入反应路径简并数，改为输入转动对称数

v1.1
* 解决了makefile与文件夹重名的问题
* 修正了对转动对称数的描述


v1.0
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
## 二、安装使用方法
iOS尚不支持
### 1. Windows
#### 方法1：下载二进制文件
下载Windows的可执行文件```QD_1_1_win64.exe```，直接运行。
#### 方法2：使用GCC编译器编译源码
下载[mingw64-gcc](https://sourceforge.net/projects/mingw-w64/files/latest/download)。

安装gcc 4.6以上版本。

将```mingw64/bin```目录加入PATH：[方法教程](https://blog.csdn.net/Flood_Dragon/article/details/12363705)

下载```source```文件夹中的源码，将```make_file```文件夹中```makefile.win64.mk```复制到```source```文件夹，并修改文件名为```makefile```。

在当前目录下启动命令提示符或PowerShell，运行```mingw32-make```，编译完成将在当前目录下生成```QD.exe```二进制文件。

运行```QD.exe```即可。
#### 方法3：使用Microsoft Visual Studio编译源码
下载[Microsoft Visual Studio 2019](https://visualstudio.microsoft.com/zh-hans/downloads/)。

下载```source```文件夹中的源码，通过文本编辑器（如：[Visual Studio Code](https://code.visualstudio.com/), [Vim](https://www.vim.org/), [Notepad++](https://notepad-plus-plus.org/)等）将源码的文字编码由UTF-8改为GBK。

使用Microsoft Visual Studio新建空项目，将源码的.cpp文件和.h文件分别导入空项目，然后编译整个项目，编译完成将生成二进制文件。
### 2. Linux
#### 方法1：下载二进制文件
下载Linux的二进制文件```QD_1_1_linux```

使用```chmod +x QD_1_1_linux```赋予可执行权限

使用```./QD_1_1_linux```运行软件
#### 方法2：使用GCC编译器编译源码
要求安装有GCC 4.6以上版本，可以在终端中通过```gcc -v```来查看系统自带GCC编译器版本。

下载```source```文件夹中的源码，将```make_file```文件夹中```makefile.linux.mk```复制到```source```文件夹，并修改文件名为```makefile```。

在当前目录下启动终端，运行```make```，编译完成将在当前目录下生成```QD```二进制文件。

运行```./QD```即可。
### 3. MacOS
#### 方法：使用Clang编译器编译源码
要求安装有Clang 3.1以上版本，可以通过```clang --version```来查看系统自带Clang编译器版本。

下载```source```文件夹中的源码，将```make_file```文件夹中```makefile.macos.mk```复制到```source```文件夹，并修改文件名为```makefile```。

在当前目录下启动终端，运行```make```，编译完成将在当前目录下生成```QD```二进制文件。

运行```./QD```即可。
### 4. Android
#### 方法：借助C4droid使用源码
下载安装```C4droid 5.96```及以上版本，```GCC for C4droid 6.1.0```及以上版本，```SDL plugin for C4droid 2.0.4```及以上版本。

修改C4droid默认编译器为g++编译器，修改g++编译器参数```-std=c++11```或更高C++语言标准。

下载```QD_1_1_android.cpp```文件，用C4droid打开，编译运行。
## 三、使用手册
<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>
### 1. 量子隧穿效应透过系数计算
#### 计算原理
##### Wigner方法计算量子隧穿效应透射系数
$$\kapa=1+{(h/k_B/T*\nu*30000000000.0)}^2/24.0$$

#### 输入参数
#### 计算实例
### 2. 单分子反应速率常数计算
#### 计算原理
#### 输入参数
#### 计算实例
### 3. 双分子反应速率常数计算
#### 计算原理
#### 输入参数
#### 计算实例
## 四、量子动力学计算库
QD使用了量子动力学计算库，计算库内提供了量子隧穿效应、过渡态理论的多个函数接口。

### 最新版本v1.1
### 获取和使用量子动力学计算库
下载源码中的```QDynamic.h```和```QDynamic.cpp```文件，使用时将其加入头文件中：
```c++
#include"QDynamic.h"
```

具体函数请参看```QDynamic.h```文件
