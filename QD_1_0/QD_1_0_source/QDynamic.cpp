/******************************************************
*                                                     *
*   基于量子过渡态理论的动力学计算库v1.0                 *
*   Author Qiaosong Lin, Wuhan University, 2019       *
*                                                     *
*******************************************************/

/******************************************************
*                                                     *
*   最新修改记录:                                      *
*   2019-08-09    修改了部分bug                        *
*                                                     *
*******************************************************/

/******************************************************
*                                                     *
*   版本更新记录:                                      *
*   2019-08-08    v1.0                                *
*   >支持基元反应量子隧穿效应透过系数的计算              *
*       >基于Wigner方法的计算                          *
*       >基于近似的Skodje-Truhlar方法的计算            *
*       >基于完整的Skodje-Truhlar方法的计算            *
*   >支持单分子反应的量子过渡态理论动力学计算            *
*       >基于自由能的计算                              *
*       >基于配分函数的计算                            *
*   >支持双分子反应的量子过渡态理论动力学计算            *
*       >基于自由能的计算                              *
*       >基于配分函数的计算                            *
*                                                     *
*******************************************************/




#include<iostream>
#include<cmath>
#include"QDynamic.h"

//不含势垒和V的构造函数
QTunnel::QTunnel(double temperature,double virtualFreq)
:Temperature(temperature),VirtualFreq(virtualFreq)  //初始化列表
{
    //函数体
    Barrier=-1.0;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算量子隧穿效应透过系数输入信息*"<<std::endl;
    std::cout<<"温度="<<Temperature<<" K"<<std::endl;
    std::cout<<"虚频="<<-VirtualFreq<<" cm-1"<<std::endl;
    std::cout<<"隧穿势垒=nan"<<std::endl;
    std::cout<<"放热方向校正系数=nan"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
    Alpha=Calculate_Alpha();
    Beta=Calculate_Beta();
}


//含势垒和V的构造函数
QTunnel::QTunnel(double temperature,double virtualFreq,double barrier,double v)
:Temperature(temperature),VirtualFreq(virtualFreq),Barrier(barrier),V(v) //初始化列表
{
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算量子隧穿效应透过系数输入信息*"<<std::endl;
    std::cout<<"温度="<<Temperature<<" K"<<std::endl;
    std::cout<<"虚频="<<-VirtualFreq<<" cm-1"<<std::endl;
    std::cout<<"隧穿势垒="<<Barrier<<" kJ/mol"<<std::endl;
    std::cout<<"放热方向校正系数="<<V<<" kJ/mol"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
    Alpha=Calculate_Alpha();
    Beta=Calculate_Beta();
}


//计算量子隧穿效应透过系数κ函数
double QTunnel::Calculate_Kapa(const bool option) const
{
    switch(option)
    {
        case false: return Wigner();
        case true: return Skodje_Truhlar();
    }
}


//Wigner方法函数
double QTunnel::Wigner() const
{
    return 1+pow((h/kB/Temperature*VirtualFreq*30000000000.0),2)/24.0;
}



//Skodje-Truhlar方法函数
double QTunnel::Skodje_Truhlar() const
{
    double approximate=Beta*pi/Alpha/sin(Beta*pi/Alpha);
    if(Barrier==-1.0) //近似的Skodje-Truhlar方法
    {
        if(Beta>Alpha)
        {
            std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
            std::cout<<"!!"<<std::endl;
            std::cout<<"!!\t未提供势垒和放热方向校正系数的数值，通过近似的Skodje-Truhlar方法计算!!"<<std::endl;
            std::cout<<"!!"<<std::endl;
            std::cout<<"!!\tAlpha<Beta，可能造成重大误差\t\t\t\t!!"<<std::endl;
            std::cout<<"!!"<<std::endl;
            std::cout<<"!!\t建议补充势垒和放热方向校正系数的数值，进行完整的Skodje-Truhlar方法计算!!"<<std::endl;
            std::cout<<"!!"<<std::endl;
            std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
        }
        return approximate; 
    }
    else //完整的Skodje-Truhlar方法
    {    
        if(Alpha>Beta)
            return approximate-Beta/(Alpha-Beta)*exp((Beta-Alpha)*(Barrier-V)*1000.0/NA);
        else
            return Beta/(Beta-Alpha)*(exp((Beta-Alpha)*(Barrier-V)*1000.0/NA)-1.0);
    }
}



//含自由能不含配分函数的构造函数
Dynamic::Dynamic(double temperature,double sigma,double kapa,double freeEnergy)
:Temperature(temperature),Sigma(sigma),Kapa(kapa),FreeEnergy(freeEnergy) {}


//含配分函数不含自由能的构造函数
Dynamic::Dynamic(double temperature,double sigma,double kapa,double barrier,double partitionFunction_TS)
:Temperature(temperature),Sigma(sigma),Kapa(kapa),Barrier(barrier),PartitionFunction_TS(partitionFunction_TS) {}


//含配分函数和自由能的构造函数
Dynamic::Dynamic(double temperature,double sigma,double kapa,double freeEnergy,double barrier,double partitionFunction_TS)
:Temperature(temperature),Sigma(sigma),Kapa(kapa),FreeEnergy(freeEnergy),Barrier(barrier),PartitionFunction_TS(partitionFunction_TS) {}


//含自由能不含配分函数的构造函数
Dynamic_Single::Dynamic_Single(double temperature,double sigma,double kapa,double freeEnergy)
:Dynamic(temperature,sigma,kapa,freeEnergy) 
{
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算单分子反应速率常数输入信息*"<<std::endl;
    std::cout<<"温度="<<Temperature<<" K"<<std::endl;
    std::cout<<"反应路径简并数="<<Sigma<<std::endl;
    std::cout<<"量子隧穿效应透过系数="<<Kapa<<std::endl;
    std::cout<<"生成过渡态的Gibbis自由能(含零点能ZPE)="<<FreeEnergy<<" kJ/mol"<<std::endl;
    std::cout<<"生成过渡态的能垒(含零点能ZPE)=nan"<<std::endl;
    std::cout<<"过渡态的配分函数(不含零点能V=0)=nan"<<std::endl;
    std::cout<<"反应物A的配分函数(不含零点能V=0)=nan"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}



//含配分函数不含自由能的构造函数
Dynamic_Single::Dynamic_Single(double temperature,double sigma,double kapa,double barrier,double partitionFunction_TS,double partitionFunction_A)
:Dynamic(temperature,sigma,kapa,barrier,partitionFunction_TS),PartitionFunction_A(partitionFunction_A)
{
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算单分子反应速率常数输入信息*"<<std::endl;
    std::cout<<"温度="<<Temperature<<" K"<<std::endl;
    std::cout<<"反应路径简并数="<<Sigma<<std::endl;
    std::cout<<"量子隧穿效应透过系数="<<Kapa<<std::endl;
    std::cout<<"生成过渡态的Gibbis自由能(含零点能ZPE)=nan"<<std::endl;
    std::cout<<"生成过渡态的能垒(含零点能ZPE)="<<Barrier<<" kJ/mol"<<std::endl;
    std::cout<<"过渡态的配分函数(不含零点能V=0)="<<PartitionFunction_TS<<" mol-1"<<std::endl;
    std::cout<<"反应物A的配分函数(不含零点能V=0)="<<PartitionFunction_A<<" mol-1"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}



//含配分函数和自由能的构造函数
Dynamic_Single::Dynamic_Single(double temperature,double sigma,double kapa,double freeEnergy,double barrier,double partitionFunction_TS,double partitionFunction_A)
:Dynamic(temperature,sigma,kapa,freeEnergy,barrier,partitionFunction_TS),PartitionFunction_A(partitionFunction_A)
{
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算单分子反应速率常数输入信息*"<<std::endl;
    std::cout<<"温度="<<Temperature<<" K"<<std::endl;
    std::cout<<"反应路径简并数="<<Sigma<<std::endl;
    std::cout<<"量子隧穿效应透过系数="<<Kapa<<std::endl;
    std::cout<<"生成过渡态的Gibbis自由能(含零点能ZPE)="<<FreeEnergy<<" kJ/mol"<<std::endl;
    std::cout<<"生成过渡态的能垒(含零点能ZPE)="<<Barrier<<" kJ/mol"<<std::endl;
    std::cout<<"过渡态的配分函数(不含零点能V=0)="<<PartitionFunction_TS<<" mol-1"<<std::endl;
    std::cout<<"反应物A的配分函数(不含零点能V=0)="<<PartitionFunction_A<<" mol-1"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}


//计算单分子反应半衰期函数
double Dynamic_Single::Calculate_Halftime(const bool option) const {return log(2)/Calculate_k(option);}



//含自由能不含配分函数的构造函数
Dynamic_Double::Dynamic_Double(double pressure,double temperature,double sigma,double kapa,double freeEnergy)
:Dynamic(temperature,sigma,kapa,freeEnergy),Pressure(pressure)
{
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算双分子反应速率常数输入信息*"<<std::endl;
    std::cout<<"压强="<<Pressure<<" bar"<<std::endl;
    std::cout<<"温度="<<Temperature<<" K"<<std::endl;
    std::cout<<"反应路径简并数="<<Sigma<<std::endl;
    std::cout<<"量子隧穿效应透过系数="<<Kapa<<std::endl;
    std::cout<<"生成过渡态的Gibbis自由能(含零点能ZPE)="<<FreeEnergy<<" kJ/mol"<<std::endl;
    std::cout<<"生成过渡态的能垒(含零点能ZPE)=nan"<<std::endl;
    std::cout<<"过渡态的配分函数(不含零点能V=0)=nan"<<std::endl;
    std::cout<<"反应物A的配分函数(不含零点能V=0)=nan"<<std::endl;
    std::cout<<"反应物B的配分函数(不含零点能V=0)=nan"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}


//含配分函数不含自由能的构造函数
Dynamic_Double::Dynamic_Double(double pressure,double temperature,double sigma,double kapa,double barrier,double partitionFunction_TS,double partitionFunction_A,double partitionFunction_B)
:Dynamic(temperature,sigma,kapa,barrier,partitionFunction_TS),Pressure(pressure),PartitionFunction_A(partitionFunction_A),PartitionFunction_B(partitionFunction_B)
{
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算双分子反应速率常数输入信息*"<<std::endl;
    std::cout<<"压强="<<Pressure<<" bar"<<std::endl;
    std::cout<<"温度="<<Temperature<<" K"<<std::endl;
    std::cout<<"反应路径简并数="<<Sigma<<std::endl;
    std::cout<<"量子隧穿效应透过系数="<<Kapa<<std::endl;
    std::cout<<"生成过渡态的Gibbis自由能(含零点能ZPE)=nan"<<std::endl;
    std::cout<<"生成过渡态的能垒(含零点能ZPE)="<<Barrier<<" kJ/mol"<<std::endl;
    std::cout<<"过渡态的配分函数(不含零点能V=0)="<<PartitionFunction_TS<<" mol-1"<<std::endl;
    std::cout<<"反应物A的配分函数(不含零点能V=0)="<<PartitionFunction_A<<" mol-1"<<std::endl;
    std::cout<<"反应物B的配分函数(不含零点能V=0)="<<PartitionFunction_B<<" mol-1"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}



//含自由能和配分函数的构造函数
Dynamic_Double::Dynamic_Double(double pressure,double temperature,double sigma,double kapa,double freeEnergy,double barrier,double partitionFunction_TS,double partitionFunction_A,double partitionFunction_B)
:Dynamic(temperature,sigma,kapa,freeEnergy,barrier,partitionFunction_TS),Pressure(pressure),PartitionFunction_A(partitionFunction_A),PartitionFunction_B(partitionFunction_B)
{
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算双分子反应速率常数输入信息*"<<std::endl;
    std::cout<<"压强="<<Pressure<<" bar"<<std::endl;
    std::cout<<"温度="<<Temperature<<" K"<<std::endl;
    std::cout<<"反应路径简并数="<<Sigma<<std::endl;
    std::cout<<"量子隧穿效应透过系数="<<Kapa<<std::endl;
    std::cout<<"生成过渡态的Gibbis自由能(含零点能ZPE)="<<FreeEnergy<<" kJ/mol"<<std::endl;
    std::cout<<"生成过渡态的能垒(含零点能ZPE)="<<Barrier<<" kJ/mol"<<std::endl;
    std::cout<<"过渡态的配分函数(不含零点能V=0)="<<PartitionFunction_TS<<" mol-1"<<std::endl;
    std::cout<<"反应物A的配分函数(不含零点能V=0)="<<PartitionFunction_A<<" mol-1"<<std::endl;
    std::cout<<"反应物B的配分函数(不含零点能V=0)="<<PartitionFunction_B<<" mol-1"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}
