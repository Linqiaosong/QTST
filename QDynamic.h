/******************************************************
*                                                     *
*   基于量子过渡态理论的动力学计算库v1.3                 *
*   Author Qiaosong Lin, Wuhan University, 2019       *
*                                                     *
*******************************************************/

/******************************************************
*                                                     *
*   最新修改记录:                                      *
*   2019-08-11    更新了构造函数提示的输入信息           *
*                                                     *
*******************************************************/

/******************************************************
*                                                     *
*   版本更新记录:                                      *
*   2019-08-08    v1.2                                *
*   >更新了构造函数提示的输入信息                        *
*   2019-08-08    v1.1                                *
*   >更新了构造函数提示的输入信息                        *
*   2019-08-08    v1.0                                *
*   >支持基元反应量子隧穿效应透射系数的计算              *
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


#ifndef QDYNAMIC_H
#define QDYNAMIC_H

#include<cmath>




//***********************************************基本常数***********************************************
constexpr double kB=1.3806503E-23;  // 玻尔兹曼常数，单位J/K
constexpr double h=6.6260696E-34;   // 普朗克常数，单位J*s
constexpr double NA=6.02214179E+23; // 阿伏伽德罗常数，单位mol-1
constexpr double R=kB*NA;           // 普适气体常数，单位 J/(mol*K)
constexpr double pi=3.14159265358979;   // 圆周率
constexpr double c=299792458;       // 光速，单位m/s
//*****************************************************************************************************





//*********************************量子隧道效应类：量子隧穿效应透射系数计算*********************************
class QTunnel
{
public:



    /************************************************************/
    //                         构造函数
    //
    // 不含势垒和V的构造函数
    // 初始化类对象，并在屏幕上显示输入信息
    QTunnel(double temperature,double virtualFreq);
    // 含势垒和V的构造函数
    // 初始化类对象，并在屏幕上显示输入信息
    QTunnel(double temperature,double virtualFreq,double barrier,double v);
    /************************************************************/

    // 使用默认的拷贝构造函数


    /************************************************************/
    //                   计算参数α函数（内联）
    //
    // 常成员函数，不修改类对象
    // 返回值为参数α值（双精度浮点数）
    inline double Calculate_Alpha() const {return 2*pi/h/(VirtualFreq*30000000000);}
    /************************************************************/




    /************************************************************/
    //                    计算参数β函数（内连）
    //
    // 常成员函数，不修改类对象
    // 返回值为参数β值（双精度浮点数）
    inline double Calculate_Beta() const {return 1/(kB*Temperature);}
    /************************************************************/



    /************************************************************/
    //                 计算量子隧穿效应透射系数κ函数
    //
    // 常成员函数，不修改类对象
    // 返回值为透射系数κ值（双精度浮点数）
    // 参数：bool型 option
    // option=false 使用Wigner方法计算
    // option=true  使用Skodje-Truhlar方法计算
    double Calculate_Kapa(const bool option) const;
    /************************************************************/




    /************************************************************/
    //                      Wigner方法函数
    //
    // 常成员函数，不修改类对象
    // 返回值为Wigner方法计算的透射系数κ值（双精度浮点数）
    double Wigner() const;
    /************************************************************/




    /************************************************************/
    //                    Skodje-Truhlar方法函数
    //
    // 常成员函数，不修改类对象
    // 调用不含势垒和V的构造函数时使用近似的Skodje-Truhlar方法计算
    // 调用含势垒和V的构造函数时使用完整的Skodje-Truhlar方法计算
    // 返回值为Skodje-Truhlar方法计算的透射系数κ值（双精度浮点数）
    double Skodje_Truhlar() const;
    /************************************************************/


    // 使用默认的析构函数



private:
    double Alpha;   // 量子隧穿效应参数α
    double Beta;    // 量子隧穿效应参数β
    double Temperature; // 温度
    double VirtualFreq; // 过渡态不对称伸缩振动虚频
    double Barrier; // 生成过渡态的能垒
    double V;   // 放热反应时数值为0；吸热反应时为产物U(T=0)减反应物U(T=0)
};
//***********************************************************************************************



//*********************************动力学基类：动力学的基本参数计算*********************************
class Dynamic
{
public:

    /************************************************************/
    //                       构造函数
    //
    // 含自由能不含配分函数的构造函数
    // 初始化类对象，并在屏幕上显示输入信息
    Dynamic(double temperature,double sigma,double kapa,double freeEnergy);
    // 含配分函数不含自由能的构造函数
    // 初始化类对象，并在屏幕上显示输入信息
    Dynamic(double temperature,double sigma,double kapa,double barrier,double partitionFunction_TS);
    // 含配分函数和自由能的构造函数
    // 初始化类对象，并在屏幕上显示输入信息
    Dynamic(double temperature,double sigma,double kapa,double freeEnergy,double barrier,double partitionFunction_TS);
    /************************************************************/


    // 使用默认的拷贝构造函数


    /************************************************************/
    //                     计算速率常数函数
    //
    // 常成员函数，不修改类对象
    // 纯虚函数，用来实现多态，需要在派生类重重构
    // 需要在调用派生类的版本
    // 基类版本不可调用
    virtual double Calculate_k(const bool option) const=0;
    /************************************************************/



    /************************************************************/
    //                计算量子隧穿贡献率函数（内联）
    //
    // 常成员函数，不修改类对象
    // 返回值为量子隧穿贡献率（双精度浮点数）
    inline double ContributionRate_Kapa() const {return (Kapa-1)/Kapa;}
    /************************************************************/


    // 使用默认的析构函数


protected:
    double Temperature; // 温度
    double Sigma;   // 对称因子
    double Kapa;    // 量子隧穿透射系数
    double FreeEnergy;  // 生成过渡态的Gibbs自由能
    double Barrier;     // 生成过渡态的能垒
    double PartitionFunction_TS;    // 过渡态的配分函数
};
//******************************************************************************************************





//******************************单分子反应类：继承基类，增加了反应物A的配分函数*****************************
class Dynamic_Single: public Dynamic
{
public:


    /************************************************************/
    //                     构造函数
    //
    // 含自由能不含配分函数的构造函数
    // 初始化类对象，并调用基类构造函数，并在屏幕上显示输入信息
    Dynamic_Single(double temperature,double sigma,double kapa,double freeEnergy);
    // 含配分函数不含自由能的构造函数
    // 初始化类对象，并调用基类构造函数，并在屏幕上显示输入信息
    Dynamic_Single(double temperature,double sigma,double kapa,double barrier,double partitionFunction_TS,double partitionFunction_A);
    // 含自由能和配分函数的构造函数
    // 初始化类对象，并调用基类构造函数，并在屏幕上显示输入信息
    Dynamic_Single(double temperature,double sigma,double kapa,double freeEnergy,double barrier,double partitionFunction_TS,double partitionFunction_A);
    /************************************************************/
    

    // 使用默认的拷贝构造函数



    /************************************************************/
    //             计算单分子反应速率常数（虚函数）
    //
    // 重构基类纯虚函数
    // 常成员函数，不修改类对象
    // 返回值为单分子反应的速率常数（双精度浮点数）
    // 参数：bool型 option
    // option=true  使用自由能方法计算
    // option=false 只用配分函数计算
    virtual double Calculate_k(const bool option) const
    {
        switch(option)
        {
            case false: return Kapa*Sigma*kB*Temperature/h*exp(-FreeEnergy*1000.0/(R*Temperature));
            case true: return Kapa*Sigma*kB*Temperature/h*PartitionFunction_TS/PartitionFunction_A*exp(-Barrier*1000.0/(R*Temperature));
        }
    }
    /************************************************************/




    /************************************************************/
    //               计算单分子反应半衰期函数
    //
    // 常成员函数，不修改类对象
    // 返回值为反应物A的半衰期（双精度浮点数）
    // 参数：bool型 option
    // option=true  使用自由能方法计算
    // option=false 只用配分函数计算
    double Calculate_Halftime(const bool option) const;
    /************************************************************/


    // 使用默认的析构函数



private:
    double PartitionFunction_A;     // 单分子反应反应物A的配分函数
};
//********************************************************************************************




//**********双分子反应类：继承基类，增加了双分子反应的压力和反应物A、B的配分函数********************
class Dynamic_Double: public Dynamic
{
public:


    /************************************************************/
    //                      构造函数
    //
    // 含自由能不含配分函数的构造函数
    // 初始化类对象，调用基类构造函数，并在屏幕上显示输入信息
    Dynamic_Double(double pressure,double temperature,double sigma,double kapa,double freeEnergy);
    // 含配分函数不含自由能的构造函数
    // 初始化类对象，调用基类构造函数，并在屏幕上显示输入信息
    Dynamic_Double(double pressure,double temperature,double sigma,double kapa,double barrier,double partitionFunction_TS,double partitionFunction_A,double partitionFunction_B);
    // 含自由能和配分函数的构造函数
    // 初始化类对象，调用基类构造函数，并在屏幕上显示输入信息
    Dynamic_Double(double pressure,double temperature,double sigma,double kapa,double freeEnergy,double barrier,double partitionFunction_TS,double partitionFunction_A,double partitionFunction_B);    
    /************************************************************/


    //使用默认拷贝构造函数
    

    /************************************************************/
    //              计算双分子反应速率常数函数（虚函数）
    // 
    // 重构基类纯虚函数
    // 常成员函数，不修改类对象
    // 返回值为双分子反应的速率常数（双精度浮点数）
    // 参数：bool型 option
    // option=true  使用自由能方法计算
    // option=false 只用配分函数计算
    virtual double Calculate_k(const bool option) const
    {
        switch(option)
        {
            case false: return Kapa*Sigma*kB*Temperature/h*(kB*Temperature/(Pressure*100000.0))*1000000.0*exp(-FreeEnergy*1000.0/(R*Temperature))*NA*0.001;
            case true: return Kapa*Sigma*kB*Temperature/h*(kB*Temperature/(Pressure*100000.0))*1000000.0*PartitionFunction_TS/PartitionFunction_A/PartitionFunction_B*exp(-Barrier*1000.0/(R*Temperature))*NA*0.001;
        }
    }
    /************************************************************/




    //使用默认析构函数


private:
    double Pressure;    // 双分子反应压强
    double PartitionFunction_A; // 双分子反应反应物A的配分函数
    double PartitionFunction_B; // 双分子反应反应物B的配分函数
};
//********************************************************************************************

#endif  