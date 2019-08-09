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

//***********************************************基本常数***********************************************
constexpr double kB=1.3806503E-23;  // 玻尔兹曼常数，单位J/K
constexpr double h=6.6260696E-34;   // 普朗克常数，单位J*s
constexpr double NA=6.02214179E+23; // 阿伏伽德罗常数，单位mol-1
constexpr double R=kB*NA;           // 普适气体常数，单位 J/(mol*K)
constexpr double pi=3.14159265358979;   // 圆周率
constexpr double c=299792458;       // 光速，单位m/s
//*****************************************************************************************************





//*********************************量子隧道效应类：量子隧穿效应透过系数计算*********************************
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
    //                 计算量子隧穿效应透过系数κ函数
    //
    // 常成员函数，不修改类对象
    // 返回值为透过系数κ值（双精度浮点数）
    // 参数：bool型 option
    // option=false 使用Wigner方法计算
    // option=true  使用Skodje-Truhlar方法计算
    double Calculate_Kapa(const bool option) const;
    /************************************************************/




    /************************************************************/
    //                      Wigner方法函数
    //
    // 常成员函数，不修改类对象
    // 返回值为Wigner方法计算的透过系数κ值（双精度浮点数）
    double Wigner() const;
    /************************************************************/




    /************************************************************/
    //                    Skodje-Truhlar方法函数
    //
    // 常成员函数，不修改类对象
    // 调用不含势垒和V的构造函数时使用近似的Skodje-Truhlar方法计算
    // 调用含势垒和V的构造函数时使用完整的Skodje-Truhlar方法计算
    // 返回值为Skodje-Truhlar方法计算的透过系数κ值（双精度浮点数）
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
    double Kapa;    // 量子隧穿透过系数
    double FreeEnergy;  // 生成过渡态的Gibbis自由能
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























/************************************************************/
//
//                      单位换算函数（内联）
// 常表达式函数
// Hartree to kJ/mol
inline constexpr double EHartree2kJ(double hartree) {return hartree*2625.49962;}
// eV to kJ/mol
inline constexpr double EeV2kJ(double eV) {return eV*96.485;}
// kcal to kJ/mol
inline constexpr double Ekcal2kJ(double kcal) {return kcal*4.1840;}
// 摄氏度 to 开尔文K
inline constexpr double TC2K(double C) {return C+273.15;}
// 华氏度 to 开尔文K
inline constexpr double TF2K(double F) {return (F-32)/1.8+273.15;} 
// Hz to cm-1
inline constexpr double Fhz2cm(double hz) {return hz/(100*c);}
// Pa to bar
inline constexpr double Ppa2bar(double pa) {return pa/1.0e+5;}
// atm to bar
inline constexpr double Patm2bar(double atm) {return atm/101325;}
/************************************************************/

//wigner_kapa方法计算量子隧穿效应透过系数模块
void calculate_wigner_kapa(const QTunnel& QTn);


//skodje_truhlar方法计算量子隧穿效应透过系数模块
void calculate_skodje_truhlar_kapa(const QTunnel& QTn);


//全方法计算量子隧穿效应透过系数模块
void calculate_all_kapa(const QTunnel& QTn);


//自由能方法计算单分子反应速率常数模块
void calculate_single_freeEnergy_k(const Dynamic_Single& DSn);


//配分函数方法计算单分子反应速率常数模块
void calculate_single_Q_k(const Dynamic_Single& DSn);

//全方法计算单分子反应速率常数模块
void calculate_single_all_k(const Dynamic_Single& DSn);


//自由能方法计算双分子反应速率常数模块
void calculate_double_freeEnergy_k(const Dynamic_Double& DDn);


//配分函数方法计算双分子反应速率常数模块
void calculate_double_Q_k(const Dynamic_Double& DDn);


//全方法计算双分子反应速率常数模块
void calculate_double_all_k(const Dynamic_Double& DDn);

























int main()
{
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t  QD v1.0 量子动力学计算软件"<<std::endl;
    std::cout<<"\t\t Qiaosong Lin"<<std::endl;
    std::cout<<"\t    Wuhan University, 2019"<<std::endl;
    std::cout<<"帮助文档：https://github.com/Linqiaosong/QD"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*请选择功能*"<<std::endl;
    std::cout<<"[1] 计算量子隧穿效应透过系数"<<std::endl;
    std::cout<<"[2] 计算单分子反应速率常数"<<std::endl;
    std::cout<<"[3] 计算双分子反应速率常数"<<std::endl;
    int option;
    std::cin>>option;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
    switch(option)
    {
        case 1:
        {
            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请选择计算方法*"<<std::endl;
            std::cout<<"[1] Wigner方法"<<std::endl;
            std::cout<<"[2] Skodje-Truhlar方法"<<std::endl;
            std::cout<<"[3] 全部计算"<<std::endl;
            int option1;
            std::cin>>option1;
            std::cout<<"====================================================="<<std::endl;
            std::cout<<std::endl;
            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入信息*"<<std::endl;
            std::cout<<"科学计数法均按如下表示："<<std::endl;
            std::cout<<"\t2.5*10^-1表示为：2.5e-1"<<std::endl;


            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入温度*"<<std::endl;
            std::cout<<"温度数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
            std::cout<<"支持单位：\n[1] 热力学温度(K)\n[2] 摄氏度(C)\n[3] 华氏度(F)"<<std::endl;
            std::cout<<"例：300.0 K输入：300.0 1"<<std::endl;
            double T;
            int unitT;
            std::cin>>T>>unitT;
            switch(unitT)
            {
                case 1: break;
                case 2: T=TC2K(T);    break;
                case 3: T=TF2K(T);    break;
                default: return 1;
            }


            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入虚频*"<<std::endl;
            std::cout<<"虚频数值与单位之间用一个空格隔开，带负号，单位用标号表示"<<std::endl;
            std::cout<<"支持单位：\n[1] 光谱学单位——波数(cm-1)\n[2] 波的频率单位——赫兹(Hz)"<<std::endl;
            std::cout<<"例：-1000.0 cm-1输入：-1000.0 1"<<std::endl;
            double f;
            int unitf;
            std::cin>>f>>unitf;
            switch(unitf)
            {
                case 1: break;
                case 2: f=Fhz2cm(f); break;
                default: return 1;
            }
            std::cout<<"====================================================="<<std::endl;

            if(option1==1)
            {
                QTunnel QT(T,-f);
                calculate_wigner_kapa(QT);
            }
            else if(option1==2||option1==3)
            {
                std::cout<<"\t\t*请输入隧穿势垒（可选）*"<<std::endl;
                std::cout<<"\t*缺省隧穿势垒数据请输入-1 -1 *"<<std::endl;
                std::cout<<std::endl;
                std::cout<<"注意，缺省隧穿势垒将自动缺省放热方向校正系数，并进行近似的Skodje-Truhlar方法计算"<<std::endl;
                std::cout<<"若缺省隧穿势垒和放热方向校正系数的数值，Skodje-Truhlar方法将通过近似方法计算，在温度小于250K时可能造成重大误差"<<std::endl;
                std::cout<<std::endl;
                std::cout<<"隧穿势垒为过渡态在T=0K的内能和反应物在T=0K的内能之差"<<std::endl;
                std::cout<<"能量数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
                std::cout<<"支持单位：\n[1] kJ/mol\n[2] kcal/mol\n[3] eV/粒子\n[4] Hartree/粒子"<<std::endl;
                std::cout<<"例：30.0 kJ/mol输入：30.0 1"<<std::endl;
                double E;
                int unitE;
                std::cin>>E>>unitE;
                switch(unitE)
                {
                    case -1: break;
                    case 1: break;
                    case 2: E=Ekcal2kJ(E); break;
                    case 3: E=EeV2kJ(E); break;
                    case 4: E=EHartree2kJ(E); break;
                    default: return 1;
                }
                std::cout<<"====================================================="<<std::endl;
                if(E==-1.0&&unitE==-1)
                {
                    QTunnel QT(T,-f);
                    if(option1==2)
                        calculate_skodje_truhlar_kapa(QT);
                    else
                        calculate_all_kapa(QT);
                }
                else
                {
                    std::cout<<"\t\t*请输入放热方向校正系数*"<<std::endl;
                    std::cout<<"\t*正反应为放热反应时请输入0 1*"<<std::endl;
                    std::cout<<"\t*正反应为吸热反应时校正系数为产物在T=0K的内能和反应物在T=0K的内能之差*"<<std::endl;
                    std::cout<<"能量数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
                    std::cout<<"支持单位：\n[1] kJ/mol\n[2] kcal/mol\n[3] eV/粒子\n[4] Hartree/粒子"<<std::endl;
                    std::cout<<"例：30.0 kJ/mol输入：30.0 1"<<std::endl;
                    double V;
                    int unitV;
                    std::cin>>V>>unitV;
                    switch(unitV)
                    {
                        case 1: break;
                        case 2: E=Ekcal2kJ(E); break;
                        case 3: E=EeV2kJ(E); break;
                        case 4: E=EHartree2kJ(E); break;
                        default: return 1;
                    }
                    std::cout<<"====================================================="<<std::endl;


                    QTunnel QT(T,-f,E,V);
                    if(option1==2)
                        calculate_skodje_truhlar_kapa(QT);
                    else
                        calculate_all_kapa(QT);
                }
            }
            else return 1;
            break;
        }
        case 2:
        {
            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请选择计算方法*"<<std::endl;
            std::cout<<"[1] Gibbis自由能方法"<<std::endl;
            std::cout<<"[2] 配分函数方法"<<std::endl;
            std::cout<<"[3] 全部计算"<<std::endl;
            int option1;
            std::cin>>option1;
            std::cout<<"====================================================="<<std::endl;
            std::cout<<std::endl;
            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入信息*"<<std::endl;
            std::cout<<"科学计数法均按如下表示："<<std::endl;
            std::cout<<"\t2.5*10^-1表示为：2.5e-1"<<std::endl;


            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入温度*"<<std::endl;
            std::cout<<"温度数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
            std::cout<<"支持单位：\n[1] 热力学温度(K)\n[2] 摄氏度(C)\n[3] 华氏度(F)"<<std::endl;
            std::cout<<"例：300.0 K输入：300.0 1"<<std::endl;
            double T;
            int unitT;
            std::cin>>T>>unitT;
            switch(unitT)
            {
                case 1: break;
                case 2: T=TC2K(T);    break;
                case 3: T=TF2K(T);    break;
                default: return 1;
            }




            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入反应路径简并数*"<<std::endl;
            std::cout<<"\t*如果计算自由能或配分函数时已经进行了对称简并校正，则此处填1*"<<std::endl;
            std::cout<<"反应路径简并数=反应物转动对称数/过渡态转动对称数"<<std::endl;
            double sigma;
            std::cin>>sigma;


            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入量子隧穿效应透过系数*"<<std::endl;
            std::cout<<"\t\t*若忽略量子隧穿效应，请输入1*"<<std::endl;
            double kapa;
            std::cin>>kapa;


            if(option1==2||option1==3)
            {
                std::cout<<"====================================================="<<std::endl;
                std::cout<<"\t\t*请输入生成过渡态的能垒*"<<std::endl;
                std::cout<<"\t*能垒为反应温度下过渡态内能（含零点能ZPE）和反应物内能（含零点能ZPE）之差*"<<std::endl;
                std::cout<<std::endl;
                std::cout<<"能量数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
                std::cout<<"支持单位：\n[1] kJ/mol\n[2] kcal/mol\n[3] eV/粒子\n[4] Hartree/粒子"<<std::endl;
                std::cout<<"例：30.0 kJ/mol输入：30.0 1"<<std::endl;
                double E;
                int unitE;
                std::cin>>E>>unitE;
                switch(unitE)
                {
                    case 1: break;
                    case 2: E=Ekcal2kJ(E); break;
                    case 3: E=EeV2kJ(E); break;
                    case 4: E=EHartree2kJ(E); break;
                    default: return 1;
                }



                std::cout<<"====================================================="<<std::endl;
                std::cout<<"\t\t*请输入过渡态的配分函数*"<<std::endl;
                std::cout<<"\t配分函数以mol-1为单位，不含零点能ZPE成分，即Q(V=0)"<<std::endl;
                double Q_TS;
                std::cin>>Q_TS;


                std::cout<<"====================================================="<<std::endl;
                std::cout<<"\t\t*请输入反应物A的配分函数*"<<std::endl;
                std::cout<<"\t配分函数以mol-1为单位，不含零点能ZPE成分，即Q(V=0)"<<std::endl;
                double Q_A;
                std::cin>>Q_A;

                if(option1==2)
                {
                    Dynamic_Single DS(T,sigma,kapa,E,Q_TS,Q_A);
                    calculate_single_Q_k(DS);
                }
                else
                {
                    std::cout<<"====================================================="<<std::endl;
                    std::cout<<"\t\t*请输入生成过渡态的Gibbis自由能*"<<std::endl;
                    std::cout<<"\t*Gibbis自由能为标准压力、反应温度下过渡态Gibbis自由能（含零点能ZPE）和反应物Gibbis自由能（含零点能ZPE）之差*"<<std::endl;
                    std::cout<<std::endl;
                    std::cout<<"能量数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
                    std::cout<<"支持单位：\n[1] kJ/mol\n[2] kcal/mol\n[3] eV/粒子\n[4] Hartree/粒子"<<std::endl;
                    std::cout<<"例：30.0 kJ/mol输入：30.0 1"<<std::endl;
                    double G;
                    int unitG;
                    std::cin>>G>>unitG;
                    switch(unitG)
                    {
                        case 1: break;
                        case 2: G=Ekcal2kJ(G); break;
                        case 3: G=EeV2kJ(G); break;
                        case 4: G=EHartree2kJ(G); break;
                        default: return 1;
                    }

                    Dynamic_Single DS(T,sigma,kapa,G,E,Q_TS,Q_A);
                    calculate_single_all_k(DS);
                }
            }
            else if(option1==1)
            {
                std::cout<<"====================================================="<<std::endl;
                std::cout<<"\t\t*请输入生成过渡态的Gibbis自由能*"<<std::endl;
                std::cout<<"\t*Gibbis自由能为标准压力、反应温度下过渡态Gibbis自由能（含零点能ZPE）和反应物Gibbis自由能（含零点能ZPE）之差*"<<std::endl;
                std::cout<<std::endl;
                std::cout<<"能量数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
                std::cout<<"支持单位：\n[1] kJ/mol\n[2] kcal/mol\n[3] eV/粒子\n[4] Hartree/粒子"<<std::endl;
                std::cout<<"例：30.0 kJ/mol输入：30.0 1"<<std::endl;
                double G;
                int unitG;
                std::cin>>G>>unitG;
                switch(unitG)
                {
                    case 1: break;
                    case 2: G=Ekcal2kJ(G); break;
                    case 3: G=EeV2kJ(G); break;
                    case 4: G=EHartree2kJ(G); break;
                    default: return 1;
                }
                Dynamic_Single DS(T,sigma,kapa,G);
                calculate_single_freeEnergy_k(DS);
            }
            break;
        }
        case 3:
        {
            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请选择计算方法*"<<std::endl;
            std::cout<<"[1] Gibbis自由能方法"<<std::endl;
            std::cout<<"[2] 配分函数方法"<<std::endl;
            std::cout<<"[3] 全部计算"<<std::endl;
            int option1;
            std::cin>>option1;
            std::cout<<"====================================================="<<std::endl;
            std::cout<<std::endl;
            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入信息*"<<std::endl;
            std::cout<<"科学计数法均按如下表示："<<std::endl;
            std::cout<<"\t2.5*10^-1表示为：2.5e-1"<<std::endl;



            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入压强*"<<std::endl;
            std::cout<<"压强数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
            std::cout<<"支持单位：\n[1] 标准压强(bar)\n[2] 帕斯卡(Pa)\n[3] 大气压(atm)"<<std::endl;
            std::cout<<"例：101325 Pa输入：101325 2"<<std::endl;
            double p;
            int unitp;
            std::cin>>p>>unitp;
            switch(unitp)
            {
                case 1: break;
                case 2: p=Ppa2bar(p);    break;
                case 3: p=Patm2bar(p);    break;
                default: return 1;
            }


            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入温度*"<<std::endl;
            std::cout<<"温度数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
            std::cout<<"支持单位：\n[1] 热力学温度(K)\n[2] 摄氏度(C)\n[3] 华氏度(F)"<<std::endl;
            std::cout<<"例：300.0 K输入：300.0 1"<<std::endl;
            double T;
            int unitT;
            std::cin>>T>>unitT;
            switch(unitT)
            {
                case 1: break;
                case 2: T=TC2K(T);    break;
                case 3: T=TF2K(T);    break;
                default: return 1;
            }




            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入反应路径简并数*"<<std::endl;
            std::cout<<"\t*如果计算自由能或配分函数时已经进行了对称简并校正，则此处填1*"<<std::endl;
            std::cout<<"反应路径简并数=反应物转动对称数乘积/过渡态转动对称数"<<std::endl;
            double sigma;
            std::cin>>sigma;


            std::cout<<"====================================================="<<std::endl;
            std::cout<<"\t\t*请输入量子隧穿效应透过系数*"<<std::endl;
            std::cout<<"\t\t*若忽略量子隧穿效应，请输入1*"<<std::endl;
            double kapa;
            std::cin>>kapa;


            if(option1==2||option1==3)
            {
                std::cout<<"====================================================="<<std::endl;
                std::cout<<"\t\t*请输入生成过渡态的能垒*"<<std::endl;
                std::cout<<"\t*能垒为反应温度下过渡态内能（含零点能ZPE）和所有反应物内能（含零点能ZPE）之差*"<<std::endl;
                std::cout<<std::endl;
                std::cout<<"能量数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
                std::cout<<"支持单位：\n[1] kJ/mol\n[2] kcal/mol\n[3] eV/粒子\n[4] Hartree/粒子"<<std::endl;
                std::cout<<"例：30.0 kJ/mol输入：30.0 1"<<std::endl;
                double E;
                int unitE;
                std::cin>>E>>unitE;
                switch(unitE)
                {
                    case 1: break;
                    case 2: E=Ekcal2kJ(E); break;
                    case 3: E=EeV2kJ(E); break;
                    case 4: E=EHartree2kJ(E); break;
                    default: return 1;
                }



                std::cout<<"====================================================="<<std::endl;
                std::cout<<"\t\t*请输入过渡态的配分函数*"<<std::endl;
                std::cout<<"\t配分函数以mol-1为单位，不含零点能ZPE成分，即Q(V=0)"<<std::endl;
                double Q_TS;
                std::cin>>Q_TS;


                std::cout<<"====================================================="<<std::endl;
                std::cout<<"\t\t*请输入反应物A的配分函数*"<<std::endl;
                std::cout<<"\t配分函数以mol-1为单位，不含零点能ZPE成分，即Q(V=0)"<<std::endl;
                double Q_A;
                std::cin>>Q_A;


                std::cout<<"====================================================="<<std::endl;
                std::cout<<"\t\t*请输入反应物B的配分函数*"<<std::endl;
                std::cout<<"\t配分函数以mol-1为单位，不含零点能ZPE成分，即Q(V=0)"<<std::endl;
                double Q_B;
                std::cin>>Q_B;


                if(option1==2)
                {
                    Dynamic_Double DD(p,T,sigma,kapa,E,Q_TS,Q_A,Q_B);
                    calculate_double_Q_k(DD);
                }
                else
                {
                    std::cout<<"====================================================="<<std::endl;
                    std::cout<<"\t\t*请输入生成过渡态的Gibbis自由能*"<<std::endl;
                    std::cout<<"\t*Gibbis自由能为标准压力、反应温度下过渡态Gibbis自由能（含零点能ZPE）和反应物Gibbis自由能（含零点能ZPE）之差*"<<std::endl;
                    std::cout<<std::endl;
                    std::cout<<"能量数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
                    std::cout<<"支持单位：\n[1] kJ/mol\n[2] kcal/mol\n[3] eV/粒子\n[4] Hartree/粒子"<<std::endl;
                    std::cout<<"例：30.0 kJ/mol输入：30.0 1"<<std::endl;
                    double G;
                    int unitG;
                    std::cin>>G>>unitG;
                    switch(unitG)
                    {
                        case 1: break;
                        case 2: G=Ekcal2kJ(G); break;
                        case 3: G=EeV2kJ(G); break;
                        case 4: G=EHartree2kJ(G); break;
                        default: return 1;
                    }

                    Dynamic_Double DD(p,T,sigma,kapa,G,E,Q_TS,Q_A,Q_B);
                    calculate_double_all_k(DD);
                }
            }
            else if(option1==1)
            {
                std::cout<<"====================================================="<<std::endl;
                std::cout<<"\t\t*请输入生成过渡态的Gibbis自由能*"<<std::endl;
                std::cout<<"\t*Gibbis自由能为标准压力、反应温度下过渡态Gibbis自由能（含零点能ZPE）和反应物Gibbis自由能（含零点能ZPE）之差*"<<std::endl;
                std::cout<<std::endl;
                std::cout<<"能量数值与单位之间用一个空格隔开，单位用标号表示"<<std::endl;
                std::cout<<"支持单位：\n[1] kJ/mol\n[2] kcal/mol\n[3] eV/粒子\n[4] Hartree/粒子"<<std::endl;
                std::cout<<"例：30.0 kJ/mol输入：30.0 1"<<std::endl;
                double G;
                int unitG;
                std::cin>>G>>unitG;
                switch(unitG)
                {
                    case 1: break;
                    case 2: G=Ekcal2kJ(G); break;
                    case 3: G=EeV2kJ(G); break;
                    case 4: G=EHartree2kJ(G); break;
                    default: return 1;
                }

                Dynamic_Double DD(p,T,sigma,kapa,G);
                calculate_double_freeEnergy_k(DD);
            }
            break;
        }
        default: return 1;
    }
    std::cout<<"按任意键退出..."<<std::endl;
    getchar(); 
    getchar();
    return 0;
}




















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















//wigner_kapa方法计算量子隧穿效应透过系数模块
void calculate_wigner_kapa(const QTunnel& QTn)
{
    double kapa;
    kapa=QTn.Calculate_Kapa(false);
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算量子隧穿效应透过系数结果*"<<std::endl;
    std::cout<<"\t\t*透过系数通过Wigner方法计算*"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"注意：Wigner方法计算精度较低，在非高温情况下通常会造成重大误差"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"透过系数Kapa="<<kapa<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}



//skodje_truhlar方法计算量子隧穿效应透过系数模块
void calculate_skodje_truhlar_kapa(const QTunnel& QTn)
{
    double kapa,alpha,beta;
    alpha=QTn.Calculate_Alpha();
    beta=QTn.Calculate_Beta();
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算量子隧穿效应透过系数结果*"<<std::endl;
    std::cout<<"\t\t*透过系数通过Skodje-Truhlar方法计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"Alpha="<<alpha<<std::endl;
    std::cout<<"Beta="<<beta<<std::endl;
    kapa=QTn.Calculate_Kapa(true);
    std::cout<<"透过系数Kapa="<<kapa<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}




//全方法计算量子隧穿效应透过系数模块
void calculate_all_kapa(const QTunnel& QTn)
{
    double kapa1,kapa2,alpha,beta;
    alpha=QTn.Calculate_Alpha();
    beta=QTn.Calculate_Beta();
    kapa1=QTn.Calculate_Kapa(false);
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算量子隧穿效应透过系数结果*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*Wigner方法计算结果*"<<std::endl;
    std::cout<<"透过系数Kapa="<<kapa1<<std::endl;
    std::cout<<std::endl;
    std::cout<<"注意：Wigner方法计算精度较低，在非高温情况下通常会造成重大误差"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*Skodje-Truhlar方法计算结果*"<<std::endl;
    std::cout<<"Alpha="<<alpha<<std::endl;
    std::cout<<"Beta="<<beta<<std::endl;
    kapa2=QTn.Calculate_Kapa(true);
    std::cout<<"透过系数Kapa="<<kapa2<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}




//自由能方法计算单分子反应速率常数模块
void calculate_single_freeEnergy_k(const Dynamic_Single& DSn)
{
    double k, rate_kapa, t;
    k=DSn.Calculate_k(false);
    rate_kapa=DSn.ContributionRate_Kapa();
    t=DSn.Calculate_Halftime(false);
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算单分子反应速率常数结果*"<<std::endl;
    std::cout<<"\t\t*速率常数通过自由能计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A->TS->P"<<std::endl;
    std::cout<<"速率常数="<<k<<" s-1"<<std::endl;
    std::cout<<"反应物A半衰期="<<t<<" s"<<std::endl;
    std::cout<<"量子隧穿效应贡献率="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST)/s-1="<<k<<"*[A]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}




//配分函数方法计算单分子反应速率常数模块
void calculate_single_Q_k(const Dynamic_Single& DSn)
{
    double k, rate_kapa, t;
    k=DSn.Calculate_k(true);
    rate_kapa=DSn.ContributionRate_Kapa();
    t=DSn.Calculate_Halftime(true);
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算单分子反应速率常数结果*"<<std::endl;
    std::cout<<"\t\t*速率常数通过配分函数计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A->TS->P"<<std::endl;
    std::cout<<"速率常数="<<k<<" s-1"<<std::endl;
    std::cout<<"反应物A半衰期="<<t<<" s"<<std::endl;
    std::cout<<"量子隧穿效应贡献率="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST)/s-1="<<k<<"*[A]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}




//全方法计算单分子反应速率常数模块
void calculate_single_all_k(const Dynamic_Single& DSn)
{
    double k1,k2,rate_kapa,t1,t2;
    k1=DSn.Calculate_k(false);
    t1=DSn.Calculate_Halftime(false);
    k2=DSn.Calculate_k(true);
    rate_kapa=DSn.ContributionRate_Kapa();
    t2=DSn.Calculate_Halftime(true);
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算单分子反应速率常数结果*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*速率常数通过自由能计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A->TS->P"<<std::endl;
    std::cout<<"速率常数(G)="<<k1<<" s-1"<<std::endl;
    std::cout<<"反应物A半衰期(G)="<<t1<<" s"<<std::endl;
    std::cout<<"量子隧穿效应贡献率="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST,G)/s-1="<<k1<<"*[A]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*速率常数通过配分函数计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A->TS->P"<<std::endl;
    std::cout<<"速率常数(Q)="<<k2<<" s-1"<<std::endl;
    std::cout<<"反应物A半衰期(Q)="<<t2<<" s"<<std::endl;
    std::cout<<"量子隧穿效应贡献率="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST,Q)/s-1="<<k2<<"*[A]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*平均值*"<<std::endl;
    std::cout<<"\t*仅供参考，当自由能结果与配分函数结果差异较大时，请谨慎考虑*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A->TS->P"<<std::endl;
    std::cout<<"速率常数="<<(k1+k2)/2<<" s-1"<<std::endl;
    std::cout<<"反应物A半衰期="<<(t1+t2)/2<<" s"<<std::endl;
    std::cout<<"量子隧穿效应贡献率="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST)/s-1="<<(k1+k2)/2<<"*[A]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}




//自由能方法计算双分子反应速率常数模块
void calculate_double_freeEnergy_k(const Dynamic_Double& DDn)
{
    double k, rate_kapa;
    k=DDn.Calculate_k(false);
    rate_kapa=DDn.ContributionRate_Kapa();
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算双分子反应速率常数结果*"<<std::endl;
    std::cout<<"\t\t*速率常数通过自由能计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A+B->TS->P"<<std::endl;
    std::cout<<"速率常数="<<k<<" s-1*(mol/L)-1"<<std::endl;
    std::cout<<"量子隧穿效应贡献率="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST)/s-1*(mol/L)-1="<<k<<"*[A]*[B]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}




//配分函数方法计算双分子反应速率常数模块
void calculate_double_Q_k(const Dynamic_Double& DDn)
{
    double k, rate_kapa;
    k=DDn.Calculate_k(true);
    rate_kapa=DDn.ContributionRate_Kapa();
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算双分子反应速率常数结果*"<<std::endl;
    std::cout<<"\t\t*速率常数通过配分函数计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A+B->TS->P"<<std::endl;
    std::cout<<"速率常数="<<k<<" s-1*(mol/L)-1"<<std::endl;
    std::cout<<"量子隧穿效应贡献率="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST)/s-1*(mol/L)-1="<<k<<"*[A]*[B]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}





//全方法计算双分子反应速率常数模块
void calculate_double_all_k(const Dynamic_Double& DDn)
{
    double k1,k2,rate_kapa;
    k1=DDn.Calculate_k(false);
    k2=DDn.Calculate_k(true);
    rate_kapa=DDn.ContributionRate_Kapa();
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算双分子反应速率常数结果*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*速率常数通过自由能计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A+B->TS->P"<<std::endl;
    std::cout<<"速率常数(G)="<<k1<<" s-1*(mol/L)-1"<<std::endl;
    std::cout<<"量子隧穿效应贡献率="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST,G)/s-1*(mol/L)-1="<<k1<<"*[A]*[B]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*速率常数通过配分函数计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A+B->TS->P"<<std::endl;
    std::cout<<"速率常数(Q)="<<k2<<" s-1*(mol/L)-1"<<std::endl;
    std::cout<<"量子隧穿效应贡献率="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST,Q)/s-1*(mol/L)-1="<<k2<<"*[A]*[B]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*平均值*"<<std::endl;
    std::cout<<"\t*仅供参考，当自由能结果与配分函数结果差异较大时，请谨慎考虑*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A+B->TS->P"<<std::endl;
    std::cout<<"速率常数="<<(k1+k2)/2<<" s-1*(mol/L)-1"<<std::endl;
    std::cout<<"量子隧穿效应贡献率="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST)/s-1*(mol/L)-1="<<(k1+k2)/2<<"*[A]*[B]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}

