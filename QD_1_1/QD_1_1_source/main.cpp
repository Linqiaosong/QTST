//Author Qiaosong Lin, Wuhan University, 2019 
//2019-08-09    完善了UI



#include<iostream>
#include<cmath>
#include"QDynamic.h"
#include"QDApplication.h"

int main()
{
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t  QD v1.1 量子动力学计算软件"<<std::endl;
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