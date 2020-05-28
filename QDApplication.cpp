//Author Qiaosong Lin, Wuhan University, 2019 
//2020-05-28    更新了提示信息 


#include<iostream>
#include"QDynamic.h"
#include"QDApplication.h"


//wigner_κ方法计算量子隧穿效应透射系数模块
void calculate_wigner_kapa(const QTunnel& QTn)
{
    double kapa;
    kapa=QTn.Calculate_Kapa(false);
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算量子隧穿效应透射系数结果*"<<std::endl;
    std::cout<<"\t\t*透射系数通过Wigner方法计算*"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"注意：Wigner方法计算精度较低，在非高温情况下通常会造成重大误差"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"透射系数κ="<<kapa<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}



//skodje_truhlar方法计算量子隧穿效应透射系数模块
void calculate_skodje_truhlar_kapa(const QTunnel& QTn)
{
    double kapa,alpha,beta;
    alpha=QTn.Calculate_Alpha();
    beta=QTn.Calculate_Beta();
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算量子隧穿效应透射系数结果*"<<std::endl;
    std::cout<<"\t\t*透射系数通过Skodje-Truhlar方法计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"α="<<alpha<<" J-1"<<std::endl;
    std::cout<<"β="<<beta<<" J-1"<<std::endl;
    kapa=QTn.Calculate_Kapa(true);
    std::cout<<"透射系数κ="<<kapa<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}




//全方法计算量子隧穿效应透射系数模块
void calculate_all_kapa(const QTunnel& QTn)
{
    double kapa1,kapa2,alpha,beta;
    alpha=QTn.Calculate_Alpha();
    beta=QTn.Calculate_Beta();
    kapa1=QTn.Calculate_Kapa(false);
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算量子隧穿效应透射系数结果*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*Wigner方法计算结果*"<<std::endl;
    std::cout<<"透射系数κ="<<kapa1<<std::endl;
    std::cout<<std::endl;
    std::cout<<"注意：Wigner方法计算精度较低，在非高温情况下通常会造成重大误差"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*Skodje-Truhlar方法计算结果*"<<std::endl;
    std::cout<<"α="<<alpha<<" J-1"<<std::endl;
    std::cout<<"β="<<beta<<" J-1"<<std::endl;
    kapa2=QTn.Calculate_Kapa(true);
    std::cout<<"透射系数κ="<<kapa2<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}




//Gibbs自由能方法计算单分子反应速率常数模块
void calculate_single_freeEnergy_k(const Dynamic_Single& DSn)
{
    double k, rate_kapa, t;
    k=DSn.Calculate_k(false);
    rate_kapa=DSn.ContributionRate_Kapa();
    t=DSn.Calculate_Halftime(false);
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算单分子反应速率常数结果*"<<std::endl;
    std::cout<<"\t\t*速率常数通过Gibbs自由能计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A->TS->P"<<std::endl;
    std::cout<<"速率常数k(TST)="<<k<<" s-1"<<std::endl;
    std::cout<<"反应物A半衰期τ="<<t<<" s"<<std::endl;
    std::cout<<"量子隧穿效应贡献率η="<<rate_kapa<<std::endl;
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
    std::cout<<"速率常数k(TST)="<<k<<" s-1"<<std::endl;
    std::cout<<"反应物A半衰期τ="<<t<<" s"<<std::endl;
    std::cout<<"量子隧穿效应贡献率η="<<rate_kapa<<std::endl;
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
    std::cout<<"\t\t*速率常数通过Gibbs自由能计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A->TS->P"<<std::endl;
    std::cout<<"速率常数k(TST,G)="<<k1<<" s-1"<<std::endl;
    std::cout<<"反应物A半衰期τ(G)="<<t1<<" s"<<std::endl;
    std::cout<<"量子隧穿效应贡献率η="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST,G)/s-1="<<k1<<"*[A]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*速率常数通过配分函数计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A->TS->P"<<std::endl;
    std::cout<<"速率常数k(TST,Q)="<<k2<<" s-1"<<std::endl;
    std::cout<<"反应物A半衰期τ(Q)="<<t2<<" s"<<std::endl;
    std::cout<<"量子隧穿效应贡献率η="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST,Q)/s-1="<<k2<<"*[A]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*平均值*"<<std::endl;
    std::cout<<"\t*仅供参考，当Gibbs自由能结果与配分函数结果差异较大时，请谨慎考虑*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A->TS->P"<<std::endl;
    std::cout<<"速率常数k(TST)="<<(k1+k2)/2<<" s-1"<<std::endl;
    std::cout<<"反应物A半衰期τ="<<(t1+t2)/2<<" s"<<std::endl;
    std::cout<<"量子隧穿效应贡献率η="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST)/s-1="<<(k1+k2)/2<<"*[A]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}




//Gibbs自由能方法计算双分子反应速率常数模块
void calculate_double_freeEnergy_k(const Dynamic_Double& DDn)
{
    double k, rate_kapa;
    k=DDn.Calculate_k(false);
    rate_kapa=DDn.ContributionRate_Kapa();
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*计算双分子反应速率常数结果*"<<std::endl;
    std::cout<<"\t\t*速率常数通过Gibbs自由能计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A+B->TS->P"<<std::endl;
    std::cout<<"速率常数k(TST)="<<k<<" s-1*(mol/L)-1"<<std::endl;
    std::cout<<"量子隧穿效应贡献率η="<<rate_kapa<<std::endl;
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
    std::cout<<"速率常数k(TST)="<<k<<" s-1*(mol/L)-1"<<std::endl;
    std::cout<<"量子隧穿效应贡献率η="<<rate_kapa<<std::endl;
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
    std::cout<<"\t\t*速率常数通过Gibbs自由能计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A+B->TS->P"<<std::endl;
    std::cout<<"速率常数k(TST,G)="<<k1<<" s-1*(mol/L)-1"<<std::endl;
    std::cout<<"量子隧穿效应贡献率η="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST,G)/s-1*(mol/L)-1="<<k1<<"*[A]*[B]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*速率常数通过配分函数计算*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A+B->TS->P"<<std::endl;
    std::cout<<"速率常数k(TST,Q)="<<k2<<" s-1*(mol/L)-1"<<std::endl;
    std::cout<<"量子隧穿效应贡献率η="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST,Q)/s-1*(mol/L)-1="<<k2<<"*[A]*[B]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"\t\t*平均值*"<<std::endl;
    std::cout<<"\t*仅供参考，当Gibbs自由能结果与配分函数结果差异较大时，请谨慎考虑*"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<"反应机理：A+B->TS->P"<<std::endl;
    std::cout<<"速率常数k(TST)="<<(k1+k2)/2<<" s-1*(mol/L)-1"<<std::endl;
    std::cout<<"量子隧穿效应贡献率η="<<rate_kapa<<std::endl;
    std::cout<<"反应动力学方程：r(TST)/s-1*(mol/L)-1="<<(k1+k2)/2<<"*[A]*[B]"<<std::endl;
    std::cout<<"====================================================="<<std::endl;
    std::cout<<std::endl;
}

