//Author Qiaosong Lin, Wuhan University, 2019 
//2019-08-09    修改了部分bug 


#ifndef QDAPPLICATION_H
#define QDAPPLICATION_H

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

#endif  