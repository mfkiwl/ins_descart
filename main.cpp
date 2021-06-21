#include <iostream>
#include "INS.h"
#include <fstream>
#include <iomanip>

#define pi  3.1415926
#define arg_deg  (pi/180)
#define arg_min (arg_deg/60)
#define arg_sec (arg_min/60)
#define deg_per_hour (arg_deg/3600)
#define Re (6378160)
#define var (0.1*deg_per_hour)

int main() {
    cout << setprecision(8);
    Vector3d posi_ini, atti_ini, vec_ini;
    posi_ini << 0.5934119457, 1.8849555922, 2.0000000000;
    atti_ini << 0, 0, 0.3491;
    vec_ini << 0, 0, 0;

    double tmp[18] = {arg_min, arg_min, arg_min, 0.5, 0.5, 0.5, 30.0, 30.0, 30.0, var, var, var, var, var, var, 1.e-3,
                      1.e-3, 1.e-3};
    Cov cov_ini = Cov::Zero();
    for (int i = 0; i < 18; ++i) {
        cov_ini(i, i) = tmp[i] * tmp[i];
    }

    INS ins(posi_ini, vec_ini, atti_ini, Err::Zero(), cov_ini);

    ifstream gps_fin, acc_fin, omega_fin, posi_fin, velc_fin;
    gps_fin.open("../gps");
    acc_fin.open("../acc");
    omega_fin.open("../omega");
    posi_fin.open("../posi");
    velc_fin.open("../vel");

//    auto R01 = AngleAxisd(0.3, Vector3d(1, 0, 0)).matrix();
//    auto R11 = AngleAxisd(0.2, Vector3d(0, 1, 0)).matrix();
//    auto R21 = AngleAxisd(0.1, Vector3d(0, 0, 1)).matrix();
//    auto R=R01*R11*R21;
//    auto vec=R.eulerAngles(1,0,2);
//    auto R0=AngleAxisd(vec(0), Vector3d(1, 0, 0)).matrix();
//    auto R1=AngleAxisd(vec(1), Vector3d(0, 1, 0)).matrix();
//    auto R2=AngleAxisd(vec(2), Vector3d(0, 0, 1)).matrix();
//    cout<<(R-R1*R0*R2)<<endl;
//    return 0;

    double tmp1, tmp2, tmp3;
    Vector3d gps, acc, omega, posi, atti, velc, posi_, atti_, velc_;
    int cnt = 0;
    for (int i = 0; i < 500; ++i) {
        acc_fin >> tmp1 >> tmp2 >> tmp3;
        acc << tmp1, tmp2, tmp3;
        omega_fin >> tmp1 >> tmp2 >> tmp3;
        omega << tmp1, tmp2, tmp3;
        posi = ins.inertialUpdate(acc, omega);
        if (cnt == 9) {
            gps_fin >> tmp1 >> tmp2 >> tmp3;
            gps << tmp1, tmp2, tmp3;
            ins.updateWithGps(gps);
            cnt = 0;
        } else {
            cnt++;
        }
        posi_fin >> tmp1 >> tmp2 >> tmp3;
        posi_ << tmp1, tmp2, tmp3;
        velc_fin >> tmp1 >> tmp2 >> tmp3;
        velc_ << tmp1, tmp2, tmp3;
        if (i % 100 == 0) cout << "culculated position: " << posi.transpose() << "\treal position:  " << ins.sph2iner(posi_).transpose() << endl;
    }
}
