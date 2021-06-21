#include "INS.h"
#include "Eigen/Geometry"
#include <cmath>

Matrix3d euler2Matrix(Vector3d atti) {
    auto R2 = AngleAxisd(atti(2), Vector3d(0, 0, 1)).matrix();
    auto R1 = AngleAxisd(atti(1), Vector3d(0, 1, 0)).matrix();
    auto R0 = AngleAxisd(atti(0), Vector3d(1, 0, 0)).matrix();
    return R1 * R0 * R2;
}

INS::INS(const Spherical &sph, const Vel &vel, const Atti &atti, const Err &err, const Cov &cov) {
    posi_i = sph2iner(sph);

    Matrix3d Rnb = euler2Matrix(atti);
    double latitude = sph(0);
    double longtitude = sph(1);
    Matrix3d Rin = AngleAxisd(-M_PI_2 + longtitude, Vector3d(0, 0, 1)).matrix()
                   * AngleAxisd(-M_PI_2 + latitude, Vector3d(1, 0, 0)).matrix();
    Rib = Rin * Rnb;

    Vi = Rib * vel;

    this->cov = cov;
    this->err = err;

    sys_noise_cov = SysNoiseCov::Zero();
    for (int i = 0; i < SYS_NOISE_DIM; ++i) {
        sys_noise_cov(i, i) = sys_noise_std[i] * sys_noise_std[i];
    }
    measure_noise_cov = MeasureNoiseCov::Zero();
    for (int i = 0; i < MEASURE_DIM; ++i) {
        measure_noise_cov(i, i) = measure_noise_std[i] * measure_noise_std[i];
    }

    sys_noise_2_err.block(9, 0, 3, 3) = Matrix3d::Identity();
    sys_noise_2_err.block(12, 3, 3, 3) = Matrix3d::Identity();
    sys_noise_2_err.block(15, 6, 3, 3) = Matrix3d::Identity();
    state2measure.block(0, 6, 3, 3) = Matrix3d::Identity();

    kf_ptr = make_shared<KF>(sys_noise_cov, measure_noise_cov, ts);

}

Posi INS::sph2iner(const Spherical &sph) const {
    double latitude = sph(0);
    double longtitude = sph(1);
    double height = sph(2);
    double x = (height + Re) * cos(latitude) * cos(longtitude);
    double y = (height + Re) * cos(latitude) * sin(longtitude);
    double z = (height + Re) * sin(latitude);
    return Posi(x, y, z);
}

inline Quaterniond q_add(Quaterniond q1, Quaterniond q2) {
    return Quaterniond(q1.w() + q2.w(), q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z());
}

inline Quaterniond q_scalar_multiplication(Quaterniond q, double scalar) {
    return Quaterniond(q.w() * scalar, q.x() * scalar, q.y() * scalar, q.z() * scalar);
}

void INS::updateWithGps(const Spherical &gps_sph) {
    Posi gps_posi = sph2iner(gps_sph);
    //err=real-cur
    //real=err+cur
    kf_ptr->updateWithMeasure(state_trans, sys_noise_2_err, state2measure, gps_posi - posi_i, cov, err);
//    cout<<"gps\t"<<gps_posi.transpose()<<"\tposi"<<posi_i.transpose()<<endl;
    Atti atti = Rib.eulerAngles(1, 0, 2);
    atti += err.block(0, 0, 3, 1);
//    cout<<err.block(0,0,9,1).transpose()<<endl;
    Rib = euler2Matrix(atti);
    Vi += err.block(3, 0, 3, 1);
    posi_i += err.block(6, 0, 3, 1);
    err.block(0, 0, 9, 1) = Matrix<double, 9, 1>::Zero();
}

Posi INS::inertialUpdate(const Acc &acc_b, Gyro omega_b_ib) {
    Quaterniond q(Rib);
    Quaterniond delta_q(0, omega_b_ib(0), omega_b_ib(1), omega_b_ib(2));
    q = q_add(q, q_scalar_multiplication(q * delta_q, 0.5 * ts));
    q.normalize();
    Rib = q.toRotationMatrix();

    double r = posi_i.norm();
    Acc gracity_i = posi_i / r * (-g0);
    auto next_Vi = Vi + (Rib * acc_b + gracity_i) * ts;
//    cout << "acc\t" << (Rib * acc_b).transpose() << "\t" << gracity_i.transpose() << "norm\t"
//         << (Rib * acc_b + gracity_i).norm() << endl;
    posi_i += (Vi + next_Vi) * ts / 2;
    Vi = next_Vi;

    setStateTrans(Rib * acc_b);
    kf_ptr->updateWithoutMeasure(state_trans, sys_noise_2_err, cov, err);
    return posi_i;
}

Matrix3d vec2AntiSymmetric(Vector3d vec) {
    Matrix3d tmp;
    tmp << 0, -vec(2), vec(1)/**/, vec(2), 0, -vec(0)/**/, -vec(1), vec(0), 0;
    return tmp;
}

void INS::setStateTrans(const Acc &acc_i) {
    state_trans.block(12, 12, 3, 3) = Matrix3d::Identity() / (-Tg);//陀螺仪Markov过程
    state_trans.block(15, 15, 3, 3) = Matrix3d::Identity() / (-Ta);//加速度计Markov过程

    state_trans.block(0, 9, 3, 3) = -Rib;//陀螺仪有色噪声对姿态误差的影响
    state_trans.block(0, 12, 3, 3) = -Rib;//陀螺仪Markov过程对姿态误差的影响
    state_trans.block(3, 15, 3, 3) = Rib;//加速度计扰动对速度误差的影响

    //TODO 这里有没有负号？？？？？？？？
    state_trans.block(3, 0, 3, 3) = (vec2AntiSymmetric(acc_i));//姿态扰动对速度误差的影响
    state_trans.block(6, 3, 3, 3) = Matrix3d::Identity();//速度扰动对位置误差的影响
}

