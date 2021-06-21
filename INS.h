#include <tuple>
#include "kalman_filter.hpp"
#include <memory>

using namespace std;
using namespace Eigen;
typedef Vector3d Vel;
typedef Vector3d Acc;
typedef Vector3d Posi;
typedef Vector3d Atti;
typedef Vector3d Gyro;
typedef Vector3d Spherical;
#define ERR_DIM (18)
#define SYS_NOISE_DIM (9)
#define MEASURE_DIM (3)
typedef Matrix<double, ERR_DIM, ERR_DIM> Cov;
typedef Matrix<double, ERR_DIM, 1> Err;//atti vec posi gyro_color gyro_markov acc_markov
typedef Matrix<double, MEASURE_DIM, 1> Measure;

typedef Matrix<double, ERR_DIM, ERR_DIM> StateTrans;                //time variant
typedef Matrix<double, ERR_DIM, SYS_NOISE_DIM> SysNoise2Err;        //time invariant
typedef Matrix<double, MEASURE_DIM, ERR_DIM> State2Measure;         //time invariant

typedef Matrix<double, SYS_NOISE_DIM, SYS_NOISE_DIM> SysNoiseCov;
typedef Matrix<double, MEASURE_DIM, MEASURE_DIM> MeasureNoiseCov;

class INS {
private:
//inertia coordinate params
    Matrix3d Rib;
    Vel Vi;
    Posi posi_i;

public:
//kalman
    SysNoise2Err sys_noise_2_err;
    StateTrans state_trans;
    State2Measure state2measure;
    Err err;
    Cov cov;
    SysNoiseCov sys_noise_cov;
    MeasureNoiseCov measure_noise_cov;
    typedef KalmanFilter<ERR_DIM, SYS_NOISE_DIM, MEASURE_DIM> KF;
    shared_ptr<KF> kf_ptr;

private:
//navi 2 coor;
    void setStateTrans(const Acc &acc_i);

private:
//earth params
    const double Re = 6378160;//%地球半径(长半轴)
    const double f = 1 / 298.3;//地球扁率
    const double e = sqrt(2 * f - f * f);//偏心率
    const double e2 = e * e;
    const double Rp = (1 - f) * Re;//短半轴
    const double ep = sqrt(Re * Re - Rp * Rp) / Rp;//第二偏心率,  此处有改动，加号改为减号
//    const double ep2 = ep * ep;
//    const double wie = 7.2921151467e-5;//地球自转角速率
    const double g0 = 9.7803267714;//重力加速度
    const double ts = 0.1;//sample time

private:
    const double Tg = 3600;//陀螺仪Markov过程相关时间
    const double Ta = 1800;//加速度计Markov过程相关时间
    const double measure_noise_std[3] = {1e-5, 1e-5, 10.3986};
    const double sys_noise_std[9] = {1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 9.780e-4, 9.780e-4, 9.780e-4};
public:
    void updateWithGps(const Spherical& sph);

    INS(const Spherical& sph, const Vel &vel, const Atti& atti, const Err& err, const Cov& cov) ;

    Posi inertialUpdate(const Acc& acc_n, Gyro omega_b_ib);

    Posi sph2iner(const Spherical& sph) const;
};