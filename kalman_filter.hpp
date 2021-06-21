//
// Created by DELL on 2021/4/20.
//

#ifndef INS_KALMANFILTER_H
#define INS_KALMANFILTER_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include <iostream>

using namespace std;

using namespace Eigen;

template<unsigned state_dim, unsigned sys_noise_dim, unsigned measure_dim>
class KalmanFilter {
private:

public:
    typedef Matrix<double, state_dim, state_dim> StateTrans;
    typedef Matrix<double, state_dim, sys_noise_dim> SysNoise2State;
    typedef Matrix<double, measure_dim, state_dim> State2Measure;

    typedef Matrix<double, measure_dim, 1> Measure;
    typedef Matrix<double, state_dim, state_dim> Cov;
    typedef Matrix<double, state_dim, 1> State;

    typedef Matrix<double, sys_noise_dim, sys_noise_dim> SysNoiseCov;
    typedef Matrix<double, measure_dim, measure_dim> MeasureNoiseCov;
    SysNoiseCov sys_noise_cov;
    MeasureNoiseCov measure_noise_cov;

    double ts;

    KalmanFilter(const SysNoiseCov &sys_noise_cov_, const MeasureNoiseCov &measure_noise_cov_, const double &ts_) {
        sys_noise_cov = sys_noise_cov_;
        measure_noise_cov = measure_noise_cov_;
        ts = ts_;
    }

    void updateWithMeasure(
            const StateTrans &state_trans,
            const SysNoise2State &sys_noise_trans,
            const State2Measure &state2measure,
            const Measure &measure,
            Cov &cov,
            State &x) {

        StateTrans E = StateTrans::Identity();
        StateTrans state_trans_discrete = E + ts * state_trans;
        SysNoise2State sys_noise_trans_discrete = (E + 0.5 * ts * state_trans) * sys_noise_trans * ts;
        State x_predict = state_trans_discrete * x;
        Cov cov_predict = state_trans_discrete * cov * state_trans_discrete.transpose()
                          + sys_noise_trans_discrete * sys_noise_cov * sys_noise_trans_discrete.transpose();
//        cout<<"state_trans_discrete"<<state_trans_discrete.norm()<<endl;
//        cout<<"cov_predict"<<cov_predict.norm()<<endl;
        MatrixXd kalman_gain = cov_predict * state2measure.transpose() *
                               (state2measure * cov_predict * state2measure.transpose() + measure_noise_cov).inverse();
//        cout<<"kalman_gain"<<kalman_gain.norm()<<endl;
        cov = (Cov::Identity() - kalman_gain * state2measure) * cov_predict;
        cov = (cov + cov.transpose()) / 2;
//        cout << "covirant norm \t" << cov.norm() << endl;
        x = x + kalman_gain * (measure - state2measure * x_predict);
    }

    void updateWithoutMeasure(
            const StateTrans &state_trans,
            const SysNoise2State &sys_noise_trans,
            Cov &cov,
            State &x) {
        StateTrans E = StateTrans::Identity();
        StateTrans state_trans_discrete = E + ts * state_trans;
//        cout<<"state_trans_discrete\n";
//        print(state_trans_discrete);
        SysNoise2State sys_noise_trans_discrete = (E + 0.5 * ts * state_trans) * sys_noise_trans * ts;
        x = state_trans_discrete * x;
//        cout<<"cov init\n";
//        print(cov);
//        cout<<"old cov norm\t"<<cov.norm()<<endl;
        cov = state_trans_discrete * cov * state_trans_discrete.transpose()
              + sys_noise_trans_discrete * sys_noise_cov * sys_noise_trans_discrete.transpose();
//        cout<<"cov end\n";
//        print(cov);
//        cout<<"new cov norm\t"<<cov.norm()<<endl;
    }

    void print(Cov cov){
        for (int i = 0; i < state_dim; ++i) {
            for(int j=0;j<state_dim;++j){
                if(cov(i,j)<0.001){
                    cov(i,j)=0;
                }
            }
        }
        cout<<cov<<endl;
    }
};


#endif //INS_KALMANFILTER_H
