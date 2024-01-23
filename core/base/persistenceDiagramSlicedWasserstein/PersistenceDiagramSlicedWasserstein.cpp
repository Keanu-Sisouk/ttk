#include "Geometry.h"
#include "PersistenceDiagramUtils.h"
#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <random>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif // TTK_ENABLE_EIGEN

#include <PersistenceDiagramSlicedWasserstein.h>

using namespace ttk;

double PersistenceDiagramSlicedWasserstein::execute(
    const DiagramType &diag1, 
    const DiagramType &diag2, 
    int &sampleNumber,
    int maxSampleNb) {

    double dist = 0.;

    std::vector<std::array<double, 2>> proj1;
    std::vector<std::array<double, 2>> proj2;

    augmentDiagram(diag1, proj2);
    augmentDiagram(diag2, proj1);

    if(useQuasiMC){
        dist = quasiMonteCarlo(diag1, diag2, proj1, proj2,sampleNumber, maxSampleNb);
    } else {
        dist = classicMonteCarlo(diag1, diag2, proj1, proj2, sampleNumber, maxSampleNb);
    }

    return dist;
}

void PersistenceDiagramSlicedWasserstein::projectionOnThetaLine(
  const DiagramType &diag,
  const std::vector<std::array<double, 2>> &proj,
  std::vector<std::array<double, 2>> &projOnTheta,
  std::vector<int> &originIndices,
  std::vector<double> &scalarProd,
  double theta,
  bool vertical){

    std::array<double, 2> vecUnit{cos(theta), sin(theta)};
    size_t counter = 0;
    for(size_t j = 0; j < diag.size(); ++j){
        auto &pair = diag[j];
        double birth = pair.birth.sfValue;
        double death = pair.death.sfValue;
        double scal = birth*vecUnit[0] + death*vecUnit[1];
        std::array<double, 2> temp{scal*vecUnit[0], 
                scal*vecUnit[1]};
        // projOnTheta.emplace_back(temp);
        projOnTheta[j] = temp;
        // scalarProd.emplace_back(scal);
        scalarProd[j] = scal;
        counter += 1;
    }

    for(size_t j = 0; j < proj.size(); ++j){
        auto &pair = proj[j];
        double birth = pair[0];
        double death = pair[1];
        double scal = birth*vecUnit[0] + death*vecUnit[1];

        std::array<double, 2> temp{scal*vecUnit[0], 
                scal*vecUnit[1]};
        // projOnTheta.emplace_back(temp);
        projOnTheta[counter + j] = temp;
        // scalarProd.emplace_back(scal);
        scalarProd[counter + j] = scal;
    }   


    // originIndices.resize(projOnTheta.size());
    std::iota(originIndices.begin(), originIndices.end(), 0);

    if(vertical != true){

        std::sort(originIndices.begin(), originIndices.end(), 
            [&](int i , int j){
                return (projOnTheta[i][0] < projOnTheta[j][0]);
            });


        // std::sort(projOnTheta.begin(), projOnTheta.end(), 
        //     [](std::array<double, 2> p1, std::array<double, 2> p2) 
        //     { 
        //         return (p1[0] < p2[0]);
        //     });
    } else {

        std::sort(originIndices.begin(), originIndices.end(), 
            [&](int i , int j){
                return (projOnTheta[i][1] < projOnTheta[j][1]);
            });

        // std::sort(projOnTheta.begin(), projOnTheta.end(), 
        //     [](std::array<double, 2> p1, std::array<double, 2> p2) 
        //     { 
        //         return (p1[1] < p2[1]);
        //     });       
    }


}

void PersistenceDiagramSlicedWasserstein::augmentDiagram(
  const DiagramType &diag,
  std::vector<std::array<double, 2>> &proj) {

    for(size_t j = 0; j < diag.size() ; ++j){
        auto &pair = diag[j];
        const double birth = pair.birth.sfValue;
        const double death = pair.death.sfValue;
        const double pointProjected = (birth + death)/2.;
        const std::array<double, 2> temp{pointProjected,pointProjected};
        proj.emplace_back(temp);
    }
}

void PersistenceDiagramSlicedWasserstein::slicedTransport(
    std::vector<ttk::MatchingType> &matchings,
    const DiagramType &diag1,
    const DiagramType &diag2,
    int sampleNumber) {

    const double gradStep = 2.;

    std::vector<double> thetaList(sampleNumber);
    for(int p = 0; p < sampleNumber ; ++p){
        const double theta = static_cast<double>(p) * 
            M_PI / static_cast<double>(sampleNumber);
        thetaList[p] = theta;
    }

    std::vector<std::array<double, 2>> proj1;
    std::vector<std::array<double, 2>> proj2;

    augmentDiagram(diag1, proj2);
    augmentDiagram(diag2, proj1);

    std::vector<std::array<double, 2>> limitMeasure;
    limitMeasure.resize(diag1.size() + proj1.size());

    for(size_t i = 0; i < diag1.size(); ++i){
        const auto &pair = diag1[i];
        const double birth = pair.birth.sfValue;
        const double death = pair.death.sfValue;

        auto &t = limitMeasure[i];
        t[0] = birth;
        t[1] = death;
    }

    for(size_t i = 0; i < proj1.size(); ++i){
        const auto &pair = proj1[i];

        auto &t = limitMeasure[diag1.size()+i];
        t[0] = pair[0];
        t[1] = pair[1]; 

    }


    int epoch = 0;
    double gradNorm = 10.;
    bool vertical = false;

#ifdef TTK_ENABLE_EIGEN
    while(epoch < EPOCH_MAX && gradNorm > 1e-5 ){

        int m = limitMeasure.size();
        Eigen::MatrixXd dummy(m, 2);
        for(int j = 0; j < m; ++j) {
            for(int k = 0; k < 2; ++k) {
                dummy(j, k) = 0.;
            }
        }

        for(size_t p = 0; p < thetaList.size(); ++p){

            const double theta = thetaList[p];
            std::vector<std::array<double, 2>> projOnTheta1;
            std::vector<std::array<double, 2>> projOnTheta2;
            std::vector<int> originIndices;
            std::vector<double> scalarProd1;
            projectionOnThetaLine(limitMeasure ,projOnTheta1, theta);
            projectionOnThetaLine(diag2, proj2, projOnTheta2, originIndices, scalarProd1,theta, vertical);


            std::vector<int> measureProjIndices(limitMeasure.size());
            std::iota(measureProjIndices.begin(), measureProjIndices.end(), 0);

            std::sort(measureProjIndices.begin(), measureProjIndices.end(), 
            [&](int i , int j){
                return (projOnTheta1[i][0] < projOnTheta1[j][0]);
            });


            for(size_t k = 0; k < projOnTheta1.size(); ++k){
                auto &p1 = projOnTheta1[measureProjIndices[k]];
                auto &p2 = projOnTheta2[k];
                const double diffX = p1[0] - p2[0];
                const double diffY = p1[1] - p2[1];

                dummy(measureProjIndices[k],0) += diffX/static_cast<double>(sampleNumber);
                dummy(measureProjIndices[k],1) += diffY/static_cast<double>(sampleNumber);

            }
        }


        for(size_t i = 0; i < limitMeasure.size(); ++i){
            auto &p1 = limitMeasure[i];
            const auto &p2_1 = dummy(i,0);
            const auto &p2_2 = dummy(i,1);

            p1[0] -= gradStep*p2_1;
            p1[1] -= gradStep*p2_2;
        }
        epoch +=1;
        gradNorm = dummy.norm();
    }

    std::cout << "GRAD NORM: " << gradNorm << std::endl;

    getMatchings(matchings, diag1, diag2, proj1, limitMeasure);

#endif
    
    
}


void PersistenceDiagramSlicedWasserstein::projectionOnThetaLine(
    const std::vector<std::array<double, 2>> &limitMeasure,
    std::vector<std::array<double, 2>> &projOnTheta,
    double theta){

    std::array<double, 2> vecUnit{cos(theta), sin(theta)};


    for(size_t j = 0; j < limitMeasure.size(); ++j){
        auto &pair = limitMeasure[j];
        double birth = pair[0];
        double death = pair[1];
        std::array<double, 2> temp{(birth*vecUnit[0] + death*vecUnit[1])*vecUnit[0], 
                (birth*vecUnit[0] + death*vecUnit[1])*vecUnit[1]};
        projOnTheta.emplace_back(temp);
    }

}

void PersistenceDiagramSlicedWasserstein::getMatchings(
    std::vector<MatchingType> &matchings,
    const ttk::DiagramType &diag1,
    const ttk::DiagramType &diag2,
    const std::vector<std::array<double, 2>> &proj,
    const std::vector<std::array<double, 2>> &limitMeasure){

    matchings.resize(limitMeasure.size());
    std::vector<int> checker(diag2.size(), 0);

    for(size_t i = 0 ; i < diag1.size() ; ++i){
        auto p = diag1[i];
        const double birth = p.birth.sfValue;
        const double death = p.death.sfValue;
        auto &pLimit = limitMeasure[i];
        auto &matching = matchings[i];
        if(abs(pLimit[0] - pLimit[1]) < 1e-1){
            std::get<0>(matching) = i;
            std::get<1>(matching) = -1;
            std::get<2>(matching) = pow(birth - pLimit[0], 2) + pow(death - pLimit[0],2);
        } else {
            for(size_t j = 0; j < diag2.size(); ++j){
                auto &pTrue = diag2[j];
                const double birth_limit = pTrue.birth.sfValue;
                const double death_limit = pTrue.birth.sfValue;

                if((pow(birth_limit - pLimit[0], 2) + pow(death_limit - pLimit[1],2) < 5e-1) && (checker[j] != 1)){
                    
                    std::get<0>(matching) = i;
                    std::get<1>(matching) = j;
                    std::get<2>(matching) = pow(birth - birth_limit, 2) + pow(death - death_limit,2);
                    checker[j] = 1;
                }
            }
        }
    }


    for(size_t i = 0; i < proj.size() ; ++i){
        auto p = proj[i];
        const double birth = p[0];
        const double death = p[1];
        auto &pLimit = limitMeasure[diag1.size() + i];
        auto &matching = matchings[diag1.size() + i];
        if(abs(pLimit[0] - pLimit[1]) < 1e-1){
            std::get<0>(matching) = -1;
            std::get<1>(matching) = -1;
            std::get<2>(matching) = 0.;            
        } else {
            for(size_t j = 0; j < diag2.size(); ++j){
                auto &pTrue = diag2[j];
                const double birth_limit = pTrue.birth.sfValue;
                const double death_limit = pTrue.birth.sfValue;

                if((pow(birth_limit - pLimit[0], 2) + pow(death_limit - pLimit[1],2) < 5e-1) && (checker[j] != 1)){
                    
                    std::get<0>(matching) = -1;
                    std::get<1>(matching) = j;
                    std::get<2>(matching) = pow(birth - birth_limit, 2) + pow(death - death_limit,2);
                    checker[j] = 1;
                }
            }           
        }
    }

    int test = 0;
    for(size_t j = 0; j < checker.size(); ++j){
        if(checker[j] == 0){
            std::cout << "WHICH INDEX: " << j << std::endl;
            test += 1;
        }
    }
    std::cout << "IF ALL CHECKED: " << test << std::endl;
}


double PersistenceDiagramSlicedWasserstein::computeNormGradient(
    const ttk::DiagramType &diag1,
    const ttk::DiagramType &diag2,
    const std::vector<std::array<double, 2>> &proj1,
    const std::vector<std::array<double, 2>> &proj2,
    const std::vector<int> &originIndices1,
    const std::vector<int> &originIndices2,
    const std::vector<double> &scalarProd1,
    const std::vector<double> &scalarProd2,
    const double &theta){


    std::array<double, 2> tempArray{0.,0.};
    std::array<double, 2> vecUnit{cos(theta), sin(theta)};

    double result = 0.;

    int overallSize = originIndices1.size();
    int size1 = diag1.size();
    int size2 = diag2.size();

    std::array<double, 2> temp1{0.,0.};
    std::array<double, 2> temp2{0.,0.};

    double scal1;
    double scal2;

    for(int i = 0; i < overallSize; ++i){
        int index1 = originIndices1[i];
        int index2 = originIndices2[i];

        scal1 = scalarProd1[index1];
        scal2 = scalarProd2[index2];

        if(index1 < size1){
            const auto &t = diag1[index1];
            temp1[0] = t.birth.sfValue;
            temp1[1] = t.birth.sfValue;
        } else {
            const auto &t = proj1[index1 - size1];
            temp1[0] = t[0];
            temp1[1] = t[1];
        }

        if(index2 < size2){
            const auto &t = diag2[index2];
            temp2[0] = t.birth.sfValue;
            temp2[1] = t.birth.sfValue;
        } else {
            const auto &t = proj2[index2 - size2];
            temp2[0] = t[0];
            temp2[1] = t[1];
        }

        tempArray[0]+= 2*Geometry::powInt(scal1 - scal2, 2)*vecUnit[0] + 2*(scal1 - scal2)*(temp1[0] - temp2[0]);
        tempArray[1]+= 2*Geometry::powInt(scal1 - scal2, 2)*vecUnit[1] + 2*(scal1 - scal2)*(temp1[1] - temp2[1]);
    }

    // result = Geometry::pow(Geometry::pow(tempArray[0],2) + Geometry::pow(tempArray[1],2), 1./2);
    result = abs(-sin(theta)*tempArray[0] + cos(theta)*tempArray[1]);
    return result;
}


double PersistenceDiagramSlicedWasserstein::classicMonteCarlo(
    const ttk::DiagramType &diag1,
    const ttk::DiagramType &diag2,
    const std::vector<std::array<double, 2>> &proj1,
    const std::vector<std::array<double, 2>> &proj2,
    int &sampleNumber,
    int maxSampleNb){

    
    int sizeDiag1 = diag1.size() + proj1.size();
    int sizeDiag2 = diag2.size() + proj2.size();


    double temp = 0.;
    bool cond = false;
    int n = 0;
    double ortho_angle = M_PI;
    double mean = 0.;
    double mean_sq = 0.;
    double tresh = 0.01;

    double prev_mean = 0.;
    double prev_sd = 0.;

    double tot_variation_f = 0.;

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0., ortho_angle);
    int total_number = 0;
    // bool vertical = false;

    std::vector<double> buffer_angle(nbPoints);


    std::vector<double> temp_mean_array(buffer_angle.size(), 0.);
    std::vector<double> temp_sd_array(buffer_angle.size(), 0.);
    std::vector<double> temp_vf_array(buffer_angle.size(), 0.);


    while (cond == false && n < maxSampleNb) {
        // int nb_sample = static_cast<int>(std::pow(2., static_cast<double>(n)));


        // std::cout << "BUFFER ANGLE SIZE:" << buffer_angle.size() << std::endl;
        // double angle = dis(gen);
        // buffer_angle.emplace_back(angle);
        // buffer_angle.emplace_back(angle + ortho_angle);    

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP   
        for(int i = 0; i < this->threadNumber_; ++i){
            double angle = dis(gen);
            buffer_angle[i] = angle;
        }


        // for(int i = 1; i < nb_sample ; i += 2){
        //     double angle = static_cast<double>(i) * ortho_angle / static_cast<double>(nb_sample);
        //     buffer_angle.emplace_back(angle);
        //     buffer_angle.emplace_back(angle + ortho_angle);
        // }

        total_number += static_cast<int>(buffer_angle.size());




        std::vector<std::vector<std::array<double, 2>>> projOnTheta1(nbPoints);
        std::vector<std::vector<std::array<double, 2>>> projOnTheta2(nbPoints);
        // plus rapide en sequentielle ou OPENMP ici
        for(int t = 0; t < this->threadNumber_; t++){
            projOnTheta1[t].resize(sizeDiag1);
            projOnTheta2[t].resize(sizeDiag2);
        }


#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        // total_number += static_cast<int>(buffer_angle.size());
        for(size_t t = 0; t < buffer_angle.size(); ++t){
            const double theta = buffer_angle[t];
            // std::vector<std::array<double, 2>> projOnTheta1;
            // std::vector<std::array<double, 2>> projOnTheta2;
            std::vector<int> originIndices1;
            std::vector<int> originIndices2;
            std::vector<double> scalarProd1;
            std::vector<double> scalarProd2;

            projectionOnThetaLine(diag1, proj1, projOnTheta1[t], originIndices1, scalarProd1, theta, false);
            projectionOnThetaLine(diag2, proj2, projOnTheta2[t], originIndices2, scalarProd2, theta, false);

            double distOneLine = 0.;

            for(size_t k = 0; k < projOnTheta1.size(); ++k){
                auto &index1 = originIndices1[k];
                auto &index2 = originIndices2[k];
                auto &p1 = projOnTheta1[t][index1];
                auto &p2 = projOnTheta2[t][index2];
                const double diffX = p1[0] - p2[0];
                const double diffY = p1[1] - p2[1];
                distOneLine += diffX*diffX + diffY*diffY;
            }

            temp_mean_array[t] = distOneLine;
            temp_sd_array[t] = std::pow(distOneLine, 2.);
        }
        
        mean +=  std::accumulate(temp_mean_array.begin(), temp_mean_array.end(), 0.);
        // std::cout << "MEAN: " << mean << std::endl;
        mean_sq += std::accumulate(temp_sd_array.begin(), temp_sd_array.end(), 0.);
        // tot_variation_f += std::accumulate(temp_vf_array.begin(), temp_vf_array.end(), 0.);
    
        // double number_temp = 2. * static_cast<double>(nb_sample);
        double number_temp = static_cast<double>(total_number);
        double current_mean = mean/(number_temp);
 
        double current_sd = std::pow((number_temp/(number_temp - 1.)) * (mean_sq/number_temp - std::pow(current_mean, 2.)), 0.5);

        temp = mean/number_temp;
        if( current_sd * 1.96 / std::pow(number_temp, 0.5) < tresh){
            cond = true;
        }

        // if( current_sd * 1.96 / std::pow(number_temp, 0.5) < tresh ||  (abs(current_mean - prev_mean) < 0.01 && abs(current_sd - prev_sd) < 0.01)){
        //     dist = mean/number_temp;
        //     cond = true;
        // }

        n = n + nbPoints;

        
    }
    // this->setNbOfProjused(n);
    sampleNumber += n;
    // std::cout << "WHAT " << temp << std::endl;
    return temp;

}

double PersistenceDiagramSlicedWasserstein::quasiMonteCarlo(
    const ttk::DiagramType &diag1,
    const ttk::DiagramType &diag2,
    const std::vector<std::array<double, 2>> &proj1,
    const std::vector<std::array<double, 2>> &proj2,
    int &sampleNumber,
    int maxSampleNb){

    int oriSize1 = diag1.size();
    int oriSize2 = diag2.size();

    int sizeDiag1 = diag1.size() + proj1.size();
    int sizeDiag2 = diag2.size() + proj2.size();

    double temp = 0.;
    bool cond = false;
    int n = 0;
    double ortho_angle = M_PI;
    double mean = 0.;
    double mean_sq = 0.;
    double tresh = 0.01;

    double prev_mean = 0.;
    // double prev_sd = 0.;

    int counter = 0;
    int anti_counter = 0;

    double tot_variation_f = 0.;

    int total_number = 0;
    // bool vertical = false;
    std::vector<double> buffer_angle(this->threadNumber_);

    std::vector<double> temp_mean_array(buffer_angle.size());
    // std::vector<double> temp_sd_array(buffer_angle.size(), 0.);
    std::vector<double> temp_vf_array(buffer_angle.size());

    std::cout << "THREAD NUMBER " << this->threadNumber_ << std::endl;
    while (cond == false && n < maxSampleNb) {
        
        // int nb_sample = static_cast<int>(std::pow(2., static_cast<double>(n)));

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP        
        for(int i = 0; i < this->threadNumber_; ++i){
            double angle=0, bk=1.0/2;
            int m = n + i;
            while (m > 0) {
                angle += (m % 2)*bk;
                m /= 2;
                bk /= 2;
            }
            angle = ortho_angle*angle;
            // angle = M_PI * 0.25 + angle * ortho_angle;
            // buffer_angle.emplace_back(angle);
            // buffer_angle.emplace_back(angle + ortho_angle);
            buffer_angle[i] = angle;
        }

        total_number += this->threadNumber_;
        // total_number += static_cast<int>(buffer_angle.size());

        std::vector<std::vector<std::array<double, 2>>> projOnTheta1(this->threadNumber_);
        std::vector<std::vector<std::array<double, 2>>> projOnTheta2(this->threadNumber_);

        std::vector<std::vector<int>> originIndices1(this->threadNumber_);
        std::vector<std::vector<int>> originIndices2(this->threadNumber_);

        std::vector<std::vector<double>> scalarProd1(this->threadNumber_);
        std::vector<std::vector<double>> scalarProd2(this->threadNumber_);
        // plus rapide en sequentielle ou OPENMP ici

// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(nbPoints)
// #endif // TTK_ENABLE_OPENMP
        for(int t = 0; t < this->threadNumber_; t++){
            projOnTheta1[t].resize(sizeDiag1);
            projOnTheta2[t].resize(sizeDiag2);
            originIndices1[t].resize(sizeDiag1);
            originIndices2[t].resize(sizeDiag2);
            scalarProd1[t].resize(sizeDiag1);
            scalarProd2[t].resize(sizeDiag2);
        }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        // total_number += static_cast<int>(buffer_angle.size());
        for(int t = 0; t < this->threadNumber_; ++t){
            const double theta = buffer_angle[t];


            bool vert = false;

            if(theta > 0.85*ortho_angle*0.5 && theta < 1.15*ortho_angle*0.5){
                vert = true;
            }

            projectionOnThetaLine(diag1, proj1, projOnTheta1[t], originIndices1[t], scalarProd1[t], theta, vert);
            projectionOnThetaLine(diag2, proj2, projOnTheta2[t], originIndices2[t], scalarProd2[t], theta, vert);

            double distOneLine = 0.;

            for(int k = 0; k < sizeDiag1; ++k){
                double diffX = 0.;
                double diffY = 0.;
                auto &index1 = originIndices1[t][k];
                auto &index2 = originIndices2[t][k];
                auto &p1 = projOnTheta1[t][index1];
                auto &p2 = projOnTheta2[t][index2];
                diffX = p1[0] - p2[0];
                diffY = p1[1] - p2[1];
                distOneLine += diffX*diffX + diffY*diffY;
            }

            // double vf = computeNormGradient(diag1, diag2, proj1, proj2, originIndices1[t], originIndices2[t], scalarProd1[t], scalarProd2[t], theta);
            // double vf = 42.;
            temp_mean_array[t] = distOneLine;
            // temp_sd_array[t] = std::pow(distOneLine, 2.);
            // temp_vf_array[t] = vf;
        }
        
        mean +=  std::accumulate(temp_mean_array.begin(), temp_mean_array.end(), 0.);
        // mean_sq += std::accumulate(temp_sd_array.begin(), temp_sd_array.end(), 0.);
        // tot_variation_f += std::accumulate(temp_vf_array.begin(), temp_vf_array.end(), 0.);
    
        double number_temp = static_cast<double>(total_number);

        temp = mean/number_temp;

        // if (tot_variation_f/(number_temp*number_temp) < tresh){
        //     cond = true;
        // }

        // if (tot_variation_f/(number_temp*number_temp) < tresh ||  abs(temp - prev_mean) < tresh*0.5){
        //     cond = true;
        // }

        // if(temp == 0){
        //     anti_counter +=1;
        // }

        // if(anti_counter > 10){
        //     cond = true;
        // }

        // if (temp > prev_mean){
        //     if( 1 - prev_mean/temp < tresh*0.1){
        //         counter+=1;
        //     }
        // } else {
        //     if (1 - temp/prev_mean < tresh*0.1){
        //         counter+=1;
        //     }
        // }

        // // std::cout << "COUNTER: " << counter << std::endl;

        // if(counter > 100){
        //     cond = true;
        // }

        // if (abs(temp - prev_mean) < tresh){
        //     cond = true;
        // }

        prev_mean = temp;

        n = n + this->threadNumber_;
    }

    // this->setNbOfProjused(n);
    sampleNumber += n;
    return temp;

}