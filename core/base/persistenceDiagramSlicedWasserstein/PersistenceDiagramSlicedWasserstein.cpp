#include "PersistenceDiagramUtils.h"
#include <algorithm>
#include <array>
#include <limits>
#include <numeric>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif // TTK_ENABLE_EIGEN

#include <PersistenceDiagramSlicedWasserstein.h>

using namespace ttk;

double PersistenceDiagramSlicedWasserstein::execute(
    const DiagramType &diag1, 
    const DiagramType &diag2, 
    int sampleNumber) {

    double dist = 0.;

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

    for(size_t p = 0; p < thetaList.size(); ++p){

        const double theta = thetaList[p];
        std::vector<std::array<double, 2>> projOnTheta1;
        std::vector<std::array<double, 2>> projOnTheta2;

        projectionOnThetaLine(diag1, proj1, projOnTheta1, theta);
        projectionOnThetaLine(diag2, proj2, projOnTheta2, theta);

        double distOneLine = 0.;

        for(size_t k = 0; k < projOnTheta1.size(); ++k){
            auto &p1 = projOnTheta1[k];
            auto &p2 = projOnTheta2[k];
            const double diffX = p1[0] - p2[0];
            const double diffY = p1[1] - p2[1];
            distOneLine += diffX*diffX + diffY*diffY;
        }

        dist += distOneLine / static_cast<double>(sampleNumber);
    }

    std::vector<ttk::MatchingType> matchings;
    slicedTransport(matchings, diag1, diag2, sampleNumber);


    return dist;
}

void PersistenceDiagramSlicedWasserstein::projectionOnThetaLine(
  const DiagramType &diag,
  const std::vector<std::array<double, 2>> &proj,
  std::vector<std::array<double, 2>> &projOnTheta,
  double theta){

    std::array<double, 2> vecUnit{cos(theta), sin(theta)};

    for(size_t j = 0; j < diag.size(); ++j){
        auto &pair = diag[j];
        double birth = pair.birth.sfValue;
        double death = pair.death.sfValue;
        std::array<double, 2> temp{(birth*vecUnit[0] + death*vecUnit[1])*vecUnit[0], 
                (birth*vecUnit[0] + death*vecUnit[1])*vecUnit[1]};
        projOnTheta.emplace_back(temp);
    }

    for(size_t j = 0; j < proj.size(); ++j){
        auto &pair = proj[j];
        double birth = pair[0];
        double death = pair[1];
        std::array<double, 2> temp{(birth*vecUnit[0] + death*vecUnit[1])*vecUnit[0], 
                (birth*vecUnit[0] + death*vecUnit[1])*vecUnit[1]};
        projOnTheta.emplace_back(temp);
    }   

    std::sort(projOnTheta.begin(), projOnTheta.end(), 
        [](std::array<double, 2> p1, std::array<double, 2> p2) 
        { 
            return (p1[0] < p2[0]);
        });


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


#ifdef TTK_ENABLE_EIGEN
    while(epoch < EPOCH_MAX && gradNorm > 1e-2 ){

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

            projectionOnThetaLine(limitMeasure ,projOnTheta1, theta);
            projectionOnThetaLine(diag2, proj2, projOnTheta2, theta);


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
        if(abs(pLimit[0] - pLimit[1]) < 1e-2){
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
        if(abs(pLimit[0] - pLimit[1]) < 1e-2){
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
            test += 1;
        }
    }
    std::cout << " IF ALL CHECKED: " << test << std::endl;
}