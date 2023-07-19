#include "PersistenceDiagramUtils.h"
#include <algorithm>
#include <array>
#include <limits>

#include <PersistenceDiagramSlicedWasserstein.h>

using namespace ttk;

double PersistenceDiagramSlicedWasserstein::execute(
    const DiagramType &diag1, 
    const DiagramType &diag2, 
    int sampleNumber) {

    double dist = 0.;

    std::vector<double> thetaList(sampleNumber);
    for(int i = 0; i < sampleNumber ; ++i){
        const double theta = static_cast<double>(i) * 
            M_PI/ static_cast<double>(sampleNumber);
        thetaList[i] = theta;
    }

    std::vector<std::array<double, 2>> proj1;
    std::vector<std::array<double, 2>> proj2;

    augmentDiagram(diag1, proj2);
    augmentDiagram(diag2, proj1);

    for(size_t j = 0; j < thetaList.size(); ++j){

        double theta = thetaList[j];
        std::vector<std::array<double, 2>> projOnTheta1;
        std::vector<std::array<double, 2>> projOnTheta2;

        projectionOnThetaLine(diag1, proj1, projOnTheta1, theta);
        projectionOnThetaLine(diag2, proj2, projOnTheta2, theta);

        double distOneLine = 0.;
        
        for(size_t k = 0; k < projOnTheta1.size(); ++k){
            auto &p1 = projOnTheta1[k];
            auto &p2 = projOnTheta2[k];
            double diffX = p1[0] - p2[0];
            double diffY = p1[1] - p2[1];
            distOneLine += diffX*diffX + diffY*diffY;
        }

        dist += distOneLine / static_cast<double>(sampleNumber);

    }

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
        { return sqrt(p1[0]*p1[0]+p1[1]*p1[1]) < sqrt(p2[0]*p2[0]+ p2[1]*p2[0]);});
    
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