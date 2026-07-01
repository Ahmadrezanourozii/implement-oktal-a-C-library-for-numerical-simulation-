#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "oktal/data/GridVector.hpp"
#include "oktal/io/VtkExport.hpp"
#include "oktal/lbm/D3Q19.hpp"
#include "oktal/lbm/LbmKernels.hpp"
#include "oktal/lbm/TaylorGreen.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"

namespace fs = std::filesystem;
using namespace oktal;
using namespace oktal::lbm;

// NOLINTNEXTLINE(bugprone-exception-escape, readability-function-cognitive-complexity)
int main(int argc, char **argv) {
  if (argc != 3) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    std::cerr << "Usage: " << argv[0] << " <refinement-level> <output-directory>\n";
    return 1;
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const size_t refinementLevel = std::stoul(argv[1]);
  if (refinementLevel < 5) {
    std::cerr << "Refinement level must be >= 5\n";
    return 1;
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const fs::path outputDir = argv[2];
  if (!fs::exists(outputDir)) {
    fs::create_directories(outputDir);
  }

  std::cout << "Starting TGV simulation (Level " << refinementLevel << ")\n";

  const TaylorGreen tgv(refinementLevel);
  auto octree = CellOctree::createUniformGrid(tgv.geometry(), refinementLevel);
  auto cells = CellGrid::create(octree)
                   .levels({size_t(refinementLevel)})
                   .neighborhood(D3Q19::CS_NO_CENTER)
                   .periodicityMapper(Torus({true, true, true}))
                   .build();

  D3Q19Lattice pdfs(cells);
  D3Q19Lattice pdfsTmp(cells);
  GridVector<double, 1> rho(cells);
  GridVector<double, 3> u(cells);
  GridVector<double, 1> uError(cells);

  auto rhoView = rho.view();
  auto uView = u.view();
  
  for (auto cell : cells) {
    const size_t i = cell;
    const Vec3D x = cell.center();
    rhoView[i] = tgv.rho(x, 0.0);
    const Vec3D startU = tgv.u(x, 0.0);
    uView[i, 0] = startU[0];
    uView[i, 1] = startU[1];
    uView[i, 2] = startU[2];
  }

  InitializePdfs initPdfs;
  std::ranges::for_each(cells, [&](auto cell) {
    initPdfs(cell, pdfs.view(), rho.const_view(), u.const_view());
  });

  Collide collide(tgv.omega());
  Stream stream;
  ComputeMacroscopicQuantities computeMacroscopic;

  const size_t numSteps = tgv.numberOfTimesteps();
  const double dt = tgv.dt();

  std::cout << "Total Steps: " << numSteps << ", Omega: " << tgv.omega() << "\n";
  const auto startTime = std::chrono::steady_clock::now();

  for (size_t t = 0; t <= numSteps; ++t) {
    const double physTime = double(t) * dt;

    if (t % 50 == 0 || t == numSteps) {
      double l2ErrorNum = 0.0;
      double l2ErrorDenom = 0.0;
      auto uErrorView = uError.view();

      for (auto cell : cells) {
        const size_t i = cell;
        const Vec3D x = cell.center();
        const Vec3D anaU = tgv.u(x, physTime);
        const Vec3D simU{uView[i, 0], uView[i, 1], uView[i, 2]};
        const Vec3D diff = simU - anaU;
        
        const double diffMag = diff.magnitude();
        uErrorView[i] = diffMag;

        l2ErrorNum += diff.sqrMagnitude(); 
        l2ErrorDenom += anaU.sqrMagnitude();
      }
      
      const double relL2 = std::sqrt(l2ErrorNum / l2ErrorDenom);

      const auto now = std::chrono::steady_clock::now();
      const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();
      const double rate = (t > 0) ? double(t) / (double(elapsed) + 1e-9) : 0.0; // steps per second
      const double remaining = (rate > 0) ? double(numSteps - t) / rate : 0.0;

      std::cout << "Step " << t << "/" << numSteps 
                << " RelL2: " << relL2 
                << " ETA: " << remaining << "s\n";

      io::vtk::exportCellGrid(cells, outputDir / ("step" + std::to_string(t) + ".vtkhdf")) 
                              .writeGridVector<double, 1>("rho", rho.const_view())
                              .writeGridVector<double, 3>("u", u.const_view())
                              .writeGridVector<double, 1>("uError", uError.const_view());
                              
      if (t == numSteps) {
          double numX = 0;
          double denX = 0;
          double numY = 0; 
          double denY = 0;
          for (auto cell : cells) {
              const size_t i = cell;
              const Vec3D x = cell.center();
              const Vec3D anaU = tgv.u(x, physTime);
              const double dx = uView[i, 0] - anaU[0];
              const double dy = uView[i, 1] - anaU[1];
              numX += dx*dx; 
              denX += anaU[0]*anaU[0];
              numY += dy*dy; 
              denY += anaU[1]*anaU[1];
          }
           const double Ex = std::sqrt(numX / denX);
           const double Ey = std::sqrt(numY / denY);
           
           std::ofstream errFile(outputDir / "errors.txt");
           errFile << Ex << " " << Ey << "\n";
      }
    }
    
    if (t < numSteps) {
        std::ranges::for_each(cells, [&](auto cell) { 
            collide(cell, pdfs.view(), rho.const_view(), u.const_view()); 
        });
        
        std::ranges::for_each(cells, [&](auto cell) { 
            stream(cell, pdfsTmp.view(), pdfs.const_view()); 
        });
        
        std::swap(pdfs, pdfsTmp);

        std::ranges::for_each(cells, [&](auto cell) { 
             computeMacroscopic(cell, pdfs.const_view(), rho.view(), u.view()); 
        });
    }
  }

  const auto end = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed_seconds = end - startTime;
  std::cout << "Simulation completed in " << elapsed_seconds.count() << "s\n";

  return 0;
}
