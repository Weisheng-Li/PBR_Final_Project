#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <dirt/common.h>

/* FGD structure to linearly interpolate the 4D datafile.
 */
struct FGD {
    FGD() {}
    FGD(const std::string path) {
        std::ifstream in(path);
        // Warn if not loaded
        if(in.bad() || in.fail()) {
            std::cout << "Unable to load FGD file." << std::endl;
        } else {
            std::cout << "Loading FGD file: " << path << std::endl;
        }
        this->load(in);
    }

    void load(std::ifstream& in) {
        // Read sizes
        in.read((char*)&Nt, sizeof(int));
        in.read((char*)&Na, sizeof(int));
        in.read((char*)&Nn, sizeof(int));
        in.read((char*)&Nk, sizeof(int));
        std::cout << "Loading FGD texture of dimension " 
        << Nt << "x" << Na << "x" << Nn << "x" << Nk << std::endl;
        sizes[0] = Nt;
        sizes[1] = Na;
        sizes[2] = Nn;
        sizes[3] = Nk;

        // Read data range (min / max)
        float mM[8];
        in.read((char*)mM, 8*sizeof(float));
        tm = mM[0]; tM = mM[1];
        am = mM[2]; aM = mM[3];
        nm = mM[4]; nM = mM[5];
        km = mM[6]; kM = mM[7];
        std::cout << "Range [" 
            << tm << "," << tM << ";" 
            << am << "," << aM << ";"
            << nm << "," << nM << ";"
            << km << "," << kM << ";"
        << "]" << std::endl;

        // Read data
        const int size = Nt*Na*Nn*Nk;
        buff.assign(size, 0.0f);
        in.read((char*)buff.data(), size*sizeof(float));
    }

    inline float operator() (float t, float a, float n, float k) const {

       // Floating point index
       float ta = Nt * (t - tm) / (tM -tm);
       float aa = Na * (a - am) / (aM -am);
       float na = Nn * (n - nm) / (nM -nm);
       float ka = Nk * (k - km) / (kM -km);

       // Integer index
       int ti = floor(ta);
       int ai = floor(aa);
       int ni = floor(na);
       int ki = floor(ka);

       // Ensure the indexes stays in the limits
       ti = clamp<int>(ti, 0, Nt-1);
       ai = clamp<int>(ai, 0, Na-1);
       ni = clamp<int>(ni, 0, Nn-1);
       ki = clamp<int>(ki, 0, Nk-1);

       //*
       // Clamp the interpolation weights
       float alphas[4] = {ta-ti, aa-ai, na-ni, ka-ki};
       for(int i=0; i<4; ++i) {
          alphas[i] = clamp(alphas[i], 0.0f, 1.0f);
       }

       // Index of the middle point
       const int indices[4] = {ti, ai, ni, ki};
       const int index = ki + Nk*(ni + Nn*(ai + Na*ti));

       // Result vector
       float v = 0.0f;

       // For every possible combinaison of index shift per dimension,
       // fetch the value in memory and do linear interpolation.
       // We fetch using shift of 0 and 1.
       //
       //     v(i+di, j+di, k+dk, l+dl),  where dk in [0,1]
       //
       const unsigned int D = pow(2, 4);
       for(unsigned int d=0; d<D; ++d) {

          float alpha = 1.0; // Global alpha
          int   cid_s = 0;   // Id shift

          // Evaluate the weight of the sample d which correspond to
          // one for the shifted configuration:
          // The weight is the product of the weights per dimension.
          //
          for(int i=0; i<4; ++i) {
             bool  bitset = ((1 << i) & d);
             float calpha =  (bitset) ? alphas[i] : 1.0-alphas[i];

             // Correct the shift to none if we go out of the grid
             if(indices[i] + 1 >= sizes[i]) {
                bitset = false;
             }

             alpha *= calpha;
             cid_s = cid_s*sizes[i] + ((bitset) ? 1 : 0);
          }

          v += alpha * buff[index+cid_s];
       }
       return v;
    }

    std::vector<float> buff;
    int sizes[4];
    int Nt, Na, Nn, Nk;
    float tm, tM, aM, am, nm, nM, km, kM;
};

