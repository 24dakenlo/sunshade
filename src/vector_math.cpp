/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <stdint.h>
#include <cstdint>
#include <string>
#include <iostream>
#include <math.h>
#include "vector_math.h"

// ノルム（ベクトルの長さ）を求める関数
double norm(VECT a) {
    return(sqrt(a.x*a.x + a.y*a.y + a.z*a.z));
  }

// 内積を求める関数
double scprod(VECT a, VECT b) {
    return(a.x*b.x + a.y*b.y + a.z*b.z);
  }

// ベクトルを成分値で定義する関数
VECT vectdef(double x, double y, double z) {
    VECT ans;
    ans.x=x; ans.y=y; ans.z=z; return(ans);
  }
