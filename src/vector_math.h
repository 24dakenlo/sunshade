/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   vector_math.h
 * Author: nishida
 *
 * Created on 2017/06/19, 20:13
 */

#ifndef VECTOR_MATH_H
#define VECTOR_MATH_H

#ifdef __cplusplus
extern "C" {
#endif


// 3次元ベクトル構造体を定義する。
struct VECT {
  double x, y, z, t;   // 構造体を定義。(x, y, z)は座標。tは予備変数 (直線上の点として表す時の媒介変数等)

// 2項演算の定義。tについては何もしない。
// 注意! 演算の優先順位（外積が加減より先, など）は保証されない。適宜, ()で優先順位を明示するべし。
  VECT operator+(VECT p) {   // 加法演算子を定義
    VECT sum;
    sum.x = x + p.x;
    sum.y = y + p.y;
    sum.z = z + p.z;
    return sum;
  }

  VECT operator-(VECT p) {    // 減法演算子を定義
    VECT sub;
    sub.x = x - p.x;
    sub.y = y - p.y;
    sub.z = z - p.z;
    return sub;
  }

  VECT operator*(double a) {    // スカラー倍演算子を定義 ... スカラーを後に書く必要あり!
    VECT scl;
    scl.x = a * x;
    scl.y = a * y;
    scl.z = a * z;
    return scl;
  }

  VECT operator/(double a) {    // スカラー割演算子を定義 ... スカラーを後に書く必要あり!(当然だが)
    VECT scl;
    scl.x = x / a;
    scl.y = y / a;
    scl.z = z / a;
    return scl;
  }

  VECT operator^(VECT p) {    // 外積を定義
    VECT vec;
    vec.x = y*p.z - z*p.y;
    vec.y = z*p.x - x*p.z;
    vec.z = x*p.y - y*p.x;
    return vec;
  }
};

double norm(VECT a);
double scprod(VECT a, VECT b);
VECT vectdef(double x, double y, double z);



#ifdef __cplusplus
}
#endif

#endif /* VECTOR_MATH_H */
